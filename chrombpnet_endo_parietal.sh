#!/bin/bash
#SBATCH --time=0-2:00:00
#SBATCH --account=def-jbl
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --job-name=slurm_slurm
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca

module --force purge
module load StdEnv/2023
module load apptainer/1.4.5
module load bedtools
module load python/3.11.5
module load scipy-stack
module load samtools/1.18
module load picard/3.1.0

source ~/envs/chrombpnet/bin/activate

#apptainer build $SCRATCH/chrombpnet.sif docker-archive:$SCRATCH/chrombpnet.tar

BASE_SC2=/lustre07/scratch/sbernarr/sc2types_cbpn
BASE_SUBSET=/lustre07/scratch/sbernarr/subset_bam/output
REF=$BASE_SC2/ref

#todo : one time troublesome downloads:#-------------------------------------------------------------------------
if [ ! -f $REF/hg38.standard.chrom.sizes ]; then
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes \
        -O $REF/hg38.chrom.sizes
    grep -E "^chr([0-9]+|X|Y)\s" $REF/hg38.chrom.sizes > $REF/hg38.standard.chrom.sizes
fi

if [ ! -f $REF/hg38-blacklist.v2.bed.gz ]; then
    wget -q https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz \
        -O $REF/hg38-blacklist.v2.bed.gz
fi
#-------------------------------------------------------------------------

# pre-processing : merge - sort - index
samtools merge -f $BASE_SC2/parietal_merged_unsorted.bam \
    $BASE_SUBSET/EBD21_crispri_1a_ex_endo_parietal.bam \
    $BASE_SUBSET/EBD21_crispri_1b_ex_endo_parietal.bam \
    $BASE_SUBSET/EBD21_crispri_2a_ex_endo_parietal.bam \
    $BASE_SUBSET/EBD21_crispri_2b_ex_endo_parietal.bam

samtools sort -@4 $BASE_SC2/parietal_merged_unsorted.bam \
    -o $BASE_SC2/parietal_merged.bam
samtools index $BASE_SC2/bam/parietal_merged.bam

samtools view -@ 4 -F 1804 -f 2 -q 30 \
    $BASE_SC2/parietal_merged.bam \
    -o $BASE_SC2/parietal_filtered.bam
samtools index $BASE_SC2/bam/parietal_filtered.bam #todo filter again ?

 #todo; i dont understand this segment as well as i could
picard MarkDuplicates \
    I=$BASE_SC2/parietal_filtered.bam \
    O=$BASE_SC2/parietal_nodup.bam \
    M=$BASE_SC2/parietal_dup_metrics.txt \
    REMOVE_DUPLICATES=true


#shift reads bc of the tn5 correction ? +4bp +strand, -5bp -strand.
alignmentSieve \
    --numberOfProcessors 4 \
    --ATACshift \
    --bam $BASE_SC2/parietal_nodup.bam \
    -o $BASE_SC2/parietal_shifted.bam


samtools sort -@4 $BASE_SC2/parietal_shifted.bam \
    -o $BASE_SC2/parietal_shifted_sorted.bam
samtools index $BASE_SC2/parietal_shifted_sorted.bam

#calling peaks
macs2 callpeak \
    -t $BASE_SC2/parietal_shifted_sorted.bam \
    -f BAMPE \
    -n parietal \
    --outdir $BASE_SC2/peaks/parietal/ \
    -g hs \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    --keep-dup all \
    -B \
    --SPMR \
    --call-summits

#wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
bedtools intersect \
    -v \
    -a $BASE_SC2/peaks/parietal/parietal_peaks.narrowPeak \
    -b $REF/hg38-blacklist.v2.bed.gz \
    > $BASE_SC2/peaks/parietal/parietal_peaks_filtered.narrowPeak


#todo before chrombpnet main pipeline command (ie finir le pre-processing)
# 1. define train, validation and test chromosome splits (standard encode cross validation)
chrombpnet prep splits \
    -op $BASE_SC2/splits/endoparietal_fold \
    -c $BASE_SC2/hg38.standard.chrom.sizes \
    -tcr chr1 \
    -vcr chr8 chr10

# 2.Generate non-peaks (background regions) bed file is coming from the macs2 output
# args here : -k9, 9nr sorts by -log10(qvalue) best peaks first, 2-10 = summit position ?
# ± 250 gives 500 window and 4-9 keep name and score.
sort -k9,9nr \
    $BASE_SC2/peaks/parietal/parietal_peaks_filtered.narrowPeak \
    | awk 'BEGIN{OFS="\t"} {print $1, $2+$10-250, $2+$10+250, $4, $9}' \
    > $BASE_SC2/peaks/parietal/parietal_peaks_final.bed


#CHROMBPNET COMMAND , MAIN PIPELINE-------------------------------------------------------------------------
chrombpnet prep nonpeaks \
    -g $REF/hg38.fa \
    -p $BASE_SC2/peaks/parietal/parietal_peaks_no_blacklist.bed \
    -c $REF/hg38.standard.chrom.sizes \
    -fl $BASE_SC2/splits/endoparietal_fold.json \
    -br $REF/blacklist/hg38-blacklist.v2.bed.gz \
    -o $BASE_SC2/output/parietal/

