#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-jbl
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --job-name=cbbpnet_endo+parietal
#SBATCH --output=cbpnet_ep_%j.out
#SBATCH --error=cbpnet_ep_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca

set -e
set -o pipefail

module --force purge
module load StdEnv/2023
module load apptainer/1.4.5
module load bedtools
module load python/3.11.5
module load scipy-stack
module load samtools/1.18
module load picard/3.1.0
java -jar $EBROOTPICARD/picard.jar

source ~/envs/chrombpnet/bin/activate

#apptainer build $SCRATCH/chrombpnet.sif docker-archive:$SCRATCH/chrombpnet.tar

#variables
BASE_SC2=/lustre07/scratch/sbernarr/sc2types_cbpn
BAM_FILES=/lustre07/scratch/sbernarr/subset_bam/bam_files
REF=$BASE_SC2/ref

#downloads of references if needed
if [ ! -f $REF/hg38.standard.chrom.sizes ]; then
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes \
        -O $REF/hg38.chrom.sizes
    grep -E "^chr([0-9]+|X|Y)\s" $REF/hg38.chrom.sizes > $REF/hg38.standard.chrom.sizes
fi

if [ ! -f $REF/hg38-blacklist.v2.bed.gz ]; then
    wget -q https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz \
        -O $REF/hg38-blacklist.v2.bed.gz
fi

# 1st part of process, merge all the subdivded bam files together.
if [ ! -f $BASE_SC2/bam/parietal_merged.bam ]; then
    echo "$(date) -- merging and sorting BAMs"
    samtools merge -@ 8 - \
        $BAM_FILES/EBD21_crispri_1a_ex_endo_parietal.bam \
        $BAM_FILES/EBD21_crispri_1b_ex_endo_parietal.bam \
        $BAM_FILES/EBD21_crispri_2a_ex_endo_parietal.bam \
        $BAM_FILES/EBD21_crispri_2b_ex_endo_parietal.bam \
        | samtools sort -@ 8 \
            -T $BASE_SC2/bam/parietal_sort_tmp \
            -o $BASE_SC2/bam/parietal_merged.bam
    samtools index $BASE_SC2/bam/parietal_merged.bam
else
    echo "$(date) parietal_merged.bam already exists, skipping merg n sort"
fi

# 2nd part of process, we're filtering hihi
if [ ! -f $BASE_SC2/bam/parietal_filtered.bam ];then
  echo "filtering bams !"
  samtools view -@ 8 -F 1804 -f 2 -q 30 \
      $BASE_SC2/bam/parietal_merged.bam \
      -o $BASE_SC2/bam/parietal_filtered.bam
  samtools index $BASE_SC2/bam/parietal_filtered.bam
else
    echo "$(date) parietal_filtered.bam already exists, skipping filterin"
fi

#Step 3 of pre-pro, removing duplicates
if [ ! -f $BASE_SC2/bam/parietal_nodup.bam ]; then
    echo "$(date) -- removing duplicates"
    java -Xmx32g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I=$BASE_SC2/bam/parietal_filtered.bam \
        O=$BASE_SC2/bam/parietal_nodup.bam \
        M=$BASE_SC2/bam/parietal_dup_metrics.txt \
        REMOVE_DUPLICATES=true
    samtools index $BASE_SC2/bam/parietal_nodup.bam
else
    echo "$(date) parietal_nodup.bam already exists, skippin duplicates removal"
fi

# 4e étape processing : shift reads bc of the tn5 correction ? +4bp +strand, -5bp -strand.
#todo here in running 17/06 2pm
#todo im not sure i understood this step 100%, DOCUMENTATION
if [ ! -f $BASE_SC2/bam/parietal_shifted_sorted.bam ]; then
    echo "$(date) shifting reads"
    alignmentSieve \
        --numberOfProcessors 8 \
        --ATACshift \
        --bam $BASE_SC2/bam/parietal_nodup.bam \
        -o $BASE_SC2/bam/parietal_shifted.bam
    samtools sort -@ 8 \
        -T $BASE_SC2/bam/parietal_shift_tmp \
        $BASE_SC2/bam/parietal_shifted.bam \
        -o $BASE_SC2/bam/parietal_shifted_sorted.bam
    samtools index $BASE_SC2/bam/parietal_shifted_sorted.bam
    rm $BASE_SC2/bam/parietal_shifted.bam
else
    echo "$(date) parietal_shifted_sorted.bam already exists, skipping shift"
fi

#5th step, calling the peaks #todo understand this step better + documentation
if [ ! -f $BASE_SC2/peaks/parietal/parietal_peaks.narrowPeak ]; then
    echo "$(date) calling peaks"
    macs2 callpeak \
        -t $BASE_SC2/bam/parietal_shifted_sorted.bam \
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
else
    echo "$(date) peaks exist, skipping"
fi

#todo docu on this step, needed for quanti qlté via reports that outputs ? columns w q values
#wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
if [ ! -f $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak ]; then
    echo "$(date) filtering blacklist from narrowPeak"
    bedtools intersect \
        -v \
        -a $BASE_SC2/peaks/parietal/parietal_peaks.narrowPeak \
        -b $REF/hg38-blacklist.v2.bed.gz \
        > $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak
else
    echo "$(date) blacklisted narrowPeak exists, skipping"
fi

#todo idem here, more docu pls im dense
# step 7, summit centered bed :
# Generate non-peaks (background regions) bed file is coming from the macs2 output
  ## args here : -k9, 9nr sorts by -log10(qvalue) best peaks first, 2-10 = summit position ?
  ## ± 250 gives 500 window and 4-9 keep name and score.
if [ ! -f $BASE_SC2/peaks/parietal/parietal_peaks_no_blacklist.bed ]; then
    echo "$(date) -- creating summit-centered BED"
    sort -k9,9nr \
        $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak \
        | awk 'BEGIN{OFS="\t"} {
            start = $2 + $10 - 250
            end   = $2 + $10 + 250
            if (start < 0) next
            print $1, start, end, $4, $9
        }' \
        > $BASE_SC2/peaks/parietal/parietal_peaks_no_blacklist.bed
else
    echo "$(date) BED file exists, skipping"
fi


# Step 8 : define train, validation and test chromosome splits (standard encode cross validation)
if [ ! -f $BASE_SC2/splits/endoparietal_fold.json ]; then
    echo "$(date) -- creating chromosome splits"
    chrombpnet prep splits \
        -op $BASE_SC2/splits/endoparietal_fold \
        -c $REF/hg38.standard.chrom.sizes \
        -tcr chr1 \
        -vcr chr8 chr10
else
    echo "$(date) -- splits exist, skipping"
fi


#Step 9 and final of preprocessing : generate nonpeak background regions
if [ ! -f $BASE_SC2/output/parietal/parietal_nonpeaks.bed ]; then
    echo "$(date) -- generating non-peak regions"
    chrombpnet prep nonpeaks \
        -g $REF/hg38.fa \
        -p $BASE_SC2/peaks/parietal/parietal_peaks_no_blacklist.bed \
        -c $REF/hg38.standard.chrom.sizes \
        -fl $BASE_SC2/splits/endoparietal_fold.json \
        -br $REF/hg38-blacklist.v2.bed.gz \
        -o $BASE_SC2/output/parietal/
else
    echo "$(date) -- non-peaks exist, skipping"
fi



#chrombpnet main command
#apptainer exec --nv $SCRATCH/chrombpnet.sif \
#    chrombpnet pipeline \
#        -ibam $BASE_SC2/bam/parietal_shifted_sorted.bam \
#        -d "ATAC" \
#        -g $REF/hg38.fa \
#        -c $REF/hg38.standard.chrom.sizes \
#        -p $BASE_SC2/peaks/parietal/parietal_peaks_no_blacklist.bed \
#        -n $BASE_SC2/output/parietal/parietal_nonpeaks.bed \
#        -fl $BASE_SC2/splits/endoparietal_fold.json \
#        -b $BASE_SC2/bias_model/ATAC_bias.h5 \
#        -o $BASE_SC2/chrombpnet_model/parietal/
