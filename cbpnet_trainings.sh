#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --account=def-jbl_gpu
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --job-name=train_parietal
#SBATCH --output=training_%j.out
#SBATCH --error=training_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca

set -e
set -o pipefail

module --force purge
module load StdEnv/2023
module load apptainer/1.4.5

BASE_SC2=/lustre07/scratch/sbernarr/sc2types_cbpn
REF=$BASE_SC2/ref


echo "$(date) sanity check no files are missing"
for f in \
    $BASE_SC2/bam/parietal_shifted_sorted.bam \
    $REF/hg38.fa \
    $REF/hg38.standard.chrom.sizes \
    $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak \
    $BASE_SC2/splits/endoparietal_fold.json \
    $BASE_SC2/bias_model/ATAC_bias.h5; do
    if [ ! -f "$f" ]; then
        echo "missing the followeing output file : $f"
        exit 1
    fi
done

# Confirm nonpeaks output from prep nonpeaks exists
NONPEAKS=$(ls $BASE_SC2/output/parietal/*negatives.bed 2>/dev/null | head -1)
if [ -z "$NONPEAKS" ]; then
    echo "cannot find non peaks files $BASE_SC2/output/parietal/"
    exit 1
fi

echo "$(date) sanity checks for outputs"
echo "  BAM:      $BASE_SC2/bam/parietal_shifted_sorted.bam"
echo "  peaks:    $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak"
echo "  nonpeaks: $NONPEAKS"
echo "  fold:     $BASE_SC2/splits/endoparietal_fold.json"
echo "  bias:     $BASE_SC2/bias_model/ATAC_bias.h5"
echo ""

echo "$(date)  so far so good, starting pipeline : bias model"
if [ ! -f $BASE_SC2/bias_model/parietal/models/bias.h5 ]; then
    apptainer exec --nv $SCRATCH/chrombpnet.sif \
        chrombpnet bias pipeline \
            -ibam $BASE_SC2/bam/parietal_shifted_sorted.bam \
            -d "ATAC" \
            -g $REF/hg38.fa \
            -c $REF/hg38.standard.chrom.sizes \
            -p $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak \
            -n $NONPEAKS \
            -fl $BASE_SC2/splits/endoparietal_fold.json \
            -b 0.5 \
            -o $BASE_SC2/bias_model/parietal/models \
            -fp parietal
else
    echo "$(date) bias mdel was previously established"
fi

#echo "$(date) so far so good, starting pipeline : main model"
#apptainer exec --nv $SCRATCH/chrombpnet.sif \
#    chrombpnet pipeline \
#        -ibam $BASE_SC2/bam/parietal_shifted_sorted.bam \
#        -d "ATAC" \
#        -g $REF/hg38.fa \
#        -c $REF/hg38.standard.chrom.sizes \
#        -p $BASE_SC2/peaks/parietal/parietal_peaks_blacklisted.narrowPeak \
#        -n $NONPEAKS \
#        -fl $BASE_SC2/splits/endoparietal_fold.json \
#        -b $BASE_SC2/bias_model/parietal/models/bias.h5 \
#        -o $BASE_SC2/chrombpnet_model/parietal/
#
#echo "$(date)  full script ran"