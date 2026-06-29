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


#todo, the ones im missing :
echo "  fold:     $BASE_SC2/splits/endoparietal_fold.json"
echo "  bias:     $BASE_SC2/bias_model/ATAC_bias.h5"

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