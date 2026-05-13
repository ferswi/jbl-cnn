#!/bin/bash
#SBATCH --time=0-10:00:00
#SBATCH --account=def-jbl
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=samtools_sort
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

module load samtools
module load bedtools

for fold in 0 1 2 3 4; do
    OUTPUT="fold_${fold}_blacklist"
    FOLD_FILE="/lustre07/scratch/sbernarr/data/splits/fold_${fold}.json"
    echo "Running fold $fold: output=$OUTPUT, fold=$FOLD_FILE"
    chrombpnet prep nonpeaks \
        -g /lustre07/scratch/sbernarr/data/downloads/hg38.fa \
        -o "$OUTPUT" \
        -p /lustre07/scratch/sbernarr/data/peaks_no_blacklist.bed \
        -c /lustre07/scratch/sbernarr/data/downloads/hg38.chrom.sizes \
        -fl "$FOLD_FILE"
done

# prep splits -c ~/scratch/data/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ~/scratch/data/splits/fold_0

#samtools index ~/scratch/data/downloads/merged.bam

#samtools sort -@8 ~/scratch/data/downloads/merged_unsorted.bam \
#              -o ~/scratch/data/downloads/merged.bam