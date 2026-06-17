#!/bin/bash
#SBATCH --job-name=bam_to_bedgraph
#SBATCH --account=def-jbl
#SBATCH --output=logs/bedgraph_%j.out
#SBATCH --error=logs/bedgraph_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# creation de bedgraphs normalises (CPM) à partir des bams stratifiés par cell type.

module load samtools/1.20
module load bedtools/2.31.0
module load ucsc_genome_toolkit

set -euo pipefail

BAM_DIR="/lustre07/scratch/sbernarr/subset_bam/bam_files"
OUT_DIR="/lustre07/scratch/sbernarr/subset_bam/output_bedgraphs"
GENOME="/lustre07/scratch/sbernarr/hg38.chrom.sizes"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
BIGWIG=true
FORCE=false


usage() {
  echo "Usage: sbatch $0 --bam-dir <dir> --out-dir <dir> --genome <chrom.sizes>"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir)  BAM_DIR="$2";  shift 2 ;;
    --out-dir)  OUT_DIR="$2";  shift 2 ;;
    --genome)   GENOME="$2";   shift 2 ;;
    --threads)  THREADS="$2";  shift 2 ;;
    --no-bigwig) MAKE_BIGWIG=false; shift ;;
    --force)     FORCE=true; shift ;;
    -h|--help)  usage ;;
  esac
done

mkdir -p "$OUT_DIR" logs

LOGFILE="$OUT_DIR/bedgraph_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

shopt -s nullglob
BAMS=("$BAM_DIR"/*.bam)
if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "Aucun fichier .bam trouvé dans $BAM_DIR"
  exit 1
fi
FAIL=0

for BAM in "$BAM_DIR"/*.bam; do
  SAMPLE=$(basename "$BAM" .bam)
  echo "traitement: $SAMPLE"

  SORTED_BAM="$OUT_DIR/${SAMPLE}.sorted.bam"
  BAI="${SORTED_BAM}.bai"
  BEDGRAPH="$OUT_DIR/${SAMPLE}.CPM.bedgraph"
  BIGWIG_FILE="$OUT_DIR/${SAMPLE}.CPM.bw"

  if [[ "$FORCE" == false && -s "$SORTED_BAM" && -s "$BAI" ]]; then
    echo "$SAMPLE :bam trié avec index found, skip it !"
  else
    samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM"\
    samtools index "$SORTED_BAM"\
      && echo "$SAMPLE done sort and index"
  fi

  MAPPED=$(samtools view -c -F 4 "$SORTED_BAM")

  if [[ "$MAPPED" -eq 0 ]]; then
    echo "$SAMPLE est vide"
    FAIL=$((FAIL+1))
    continue
  fi

  if [[ "$FORCE" == false && -s "$BEDGRAPH" ]]; then
    echo "$SAMPLE : bedgraph déjà présent, on saute genomecov"
  else
    SCALE=$(awk "BEGIN {printf \"%.6f\", 1000000 / $MAPPED}") #todo docu facteur de normalisation, bedtools, scale is 1e6 ?
    bedtools genomecov \
      -ibam "$SORTED_BAM" \
      -bg \
      -pc \
      -scale "$SCALE" \
      | sort -k1,1 -k2,2n \
      > "$BEDGRAPH" \
      && echo "bg généré: $SAMPLE" \
      || { echo "X, bg fail for $SAMPLE"; FAIL=$((FAIL+1)); continue; }
  fi

  if [[ "$MAKE_BIGWIG" == true ]]; then
    if [[ "$FORCE" == false && -s "$BIGWIG_FILE" ]]; then
      echo "$SAMPLE : bigwig déjà présent, on saute bedGraphToBigWig"
    else
      bedGraphToBigWig "$BEDGRAPH" "$GENOME" "$BIGWIG_FILE" \
        && echo "bw généré: $SAMPLE" \
        || { echo "X, bw fail for $SAMPLE"; FAIL=$((FAIL+1)); }
    fi
  fi
done


