#!/bin/bash
#SBATCH --job-name=bam_to_bedgraph
#SBATCH --account=def-jbl
#SBATCH --output=logs/bedgraph_%j.out
#SBATCH --error=logs/bedgraph_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# creation de bedgraphs normalises (CPM) à partir des bams stratifiés par cell type.
#todo : CPM documentation here :
module load bedtools

set -euo pipefail

BAM_DIR="/lustre07/scratch/sbernarr/subset_bam/bam_files"
OUT_DIR="/lustre07/scratch/sbernarr/subset_bam/output_bedgraphs"
GENOME="/lustre07/scratch/sbernarr/hg38.chrom.sizes"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
BIGWIG=true


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

  samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM"\
    && echo "$SAMPLE trié"
  samtools index "$SORTED_BAM"\
    && echo "$SAMPLE indexé"
  MAPPED=$(samtools view -c -F 4 "$SORTED_BAM")

  if [[ "$MAPPED" -eq 0 ]]; then
    echo "$SAMPLE est vide"
    FAIL=$((FAIL+1))
  fi
  #facteur de normalisation, bedtools, scale is 1e6 ?
  SCALE=$(awk "BEGIN {printf \"%.6f\", 1000000 / $MAPPED}")

  #générer le bedgraph w genomecov -bg, -pc (paired ends for cover calculation) -scale (see above)
  #todo should i sanity check the pc-ness of the bams ?
  BEDGRAPH="$OUT_DIR/${SAMPLE}.CPM.bedgraph"
  bedtools genomecov \
    -ibam "$SORTED_BAM"\
    -bg \
    -pc \
    -scale "$SCALE" \
    -g "$GENOME" \
    | sort -k1,1 -k2,2n \
    > "$BEDGRAPH" \
    && echo " bg généré : $(" $SAMPLE")" \
    || { echo "X, bg fail for $SAMPLE"; }
  if [[ "$BIGWIG" == true ]]; then
    BIGWIG="OUT_DIR/${SAMPLE}.CPM.bw"
    bedGraphToBigWig "$BEDGRAPH" "$GENOME" "$BIGWIG" \
    && echo "bw généré: $SAMPLE" \
    || { echo "X, bw fail for $SAMPLE"; FAIL=$((FAIL+1)); }
  fi
done