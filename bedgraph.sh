#!/bin/bash
#SBATCH --job-name=bam_to_bedgraph
#SBATCH --account=def-jbl
#SBATCH --output=logs/bedgraph_%j.out
#SBATCH --error=logs/bedgraph_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# creation de bedgraphs normalises (CPM) à partir des bams stratifiés par cell type.
# to module load before

#todo:path à bien remplir :
BAM_DIR="output"
OUT_DIR="output_bedgraph"
GENOME="hg38.chrom.sizes"
THREADS="${SLURM_CPUS_PER_TASK:-4}"


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

FAIL=0

for BAM in "$BAM_DIR"/*.bam; do
  SAMPLE=$(basename "$BAM" .bam)
  echo "traitement: $SAMPLE"

SORTED_BAM="$OUT_DIR/${SAMPLE}.sorted.bam"
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM"\
  && echo "$SAMPLE trié"
samtools index "$SORTED_BAM"\
  && echo"$SAMPLE indexé"
MAPPED=$(samtools view -c -F 4 "$SORTED_BAM")

if [[ "$MAPPED" -eq 0 ]]; then
  echo  "$SAMPLE est vide"
  FAIL=$((FAIL+1))
fi

#todo unsure of this
#facteur de normalisation, bedtools, scale is 1e6 ?
SCALE=$(awk "BEGIN {printf \"%.6f\", 1000000 / $MAPPED}")

#générer le bedgraph w genomecov -bg, -pc (paired ends for cover calculation) -scale (see above)
BEDGRAPH="$OUT_DIR/${SAMPLE}.CPM.bedgraph"
bedtools genomecov \
-ibam "$SORTED_BAM"\
-bg \
-pc \
-scale "$SCALE"
-g "$GENOME" \
| sort -k1,1 -k2,2n \
> "$BEDGRAPH" \
&& echo " bg généré : $(basename " $SAMPLE")" \
|| { echo "X, bg fail for $SAMPLE"; }
done