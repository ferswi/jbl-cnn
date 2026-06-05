#!/bin/bash
#SBATCH --time=0-6:00:00
#SBATCH --account=def-jbl
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --job-name=slurm_slurm
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca
#SBATCH --array=0-3

#"""
# THIS SCRIPT SERVES TO RUN CHROMBPNET ON THE SC REPLICATES, THEN COMPARE IT TO THE TEST RUN.
#"""

#todo FIX THESE FKIN PATHS AAAH
# input this command in the thing, thats what i have to run i think
# ./subset-bam_linux -b /net/shendure/vol10/projects/david/seq02_DL_T7_CMV/12_07_21_CR_final_run/T7_CMV/outs/possorted_genome_bam.bam \
  #-c cBC_clone_26_20220112.tsv
  # -o bam_10x_DL_clone_26_20220112.bam

#todo: prendre les possorted_bam.bam extraire les données pour les cellules correspondante via script 10x genomics 1

SUBSET_BAM=/lustre07/scratch/sbernarr/subset_bam/
INPUT_BAM=/net/shendure/vol10/projects/david/seq02_DL_T7_CMV/12_07_21_CR_final_run/T7_CMV/outs/possorted_genome_bam.bam
BARCODE_DIR=/lustre07/scratch/sbernarr/subset_bam
OUT_DIR=/lustre07/scratch/sbernarr/subset_bam
export TMPDIR=/lustre07/scratch/sbernarr/subset_bam/
mkdir -p "${OUT_DIR}" logs "${TMPDIR}"


#todo: customise for each cell type.
CELL_TYPES=("ex_endo_parietal"
    "pluripotent_epiblast")
CELL_TYPE=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}
BARCODE_TSV="${BARCODE_DIR}/${CELL_TYPE}.tsv"
OUT_BAM="${OUT_DIR}/${CELL_TYPE}.bam"

# RUN subset-bam
${SUBSET_BAM} \
    --bam        "${INPUT_BAM}" \
    --cell-barcodes "${BARCODE_TSV}" \
    --out-bam    "${OUT_BAM}" \
    --cores      "${SLURM_CPUS_PER_TASK}" \
    --log-level  info

EXIT_CODE=$?

if [[ ${EXIT_CODE} -ne 0 ]]; then
    echo "ERROR: subset-bam failed with exit code ${EXIT_CODE}" >&2
    exit ${EXIT_CODE}
fi

