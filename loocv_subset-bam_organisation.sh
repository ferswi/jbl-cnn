#!/bin/bash
#SBATCH --time=0-6:00:00
#SBATCH --account=def-jbl
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --job-name=loocv_subset_bam
#SBATCH --output=logs/%x-%A_fold%a.out
#SBATCH --error=logs/%x-%A_fold%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca
#SBATCH --array=0-3

#"""
#THIS SCRIPT SERVES TO DIVIDE THE FOUR REPLICATES, USING LOOCV (LEAVE ONE OUT CROSS VALIDATION)
#BEFORE RUNNING CBPNET USING THOSE SPLITS
#"""

REPLICATES=(
    "/lustre07/scratch/sbernarr/subset_bam/CellRanger_outs_SR_mEB_10x_scATAC_2020_rep1a.tar"
    "/lustre07/scratch/sbernarr/subset_bam/CellRanger_outs_SR_mEB_10x_scATAC_2020_rep1b.tar"
    "/lustre07/scratch/sbernarr/subset_bam/CellRanger_outs_SR_mEB_10x_scATAC_2020_rep2a.tar"
    "/lustre07/scratch/sbernarr/subset_bam/CellRanger_outs_SR_mEB_10x_scATAC_2020_rep2b.tar"
)


SUBSET_BAM=/lustre07/scratch/sbernarr/subset_bam/subset-bam_linux
TAR_DIR=/home/sbernarr/scratch/subset_bam
EXTRACT_DIR=/lustre07/scratch/sbernarr/subset_bam/extracted
BARCODE_DIR=/lustre07/scratch/sbernarr/subset_bam
OUT_DIR=/lustre07/scratch/sbernarr/subset_bam/output
export TMPDIR=/lustre07/scratch/sbernarr/subset_bam/tmp

CELL_TYPES=(
    "ex_endo_parietal"
    "pluripotent_epiblast"
)

# ------------------------------------------------------------
# ici au début de chaque job, we determine which one will be the test set and which one will be the train.

FOLD=${SLURM_ARRAY_TASK_ID}
TEST_REP="${REPLICATES[$FOLD]}"

# replicates for training
TRAIN_REPS=()
for i in 0 1 2 3; do
    if [[ $i -ne $FOLD ]]; then
        TRAIN_REPS+=("${REPLICATES[$i]}")
    fi
done
#calls so i dont lose my mind
echo "FOLD ${FOLD}: LEFT-OUT (TEST) = ${TEST_REP}"
echo "train replicates: ${TRAIN_REPS[*]}"

mkdir -p "${EXTRACT_DIR}" "${OUT_DIR}" logs "${TMPDIR}"

# random function i found that unpacks a tar, no overwriting + path building
extract_replicate() {
    local tarfile="$1"
    local basename="${tarfile%.tar}"
    local dest="${EXTRACT_DIR}/${basename}"

    if [[ ! -d "${dest}" ]]; then
        echo "[INFO] extraction ${tarfile}..."
        tar -xf "${TAR_DIR}/${tarfile}" -C "${EXTRACT_DIR}"
    else
        echo "[INFO] ${basename} alr extracted "
    fi

    echo "${dest}/outs/possorted_genome_bam.bam"
}

# now we run subset bam on the left out replicate one time per cell type, so twice.
echo "test replicates ! :D "
TEST_BAM=$(extract_replicate "${TEST_REP}")
TEST_BASE="${TEST_REP%.tar}"

for CELL_TYPE in "${CELL_TYPES[@]}"; do
    BARCODE_TSV="${BARCODE_DIR}/${CELL_TYPE}.tsv"
    OUT_BAM="${OUT_DIR}/fold${FOLD}_TEST_${TEST_BASE}_${CELL_TYPE}.bam"
    echo "test, cell type: ${CELL_TYPE}"
    ${SUBSET_BAM} \
        --bam           "${TEST_BAM}" \
        --cell-barcodes "${BARCODE_TSV}" \
        --out-bam       "${OUT_BAM}" \
        --cores         "${SLURM_CPUS_PER_TASK}" \
        --log-level     info

    EXIT_CODE=$?
    if [[ ${EXIT_CODE} -ne 0 ]]; then
        echo "ERROR: subset-bam failed on TEST replicate ${TEST_BASE}, cell type ${CELL_TYPE}" >&2
        exit ${EXIT_CODE}
    fi
done

# same thing for the training replicates
echo "train replicates !"
for TRAIN_REP in "${TRAIN_REPS[@]}"; do
    TRAIN_BAM=$(extract_replicate "${TRAIN_REP}")
    TRAIN_BASE="${TRAIN_REP%.tar}"

    for CELL_TYPE in "${CELL_TYPES[@]}"; do
        BARCODE_TSV="${BARCODE_DIR}/${CELL_TYPE}.tsv"
        OUT_BAM="${OUT_DIR}/fold${FOLD}_TRAIN_${TRAIN_BASE}_${CELL_TYPE}.bam"

        echo "[TRAIN] Replicate: ${TRAIN_BASE} | Cell type: ${CELL_TYPE}"
        ${SUBSET_BAM} \
            --bam           "${TRAIN_BAM}" \
            --cell-barcodes "${BARCODE_TSV}" \
            --out-bam       "${OUT_BAM}" \
            --cores         "${SLURM_CPUS_PER_TASK}" \
            --log-level     info

        EXIT_CODE=$?
        if [[ ${EXIT_CODE} -ne 0 ]]; then
            echo "ERROR: subset-bam failed on TRAIN replicate ${TRAIN_BASE}, cell type ${CELL_TYPE}" >&2
            exit ${EXIT_CODE}
        fi
    done
done
