#!/usr/bin/env bash
#SBATCH -J tQTL_cis
#SBATCH -p standard
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 08:00:00
#SBATCH -o %x.%A_%a.out
#SBATCH -e %x.%A_%a.err

set -euo pipefail

# --- user paths ---
BASE="/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"
PLINK_PREFIX="${BASE}/GxP-eQTL_DNA_genotypes_filtered_SNPs_noChr"
COND_TSV="${BASE}/tensorqtl_conditions.tsv"
OUT_DIR="${BASE}/tensorqtl_output_cis"
mkdir -p "${OUT_DIR}"

# --- environment (edit to match your env) ---
module purge
module load plink/2.0
source activate tensorqtl_p3.11
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH

python3 -m tensorqtl --help >/dev/null 2>&1 || { echo "tensorQTL not on PATH/venv"; exit 1; }

# --- pick the row for this array task ---
# Skip header; get the ${SLURM_ARRAY_TASK_ID}-th data line
line=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==1{next} NR==n+1{print; exit}' "${COND_TSV}") || true
if [[ -z "${line:-}" ]]; then
  echo "No line for task ${SLURM_ARRAY_TASK_ID} in ${COND_TSV}" >&2
  exit 0
fi

IFS=$'\t' read -r CONDITION PHENO_BED COV_FILE <<< "${line}"

# --- outputs ---
OUT_PREFIX="${CONDITION}_cis"

echo "[$(date)] Starting tensorQTL cis"
echo "Condition   : ${CONDITION}"
echo "Genotypes   : ${PLINK_PREFIX}.[pgen|pvar|psam]"
echo "Phenotypes  : ${PHENO_BED}"
echo "Covariates  : ${COV_FILE}"
echo "Out dir     : ${OUT_DIR}"
echo "Out prefix  : ${OUT_PREFIX}"

# --- run ---
python3 -m tensorqtl \
  "${PLINK_PREFIX}" \
  "${PHENO_BED}" \
  "${OUT_PREFIX}" \
  --covariates "${COV_FILE}" \
  --mode cis \
  --fdr 0.1 \
  -o "${OUT_DIR}"

echo "[$(date)] Done."
