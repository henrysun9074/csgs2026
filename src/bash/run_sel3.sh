#!/bin/bash
#SBATCH --job-name=csgs_coxp
#SBATCH --partition=schultzlab
#SBATCH --output=/work/hs325/csgs2026/src/logs/csgs_coxp.out
#SBATCH --error=/work/hs325/csgs2026/src/logs/csgs_coxp.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate gsAI

echo "Job started on: $(date)"
echo "Running from directory: $(pwd)"
echo "--------------------------------------------------"

DATA_ROOT="/work/hs325/csgs2026/data"

# Format: dataset_name|genotype_path|phenotype_path

DATASETS=(
"sel_n9|selectedlines/geno/sel_n9_geno.csv|selectedlines/pheno/sel_n9_pheno.csv"
"sel_all|selectedlines/geno/sel_all_geno.csv|selectedlines/pheno/sel_all_pheno.csv"
)

HPT_ITERATIONS=50

for entry in "${DATASETS[@]}"; do

IFS="|" read -r DF_NAME GENO_SUBPATH PHENO_SUBPATH <<< "$entry"

TARGET_OUTDIR="gebvs_coxp_corr/${DF_NAME}"

FULL_GENO_PATH="${DATA_ROOT}/${GENO_SUBPATH}"
FULL_PHENO_PATH="${DATA_ROOT}/${PHENO_SUBPATH}"

echo "Processing Population: ${DF_NAME}"
echo " -> Genotype Input:  ${FULL_GENO_PATH}"
echo " -> Phenotype Input: ${FULL_PHENO_PATH}"
echo " -> Output Directory: /work/hs325/csgs2026/${TARGET_OUTDIR}"
echo ""

python3 ../csgs_predCoxP_corr.py \
    --genofile "${FULL_GENO_PATH}" \
    --phenofile "${FULL_PHENO_PATH}" \
    --outdir "${TARGET_OUTDIR}" \
    --hpt_iter "${HPT_ITERATIONS}" \
    --verbose

EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: ${DF_NAME} failed with exit code ${EXIT_CODE}"
    exit $EXIT_CODE
fi

echo "Finished evaluating dataset: ${DF_NAME}"
echo "--------------------------------------------------"

done

echo "All datasets completed successfully."
echo "Finished at: $(date)"
