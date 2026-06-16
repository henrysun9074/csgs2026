#!/bin/bash
#SBATCH --job-name=csgs_21n2GB
#SBATCH --partition=common
#SBATCH --output=/work/hs325/csgs2026/src/logs/csgs_21n2GB.out
#SBATCH --error=/work/hs325/csgs2026/src/logs/csgs_21n2GB.err
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

# Base paths
DATA_ROOT="/work/hs325/csgs2026/data/21n2"
GENO_PATH="${DATA_ROOT}/sel_21n2_geno.csv"
PHENO_PATH="${DATA_ROOT}/sel_21n2_pheno.csv"
HPT_ITERATIONS=50

# Array of scripts to run
SCRIPTS=("../csgs_predGB.py" "../csgs_predCoxGB.py")

echo "Targeting 21N2 Dataset:"
echo " -> Genotype Input:  ${GENO_PATH}"
echo " -> Phenotype Input: ${PHENO_PATH}"
echo "--------------------------------------------------"

for SCRIPT in "${SCRIPTS[@]}"; do
    # Extract the script name without extension to create a distinct subfolder name
    SCRIPT_NAME=$(basename "$SCRIPT" .py)
    TARGET_OUTDIR="gebvs_21n2/${SCRIPT_NAME}"
    
    echo "Running model: ${SCRIPT_NAME}"
    echo " -> Output Directory: /work/hs325/csgs2026/${TARGET_OUTDIR}"
    echo ""

    python3 "$SCRIPT" \
        --genofile "${GENO_PATH}" \
        --phenofile "${PHENO_PATH}" \
        --outdir "${TARGET_OUTDIR}" \
        --hpt_iter "${HPT_ITERATIONS}" \
        --verbose

    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: ${SCRIPT_NAME} failed with exit code ${EXIT_CODE}"
        exit $EXIT_CODE
    fi

    echo "Finished evaluating model: ${SCRIPT_NAME}"
    echo "--------------------------------------------------"
done

echo "All models completed successfully for 21n2."
echo "Finished at: $(date)"
