#!/bin/bash
#SBATCH --job-name=csgs_gwas
#SBATCH --partition=common
#SBATCH --output=/work/hs325/csgs2026/src/logs/csgs_gwas.out
#SBATCH --error=/work/hs325/csgs2026/src/logs/csgs_gwas.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

# GWAS of survival status: 
# depends: PLINK for QC, GEMMA for association
#
# Inputs: 
# Master.ped + 21.map : raw genotypes (4192 individuals, 65893 SNPs)
# phenotype data from master_gebvs_13nov.csv
### UPDATE LATER WITH COX 

# QC applied sequentially on all 4192 individuals
# SNP call rate        > 90%  (--geno 0.10)
# Individual call rate > 95%  (--mind 0.05)
# MAF                  > 0.01 (--maf 0.01)

# run separately for each pool (selected vs wild) to avoid confounding by line
# kinship computed within each pool

# =============================================================================

set -euo pipefail
module load Plink/1.90
module load GEMMA/0.98.5

DATA_DIR="/work/hs325/csgs2026/data/gwas"
PHENO_CSV="/work/hs325/csgs2026/data/master_gebvs_13nov.csv"
OUT_DIR="/work/hs325/csgs2026/gwas_out"
GEMMA_DIR="${OUT_DIR}/gemma"
mkdir -p "${OUT_DIR}" "${GEMMA_DIR}"
cd "${OUT_DIR}"

THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${THREADS}
PLINK="plink --chr-set 10 --allow-no-sex --threads ${THREADS} --memory 28000"

declare -A POOL_LINES=(
  [selected]="20D4 20N1 21N2 22N9"
  [wild]="21DW 22DBW"
)

# drop MT
grep -v '^#' "${DATA_DIR}/21.map" > master.map

#!!!!!!!!!!!!!!!!!!!!!!
# Master.ped has the "T" allele turned into "TRUE" for
# the 20D4/20N1/21DW/21N2 rows (~813k genotype calls). Checked against
# 21DW21N2_allSNPs.ped: TRUE -> T is the only difference, so restore "T"
awk 'BEGIN { FS = OFS = "\t" } { for (i = 7; i <= NF; i++) if ($i == "TRUE") $i = "T"; print }' \
  "${DATA_DIR}/Master.ped" > master.ped

${PLINK} --ped master.ped --map master.map \
  --chr 1-10 \
  --make-bed --out master_raw

# QC
${PLINK} --bfile master_raw --geno 0.10 --make-bed --out qc1_geno
${PLINK} --bfile qc1_geno   --mind 0.05 --make-bed --out qc2_mind
${PLINK} --bfile qc2_mind   --maf 0.01  --make-bed --out master_qc

# missingness / frequency reports on the final set
${PLINK} --bfile master_qc --missing --freq --out master_qc_stats

# print QC summary
echo "=== QC summary ==="
echo "raw:   $(wc -l < master_raw.bim) SNPs, $(wc -l < master_raw.fam) individuals"
echo "geno:  $(wc -l < qc1_geno.bim) SNPs"
echo "mind:  $(wc -l < qc2_mind.fam) individuals"
echo "final: $(wc -l < master_qc.bim) SNPs, $(wc -l < master_qc.fam) individuals"

# generate phenotype + group files
awk -F',' 'NR > 1 { gsub(/"/, ""); if (!($7 in seen)) { seen[$7]; print $7, $8, $12 } }' \
  "${PHENO_CSV}" > id_group_status.txt

awk 'NR == FNR { grp[$1] = $2; st[$1] = $3; next }
     ($2 in st) { print $1, $2, grp[$2], st[$2] }' \
  id_group_status.txt master_qc.fam > group_status.txt

echo "phenotyped animals passing QC: $(wc -l < group_status.txt)"
awk '{print $3}' group_status.txt | sort | uniq -c

# do GEMMA for both selected and wild groups
for pool in selected wild; do
  lines="${POOL_LINES[$pool]}"
  echo "=== GEMMA: ${pool} (${lines}) ==="

  # phenotyped animals in this pool
  awk -v L="${lines}" 'BEGIN { n = split(L, a, " "); for (i = 1; i <= n; i++) keep[a[i]] }
                       ($3 in keep) { print $1, $2 }' group_status.txt > "keep_${pool}.txt"

  # reapply MAF within the pool
  ${PLINK} --bfile master_qc --keep "keep_${pool}.txt" --maf 0.01 \
    --make-bed --out "${pool}_qc"

  awk 'NR == FNR { st[$2] = $4; next } { $6 = st[$2]; print }' \
    group_status.txt "${pool}_qc.fam" > "${pool}_qc.fam.tmp" && mv "${pool}_qc.fam.tmp" "${pool}_qc.fam"

  awk -v L="${lines}" 'BEGIN { n = split(L, a, " ") }
       NR == FNR { grp[$2] = $3; next }
       { row = "1"; for (i = 2; i <= n; i++) row = row " " (grp[$2] == a[i] ? 1 : 0); print row }' \
    group_status.txt "${pool}_qc.fam" > "${pool}_covar.txt"

  echo "${pool}: $(wc -l < ${pool}_qc.fam) animals, $(wc -l < ${pool}_qc.bim) SNPs"

  # kinship (centred GRM) within the pool
  gemma -bfile "${pool}_qc" -gk 1 -outdir "${GEMMA_DIR}" -o "${pool}_grm"

  # LMM
  gemma -bfile "${pool}_qc" \
    -k "${GEMMA_DIR}/${pool}_grm.cXX.txt" \
    -c "${pool}_covar.txt" \
    -lmm 4 -maf 0.01 -miss 0.1 \
    -outdir "${GEMMA_DIR}" -o "${pool}_status"

  grep -i "pve estimate\|se(pve)" "${GEMMA_DIR}/${pool}_status.log.txt" || true
done

echo "Done. Results: ${GEMMA_DIR}/*_status.assoc.txt"
