#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
XO_DATA_DIR="$1"
FAM_FILE="$2"
PED_FILE="$3"
OUT_DIR="$4"
SCRIPT_DIR="$5"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Check if output directories exist
mkdir -p "${OUT_DIR}"

# Prepare header for output file
# All files have the same header lines so we'll pull the header from the first file in the directory
head -n 1 $(ls ${XO_DATA_DIR} | head -n 1) > ${XO_DATA_DIR}/all_families_total_XO_count.txt
# Combine phenotype for each family into a single column file containing only total counts
for i in $(ls ${XO_DATA_DIR}/*_xo_count.txt)
do
    cat $i | tail -n +2 >> ${XO_DATA_DIR}/all_families_total_XO_count.txt
done

# Add XO phenotypes to PLINK files
# Update FAM file
FAM_PREFIX=$(basename ${FAM_FILE} .fam)
combine_pheno_and_plink_fam.py \
    ${XO_DATA_DIR}/all_families_total_XO_count.txt \
    ${FAM_FILE} \
    ${OUT_DIR}/${FAM_PREFIX}_pheno.fam

# Update PED file
PED_PREFIX=$(basename ${PED_FILE} .ped)
combine_pheno_and_plink_ped.py \
    ${XO_DATA_DIR}/all_families_total_XO_count.txt \
    ${PED_FILE} \
    ${OUT_DIR}/${PED_PREFIX}_pheno.ped
