#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
XO_DATA_DIR="$1"
FINAL_SPLIT_PED_DIR="$2" # Contains split by family (corrected of errors) files
FOUNDERS_PED="$3" # Contains founder lines only
OUT_DIR="$4"
SCRIPT_DIR="$5"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Check if output directories exist
mkdir -p "${OUT_DIR}/gemma_analysis" "${OUT_DIR}/gemma_analysis/plink_pheno_files"

# Prepare header for output file
# All files have the same header lines so we'll pull the header from the first file in the directory
TEMP=$(ls ${XO_DATA_DIR} | head -n 1)
head -n 1 ${XO_DATA_DIR}/${TEMP} > ${XO_DATA_DIR}/all_families_pheno_xo.txt
# Combine phenotype for each family into a single column file containing only total counts
for i in $(ls ${XO_DATA_DIR}/*_pheno.txt)
do
    cat $i | tail -n +2 >> ${XO_DATA_DIR}/all_families_pheno_xo.txt
done

# Build list of final cleaned split PED files
find ${FINAL_SPLIT_PED_DIR} -name "*.ped" | sort -V > ${OUT_DIR}/split_by_family_cleaned_ped_list.txt

# Combine cleaned split PED files into a single PED file
combine_split_ped.py ${FOUNDERS_PED} ${OUT_DIR}/split_by_family_cleaned_ped_list.txt > ${OUT_DIR}/gemma_analysis/all_families.ped

# Update PED file
PED_PREFIX=$(basename ${PED_FILE} .ped)
combine_pheno_and_plink_ped.py \
    ${XO_DATA_DIR}/all_families_pheno_xo.txt \
    ${OUT_DIR}/gemma_analysis/all_families.ped \
    ${OUT_DIR}/gemma_analysis/plink_pheno_files

# Generate FAM file using plink
for i in $(find ${OUT_DIR}/gemma_analysis/plink_pheno_files -name "pheno*.ped" | sort -V)
do
    filename=$(basename ${i} .ped)
    # FAM file is the same as the first six fields in a PED file
    cut -d' ' -f 1-6 ${i} > ${OUT_DIR}/gemma_analysis/plink_pheno_files/${filename}.fam
done
