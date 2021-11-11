#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
XO_DATA_DIR="$1"
FINAL_SPLIT_PED_DIR="$2" # Contains split by family (corrected of errors) files
FOUNDERS_PED="$3" # Contains founder lines only
MAP_FILE="$4"
OUT_DIR="$5"
SCRIPT_DIR="$6"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Check if output directories exist
mkdir -p "${OUT_DIR}/gemma_analysis"

# Prepare header for output file
# All files have the same header lines so we'll pull the header from the first file in the directory
TEMP=$(ls ${XO_DATA_DIR} | head -n 1)
head -n 1 ${XO_DATA_DIR}/${TEMP} > ${XO_DATA_DIR}/all_families_pheno_xo.txt
# Combine phenotype for each family into a single column file containing only total counts
for i in $(ls ${XO_DATA_DIR}/*_pheno.txt)
do
    cat $i | tail -n +2 >> ${XO_DATA_DIR}/all_families_pheno_xo.txt
done

# Only run GEMMA analysis on families that have phenotype tables
#   successfully generated. Some families may not have phenotype tables
#   generated due to errors occurring.
# Generate list of family names where phenotype table was generated
find ${XO_DATA_DIR} -name "*_pheno.txt" | sort -V | sed -e "s,${XO_DATA_DIR}/,," -e 's,_pheno.txt,,' > ${XO_DATA_DIR}/all_families_pheno_names.txt

# Build list of final cleaned split PED files
find ${FINAL_SPLIT_PED_DIR} -name "*.ped" | sort -V > ${OUT_DIR}/split_by_family_cleaned_ped_list.txt
# Pull out PED files where a phenotype table was successfully generated
grep -f ${XO_DATA_DIR}/all_families_pheno_names.txt ${OUT_DIR}/split_by_family_cleaned_ped_list.txt > ${OUT_DIR}/split_by_family_cleaned_ped_wPheno_list.txt

# Combine cleaned split PED files into a single PED file
combine_split_ped.py ${FOUNDERS_PED} ${OUT_DIR}/split_by_family_cleaned_ped_wPheno_list.txt > ${OUT_DIR}/gemma_analysis/all_families.ped

PLINK_PREFIX=$(basename ${OUT_DIR}/gemma_analysis/all_families.ped .ped)

# Copy MAP file and rename file to match with the all_families.ped file
# Except strip out any characters in the chromosome name
cp ${MAP_FILE} ${OUT_DIR}/gemma_analysis/${PLINK_PREFIX}.map

# Generate FAM file without phenotypes from recombined PED files
plink --file ${OUT_DIR}/gemma_analysis/${PLINK_PREFIX} \
    --make-just-fam \
    --allow-extra-chr \
    --out ${OUT_DIR}/gemma_analysis/${PLINK_PREFIX}

# Generate binary version of PED file (BED)
plink --file ${OUT_DIR}/gemma_analysis/${PLINK_PREFIX} \
    --make-bed \
    --allow-extra-chr \
    --out ${OUT_DIR}/gemma_analysis/${PLINK_PREFIX}

# First rename existing FAM file
mv ${OUT_DIR}/gemma_analysis/${PLINK_PREFIX}.fam ${OUT_DIR}/gemma_analysis/temp_${PLINK_PREFIX}_no_pheno.fam

# Generate a FAM file for founders only
FOUNDERS_PED_DIR=$(dirname ${FOUNDERS_PED})
FOUNDERS_PREFIX=$(basename ${FOUNDERS_PED} .ped)
# The FAM file is just the first 6 fields of the PED file
#   This is only needed for the purposes of knowing which individuals
#   are the parents in the FAM file, so we don't need the full plink set
#   of files
cut -d' ' -f 1-6 ${FOUNDERS_PED} > ${FOUNDERS_PED_DIR}/${FOUNDERS_PREFIX}.fam

# Update FAM file with all phenotypes from phenotype table
# Note: Multiple phenotype columns in the FAM file is specifically
#   formatted for use with GEMMA
combine_pheno_and_plink_fam.py \
    ${XO_DATA_DIR}/all_families_pheno_xo.txt \
    ${OUT_DIR}/gemma_analysis/temp_${PLINK_PREFIX}_no_pheno.fam \
    ${FOUNDERS_PED_DIR}/${FOUNDERS_PREFIX}.fam \
    ${PLINK_PREFIX} \
    ${OUT_DIR}/gemma_analysis
