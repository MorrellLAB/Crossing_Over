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

# Only run GEMMA analysis on families that have phenotype tables
#   successfully generated. Some families may not have phenotype tables
#   generated due to errors occurring.
# Generate list of family names where phenotype table was generated
find ${XO_DATA_DIR} -name "*_pheno.txt" | sort -V | sed -e "s,${XO_DATA_DIR},," -e 's,_pheno.txt,,' > ${XO_DATA_DIR}/all_families_pheno_names.txt

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

# Update FAM file with all phenotypes from phenotype table
# Note: Multiple phenotype columns in the FAM file is specifically
#   formatted for use with GEMMA
combine_pheno_and_plink_fam.py \
    ${XO_DATA_DIR}/all_families_pheno_xo.txt \
    ${OUT_DIR}/gemma_analysis/temp_${PLINK_PREFIX}_no_pheno.fam \
    ${PLINK_PREFIX} \
    ${OUT_DIR}/gemma_analysis

# Update PED file
# PED_PREFIX=$(basename ${PED_FILE} .ped)
# combine_pheno_and_plink_ped.py \
#     ${XO_DATA_DIR}/all_families_pheno_xo.txt \
#     ${OUT_DIR}/gemma_analysis/all_families.ped \
#     ${OUT_DIR}/gemma_analysis/plink_pheno_files

# # Generate FAM/MAP/BED/BIM files using plink
# for i in $(find ${OUT_DIR}/gemma_analysis/plink_pheno_files -name "pheno*.ped" | sort -V)
# do
#     filename=$(basename ${i} .ped)
#     # FAM file is the same as the first six fields in a PED file
#     cut -d' ' -f 1-6 ${i} > ${OUT_DIR}/gemma_analysis/plink_pheno_files/${filename}.fam
#     # Generate MAP file and put in same directory as plink files
#     cp ${MAP_FILE} ${OUT_DIR}/gemma_analysis/plink_pheno_files/${filename}.map
#     # Generate BED/BIM files
#     plink \
#         --file ${OUT_DIR}/gemma_analysis/plink_pheno_files/${filename} \
#         --make-bed \
#         --allow-extra-chr \
#         --out ${OUT_DIR}/gemma_analysis/plink_pheno_files/${filename}
# done
