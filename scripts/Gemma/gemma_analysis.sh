#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
PLINK_PED="$1"
OUT_DIR="$2"
SCRIPT_DIR="$3"
TEMP="$4"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Compute kinship matrix
# Eventually, GEMMA installation will be automatically taken care of
#   but for development purposes, to keep the software reference consistent,
#   we created a symbolic link called "gemma" to the linux static build
#   Example command run from within the same directory as the linux static build"
#       ln -s /path/to/static_build/GEMMA-0.98.1/gemma-0.98.1-linux-static gemma

# Get PED file basename
PLINK_DIR=$(dirname ${PLINK_PED})
PLINK_BN=$(basename ${PLINK_PED} .ped)

PHENO_COLNAMES=${PLINK_DIR}/${PLINK_BN}_pheno_order_in_fam.txt
PHENO_COLNAME_ARR=($(tr ' ' '\n' < ${PLINK_DIR}/${PLINK_BN}_pheno_order_in_fam.txt))
# Add placeholder as first element of array so column numbers are easier to reference later on
PHENO_COLNAME_ARR=("placeholder" "${PHENO_COLNAME_ARR[@]}")
# Remember the first element is a placeholder, so don't count it
NUM_PHENO=$(echo $((${#PHENO_COLNAME_ARR[@]}-1)))

function run_gemma() {
    local plink_bfilename="$1"
    local out_dir="$2"
    local col_num="$3" # Which phenotype column are we processing?
    local pheno_colname="$4" # Pulled from phenotype table header
    echo "Processing phenotype: ${pheno_colname} in column ${col_num}..."
    gemma \
        -bfile ${plink_bfilename} \
        -n ${col_num} \
        -gk \
        -o ${pheno_colname}_${plink_bfilename}
    
    # Run univariate LMM and include column to indicate direction of effects
    gemma \
        -bfile ${plink_bfilename} \
        -n ${col_num} \
        -k ${out_dir}/output/${pheno_colname}_${plink_bfilename}.cXX.txt \
        -lmm 4 \
        -o ${pheno_colname}_${plink_bfilename}_lmm
}

export -f run_gemma

# Gemma only allows writing output to current working directory
cd ${PLINK_DIR}
# Each phenotype is one column in the FAM file starting at column 6
#   Remember index 0 is a placeholder, so we start with index 1
for i in $(seq 1 ${NUM_PHENO})
do
    curr_pheno_colname=${PHENO_COLNAME_ARR[$i]}
    run_gemma ${PLINK_BN} ${OUT_DIR}/gemma_analysis ${i} ${curr_pheno_colname}
done

# Run in parallel
# parallel --tmpdir ${TEMP} run_gemma {} ${OUT_DIR}/gemma_analysis/plink_pheno_files ::: ${PLINK_BN_ARR[@]}
