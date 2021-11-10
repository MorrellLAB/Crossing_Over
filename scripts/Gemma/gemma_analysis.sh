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

# Check if temp directory exists, if not, make it
mkdir -p ${TEMP}

# Get PED file basename
PLINK_DIR=$(dirname ${PLINK_PED})
PLINK_BN=$(basename ${PLINK_PED} .ped)

PHENO_COLNAMES=${PLINK_DIR}/${PLINK_BN}_pheno_order_in_fam.txt
PHENO_COLNAME_ARR=($(tr ' ' '\n' < ${PHENO_COLNAMES}))
# Add placeholder as first element of array so column numbers are easier to reference later on
PHENO_COLNAME_ARR=("placeholder" "${PHENO_COLNAME_ARR[@]}")
# Remember the first element is a placeholder, so don't count it
NUM_PHENO=$(echo $((${#PHENO_COLNAME_ARR[@]}-1)))

function run_gemma() {
    local plink_bfilename="$1"
    local out_dir="$2"
    local col_pheno="$3" # Which column and phenotype are we processing?
    col_num=$(echo ${col_pheno} | tr '..' '\t' | awk '{print $1}')
    pheno_colname=$(echo ${col_pheno} | tr '..' '\t' | awk '{print $2}')
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
SEQ_ARR=($(seq 0 ${NUM_PHENO}))
SEQ_PHENO_ARR=()
# Iteratively build array values
for num in ${SEQ_ARR[@]}
do
    if [[ ${num} != "0" ]]; then
        echo "Phenotype: ${PHENO_COLNAME_ARR[${num}]} in column ${num}."
        new_elem=$(echo "${num}..${PHENO_COLNAME_ARR[${num}]}")
        SEQ_PHENO_ARR+=("${new_elem}")
    fi
done

# Run in parallel
parallel --verbose --tmpdir ${TEMP} run_gemma ${PLINK_BN} ${OUT_DIR}/gemma_analysis {} ::: ${SEQ_PHENO_ARR[@]}
