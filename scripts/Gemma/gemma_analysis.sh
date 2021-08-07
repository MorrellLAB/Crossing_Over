#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
PED_FILE="$1"
OUT_DIR="$2"
SCRIPT_DIR="$3"
TEMP="$4"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Check if output directories exist
mkdir -p "${OUT_DIR}/gemma_analysis/plink_pheno_files"

# Get PED file basename
PED_BASENAME=$(basename ${PED_FILE} .ped)

# Compute kinship matrix
# Eventually, GEMMA installation will be automatically taken care of
#   but for development purposes, to keep the software reference consistent,
#   we created a symbolic link called "gemma" to the linux static build
#   Example command run from within the same directory as the linux static build"
#       ln -s /path/to/static_build/GEMMA-0.98.1/gemma-0.98.1-linux-static gemma

function run_gemma() {
    local plink_bfilename="$1"
    local out_dir="$2"
    gemma \
        -bfile ${plink_bfilename} \
        -gk \
        -o gw_xo_count_${plink_bfilename}
    
    # Run univariate LMM and include column to indicate direction of effects
    gemma \
        -bfile ${plink_bfilename} \
        -k ${out_dir}/output/gw_xo_count_${plink_bfilename}.cXX.txt \
        -lmm 4 \
        -o gw_xo_count_${plink_bfilename}_lmm
}

export -f run_gemma

# Gemma only allows writing output to current working directory
cd ${OUT_DIR}/gemma_analysis/plink_pheno_files
PLINK_BN_ARR=($(find . -name "pheno*.ped" | sort -V | sed -e 's,./,,' -e 's,.ped,,'))

for i in ${PLINK_BN_ARR[@]}
do
    echo "Processing plink file basename: ${i}..."
    run_gemma ${i} ${OUT_DIR}/gemma_analysis/plink_pheno_files
done

# Run in parallel
# parallel --tmpdir ${TEMP} run_gemma {} ${OUT_DIR}/gemma_analysis/plink_pheno_files ::: ${PLINK_BN_ARR[@]}

# Generate Manhattan plot
GEMMA_Manhattan_Plot_GW_XO.R ${OUT_DIR}/gemma_analysis/plink_pheno_files/output/gw_xo_count_${PLINK_BN_ARR}_lmm.assoc.txt
