#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
PLINK_PED="$1"
OUT_DIR="$2"
SCRIPT_DIR="$3"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Get PED file basename
PLINK_DIR=$(dirname ${PLINK_PED})

ARR=($(find ${PLINK_DIR}/output -name "*_lmm.assoc.txt" | sort -V))

function plot_gemma_output() {
    local filename="$1"
    local out_dir="$2"
    curr_pheno_name=$(basename ${filename} _lmm.assoc.txt)
    # Generate Manhattan plot
    echo "Generating Manhattan plot for phenotype: ${curr_pheno_name}..."
    GEMMA_Manhattan_Plot_GW_XO.R ${filename} ${curr_pheno_name} ${out_dir}/gemma_analysis
}

export -f plot_gemma_output

# Generate Manhattan plot for each phenotype
parallel --verbose plot_gemma_output {} ${OUT_DIR} ::: ${ARR[@]}
echo "Done generating plots."
