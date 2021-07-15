#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
PED_FILE="$1"
OUT_DIR="$2"
SCRIPT_DIR="$3"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/Gemma:"${PATH}"

# Check if output directories exist
mkdir -p "${OUT_DIR}"

# Get PED file basename
PED_BASENAME=$(basename ${PED_FILE} .ped)

# Compute kinship matrix
# !!! Calling on Gemma tool needs to be generalized to work on all Gemma versions and not rely on specific installation path
~/Shared/Software/GEMMA-0.98.1/gemma-0.98.1-linux-static \
    -bfile ${PED_BASENAME} \
    -gk \
    -o ${OUT_DIR}/output/gw_xo_count_${PED_BASENAME}

# Run univariate LMM and include column to indicate direction of effects
~/Shared/Software/GEMMA-0.98.1/gemma-0.98.1-linux-static \
    -bfile ${PED_BASENAME} \
    -k ${OUT_DIR}/output/gw_xo_count_${PED_BASENAME}.cXX.txt \
    -lmm 4 \
    -o ${OUT_DIR}/output/gw_xo_count_${PED_BASENAME}_lmm

# Generate Manhattan plot
GEMMA_Manhattan_Plot_GW_XO.R ${OUT_DIR}/output/gw_xo_count_${PED_BASENAME}_lmm.assoc.txt
