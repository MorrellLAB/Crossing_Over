#!/bin/bash

set -e
set -o pipefail

# Generate boxplots for each group listed in PLOT_GROUPINGS
#   See toy dataset family_plot_groupings.csv for more details on file format

# User provided input arguments
XO_DATA_DIR="$1"
OUT_DIR="${2}/rqtl2"
PEDIGREE="$3"
PLOT_GROUPINGS="$4"
SCRIPT_DIR="$5"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/crossovers:"${PATH}"

# Check if out dir exists, if not make it
mkdir -p ${OUT_DIR}

# Generate summary plots
xo_count_summary_plotting.R \
    ${XO_DATA_DIR} \
    ${OUT_DIR} \
    ${PEDIGREE} \
    ${YAML_DIR} \
    ${PLOT_GROUPINGS}
