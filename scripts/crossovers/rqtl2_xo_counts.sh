#!/bin/bash

set -e
set -o pipefail

# Driver script that runs R/qtl2 crossover count analysis on all families

# User provided input arguments
YAML_DIR=$(realpath "$1") # Containing yaml files corresponding to each family
PCENT_FP="$2"
OUT_DIR="${3}/rqtl2"
SCRIPT_DIR="$4"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/crossovers:"${PATH}"

# Generate list of control files (i.e., .yaml files)
# If data diagnostics step was run, there may be extra YAML control files
#   that don't apply to this specific step. Exclude those.
# These will automatically have a prefix of "data_diagnostics"
# Exclude YAML files that don't apply to this step, these files
#   have a file prefix of "data_diagnostics"
echo "Generating list of control files for each family."
if test -n "$(find ${YAML_DIR} -name 'data_diagnostics*' -print -quit)"
then
    echo "exists"
    echo "Files with prefix "data_diagnostics" exists, excluding them"
    find "${YAML_DIR}" -name "*.yaml" | grep -v "data_diagnostics" | sort -V > "${YAML_DIR}/all_yaml_files_list.txt"
else
    # Don't have files with prefix "data_diagnostics"
    # don't need to do anything special
    find "${YAML_DIR}" -name "*.yaml" | sort -V > "${YAML_DIR}/all_yaml_files_list.txt"
fi

# Check if out dir exists, if not make it
mkdir -p "${OUT_DIR}" \
         "${OUT_DIR}/genetic_map_plots" \
         "${OUT_DIR}/other_summaries" \
         "${OUT_DIR}/physical_map_plots"

function xo_counts() {
    local yaml_file="$1"
    local pcent_file="$2"
    local out_dir="$3"
    name=$(basename "${yaml_file}" _forqtl2.yaml)
    echo "Processing sample: ${name}..."
    # Run crossover analysis
    rqtl2_xo_counts.R \
        "${yaml_file}" \
        "${pcent_file}" \
        "${out_dir}" 2>&1 | tee "${out_dir}"/other_summaries/"${name}".log
}

export -f xo_counts

# Run program in parallel
echo "Counting crossovers..."
parallel xo_counts {} "${PCENT_FP}" "${OUT_DIR}" :::: "${YAML_DIR}/all_yaml_files_list.txt"

# Reorganize output files
echo "Cleaning up intermediate files..."
# Cleanup intermediate files
rm -rf "${OUT_DIR}"/percent_missing_interactive_plots/*_percent_missing_files
rm -rf "${OUT_DIR}"/num_xo_interactive_plots/*_num_xo_files
echo "Done."
