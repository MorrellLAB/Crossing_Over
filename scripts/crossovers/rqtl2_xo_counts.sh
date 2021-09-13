#!/bin/bash

set -e
set -o pipefail

# Driver script that runs R/qtl2 crossover count analysis on all families

# User provided input arguments
YAML_DIR=$(realpath "$1") # Containing yaml files corresponding to each family
PCENT_FP="$2"
USERDEF_ERR_PROB="$3"
OUT_DIR="${4}/rqtl2"
SCRIPT_DIR="$5"

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
    local userdef_err_prob="$3"
    local out_dir="$4"
    name=$(basename "${yaml_file}" _forqtl2.yaml)
    echo "Processing sample: ${name}..."
    # Run crossover analysis
    rqtl2_xo_counts.R \
        "${yaml_file}" \
        "${pcent_file}" \
        "${userdef_err_prob}" \
        "${out_dir}" 2>&1 | tee "${out_dir}"/other_summaries/"${name}".log
}

export -f xo_counts

# Check that error probability is a single non-negative number
if [[ "${USERDEF_ERR_PROB}" < 0 ]]; then
    echo "User defined error probability is negative, please specify a non-negative number in the config file and re-run. Exiting..."
    exit 1
else
    echo "User defined error probability is non-negative, proceeding to counting crossovers..."
fi

# Run program in parallel
echo "Counting crossovers..."
parallel xo_counts {} "${PCENT_FP}" "${USERDEF_ERR_PROB}" "${OUT_DIR}" :::: "${YAML_DIR}/all_yaml_files_list.txt"

# Reorganize output files
echo "Cleaning up intermediate files..."
# Cleanup intermediate files
rm -rf "${OUT_DIR}"/percent_missing_interactive_plots/*_percent_missing_files
rm -rf "${OUT_DIR}"/num_xo_interactive_plots/*_num_xo_files
echo "Done."
