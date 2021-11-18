#!/bin/bash

set -e
set -o pipefail

# Driver script that runs R/qtl2 crossover count analysis on all families

# User provided input arguments
YAML_DIR=$(realpath "$1") # Containing yaml files corresponding to each family
PCENT_FP="$2"
USERDEF_ERR_PROB="$3"
USERDEF_MAP_FN="$4"
USERDEF_ERROR_LOD_CUTOFF="$5"
USERDEF_PGENO_ERR_CUTOFF="$6"
USERDEF_PMISS_CUTOFF="$7"
MAXMARG_MINPROB="$8"
OUT_DIR="${9}/rqtl2"
SCRIPT_DIR="${10}"

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
         "${OUT_DIR}/physical_map_plots" \
         "${OUT_DIR}/physical_map_plots_miss" \
         "${OUT_DIR}/log_files"

function xo_counts() {
    local yaml_file="$1"
    local pcent_file="$2"
    local userdef_err_prob="$3"
    local userdef_map_fn="$4"
    local userdef_error_lod_cutoff="$5"
    local userdef_pgeno_err_cutoff="$6"
    local userdef_pmiss_cutoff="$7"
    local maxmarg_minprob="$8"
    local out_dir="$9"
    local fam_log_dir="${10}"
    name=$(basename "${yaml_file}" _forqtl2.yaml)
    printf "\n"
    echo "Processing sample: ${name}..."
    # Run crossover analysis
    rqtl2_xo_counts.R \
        "${yaml_file}" \
        "${pcent_file}" \
        "${userdef_err_prob}" \
        "${userdef_map_fn}" \
        "${userdef_error_lod_cutoff}" \
        "${userdef_pgeno_err_cutoff}" \
        "${userdef_pmiss_cutoff}" \
        "${maxmarg_minprob}" \
        "${out_dir}" \
        "${fam_log_dir}" 2>&1 | tee "${fam_log_dir}/${name}.log"
}

export -f xo_counts

# Check that error probability is a single non-negative number
if [[ "${USERDEF_ERR_PROB}" < 0 ]]; then
    echo "User defined error probability is negative, please specify a non-negative number in the config file and re-run. Exiting..."
    exit 1
else
    echo "User defined error probability is set to: ${USERDEF_ERR_PROB}"
    echo "User defined error probability is non-negative, proceeding to counting crossovers..."
fi

# Print out maxmarg cutoff being used
echo "User defined threshold for the minimum probability used to pick the most likely genotype from a set of genotype probabilities: ${MAXMARG_MINPROB}"

# Run program in parallel
printf "\n"
echo "##########################"
echo "Counting crossovers..."
parallel xo_counts {} "${PCENT_FP}" "${USERDEF_ERR_PROB}" "${USERDEF_MAP_FN}" "${USERDEF_ERROR_LOD_CUTOFF}" "${USERDEF_PGENO_ERR_CUTOFF}" "${USERDEF_PMISS_CUTOFF}" "${MAXMARG_MINPROB}" "${OUT_DIR}" "${OUT_DIR}/log_files" :::: "${YAML_DIR}/all_yaml_files_list.txt"

# Print some file number summaries to help catch errors
printf "\n"
echo "##########################"
echo "Printing some file number summaries to help detect errors/issues..."
num_yaml_files=$(wc -l "${YAML_DIR}/all_yaml_files_list.txt")
echo "Number of YAML files we started with (i.e., attempted to process): ${num_yaml_files}"
num_phys_map_plots_miss=$(find ${OUT_DIR}/physical_map_plots_miss/pmap_by_chr -name "*.pdf" | wc -l)
echo "Number of plots generated in ${OUT_DIR}/physical_map_plots_miss/pmap_by_chr directory: ${num_phys_map_plots_miss}"
num_pheno_tables=$(find ${OUT_DIR}/phenotype_tables -name "*pheno.txt" | wc -l)
echo "Number of phenotype tables output from script: ${num_pheno_tables}"
echo "If no errors occurred, the number of phenotype tables should be the same as the number of YAML files we started with."
echo "Identifying files where phenotype table did not get generated (likely due to errors)..."
find ${OUT_DIR}/phenotype_tables -name "*pheno.txt" | sed -e "s,${OUT_DIR}/phenotype_tables/,," -e "s,_pheno.txt,," | sort -V > ${OUT_DIR}/log_files/temp_pheno_tables_names.txt
sed -e "s,${YAML_DIR}/,," -e "s,_forqtl2.yaml,," "${YAML_DIR}"/all_yaml_files_list.txt | sort -V > ${OUT_DIR}/log_files/temp_yaml_names.txt
grep -vf ${OUT_DIR}/log_files/temp_pheno_tables_names.txt ${OUT_DIR}/log_files/temp_yaml_names.txt > ${OUT_DIR}/log_files/no_pheno_table_sample_names.txt
echo "List of sample names where phenotype table did not get generated is located here: ${OUT_DIR}/log_files/no_pheno_table_sample_names.txt"
echo "Log files are located at: ${OUT_DIR}/log_files"

# Reorganize output files
echo "Cleaning up intermediate files..."
# Cleanup intermediate files
rm -rf "${OUT_DIR}"/percent_missing_interactive_plots/*_percent_missing_files
rm -rf "${OUT_DIR}"/num_xo_interactive_plots/*_num_xo_files
rm -rf "${OUT_DIR}"/genotyping_error_lod_plots/*_files
echo "Done."
