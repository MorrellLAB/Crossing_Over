#!/bin/bash

set -e
set -o pipefail

# Driver script that runs R/qtl2 crossover count analysis on all families

# User provided input arguments
YAML_DIR="$1"
PCENT_FP="$2"
OUT_DIR="${3}/rqtl2"
SCRIPT_DIR="$4"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/crossovers:"${PATH}"

# Generate list of control files (i.e., .yaml files)
find "${YAML_DIR}" -name "*.yaml" | sort -V > "${YAML_DIR}/all_yaml_files_list.txt"

# Check if out dir exists, if not make it
mkdir -p "${OUT_DIR}" \
         "${OUT_DIR}/genetic_map_plots" \
         "${OUT_DIR}/genetic_map_plots/gmap_by_chr" \
         "${OUT_DIR}/genetic_map_plots/gmap_by_ind" \
         "${OUT_DIR}/other_summaries" \
         "${OUT_DIR}/phenotype_tables" \
         "${OUT_DIR}/physical_map_plots" \
         "${OUT_DIR}/physical_map_plots/pmap_by_chr" \
         "${OUT_DIR}/physical_map_plots/pmap_by_ind" \
         "${OUT_DIR}/too_few_markers" \
         "${OUT_DIR}/percent_missing_interactive_plots" \
         "${OUT_DIR}/percent_missing_data" \
         "${OUT_DIR}/num_xo_interactive_plots" \
         "${OUT_DIR}/num_xo_data"

function xo_counts() {
    local yaml_file="$1"
    local pcent_file="$2"
    local out_dir="$3"
    name=$(basename "${yaml_file}" _forqtl2.yaml)
    # Run crossover analysis
    rqtl2_xo_counts.R \
        "${yaml_file}" \
        "${pcent_file}" \
        "${out_dir}" &> "${out_dir}"/"${name}".log
}

export -f xo_counts

# Run program in parallel
parallel xo_counts {} "${PCENT_FP}" "${OUT_DIR}" :::: "${YAML_DIR}/all_yaml_files_list.txt"

# Reorganize output files
cd "${OUT_DIR}"
echo "Reorganizing..."
set -x
mv too_few_markers_per_chr*.txt "${OUT_DIR}/too_few_markers/"
mv *gmap_xo_cM_by_chr.pdf "${OUT_DIR}/genetic_map_plots/gmap_by_chr/"
mv *gmap_xo_cM_by_ind.pdf "${OUT_DIR}/genetic_map_plots/gmap_by_ind/"
mv *pheno.txt "${OUT_DIR}/phenotype_tables/"
mv *pmap_xo_Mbp_by_chr.pdf "${OUT_DIR}/physical_map_plots/pmap_by_chr/"
mv *pmap_xo_Mbp_by_ind.pdf "${OUT_DIR}/physical_map_plots/pmap_by_ind/"
mv *.log "${OUT_DIR}/other_summaries/"
mv *_duplicates_summary.txt "${OUT_DIR}/other_summaries/"
mv *_percent_missing.html "${OUT_DIR}/percent_missing_interactive_plots/"
mv *_percent_missing.txt "${OUT_DIR}/percent_missing_data/"
mv *_num_xo.html "${OUT_DIR}/num_xo_interactive_plots/"
mv *_xo_count.txt "${OUT_DIR}/num_xo_data/"
# Cleanup intermediate files
rm -rf "${OUT_DIR}/*_percent_missing_files"
rm -rf "${OUT_DIR}/*_num_xo_files"
set +x
