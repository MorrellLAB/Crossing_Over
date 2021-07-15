#!/bin/bash

set -e
set -o pipefail

# Driver script that creates one R/qtl2 YAML control file for each family

# User provided input arguments
# Full path to qtl2inputs directory
#   .yaml control files will be created in this same directory
QTL_INPUTS_DIR="$1"
SCRIPT_DIR="$2"
CROSSTYPE="$3"
SUFFIX="_fillIn_progeny_AB_geno.csv"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/data_handling:"${PATH}"

# Check if family prefixes exist, if so delete so we start from a clean slate
if [[ -f "${QTL_INPUTS_DIR}/temp_family_prefixes.txt" ]]; then
    rm "${QTL_INPUTS_DIR}/temp_family_prefixes.txt"
fi
# Generate a file containing file prefixes
for i in $(ls "${QTL_INPUTS_DIR}"/*"${SUFFIX}")
do
    basename "$i" "${SUFFIX}" >> "${QTL_INPUTS_DIR}/temp_family_prefixes.txt"
done

function make_yaml() {
    local prefix="$1"
    local qtl2inputs_dir="$2"
    local crosstype="$3"
    outfile_prefix=$(basename "${ped_file}" _Mendel_fillIn.ped)
    # Make YAML control file for each family
    # Output dir is same as qtl2inputs_dir because yaml file
    # needs to be in same directory as inputs
    make_control_file.R \
        ${prefix} \
        ${prefix}_fillIn_progeny_AB_geno.csv \
        ${prefix}_fillIn_founder_AB_geno.csv \
        ${prefix}_fillIn_gmap.csv \
        ${prefix}_fillIn_pmap.csv \
        ${qtl2inputs_dir} \
        ${crosstype}
}

export -f make_yaml

# Run program in parallel
echo "Generating .yaml Rqtl2 input files..."
cd ${OUT_DIR}
parallel make_yaml {} ${QTL_INPUTS_DIR} ${CROSSTYPE} :::: "${QTL_INPUTS_DIR}/temp_family_prefixes.txt"
echo "Done generating .yaml Rqtl2 input files."

# !!! Add check to make sure each set of Rqtl2 input files has the corresponding .yaml control file.
