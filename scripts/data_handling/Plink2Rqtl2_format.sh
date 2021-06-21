#!/bin/bash

set -e
set -o pipefail

# Driver script that runs mendel checking with Plink on each family

# User provided input arguments
PED_DIR="$1"
MAP_FILE="$2"
LOOKUP_TABLE="$3"
OUT_DIR="$4"
SCRIPT_DIR="$5"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/data_handling:"${PATH}"

# Check that out dir exists, if not make it
mkdir -p "${OUT_DIR}/qtl2inputs"

# Store ped files in array
# Pull only progeny and exclude parents (denoted by "-9" in filename)
find "${PED_DIR}" -name "*.ped" | sort -V | grep -v "\-9" > "${PED_DIR}"/ped_list.txt

function plink2rqtl2() {
    local ped_file="$1"
    local map_file="$2"
    local lookup_table="$3"
    local out_dir="$4"
    outfile_prefix=$(basename "${ped_file}" _Mendel_fillIn.ped)
    # Convert PED and MAP files to rqtl2 file formats
    Plink2Rqtl2.py \
        "${ped_file}" \
        "${map_file}" \
        "${lookup_table}" \
        "${out_dir}/${outfile_prefix}_fillIn"
}

export -f plink2rqtl2

# Run program in parallel
echo "Converting Plink to Rqtl2 input file formats..."
parallel plink2rqtl2 {} "${MAP_FILE}" "${LOOKUP_TABLE}" "${OUT_DIR}/qtl2inputs" :::: "${PED_DIR}/ped_list.txt"
echo "Done converting Plink to Rqtl2 input file formats."
