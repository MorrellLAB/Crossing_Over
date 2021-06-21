#!/bin/bash

set -e
set -o pipefail

# Driver script to run PED_parental_fill-in.py script on each family

# User provided input arguments
PED_LIST="$1"
OUT_DIR="$2"
SCRIPT_DIR="$3"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/data_handling:"${PATH}"

# Check that out dir exists, if not make it
mkdir -p "${OUT_DIR}"/split_by_family_fillIn
cd "${OUT_DIR}"

# Store ped files in array
PED_ARR=($(cat "${PED_LIST}"))

function ped_parental_fill_in() {
    local ped_file="$1"
    local out_dir="$2"
    #local script_dir="$3"
    name=$(basename "${ped_file}" _Mendel.ped)
    # Run script for parental fill in
    PED_parental_fill-in.py \
        "${ped_file}" \
        "${out_dir}"/${name}_Mendel_fillIn.ped \
        "${out_dir}"/${name}_fillIn_tracking.txt
}

export -f ped_parental_fill_in

# Run program in parallel
parallel ped_parental_fill_in {} "${OUT_DIR}"/split_by_family_fillIn ::: "${PED_ARR[@]}"
