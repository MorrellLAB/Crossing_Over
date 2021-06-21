#!/bin/bash

set -e
set -o pipefail

# User provided input arguments
PED_FILE="$1"
OUT_DIR="$2"
SCRIPT_DIR="$3"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/data_handling:"${PATH}"

# Check if output directories exist
mkdir -p "${OUT_DIR}"

# Run script
split_by_family.py \
    "${PED_FILE}" \
    "${OUT_DIR}"

# Generate list of split PED files
#   Exclude the parents themselves, which have a naming scheme similar to "-9_Mendel.ped"
find "${OUT_DIR}/split_by_family" -name "*.ped" | sort -V | grep -v "\-9" > "${OUT_DIR}/split_by_family/split_ped_files_list.txt"

# Need to add some checks here to make sure the expected files are produced
