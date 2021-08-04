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
SUFFIX="$6" # Example: "_Mendel_fillIn.ped"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/data_handling:"${PATH}"
# Load functions defined in utilities script
source "${SCRIPT_DIR}"/scripts/utils.sh
# Check that required functions exist
test_function plink2rqtl2

# Check that out dir exists, if not make it
mkdir -p "${OUT_DIR}/qtl2inputs"

# Store ped files in array
# Pull only progeny and exclude parents (denoted by "-9" in filename)
find "${PED_DIR}" -name "*.ped" | sort -V | grep -v "\-9" > "${PED_DIR}"/ped_list.txt

# Run program in parallel
echo "Converting Plink to Rqtl2 input file formats..."
parallel plink2rqtl2 {} "${MAP_FILE}" "${LOOKUP_TABLE}" "${OUT_DIR}/qtl2inputs" "${SUFFIX}" :::: "${PED_DIR}/ped_list.txt"
echo "Done converting Plink to Rqtl2 input file formats."
