#!/bin/bash

set -e
set -o pipefail

# Driver script that runs mendel checking with Plink on each family

# User provided input arguments
YAML_PREFIX="$1"
PROGENY_GENO_LIST="$2"
FOUNDER_GENO_LIST="$3"
GMAP_LIST="$4"
PMAP_LIST="$5"
CROSSTYPE="$6"
OUT_DIR="$7"
SCRIPT_DIR="$8"

# We'll assume file ends in .ped and use that as the file suffix
#   since this step is working with all PED data combined into one file
SUFFIX=".ped"

# Export path to directory that contains executable script
export PATH="${SCRIPT_DIR}"/scripts/data_handling:"${PATH}"
# # Load functions defined in utilities script
# source "${SCRIPT_DIR}"/scripts/utils.sh

# Check that out dir exists, if not make it
mkdir -p "${OUT_DIR}/data_diagnostics"

# Build qtl inputs directory path from one of files
#   Reason: the YAML file needs to be in the same directory as the inputs
QTL2INPUTS_DIR=$(dirname $(head -n 1 ${PROGENY_GENO_LIST}))

# Make one YAML control file for given list of files
# Output dir is same as qtl2inputs_dir because yaml file
# needs to be in same directory as inputs
make_control_file.R \
    ${YAML_PREFIX} \
    ${PROGENY_GENO_LIST} \
    ${FOUNDER_GENO_LIST} \
    ${GMAP_LIST} \
    ${PMAP_LIST} \
    ${QTL2INPUTS_DIR} \
    ${CROSSTYPE}

# Run data diagnostics script
