#!/bin/bash

# Stores functions shared across scripts
#   Functions here peither perform some sort of check or processes data/runs analyses

#-----------------------------------
# Functions that perform checks

# Test if required function has been successfully loaded into the environment
function test_function() {
    local fn_name="$1"
    if [[ $(type -t ${fn_name}) == function ]]
    then
        echo "${fn_name} function exists, proceeding..."
    else
        echo "${fn_name} doesn't exist."
        echo "Please check that the function exists and has been defined, exiting..."
        exit 127
    fi
}

export -f test_function

#-----------------------------------
# Functions that processes data or runs analyses

# Convert Plink PED/MAP format to R/qtl2 input file formats
function plink2rqtl2() {
    local ped_file="$1"
    local map_file="$2"
    local lookup_table="$3"
    local out_dir="$4"
    local file_suffix="$5"
    outfile_prefix=$(basename ${ped_file} ${file_suffix})
    # Convert PED and MAP files to rqtl2 file formats
    Plink2Rqtl2.py \
        "${ped_file}" \
        "${map_file}" \
        "${lookup_table}" \
        "${out_dir}/${outfile_prefix}"
}

export -f plink2rqtl2
