#!/usr/bin/env Rscript

# This script takes in multiple input arguments and outputs
# a single YAML control file for R/qtl2

# Usage:
#   ./make_control_file.R [prefix] [geno_fp] [founder_geno_fp] [gmap_fp] [pmap_fp] [out_dir]

library(qtl2)

# Take command line arguments
args <- commandArgs(trailingOnly = TRUE)

# User provided input arguments
prefix <- args[1]
geno_fp <- args[2] # Single file or comma separated list of files
founder_geno_fp <- args[3] # Single file or comma separated list of files
gmap_fp <- args[4] # Single file or comma separated list of files
pmap_fp <- args[5] # Single file or comma separated list of files
out_dir <- args[6]
cross_type <- args[7]

outfile <- paste0(out_dir, "/", prefix, "_forqtl2.yaml")

# Check if we are working with a single file or a list of files
if (grepl(".csv", geno_fp)) {
    print("Generating control file for one set (family) of files.")
    geno_files <- geno_fp
    founder_geno_files <- founder_geno_fp
    gmap_files <- gmap_fp
    pmap_files <- pmap_fp
} else if (grepl(".txt", geno_fp)) {
    print("Generating control file for multiple sets of files provided in input file lists.")
    # Read in lists of files
    geno_files <- read.delim(geno_fp, header = FALSE)$V1
    founder_geno_files <- read.delim(founder_geno_fp, header = FALSE)$V1
    gmap_files <- read.delim(gmap_fp, header = FALSE)$V1
    pmap_files <- read.delim(pmap_fp, header = FALSE)$V1
}

# Make control file
if (cross_type == "f2") {
    write_control_file(
        output_file = outfile,
        crosstype = "f2",
        geno_file = geno_files,
        founder_geno_file = founder_geno_files,
        gmap_file = gmap_files,
        pmap_file = pmap_files,
        geno_codes = list(AA=1, AB=2, BB=3),
        na.strings = "NA",
        overwrite = TRUE
    )
} else if (cross_type == "dh") {
    write_control_file(
        output_file = outfile,
        crosstype = "dh",
        geno_file = geno_files,
        founder_geno_file = founder_geno_files,
        gmap_file = gmap_files,
        pmap_file = pmap_files,
        geno_codes = list(AA=1, BB=2),
        na.strings = "NA",
        overwrite = TRUE
    )
}
