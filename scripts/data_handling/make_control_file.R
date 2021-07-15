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
geno_fp <- args[2]
founder_geno_fp <- args[3]
gmap_fp <- args[4]
pmap_fp <- args[5]
out_dir <- args[6]
cross_type <- args[7]

outfile <- paste0(out_dir, "/", prefix, "_forqtl2.yaml")

# Make control file
if (cross_type == "f2") {
    # Note: currently some hardcoded values below are tailored
    # for Kono et al. GP population
    write_control_file(
        output_file = outfile,
        crosstype = "f2",
        geno_file = geno_fp,
        founder_geno_file = founder_geno_fp,
        gmap_file = gmap_fp,
        pmap_file = pmap_fp,
        geno_codes = list(AA=1, AB=2, BA=3, BB=4),
        na.strings = "NA",
        overwrite = TRUE
    )
} else if (cross_type == "dh") {
    write_control_file(
        output_file = outfile,
        crosstype = "dh",
        geno_file = geno_fp,
        founder_geno_file = founder_geno_fp,
        gmap_file = gmap_fp,
        pmap_file = pmap_fp,
        geno_codes = list(AA=1, BB=2),
        na.strings = "NA",
        overwrite = TRUE
    )
}
