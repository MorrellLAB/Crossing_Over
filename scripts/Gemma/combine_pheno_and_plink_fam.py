#!/usr/bin/env python3
"""This script adds xo phenotypes from the phenotype table generated to
the FAM files in preparation for dowstream analysis with GEMMA. This
script handles the multi-column phenotype table and creates a FAM file
that includes every phenotype column. Currently works with PLINK 1.9 FAM
file format.

Usage: ./combine_pheno_and_plink_fam.py [xo_pheno_fp] [fam_fp] [plink_prefix] [out_dir]

Where:
1) [xo_pheno_fp] is the full filepath to the xo phenotypes table
2) [fam_fp] is the full filepath to the FAM file to update
3) [plink_prefix] is the prefix that matches the set of plink PED/MAP files
4) [out_dir] is the full filepath to our output directory
"""

import os
import sys


def read_xo_counts(xo_pheno_fp):
    """Read in xo counts file"""
    counts_dict = {}
    with open(xo_pheno_fp, "rt") as handle:
        for record in handle:
            if record.startswith("sampleID"):
                counts_dict["header_line"] = record.split()
            else:
                tmp = record.split()
                ind = tmp[0]
                counts_dict[ind] = tmp
    return(counts_dict)


def read_plink_fam(plink_fam_fp):
    """Read in plink .ped file."""
    parents_dict = {}
    fam_dict = {}
    with open(plink_fam_fp, "rt") as handle:
        for record in handle:
            if record.startswith("-9"):
                parents_dict[record.split()[1]] = record.split()
            else:
                fam_dict[record.split()[1]] = record.split()
    return(parents_dict, fam_dict)


def add_xo_pheno(fam, xo_pheno):
    """Add xo phenotype to fam dictionary. xo_pheno is the
    phenotype table stored in a dictionary."""
    # Stores updated PED with xo phenotype added
    updated_fam_dict = {}

    # Add xo phenotype
    # Iteratively add one phenotype column to ped
    for key in fam.keys():
        # Pull corresponding phenotype column
        #   Individual IDs should match
        updated_values = fam[key][0:5] + xo_pheno[key][1:]
        # Add phenotypes to fam_dict
        updated_fam_dict[key] = updated_values
    return(updated_fam_dict)


def write_to_fam(parents, fam_w_pheno, out_fp):
    """Write PED file for individuals with currently processing phenotype."""
    # Start from clean file, check if file exists
    if os.path.exists(os.path.expanduser(out_fp)):
        os.remove(os.path.expanduser(out_fp))
    # Start with parents first
    with open(os.path.expanduser(out_fp), 'a') as file:
        for counter, i in enumerate(parents):
            file.write(" ".join(parents[i]) + "\n")
    # Now, add progeny
    with open(os.path.expanduser(out_fp), 'a') as file:
        for counter, i in enumerate(fam_w_pheno):
            file.write(" ".join(fam_w_pheno[i]) + "\n")
    return


def main(xo_pheno_fp, fam_fp, plink_prefix, out_dir):
    """Driver function."""
    # Output file prefix
    # Read in files
    xo_pheno = read_xo_counts(xo_pheno_fp)
    parents, fam = read_plink_fam(fam_fp)
    # Add xo phenotype to fam dict
    updated_fam = add_xo_pheno(fam, xo_pheno)
    # Prepare parent phenotype placeholders (missing phenotype)
    #   so we have the same number of columns as the progeny
    updated_parents = {}
    # Get the number of phenotype columns
    num_pheno = len(xo_pheno['header_line'][1:])
    # Prepare list of missing phenotypes
    miss_pheno = ["NA"] * num_pheno
    for key in parents.keys():
        updated_values = parents[key][0:5] + miss_pheno
        updated_parents[key] = updated_values

    # Prepare output file path
    out_fp = os.path.expanduser(out_dir.rstrip('/')) + '/' + plink_prefix + '.fam'
    out_pheno_order_fp = os.path.expanduser(out_dir.rstrip('/')) + '/' + plink_prefix + '_pheno_order_in_fam.txt'
    # Save to output file
    write_to_fam(updated_parents, updated_fam, out_fp)
    # Save phenotype column names/order to file for reference
    with open(os.path.expanduser(out_pheno_order_fp), 'w') as file:
        file.write(" ".join(xo_pheno['header_line'][1:]) + "\n")
    print("Done adding phenotypes to FAM file.")
    return

# Print usage message if we don't have enough input arguments
if len(sys.argv) < 2:
    print(__doc__)
    exit(1)
else:
    # Run the program
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) # Run the program
