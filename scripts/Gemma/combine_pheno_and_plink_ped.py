#!/usr/bin/env python3
"""This script adds xo phenotypes from the phenotype table generated to
the PED files in preparation for dowstream analysis. This script handles
the multi-column phenotype table and creates a PED file for every phenotype column. Currently works with PLINK 1.9 PED file format.

Usage: ./combine_pheno_and_plink_ped.py [xo_pheno_fp] [ped_fp] [out_dir]

Where:
1) [xo_pheno_fp] is the full filepath to the xo phenotypes table
2) [ped_fp] is the full filepath to the recombined PED file
3) [out_dir] is the full filepath to our output directory
"""

import os
import sys
import pandas as pd

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


def dict_to_pddf(input_dict):
    """Convert dictionary to pandas data frame."""
    df = pd.DataFrame.from_dict(input_dict, orient='index')
    return(df)


def read_plink_ped(plink_ped_fp):
    """Read in plink .ped file."""
    parents_dict = {}
    ped_dict = {}
    with open(plink_ped_fp, "rt") as handle:
        for record in handle:
            if record.startswith("-9"):
                parents_dict[record.split()[1]] = record.split()
            else:
                ped_dict[record.split()[1]] = record.split()
    return(parents_dict, ped_dict)


def add_xo_pheno(ped, xo_pheno):
    """Add xo phenotype to ped dictionary. If individual
    in ped dictionary does not have a xo phenotype, save to
    separate dictionary remove from output PED. xo_pheno is the
    phenotype table stored in a dictionary."""
    # Convert xo_pheno dict to pandas data frame
    pheno_df = dict_to_pddf(xo_pheno)
    # Add xo phenotype
    # Stores updated PED with xo phenotype added
    updated_ped_dict = {}

    # Iteratively add one phenotype column to ped
    for cidx in pheno_df.columns:
        if cidx == 0:
            continue
        else:
            # Contains sampleID and one phenotype column
            curr_pheno_df = pheno_df[[0, cidx]]
            # Convert pandas dataframe to dictionary for downstream processing
            curr_pheno_dict = curr_pheno_df.set_index(0).T.to_dict(orient="list")
            for key in curr_pheno_dict.keys():
                if key == "sampleID":
                    continue
                elif key in ped.keys():
                    # Save current phenotype name (e.g., chr7_total)
                    # Add as a new value in dictionary
                    pheno_name = curr_pheno_dict['sampleID'][0]
                    # Pull out xo phenotype
                    curr_pheno_val = curr_pheno_dict[key][0]
                    # Work corresponding sampleID in PED
                    # Replace missing phenotype (-9) with new phenotype
                    new_key = key + ":" + pheno_name
                    current = ped[key]
                    updated_ped_dict[new_key] = current
                    updated_ped_dict[new_key][5] = curr_pheno_val
    return(updated_ped_dict)


def write_to_ped(parents, curr_pheno_ped, out_fp):
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
        for counter, i in enumerate(curr_pheno_ped):
            file.write(" ".join(curr_pheno_ped[i]) + "\n")
    return


def main(xo_pheno_fp, ped_fp, out_dir):
    """Driver function."""
    # Output file prefix
    # Read in files
    xo_pheno = read_xo_counts(xo_pheno_fp)
    parents, ped = read_plink_ped(ped_fp)
    # Add xo phenotype to ped dict
    updated_ped = add_xo_pheno(ped, xo_pheno)
    # Save to output file
    # One phenotype per PED file
    for pelem in xo_pheno['header_line']:
        if pelem == "sampleID":
            continue
        else:
            # Prepare output file paths
            out_suffix = pelem
            # Make sure to remove trailing slash from directory path
            out_fp = os.path.expanduser(out_dir.rstrip('/')) + '/pheno-' + out_suffix + '.ped'
            # Make placeholder dictionary for current phenotype
            tmp_pheno = {}
            for key in updated_ped.keys():
                tmp = key.split(':')
                if pelem == tmp[1]:
                    # Add individuals for current phenotype to tmp dictionary
                    tmp_pheno[tmp[0]] = updated_ped[key]
            # Write to file
            write_to_ped(parents, tmp_pheno, out_fp)
    print("Done.")
    return

main(sys.argv[1], sys.argv[2], sys.argv[3]) # Run the program
