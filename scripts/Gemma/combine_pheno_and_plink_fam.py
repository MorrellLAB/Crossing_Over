#!/usr/bin/env python3

# Currently, script handles phenotype in a single column and
# works for PLINK 1.9 FAM file format.
# In the future, will handle multiple columns of phenotypes.

import os
import sys

def read_xo_counts(xo_counts_fp):
    """Read in xo counts file"""
    counts_dict = {}
    with open(xo_counts_fp, "rt") as handle:
        for record in handle:
            if record.startswith("ind"):
                continue
            else:
                ind = record.split()[0]
                xo_count = record.split()[1]
                counts_dict[ind] = xo_count
    return counts_dict


def read_plink_fam(plink_fam_fp):
    """Read in plink .fam file."""
    parents_dict = {}
    fam_dict = {}
    with open(plink_fam_fp, "rt") as handle:
        for record in handle:
            if record.startswith("-9"):
                parents_dict[record.split()[1]] = record.split()
            else:
                fam_dict[record.split()[1]] = record.split()
    return(parents_dict, fam_dict)


def add_xo_counts(out_file, fam, xo_counts):
    """Add xo count phenotype to fam dictionary. If individual
    in fam dictionary does not have a xo phenotype, save to
    log file and remove from output .fam file."""
    # Make sure we start with clean log file
    # log_dir = os.path.dirname(out_file)
    # temp = log_dir, "/individuals_missing_xo_phenotype.log"
    # log_filename = ''.join(temp)
    # if os.path.exists(os.path.expanduser(log_filename)):
    #     os.remove(os.path.expanduser(log_filename))
    # Add xo phenotype
    updated_fam_dict = {}
    for counter, i in enumerate(fam):
        if i in xo_counts.keys():
            current = fam[i]
            current[5] = xo_counts[i]
            updated_fam_dict[i] = current
        else:
            # Leave as is
            updated_fam_dict[i] = fam[i]
            # print("Individual not in xo counts phenotype file:", i)
            # print("Saving to log file...")
            # with open(os.path.expanduser(log_filename), 'a') as handle:
            #     handle.write(" ".join(fam[i]) + "\n")
    return updated_fam_dict


def main(XO_COUNTS, PLINK_FAM, OUT_FILE):
    """Driver function."""
    # Read in files
    xo_counts = read_xo_counts(XO_COUNTS)
    parents, fam = read_plink_fam(PLINK_FAM)
    # Add xo counts phenotype to fam dict
    updated_fam = add_xo_counts(OUT_FILE, fam, xo_counts)
    # Save to output file
    # Start from clean file, check if file exists
    if os.path.exists(os.path.expanduser(OUT_FILE)):
        os.remove(os.path.expanduser(OUT_FILE))
    # Start with parents first
    with open(os.path.expanduser(OUT_FILE), 'a') as file:
        for counter, i in enumerate(parents):
            file.write(" ".join(parents[i]) + "\n")
    # Now, add progeny
    with open(os.path.expanduser(OUT_FILE), 'a') as file:
        for counter, i in enumerate(updated_fam):
            file.write(" ".join(updated_fam[i]) + "\n")
    print("Done.")
    return

main(sys.argv[1], sys.argv[2], sys.argv[3]) # Run the program


