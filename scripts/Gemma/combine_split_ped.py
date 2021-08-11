#!/usr/bin/env python3
"""This script combines cleaned split by family PED file into a single PED
file for downstream processing. The script will print final outputs to
stdout and the user will need to redirect the outputs to a file to save

Usage: ./combine_split_ped.py [founder_ped] [split_ped_list] > /path/to/user_defined/output_file.ped

Where:
1) [founder_ped] is the full filepath to the founders only PED file
2) [split_ped_list] is the full filepath to a list containing filepaths to cleaned split by family PED files
"""

import os
import sys


def parse_ped_list(split_ped_list):
    """Read in list of filepaths to split by family
    PED files."""
    ped_list = []
    with open(split_ped_list, "rt") as handle:
        for record in handle:
            ped_list.append(record.strip())
    return ped_list


def parse_ped(plink_ped_fp):
    """Read in plink .ped file and do some processing
    before printing output to stdout."""
    with open(plink_ped_fp, "rt") as handle:
        for record in handle:
            if record.startswith("-9"):
                print(record.strip())
            else:
                # Do some processing to exclude parents in each
                #   family since the parents are redundant info
                #   given how we initially split the single PED
                #   by families
                tmp = record.strip().split()
                # Given how we originally split the PED by family
                #   column 3 (Within-family ID of father) and
                #   column 4 (Within-family ID of mother) should both
                #   be 0. We can check for this here, but remember script
                #   is 0-based indexing.
                if tmp[2] == "0" and tmp[3] == "0":
                    # These are the parents in each family PED file, exclude these in the final PED
                    continue
                else:
                    # Only print progeny
                    print(record.strip())


def main(founder_ped_fp, split_ped_list):
    """Driver function."""
    # Read in file list
    ped_list = parse_ped_list(os.path.expanduser(split_ped_list))

    # Founder lines only
    parse_ped(os.path.expanduser(founder_ped_fp))
    # Progeny split into families
    for fp in ped_list:
        parse_ped(fp)
    return

# Print usage message if we don't have enough input arguments
if len(sys.argv) < 2:
    print(__doc__)
    exit(1)
else:
    # Run the program
    main(sys.argv[1], sys.argv[2])
