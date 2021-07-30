#!/usr/bin/env python3
"""This script takes in a single family formatted as a PED file and performs
checks for violations of Mendelian inheritance.

Usage: ./mendel_checking.py [PED] [crosstype] [out_dir]

Where:
1) [PED] is a single family formatted as a PED file
2) [crosstype] is the crossing scheme
3) [out_dir] is the full filepath to where the mendel reports are stored
"""

import sys
import os


def parse_ped(pedfile):
    """Parse the PED file and store it in a dictionary of
    {lineid: [PED_data]}"""
    ped_data = {}
    with open(pedfile, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            lid = tmp[1]
            # add the line to the ped_data
            ped_data[lid] = tmp
    return (ped_data)


def main(pedfile, crosstype, out_dir):
    """Driver function."""
    # Read in data
    ped_data = parse_ped(os.path.expanduser(pedfile))
    # Add calls to function here to perform Mendel checking
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
