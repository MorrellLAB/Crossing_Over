#!/usr/bin/env python3
"""This script takes in a single family formatted as a PED file and fills in
missing parental allele calls based on progeny allele calls.

Usage: ./PED_parental_fill-in.py [PED] [output_PED] [output_file]

Where:
1) [PED] is a single family formatted as a PED file
2) [output_PED] is the full path to output PED file we will save our filled-in parents to
3) [output_file] is the full path to output file used to track fill-in indices
"""

import sys
import numpy as np


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


def parent_fill_in(ped_data):
    """Fill-in missing parental allele calls based on progeny allele calls
    where possible."""
    # Extract LineIDs of parents in each family
    for key in ped_data.keys():
        if ped_data[key][4] == '1':
            dad = key
        elif ped_data[key][4] == '2':
            mom = key
    out_df = ped_data
    # Initialize lists to store index of missing calls
    # These will be saved into output file to track who we filled-in and where
    mom_idx = []
    dad_idx = []
    # Identify index in parent genotypes where we have missing calls
    for i in range(6, len(ped_data[key][6:])):
        if ped_data[mom][i] == '0':
            mom_idx.append(i)
        if ped_data[dad][i] == '0':
            dad_idx.append(i)
    mom_df = ped_data[mom][0:4] + [','.join(map(str, mom_idx))]
    dad_df = ped_data[dad][0:4] + [','.join(map(str, dad_idx))]
    idx = list(set().union(mom_idx, dad_idx))
    # Fill-in parental calls based on progeny calls
    for i in idx:
        tmp_progeny = []
        # Check progeny genotypes to determine possible parental calls
        for key in ped_data.keys():
            if (key == mom) ^ (key == dad):
                continue
            if not ped_data[key][i] in tmp_progeny:
                tmp_progeny.append(ped_data[key][i])
        # If we are missing one parent, fill in for one parent
        if (ped_data[mom][i] == '0') ^ (ped_data[dad][i] == '0'):
            # Mom missing, progeny not segregating
            if (ped_data[mom][i] == '0') & (len(tmp_progeny) == 1):
                out_df[mom][i] = tmp_progeny[0]
            # Dad missing, progeny not segregating
            elif (ped_data[dad][i] == '0') & (len(tmp_progeny) == 1):
                out_df[dad][i] = tmp_progeny[0]
            # Mom missing, progeny segregating
            elif (ped_data[mom][i] == '0') & (len(tmp_progeny) == 2):
                # Then need to check what call dad is
                dad_call_idx = tmp_progeny.index(ped_data[dad][i])
                fill_mom = np.setdiff1d(tmp_progeny, tmp_progeny[dad_call_idx])
                out_df[mom][i] = fill_mom[0]
            # Dad missing, progeny segregating
            elif (ped_data[dad][i] == '0') & (len(tmp_progeny) == 2):
                # Then need to check what call mom is
                mom_call_idx = tmp_progeny.index(ped_data[mom][i])
                fill_dad = np.setdiff1d(tmp_progeny, tmp_progeny[mom_call_idx])
                out_df[dad][i] = fill_dad[0]
            else:
                pass
        # If we are missing both parents, try to fill in both parents
        elif (ped_data[mom][i] == '0') & (ped_data[dad][i] == '0'):
            # Both missing, progeny not segregating
            if len(tmp_progeny) == 1:
                out_df[mom][i] = tmp_progeny[0]
                out_df[dad][i] = tmp_progeny[0]
            else:
                pass
        else:
            continue
    return (out_df, mom_df, dad_df)


def main(pedfile, outPED, outfile):
    """Driver function."""
    ped_data = parse_ped(pedfile)
    out_df, mom_df, dad_df = parent_fill_in(ped_data)
    # Save fill-in tracking info to output file
    header = ['FID', 'LID', 'FatherID', 'MotherID', 'FillInIdx']
    with open(outfile, 'a') as out:
        out.write('\t'.join(header) + '\n')
        out.write('\t'.join(mom_df) + '\n')
        out.write('\t'.join(dad_df) + '\n')

    # Save filled-in parental allele calls to PED file including progeny
    with open(outPED, 'a') as out_ped:
        for key in out_df.keys():
            if out_df[key][4] == '1':
                dad = key
            elif out_df[key][4] == '2':
                mom = key
        out_ped.write(' '.join(out_df[mom]) + '\n')
        out_ped.write(' '.join(out_df[dad]) + '\n')
        for key in out_df.keys():
            if (out_df[key] != out_df[mom]) & (out_df[key] != out_df[dad]):
                out_ped.write(' '.join(out_df[key]) + '\n')
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
