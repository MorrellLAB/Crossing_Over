#!/usr/bin/env python3
"""This script takes in Plink 1.9 PED and MAP files and an AB genotype lookup
table. It reformats data to R/qtl2 required input formats. Script currently
outputs in your current working directory.

Usage: ./PED_parental_fill-in.py [ped] [map] [lookup_table] [out_fp]

Where:
1) [ped] is a Plink PED file
2) [map] is a Plink MAP file with physical positions in bp
    Note: this is important because script converts physical positions from
    bp to Mbp as required by R/qtl2
3) [lookup_table] is a file used to convert genotypes into AB genotypes
4) [out_fp] is the full path to the output file directory plus the output file prefix

Note: Currently converts intercrosses only. Doesn't generate phenotype CSV
files yet, but feature coming soon. Script assumes user is running Python
v3.7 or later due to the assumption that dictionaries maintain insertion order.
"""

import sys
import os
import re
import pandas as pd


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


def parse_map(mapfile):
    """Parse the MAP file and store it in a dictionary of
    {snpid: [MAP_data]}"""
    map_data = {}
    with open(mapfile, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            snpid = tmp[1]
            # add the line to map_data
            map_data[snpid] = tmp
    return (map_data)


def parse_lookup_table(lookupfile):
    """Parse the lookup table that contains AB genotypes."""
    ab_geno = {}
    with open(lookupfile, 'r') as l:
        for line in l:
            tmp = line.strip().split(',')
            snp = tmp[0]
            # add the line to ab_geno
            ab_geno[snp] = tmp
    return(ab_geno)


def make_geno(ped_dat):
    """Use PED data to create Rqtl2's required genotype data frame.
    If founder genotypes are present, this function will return two dataframes:
    1) progeny genotypes, 2) founder genotypes. Otherwise, only progeny
    genotypes dataframe will be returned."""
    # If parental genotypes are present (Sex code = 1 or 2)
    # Create separate founders genotype file
    founder_dict = {}
    progeny_dict = {}
    for key in ped_dat:
        if (ped_dat[key][4] == '1') or (ped_dat[key][4] == '2'):
            tmp_founder = ped_dat[key][6:]
            # concatenate alleles stored in separate columns into genotypes
            tmp_founder_even = tmp_founder[0::2]  # extract even
            tmp_founder_odd = tmp_founder[1::2]  # extract odd
            founder_geno = []
            for i in range(0, len(tmp_founder_even)):
                tmp_geno = tmp_founder_even[i] + tmp_founder_odd[i]
                founder_geno.append(tmp_geno)
            # add to founder dictionary
            founder_dict[key] = founder_geno
        else:
            tmp_progeny = ped_dat[key][6:]
            # concatenate alleles stored in separate columns into genotypes
            tmp_progeny_even = tmp_progeny[0::2]  # extract even
            tmp_progeny_odd = tmp_progeny[1::2]  # extract odd
            progeny_geno = []
            for i in range(0, len(tmp_progeny_even)):
                tmp_geno = tmp_progeny_even[i] + tmp_progeny_odd[i]
                progeny_geno.append(tmp_geno)
            # add to progeny dictionary
            progeny_dict[key] = progeny_geno
    return(founder_dict, progeny_dict)


def make_dataframe(geno_dat, map_dat):
    """Reformat PED founder/progeny genotypes that are stored in dictionaries
    into dataframe where row names are sampleIDs and column names are marker
    names."""
    geno_df = pd.DataFrame.from_dict(geno_dat, orient='index',
                                     columns=list(map_dat.keys()))
    return(geno_df)


def convert_to_abgeno(geno_df, lookup_table):
    """Use lookup table to convert genotypes to AB genotype
    representation. This function assumes lookup table has alleles
    corresponding to A and B alleles in columns 2 and 3."""
    abgeno_df = geno_df
    for snp_col in abgeno_df.columns:
        if snp_col in list(lookup_table.keys()):
            tmp_rep = []
            for i in range(0, len(abgeno_df[snp_col])):
                tmp_rep1 = re.sub(lookup_table[snp_col][1], 'A',
                                  abgeno_df[snp_col][i])
                tmp_rep2 = re.sub(lookup_table[snp_col][2], 'B', tmp_rep1)
                tmp_rep.append(tmp_rep2)
            # In place replacement of column containing AB genotypes
            abgeno_df[snp_col] = tmp_rep
        else:
            continue
    return(abgeno_df)


def recode_missing(geno_df):
    """Recode missing data from '0' to 'NA'. Plink codes missing genotypes
    using '0', but Rqtl2 does not allow '0' character to represent missing
    genotypes because in R, fread() considers '0' as type boolean."""
    for snp_col in geno_df.columns:
        geno_df[snp_col] = geno_df[snp_col].replace('00', 'NA')
    return(geno_df)


def make_pos_map(map_dat):
    """Create genetic map and physical map dataframes from
    map file."""
    gmap_dict = {}
    pmap_dict = {}
    for key in list(map_dat.keys()):
        # Reformat chromosomes
        if "chr" and "H" in map_dat[key][0]:
            tmp_chr = map_dat[key][0].strip("chr").strip("H")
        elif "chr" in map_dat[key][0]:
            tmp_chr = map_dat[key][0].strip("chr")
        else:
            tmp_chr = map_dat[key][0]
        # Physical positions need to be in Mbp
        tmp_phys = int(map_dat[key][3])/1000000
        # Reorder to: marker, chr, pos
        tmp_gmap = [map_dat[key][1], tmp_chr, map_dat[key][2]]
        tmp_pmap = [map_dat[key][1], tmp_chr, tmp_phys]
        gmap_dict[key] = tmp_gmap
        pmap_dict[key] = tmp_pmap
    # Store in pandas dataframe
    gmap_df = pd.DataFrame.from_dict(gmap_dict, orient='index',
                                     columns=['marker', 'chr', 'pos'])
    pmap_df = pd.DataFrame.from_dict(pmap_dict, orient='index',
                                     columns=['marker', 'chr', 'pos'])
    # Sort dataframes by position column
    gmap_df_sorted = gmap_df.sort_values(by=['chr', 'pos'])
    pmap_df_sorted = pmap_df.sort_values(by=['chr', 'pos'])
    return(gmap_df_sorted, pmap_df_sorted)


def write_to_csv(df, outfile_name):
    """Write pandas dataframe to output CSV file. CSV file format
    is the required input file format for R/qtl2."""
    df.to_csv(outfile_name, sep=',', header=True, na_rep='NA',
              line_terminator=os.linesep)
    return


def write_map_to_csv(df, outfile_name):
    """Write genetic/physical map pandas dataframe to output CSV file,
    exclude rownames from being written."""
    df.to_csv(outfile_name, sep=',', header=True, na_rep='NA',
              index=False, line_terminator=os.linesep)
    return


def main(pedfile, mapfile, lookupfile, out_fp):
    """Driver function."""
    # Read in data
    ped_data = parse_ped(os.path.expanduser(pedfile))
    map_data = parse_map(os.path.expanduser(mapfile))
    lookup_table = parse_lookup_table(os.path.expanduser(lookupfile))

    # 1) Create required geno.csv files for founders and progeny
    # Convert Plink PED data format into genotypes
    founder_data, progeny_data = make_geno(ped_data)
    # Use lookup table to convert genotypes into AB genotypes
    # Founder conversion to dataframe
    founder_df = make_dataframe(founder_data, map_data)
    # Progeny conversion to dataframe
    progeny_df = make_dataframe(progeny_data, map_data)
    # Translate genotypes to AB genotypes
    founder_ab = convert_to_abgeno(founder_df, lookup_table)
    progeny_ab = convert_to_abgeno(progeny_df, lookup_table)
    # Recode missing genotypes from '0' to 'NA'
    founder_abr = recode_missing(founder_ab)
    progeny_abr = recode_missing(progeny_ab)
    # Save dataframe to output CSV file
    write_to_csv(founder_abr, os.path.expanduser(out_fp) + "_founder_AB_geno.csv")
    write_to_csv(progeny_abr, os.path.expanduser(out_fp) + "_progeny_AB_geno.csv")

    # 2) Create gmap.csv file with genetic map cM positions
    # 3) Create pmap.csv file with physical positions
    gmap, pmap = make_pos_map(map_data)
    # Save gmap and pmap to output CSV file
    write_map_to_csv(gmap, os.path.expanduser(out_fp) + "_gmap.csv")
    write_map_to_csv(pmap, os.path.expanduser(out_fp) + "_pmap.csv")
    return

# Print usage message if we don't have enough input arguments
if len(sys.argv) < 2:
    print(__doc__)
    exit(1)
else:
    # Run the program
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
