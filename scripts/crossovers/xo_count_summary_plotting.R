#!/usr/bin/env Rscript

library(qtl2)
library(data.table)
library(plyr)
library(ggplot2)

# Take command line arguments
args <- commandArgs(trailingOnly = TRUE)
# User provided command line arugments
xo_data_dir <- args[1]
out_dir <- args[2]
pedigree_fp <- args[3] # Full pedigree file with 3 columns: family_id, parent1, parent2
yaml_dir <- args[4]
# CSV file containing comma separated list of family_id prefixes used to group plots where
# Column 1: family_id prefixes that match with filenames
# Column 2: Plot name to use for each family
# Column 3: Output file name excluding file extension
#   Example: Genomic Prediction population cycles 1-3
#     Cycle 1 prefix,plot name: MS10,Cycle 1 Families,c1_nxo_summary_plot
#     Cycle 2 prefix,plot name: MS11,Cycle 2 Families,c2_nxo_summary_plot
#     Cycle 3 prefix,plot name: MS12,Cycle 3 Families,c3_nxo_summary_plot
#   See toy dataset family_plot_groupings.csv for more details
plot_groupings_fp <- args[5]

#----------------------------------

# Read in plot groupings
plot_groupings_df <- read.csv(file = plot_groupings_fp, header = FALSE)
# Change column names
colnames(plot_groupings_df) <- c("family_id_prefix", "plot_group_name", "output_filename")

# Read in full pedigree
pedigree <- read.csv(file = pedigree_fp, header = FALSE)
# Change column names
colnames(pedigree) <- c("family_id", "parent1", "parent2")

# Prepping filepaths
files <- list.files(path = xo_data_dir, pattern = ".txt")
fp <- paste0(xo_data_dir, "/", files)

# Read in SNP data for each family
# Match end of string only
yaml_files <- list.files(path = yaml_dir, pattern = ".yaml$")
yaml_fp <- paste0(yaml_dir, "/", yaml_files)
dcross2 <- list()
datcross2 <- list()
total_markers <- list()
for (i in 1:length(yaml_fp)) {
    yaml_prefix <- basename(yaml_fp[i])
    fam_name <- sub(pattern = "_forqtl2.yaml", replacement = "", x = yaml_prefix)
    # Read in files
    dcross2[[i]] <- read_cross2(yaml_fp[i])
    # R/qtl2 complained about pmap not being sorted correctly, so manually sort
    for (j in seq_along(dcross2[[i]]$pmap)) {
        dcross2[[i]]$pmap[[j]] <- sort(dcross2[[i]]$pmap[[j]])
    }
    # Omit markers without any genotype data or noninformative genotypes
    datcross2[[i]] <- drop_nullmarkers(dcross2[[i]])
    temp <- summary(datcross2[[i]])
    total_markers[[i]] <- data.frame(family=fam_name, tot_markers=temp$totmar)
}
# Check if physical map is now ordered correctly
check_cross2(datcross2[[1]])
check_cross2(datcross2[[10]])
tmp <- summary(datcross2[[1]])
tmp$nmar

# Read in coverage histogram file
xocount <- list()
# Using fread to read in file instead of read.table or read.delim2
# because fread is significantly faster than both
for (i in 1:length(fp)) {
    xocount[[i]] <- fread(file = fp[i], header = TRUE, sep = "\t")
}

#----------------------------------

# Restructure data into long format
xo_long <- list()
for (i in 1:length(fp)) {
    tmp_family <- strsplit(xocount[[i]]$sampleID, split = "-")
    tmp_fam_name <- tmp_family[[1]][1]
    xo_long[[i]] <- data.frame(
        ind = xocount[[i]]$sampleID,
        total_nxo = xocount[[i]]$total_xo,
        family = rep(tmp_fam_name, times = length(xocount[[i]]$total_xo)),
        all = rep("together", times = length(xocount[[i]]$total_xo))
    )
}
# Concatenate list of dataframes
all_xo_long <- ldply(xo_long, rbind)

# Generate boxplots for each grouping provided from plot_groupings_fp
MakeBoxplotsByGroup <- function(family_id_prefix, curr_xlab, output_dir, output_filename, df) {
    # Format prefix to work with grep
    #   so we only match the beginning of the string
    prefix <- paste0("^", family_id_prefix)
    # Extract counts for current grouping
    curr_group_xo_long <- df[grep(prefix, df$family), ]
    # Generate boxplot for current grouping
    curr_plot <- ggplot(curr_group_xo_long, aes(family, total_nxo)) +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(axis.title = element_text(size = 18)) +
        xlab(curr_xlab) + # Name taken from plot_groupings_fp column 2
        ylab("Number of Crossovers") +
        ylim(0, max(df$total_nxo) + 5) # scale all boxplots the same
    # Prepare output filepath
    output_fn <- paste0(output_filename, ".pdf")
    # Save to file
    ggsave(output_fn, plot = curr_plot, device = "pdf", path = output_dir, width = 10)
}

# Plotting
# For each grouping, generate boxplot
for (r in 1:nrow(plot_groupings_df)) {
    curr_group <- plot_groupings_df[r, ]
    # Make boxplot
    MakeBoxplotsByGroup(family_id_prefix = curr_group$family_id_prefix,
                        curr_xlab = curr_group$plot_group_name,
                        output_dir = out_dir,
                        output_filename = curr_group$output_filename,
                        df = all_xo_long)
}
