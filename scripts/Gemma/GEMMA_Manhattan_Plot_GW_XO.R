#!/usr/bin/env Rscript
# Make a Manhattan plot of -log(p_lrt) from the GEMMA results. These will be
# a little weird because the population is highly structured and the founders
# are pretty closely related.

library(qqman)
library(dplyr)
library(ggrepel)

# Take command-line arguments
args <- commandArgs(trailingOnly = TRUE)
xoassoc_fp <- args[1]
pheno_name <- args[2]
#suggestive_line_val <- -log10(as.numeric(args[3]))
out_dir <- args[3]

# Prepare output subdirectory
# Remove trailing slash, it could mess with building filepaths
out_dir_fp <- gsub("/$", "", out_dir)
out_subdir <- paste0(out_dir_fp, "/plots")
# Create output directory if it doesn't exist already
if (!dir.exists(file.path(out_subdir))) {
    dir.create(file.path(out_subdir))
}

# Prepare output filename
out_fn <- paste0(out_subdir, "/", pheno_name, "_manhattan.pdf")

# Read the GEMMA association results
xo_pheno <- read.table(xoassoc_fp, header = TRUE)
# qqman "manhattan" function requires chromosomes to be numeric
#   strip the non-numeric characters from chromosome column
xo_pheno$chr <- as.numeric(gsub("[^0-9]", "", xo_pheno$chr))

# Adjust P-values for multiple comparisons
#xo_pheno$p_lrt_adjusted <- p.adjust(xo_pheno$p_lrt, method = "BH")

# Extract SNPs of interest
# We'll highlight the top hit in each chromosome as the default
#snps_above_suggline <- xo_pheno[-log10(xo_pheno$p_lrt) > -log10(1e-2), "rs"]
#snps_above_suggline <- xo_pheno[-log10(xo_pheno$p_lrt) > suggestive_line_val, "rs"]
snps_of_interest <- xo_pheno %>%
    group_by(chr) %>%
    filter(-log10(p_lrt)==max(-log10(p_lrt))) %>%
    ungroup %>%
    select(rs)

# Prepare the dataset
xo_pheno_ann <- xo_pheno %>%
    # Compute chromosome size
    group_by(chr) %>%
    summarise(chr_len=max(ps)) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(xo_pheno, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, ps) %>%
    mutate(BPcum=ps+tot) %>%
    
    # Add highlight and annotation information
    mutate(is_highlight=ifelse(rs %in% snps_of_interest$rs, "yes", "no")) #%>%
    #mutate(is_annotate=ifelse(-log10(p_lrt) > 2, "yes", "no"))

# Prepare x-axis
axis_df <- xo_pheno_ann %>% group_by(chr) %>%
    summarize(center=(max(BPcum) + min(BPcum)) / 2)

# Generate the plot
ggplot(xo_pheno_ann, aes(x=BPcum, y=-log10(p_lrt))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey10", "grey60"), 22 )) +
    
    ## Add suggestive line
    #geom_hline(aes(yintercept=suggestive_line_val), colour="blue") +
    
    # custom X axis:
    scale_x_continuous( label = axis_df$chr, breaks = axis_df$center ) +
    # remove space between plot area and x axis
    scale_y_continuous(expand = c(0, 0), limits=c(0, max(-log10(xo_pheno_ann$p_lrt)+0.5))) +
    
    # Add highlighted points
    geom_point(data=subset(xo_pheno_ann, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data=subset(xo_pheno_ann, is_highlight=="yes"), 
                    aes(label=rs), size=3.8, fontface="bold",
                    box.padding=0.3,
                    # Do not repel away from left edge, only right edge
                    xlim=c(-Inf, NA),
                    # Do not repel from top or bottom edges
                    ylim=c(-Inf, Inf)) +

    # Custom the theme:
    theme_bw() +
    theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    labs(x="Chromosome")

# Save plot to file
ggsave(filename=out_fn)

#--------------------------
# Alternative manhattan plotting approach currently unused because it is less flexible
# manhattan(xo_pheno, chr="chr", bp="ps", p="p_lrt", snp="rs", chrlabs=c(unique(as.character(xo_pheno$chr))),
#           logp=TRUE,
#           alpha(c("gray10", "skyblue"), 0.7),
#           suggestiveline=-log10(1e-05),
#           genomewideline=-log10(5e-08),
#           annotatePval=-log10(1e-05),
#           annotateTop=TRUE)
