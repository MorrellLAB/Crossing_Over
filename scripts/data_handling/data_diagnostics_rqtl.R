#!/usr/bin/env Rscript

library(broman)
library(qtl2)
library(qtlcharts)
library(fst)

# Diagnostics below are considered basic diagnostics
# For additional custom diagnostics, please use the YAML control input file
#   created for the data diagnostics step to perform preferred custom diagnostic tests

PrepOutDir <- function(out_dir, subdir_name) {
  sub_dir <- paste0(out_dir, "/", subdir_name, sep = '')
  # Create output subdirectory if it doesn't exist
  if (!dir.exists(sub_dir)) {
    dir.create(file.path(sub_dir))
  }
  return(sub_dir)
}

ReadFile <- function(yaml_file) {
  # Read in files
  dat_cross2 <- read_cross2(yaml_file)
  return(dat_cross2)
}

PlotMissing <- function(dat_on, samp_name, out_dir) {
  percent_missing <- n_missing(dat_on, by = "individual", summary = "proportion")*100
  labels <- paste0(names(percent_missing), " (", round(percent_missing), "%)")
  missing_plot <- iplot(seq_along(percent_missing), percent_missing, indID=labels,
                        chartOpts=list(xlab="Individuals", ylab="Percent missing genotype data",
                                       ylim=c(0, 100)))
  # Prepare subdirectories
  out_dir_pmiss_html <- PrepOutDir(out_dir, "percent_missing_interactive_plots")
  out_dir_pmiss_txt <- PrepOutDir(out_dir, "percent_missing_data")
  htmlwidgets::saveWidget(missing_plot, file=paste0(out_dir_pmiss_html, "/", samp_name, "_percent_missing.html"))
  # Save data to file
  temp_df <- as.data.frame(percent_missing)
  out_df <- data.frame(ind=rownames(temp_df), pmissing=temp_df$percent_missing)
  write.table(
    x = out_df,
    file = paste0(out_dir_pmiss_txt, "/", samp_name, "_percent_missing.txt"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

PlotHeterozygotes <- function(dat_on, genopr, samp_name, out_dir) {
  # for (chrom in seq(length(dat_on$geno))) {
  #   # Heterozygotes are set as 2 (see make_control_file.R script)
  #   # 0 means missing
  #   phet <- rowSums(dat_on$geno[[chrom]] == 2) / rowSums(dat_on$geno[[chrom]] != 0)
  #   # Generate interactive plot per chromosome
  #   iplot(xint)
  # }
  
  # Check proportion of heterozygous genotypes
  hets_ind <- calc_het(probs = genopr, by = "individual", omit_x = TRUE)
  hets_mar <- calc_het(probs = genopr, by = "marker", omit_x = TRUE)
  # iplot(x = names(hets_ind), hets_ind)
}

# Helps us identify duplicate samples
PlotMatchingGenotypes <- function(dat_on, out_prefix, out_dir) {
  # Look for sample duplicates
  cgeno <- compare_geno(dat_on, cores=0)
  # Prepare subdirectories
  out_dir_dups <- PrepOutDir(out_dir, "other_summaries")
  capture.output(summary(cgeno), file = paste0(out_dir_dups, "/", samp_name, "_duplicates_summary.txt"))
  
  # Plot histogram of proportion of matching genotypes
  cgut <- cgeno[upper.tri(cgeno)]
  par(mar=c(5.1,0.6,0.6, 0.6))
  hist(cgut, breaks=seq(0, 1, length=201),
       main="", yaxt="n", ylab="",
       xlab="Proportion matching genotypes")
  rug(cgut)
  
  # Exclude samples with >50% missing genotypes and compare
  cgeno_sub <- cgeno[percent_missing < 50, percent_missing < 50]
  cgutsub <- cgeno_sub[upper.tri(cgeno_sub)]
  par(mar=c(5.1,0.6,0.6, 0.6))
  hist(cgutsub, breaks=seq(0, 1, length=201),
       main="", yaxt="n", ylab="",
       xlab="Proportion matching genotypes")
}

RunDataDiagnostics <- function(dat_on, out_prefix, out_dir) {
  # Set the total number of chromosomes
  n_chrom <- length(dat_on$is_x_chr)
  
  # Plot percent missing
  PlotMissing(dat_on, out_prefix, out_dir)
  
  # Look for sample duplicates
  
  
  
  
  # Calculate genotype probabilities
  # Error probability based on SNP array error rate
  # Reported in: https://www.nature.com/articles/nmeth842 (Steemers et al. 2006 Nature Methods)
  genopr <- calc_genoprob(dat_on, error_prob=0.002, map_function="c-f", cores=0)
  
  
}

Main <- function() {
  # Take command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  yaml_fp <- args[1]
  out_dir <- args[2]
  
  # Read in files
  dat <- ReadFile(yaml_fp)
  # Add code here to check if there were warnings, if so print them
  
  # Sometimes, R/qtl2 will complain about the pman not being sorted correctly
  # So, we'll add a sorting step here
  for (i in seq_along(dat$pmap)) {
    dat$pmap[[i]] <- sort(dat$pmap[[i]])
  }
  # Check if there are additional issues with the cross object
  checks <- check_cross2(dat)
  if (checks == TRUE) {
    print("Everything is correct.")
  } else {
    print("Check data.")
  }
  # Omit markers without any genotype data
  dat_omit_null <- drop_nullmarkers(dat)
  # Print cross object summary
  dat_omit_null
  
  # Extract prefix from filename
  yaml_prefix <- basename(yaml_fp)
  out_prefix <- sub(pattern = "_forqtl2.yaml", replacement = "", x = yaml_prefix)
  
  # We need a min of 3 markers per chr, otherwise rqtl2 will return error
  # Check for at least 3 markers per chr
  temp_summary <- summary(dat_omit_null)
  temp_check <- sum(temp_summary$nmar < 3)
  if (temp_check != "0") {
    print(samp_name)
    print("Not enough markers per chromosome, saving sample to log file and exiting analysis")
    # Prep subdirectories
    out_dir_too_few <- PrepOutDir(out_dir, "too_few_markers")
    # Send summary output to a file
    capture.output(summary(dat_omit_null), file = paste0(out_dir_too_few, "/", "too_few_markers_per_chr-", samp_name, ".txt"))
  } else {
    # Proceed with analysis
    RunDataDiagnostics(dat_omit_null, out_prefix, out_dir)
  }
}
