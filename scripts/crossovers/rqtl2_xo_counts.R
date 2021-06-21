#!/usr/bin/env Rscript

# Usage:
#   ./rqtl2_xo_counts.R [yaml_fp] [pcent_fp] [out_dir]

library(qtl2)
library(broman)
library(qtlcharts)
library(reshape2) # Use melt from this package to reshape multiply nested lists
library(ggplot2)
library(htmlwidgets)

readFile <- function(yaml_file) {
    # Read in files
    dat_cross2 <- read_cross2(yaml_file)
    return(dat_cross2)
}

readPericentromerePos <- function(pcent_file) {
    # Store centromere physical positions
    pcent <- read.table(file = pcent_file, sep = " ",
                        col.names = c("chr", "pstart", "pend", "type"))
    # Convert bp to Mbp
    pcent["pstartMb"] <- as.numeric(pcent$pstart)/1000000
    pcent["pendMb"] <- as.numeric(pcent$pend)/1000000
    # Strip "chr" characters
    pcent["chr"] <- gsub("chr", "", pcent$chr)
    pcent["chr"] <- gsub("H", "", pcent$chr)
    return(pcent)
}

plotMissing <- function(dat, samp_name, out_dir) {
    percent_missing <- n_missing(dat, by = "individual", summary = "proportion")*100
    labels <- paste0(names(percent_missing), " (", round(percent_missing), "%)")
    missing_plot <- iplot(seq_along(percent_missing), percent_missing, indID=labels,
          chartOpts=list(xlab="Individuals", ylab="Percent missing genotype data",
                         ylim=c(0, 100)))
    htmlwidgets::saveWidget(missing_plot, file=paste0(out_dir, "/", samp_name, "_percent_missing.html"))
    # Save data to file
    temp_df <- as.data.frame(percent_missing)
    out_df <- data.frame(ind=rownames(temp_df), pmissing=temp_df$percent_missing)
    write.table(
        x = out_df,
        file = paste0(out_dir, "/", samp_name, "_percent_missing.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
}

plotNumXO <- function(dat, btotxo, xaxis_title, samp_name, out_dir) {
    percent_missing <- n_missing(dat, by = "individual", summary = "proportion")*100
    numxo_plot <- iplot(seq_along(btotxo)[percent_missing < 19.97],
          btotxo[percent_missing < 19.97],
          chartOpts=list(
              xlab=xaxis_title, ylab="Number of crossovers",
              margin=list(left=80, top=40, right=40, bottom=40, inner=5),
              axispos=list(xtitle=25, ytitle=50, xlabel=5, ylabel=5))
          )
    htmlwidgets::saveWidget(numxo_plot, file=paste0(out_dir, "/", samp_name, "_num_xo.html"))
}

geneticMapPlotting <- function(dat, blxo, samp_name, out_dir) {
    # Reformat multiply nested lists
    lxodf_gmap <- melt(blxo)
    colnames(lxodf_gmap) <- c("xo_pos", "ind", "chr")
    physpos_chr <- melt(dat$pmap)
    physpos <- as.data.frame(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("chr", "snp", "pos"))))
    for (i in seq_along(dat$pmap)) {
        tmp_melt <- melt(dat$pmap[[i]])
        tmp_chr_pos <- data.frame(snp=rownames(tmp_melt), pos=tmp_melt[,1])
        chr <- rep(i, length(tmp_chr_pos$snp))
        new_df <- cbind(chr, tmp_chr_pos)
        physpos <- rbind(physpos, new_df)
    }
    # Reformat genetic map
    bgmap <- melt(dat$gmap)
    colnames(bgmap) <- c("pos", "chr")
    # Add column of 0 values
    bgmap['y'] <- rep(0.3, times=length(bgmap$pos))
    
    # Prep to make symbols a different color when crossovers land in the same location
    lxodf_gmap["duplicate"] <- duplicated(lxodf_gmap)
    # Change labels: TRUE = dup, FALSE = not_dup
    lxodf_gmap["duplicate"] <- gsub("TRUE", "> 1 XO", lxodf_gmap$duplicate)
    lxodf_gmap["duplicate"] <- gsub("FALSE", "1 XO", lxodf_gmap$duplicate)
    
    # Genetic map - Plotting all chr for each individual
    lpallind_gmap <- ggplot(lxodf_gmap, aes(x=xo_pos, y=as.numeric(chr)))
    lpallind_gmap + geom_hline(yintercept = 1:7, color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_color_manual(values=c("#ff9933", "#cc0000")) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25), legend.title = element_blank()) +
        scale_y_continuous(breaks=1:7) +
        facet_wrap(~ind) +
        xlab("Crossover Position (cM)") +
        ylab("Chromosome")
    # Save plot
    ggsave(filename = paste0(out_dir, "/", samp_name, "_gmap_xo_cM_by_ind.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 12, height = 10)
    
    # Genetic map - Plotting all individuals and facet wrap by chr
    lp_gmap <- ggplot(lxodf_gmap, aes(x=xo_pos, y=ind))
    lp_gmap + geom_hline(yintercept = 1:24, color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_colour_manual(values=c("#ff9933", "#cc0000")) +
        geom_point(data=bgmap, inherit.aes=FALSE, aes(x=pos, y=y),
                   shape=17, color="#ff9933", alpha=0.8) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25)) +
        facet_wrap(~chr) +
        scale_y_discrete() +
        theme(legend.title = element_blank()) +
        xlab("Crossover Position (cM)") +
        ylab("Individual")
    # Save plot
    ggsave(filename = paste0(out_dir, "/", samp_name, "_gmap_xo_cM_by_chr.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 12, height = 10)
}

# Plot using physical map positions
physicalMapPlotting <- function(dat, blxo_phys, pcent, samp_name, out_dir) {
    # Reformat multiply nested lists
    lxodf <- melt(blxo_phys)
    colnames(lxodf) <- c("xo_pos", "ind", "chr")
    # Reformat physical map
    bpmap <- melt(dat$pmap)
    colnames(bpmap) <- c("pos", "chr")
    # Add column of 0 values
    bpmap['y'] <- rep(0.3, times=length(bpmap$pos))
    
    # Prep to make symbols a different color when crossovers land in the same location
    lxodf["duplicate"] <- duplicated(lxodf)
    # Change labels: TRUE = dup, FALSE = not_dup
    lxodf["duplicate"] <- gsub("TRUE", "> 1 XO", lxodf$duplicate)
    lxodf["duplicate"] <- gsub("FALSE", "1 XO", lxodf$duplicate)
    
    
    # Physical map - Plotting all chr for each individual
    par(mfrow = c(1,1))
    lpallind <- ggplot(lxodf, aes(x=xo_pos, y=as.numeric(chr)))
    lpallind + geom_hline(yintercept = 1:7, color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_colour_manual(values=c("#00ccff", "blue")) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25), legend.title = element_blank()) +
        scale_y_continuous(breaks=1:7) +
        facet_wrap(~ind) +
        xlab("Crossover Position (Mbp)") +
        ylab("Chromosome")
    # Save plot
    ggsave(filename = paste0(out_dir, "/", samp_name, "_pmap_xo_Mbp_by_ind.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 12, height = 10)
    
    # Physical map - Plot all individuals and facet wrap by chr
    lp <- ggplot(lxodf, aes(x=xo_pos, y=ind))
    lp + geom_rect(data=pcent, inherit.aes=FALSE, 
                   mapping=aes(xmin = pstartMb, xmax = pendMb, ymin = 0, ymax = Inf), 
                   fill = "grey80", alpha = 0.7) +
        geom_hline(yintercept = 1:24, color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_colour_manual(values=c("#00ccff", "blue")) +
        geom_point(data=bpmap, inherit.aes=FALSE, aes(x=pos, y=y),
                   shape=17, color="#ff9933", alpha=0.8) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25)) +
        facet_wrap(~chr) +
        scale_y_discrete() +
        theme(legend.title = element_blank()) +
        xlab("Crossover Position (Mbp)") +
        ylab("Individual")
    # Save plot
    ggsave(filename = paste0(out_dir, "/", samp_name, "_pmap_xo_Mbp_by_chr.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 12, height = 10)
    
    return(lxodf)
}

makePhenoTable <- function(num_chr, lxodf, pcent, samp_name, out_dir) {
    pcolnames <- c("sampleID")
    # Create column names
    for (c in 1:num_chr) {
        lp <- paste("chr", c, "_LP", sep="")
        p <- paste("chr", c, "_P", sep="")
        rp <- paste("chr", c, "_RP", sep="")
        # Append to existing vector of column names
        pcolnames <- c(pcolnames, lp, p, rp)
    }
    
    pheno_df <- as.data.frame(matrix(nrow = 0, ncol = length(pcolnames),
                                     dimnames = list(NULL, c(pcolnames))))
    
    for (i in unique(lxodf$ind)) {
        tmp_row <- c(i)
        for (chrom in 1:7) {
            # Subset xo locations for current individual and current chromosome
            tmp_pchr <- lxodf[lxodf$ind == i & lxodf$chr == chrom, ]$xo_pos
            ############### Add feature ###############
            # Check if xo location is inbetween two markers that span a breakpoint
            # In this case, we can't be sure if the xo falls on one of the chr arms
            # or in the pericentromere.
            marker_pos <- blxo_phys[[chrom]][[i]]
            for (xopos in tmp_pchr) {
                print(xopos)
                # Identify the nearest 2 markers to the xo location
                distance_df <- data.frame(
                    xo_position = tmp_pchr[xopos],
                    marker_position = marker_pos,
                    distance = abs(tmp_pchr[1]-marker_pos))
                #order(distance_df$) ################## This is where I left off last
            }
            # Left of pericentromere
            lp <- sum(tmp_pchr <= pcent[pcent$chr == chrom, ]$pstartMb)
            # In pericentromere
            p <- sum(tmp_pchr > pcent[pcent$chr == chrom, ]$pstartMb & tmp_pchr < pcent[pcent$chr == chrom, ]$pendMb)
            # Right of pericentromere
            rp <- sum(tmp_pchr >= pcent[pcent$chr == chrom, ]$pendMb)
            tmp_row <- c(tmp_row, lp, p, rp)
        }
        pheno_df[i, ] <- matrix(tmp_row, ncol = length(pcolnames))
    }
    # Save phenotype table
    # write.table(x = pheno_df, 
    #             file = paste(out_dir, "/", samp_name, "_pheno.txt", sep = ""),
    #             quote = FALSE,
    #             sep = "\t",
    #             row.names = FALSE)
    return(pheno_df)
}

runXOAnalysis <- function(dat, pcent, samp_name, out_dir) {
    # Plot percent missing
    plotMissing(dat, samp_name, out_dir)
    # Look for sample duplicates
    cg <- compare_geno(dat, cores=0)
    capture.output(summary(cg), file = paste0(out_dir, "/", samp_name, "_duplicates_summary.txt"))
    
    # Calculate genotype probabilities first
    # Error probability based on SNP array error rate
    # Reported in: https://www.nature.com/articles/nmeth842 (Steemers et al. 2006 Nature Methods)
    bpr <- calc_genoprob(dat, error_prob=0.002, map_function="c-f", cores=0)
    # Identify most probable genotype at each position then count exchanges
    bm <- maxmarg(bpr, minprob=0.95, cores=0)
    
    # Crossover counts
    # Returns counts of crossovers on each chromosome (as columns) in each mouse
    bnxo <- count_xo(bm, cores=0)
    # Generate column names for per chr xo count
    bnxo_colnames <- gsub(pattern="^", replacement="chr", x=colnames(bnxo))
    colnames(bnxo) <- bnxo_colnames
    # Sum each row to get genome-wide estimates of total numbers of crossovers
    btotxo <- rowSums(bnxo)
    # Save # of crossovers data to file
    #percent_missing <- n_missing(dat, by = "individual", summary = "proportion")*100
    #temp_btotxo_pm <- btotxo[percent_missing < 19.97]
    temp_df <- as.data.frame(btotxo)
    temp_out_df <- data.frame(sampleID=rownames(temp_df), total_xo=temp_df$btotxo)
    xo_count_df <- cbind(temp_out_df, bnxo)
    rownames(xo_count_df) <- c()
    write.table(
        x = xo_count_df,
        file = paste0(out_dir, "/", samp_name, "_xo_count.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
    # Plot of # of crossovers
    plotNumXO(dat, btotxo, "F3 Barley Individuals", samp_name, out_dir)
    
    # Locate crossovers
    blxo <- locate_xo(bm, map = dat$gmap, cores = 0)
    # Try using physical map
    blxo_phys <- locate_xo(bm, map = dat$pmap, cores = 0)
    
    # Make crossover plots
    # Genetic map plot
    geneticMapPlotting(dat, blxo, samp_name, out_dir)
    # Physical map plot
    lxodf <- physicalMapPlotting(dat, blxo_phys, pcent, samp_name, out_dir)
    
    # Create and save phenotype table
    pheno_df <- makePhenoTable(num_chr = 7, lxodf, pcent, samp_name, out_dir)
    
    ############### Add feature ############### 
    # Count the number of intervals
    # If xo falls between two markers that span the boundary, set to missing
    #length(blxo_phys$`7`$`MS12_2184-001`)-1
    
    # Combine total xo counts and xo counts by chr with pheno table
    ###########################################
    
    # Save phenotype table
    write.table(x = pheno_df, 
                file = paste(out_dir, "/", samp_name, "_pheno.txt", sep = ""),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
}

main <- function() {
    # Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    # User provided command line arguments
    yaml_fp <- args[1]
    pcent_fp <- args[2]
    out_dir <- args[3]
    
    # Read in files
    dat <- readFile(yaml_fp)
    pcent <- readPericentromerePos(pcent_fp)
    
    # Do some processing
    # R/qtl2 complained about pmap not being sorted correctly, so manually sort
    for (i in seq_along(dat$pmap)) {
        dat$pmap[[i]] <- sort(dat$pmap[[i]])
    }
    checks <- check_cross2(dat)
    if (checks == TRUE) {
        print("Everything is correct.")
    } else {
        print("Check data.")
    }
    # Omit markers without any genotype data or noninformative genotypes
    dat_omit_null <- drop_nullmarkers(dat)
    
    # Extract prefix from filename
    yaml_prefix <- basename(yaml_fp)
    samp_name <- sub(pattern = "_forqtl2.yaml", replacement = "", x = yaml_prefix)
    
    # We need a min of 3 markers per chr, otherwise rqtl2 will return error
    # Check for at least 3 markers per chr
    temp_summary <- summary(dat_omit_null)
    temp_check <- sum(temp_summary$nmar < 3)
    if (temp_check != "0") {
        print(samp_name)
        print("Not enough markers per chromosome, saving sample to log file and exiting analysis")
        # Send summary output to a file
        capture.output(summary(dat_omit_null), file = paste0(out_dir, "/", "too_few_markers_per_chr-", samp_name, ".txt"))
    } else {
        # Proceed with analysis
        runXOAnalysis(dat_omit_null, pcent, samp_name, out_dir)
    }
}

# Run the program
main()
