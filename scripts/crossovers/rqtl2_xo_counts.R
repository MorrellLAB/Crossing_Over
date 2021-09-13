#!/usr/bin/env Rscript

# Usage:
#   ./rqtl2_xo_counts.R [yaml_fp] [pcent_fp] [out_dir]

suppressWarnings(suppressPackageStartupMessages(library(qtl2)))
suppressWarnings(suppressPackageStartupMessages(library(broman)))
suppressWarnings(suppressPackageStartupMessages(library(qtlcharts)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2))) # Use melt from this package to reshape multiply nested lists
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(htmlwidgets)))

PrepOutDir <- function(out_dir, subdir_name) {
    sub_dir <- paste0(out_dir, "/", subdir_name, sep = '')
    # Create output subdirectory if it doesn't exist
    if (!dir.exists(sub_dir)) {
        dir.create(file.path(sub_dir))
    }
    return(sub_dir)
}

OrderByChrom <- function(curr_list) {
    # Order chromosomes from smallest to largest for lists of lists in cross object
    curr_list <- curr_list[order(as.numeric(names(curr_list)))]
    return(curr_list)
}

ReadFile <- function(yaml_file) {
    # Read in files
    dat_cross2 <- read_cross2(yaml_file)
    # Order chromosomes from smallest to largest for lists of lists in cross object
    dat_cross2$geno <- OrderByChrom(dat_cross2$geno)
    dat_cross2$gmap <- OrderByChrom(dat_cross2$gmap)
    dat_cross2$pmap <- OrderByChrom(dat_cross2$pmap)
    dat_cross2$founder_geno <- OrderByChrom(dat_cross2$founder_geno)
    dat_cross2$is_x_chr <- OrderByChrom(dat_cross2$is_x_chr)
    return(dat_cross2)
}

ReadPericentromerePos <- function(pcent_file) {
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

PlotMissing <- function(dat, samp_name, out_dir) {
    percent_missing <- n_missing(dat, by = "individual", summary = "proportion")*100
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

PlotNumXO <- function(dat, btotxo, xaxis_title, samp_name, out_dir) {
    percent_missing <- n_missing(dat, by = "individual", summary = "proportion")*100
    numxo_plot <- iplot(seq_along(btotxo)[percent_missing < 19.97],
          btotxo[percent_missing < 19.97],
          chartOpts=list(
              xlab=xaxis_title, ylab="Number of crossovers",
              margin=list(left=80, top=40, right=40, bottom=40, inner=5),
              axispos=list(xtitle=25, ytitle=50, xlabel=5, ylabel=5))
          )
    # Prepare subdirectories
    out_dir_num_xo_html <- PrepOutDir(out_dir, "num_xo_interactive_plots")
    htmlwidgets::saveWidget(numxo_plot, file=paste0(out_dir_num_xo_html, "/", samp_name, "_num_xo.html"))
}

# Reformat marker physical positions so they are easier to work with
PullMarkerPhysPos <- function(dat) {
    physpos_chr <- melt(dat$pmap)
    # Prepare column headers
    physpos <- as.data.frame(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("chr", "snp", "pos"))))
    # Iteratively build dataframe containing markers and their physical positions
    for (i in seq_along(dat$pmap)) {
        tmp_melt <- melt(dat$pmap[[i]])
        tmp_chr_pos <- data.frame(snp=rownames(tmp_melt), pos=tmp_melt[,1])
        chr <- rep(i, length(tmp_chr_pos$snp))
        new_df <- cbind(chr, tmp_chr_pos)
        physpos <- rbind(physpos, new_df)
    }
    return(physpos)
}

GeneticMapPlotting <- function(dat, blxo, samp_name, out_dir) {
    # Reformat multiply nested lists
    lxodf_gmap <- melt(blxo)
    colnames(lxodf_gmap) <- c("xo_pos", "ind", "chr")
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
    
    # Prepare output subdirectories
    out_dir_by_ind <- paste0(out_dir, "/gmap_by_ind", sep = '')
    if (!dir.exists(out_dir_by_ind)) {
        dir.create(file.path(out_dir_by_ind)) # Create output subdirectory
    }
    out_dir_by_chr <- paste0(out_dir, "/gmap_by_chr", sep = '')
    if (!dir.exists(out_dir_by_chr)) {
        dir.create(file.path(out_dir_by_chr)) # Create output subdirectory
    }
    
    # Prepare chromosome labels
    chrom_vector <- as.vector(as.numeric(names(blxo)))
    first_chr <- chrom_vector[1] # Pull first chromosome
    last_chr <- chrom_vector[length(chrom_vector)] # Pull last chromosome
    # Genetic map - Plotting all chr for each individual
    lpallind_gmap <- ggplot(lxodf_gmap, aes(x=xo_pos, y=as.numeric(chr)))
    lpallind_gmap + geom_hline(yintercept = first_chr:last_chr, color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_color_manual(values=c("#ff9933", "#cc0000")) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25), legend.title = element_blank()) +
        scale_y_continuous(breaks=first_chr:last_chr) +
        facet_wrap(~ind) +
        xlab("Crossover Position (cM)") +
        ylab("Chromosome")
    # Save plot
    ggsave(filename = paste0(out_dir_by_ind, "/", samp_name, "_gmap_xo_cM_by_ind.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 16, height = 14) # Adjust this to scale automatically if possible
    
    # Genetic map - Plotting all individuals and facet wrap by chr
    lp_gmap <- ggplot(lxodf_gmap, aes(x=xo_pos, y=ind))
    lp_gmap + geom_hline(yintercept = 1:length(unique(lxodf_gmap$ind)), color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_colour_manual(values=c("#ff9933", "#cc0000")) +
        geom_point(data=bgmap, inherit.aes=FALSE, aes(x=pos, y=y),
                   shape=17, color="#ff9933", alpha=0.8) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25)) +
        facet_wrap(~as.numeric(chr)) +
        scale_y_discrete() +
        theme(legend.title = element_blank()) +
        xlab("Crossover Position (cM)") +
        ylab("Individual")
    # Save plot
    ggsave(filename = paste0(out_dir_by_chr, "/", samp_name, "_gmap_xo_cM_by_chr.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 26, height = 26) # Adjust this to scale automatically if possible
}

# Plot using physical map positions
PhysicalMapPlotting <- function(dat, blxo_phys, pcent, samp_name, out_dir) {
    # Reformat multiply nested lists
    lxodf <- melt(blxo_phys)
    colnames(lxodf) <- c("xo_pos", "ind", "chr")
    # Reformat physical map for plotting
    bpmap <- melt(dat$pmap)
    colnames(bpmap) <- c("pos", "chr")
    # Add column of 0 values
    bpmap['y'] <- rep(0.3, times=length(bpmap$pos))
    
    # Prep to make symbols a different color when crossovers land in the same location
    lxodf["duplicate"] <- duplicated(lxodf)
    # Change labels: TRUE = dup, FALSE = not_dup
    lxodf["duplicate"] <- gsub("TRUE", "> 1 XO", lxodf$duplicate)
    lxodf["duplicate"] <- gsub("FALSE", "1 XO", lxodf$duplicate)
    
    # Prepare output subdirectories
    out_dir_by_ind <- paste0(out_dir, "/pmap_by_ind", sep = '')
    if (!dir.exists(out_dir_by_ind)) {
        dir.create(file.path(out_dir_by_ind)) # Create output subdirectory
    }
    out_dir_by_chr <- paste0(out_dir, "/pmap_by_chr", sep = '')
    if (!dir.exists(out_dir_by_chr)) {
        dir.create(file.path(out_dir_by_chr)) # Create output subdirectory
    }
    
    # Prepare chromosome labels
    chrom_vector <- as.vector(as.numeric(names(blxo_phys)))
    first_chr <- chrom_vector[1] # Pull first chromosome
    last_chr <- chrom_vector[length(chrom_vector)] # Pull last chromosome
    # Physical map - Plotting all chr for each individual
    par(mfrow = c(1,1))
    lpallind <- ggplot(lxodf, aes(x=xo_pos, y=as.numeric(chr)))
    lpallind + geom_hline(yintercept = first_chr:last_chr, color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_colour_manual(values=c("#00ccff", "blue")) +
        theme_classic() +
        theme(axis.line = element_line(size = 0.25), legend.title = element_blank()) +
        scale_y_continuous(breaks=first_chr:last_chr) +
        facet_wrap(~ind) +
        xlab("Crossover Position (Mbp)") +
        ylab("Chromosome")
    # Save plot
    ggsave(filename = paste0(out_dir_by_ind, "/", samp_name, "_pmap_xo_Mbp_by_ind.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 16, height = 14) # Adjust this to scale automatically if possible
    
    # Physical map - Plot all individuals and facet wrap by chr
    lp <- ggplot(lxodf, aes(x=xo_pos, y=ind))
    lp + geom_rect(data=pcent, inherit.aes=FALSE, 
                   mapping=aes(xmin = pstartMb, xmax = pendMb, ymin = 0, ymax = Inf), 
                   fill = "grey80", alpha = 0.7) +
        geom_hline(yintercept = 1:length(unique(lxodf$ind)), color = "grey") +
        geom_point(shape = 4, aes(colour=factor(duplicate))) +
        scale_colour_manual(values=c("#00ccff", "blue")) +
        geom_point(data=bpmap, inherit.aes=FALSE, aes(x=pos, y=y),
                   shape=17, color="#ff9933", alpha=0.8) + # Markers shown as triangles
        theme_classic() +
        theme(axis.line = element_line(size = 0.25)) +
        facet_wrap(~as.numeric(chr)) +
        scale_y_discrete() +
        theme(legend.title = element_blank()) +
        xlab("Crossover Position (Mbp)") +
        ylab("Individual")
    # Save plot
    ggsave(filename = paste0(out_dir_by_chr, "/", samp_name, "_pmap_xo_Mbp_by_chr.pdf"),
           plot = last_plot(),
           device = "pdf",
           width = 26, height = 26) # Adjust this to scale automatically if possible
    
    return(lxodf)
}

# Find the 2 closest markers just upstream and downstream of the crossover
# Function works on a chromosome by chromosome basis, and one xo position at a time
ClosestFlankingMarkers <- function(marker_pos, curr_xo_pos, curr_chrom, curr_indv) {
    # Pull markers for current chromosome
    curr_markers <- marker_pos[marker_pos$chr == curr_chrom, ]
    # Calculate distances of markers from current xo position
    curr_markers$curr_xo_dist <- abs(curr_xo_pos - curr_markers$pos)
    # Sort by distances of markers from current xo position
    sorted_curr_markers <- curr_markers[order(curr_markers$curr_xo_dist), ]
    
    # Get closest left flanking position
    # Note: may need to deal with case where marker  has the same position as the xo position
    #   although, I'm not sure if rqtl2 actually allows this case when inferring crossovers
    if (any(sorted_curr_markers$pos < curr_xo_pos)) {
        left_of_xo <- sorted_curr_markers[sorted_curr_markers$pos < curr_xo_pos, ]
        left_pos <- left_of_xo[which(left_of_xo$curr_xo_dist == min(left_of_xo$curr_xo_dist)), ]
    } else {
        # We don't have a marker position to the left of the current xo position
        message("Edge case...")
        message("Current individual, Current chromosome, Current XO position:")
        message(c(curr_indv, as.character(curr_chrom), curr_xo_pos))
        stop("There is no marker position to the left of the current xo position.")
    }
    # Get closest right flanking position
    if (any(sorted_curr_markers$pos > curr_xo_pos)) {
        right_of_xo <- sorted_curr_markers[sorted_curr_markers$pos > curr_xo_pos, ]
        right_pos <- right_of_xo[which(right_of_xo$curr_xo_dist == min(right_of_xo$curr_xo_dist)), ]
    } else {
        # We don't have a marker position to the right of the current xo position
        message("Edge case...")
        message("Current individual, Current chromosome, Current XO position:")
        message(c(curr_indv, as.character(curr_chrom), curr_xo_pos))
        stop("There is no marker position to the right of the current xo position.")
    }
    
    # Add check for cases where two or more markers have the same physical position.
    # The user should pick the marker with the least amount of missing data as part of the data cleaning step
    if (length(left_pos$pos) > 1) {
        message("Two or more markers have the same physical position.")
        message("Current xo position:", curr_xo_pos)
        message("Current chromosome:", curr_chrom)
        message("Current individual:", curr_indv)
        message("Current left positions processing:", left_pos)
        stop("Please investigate before proceeding")
    }
    
    if (length(right_pos$pos) > 1) {
        message("Two or more markers have the same physical position. Please investigate.")
        message("Current xo position:", curr_xo_pos)
        message("Current chromosome:", curr_chrom)
        message("Current individual:", curr_indv)
        message("Current left positions processing:", left_pos)
        stop("Please investigate before proceeding")
    }
    
    # Check that current xo position falls between closest markers
    # Left position should be smaller than current xo position
    if (left_pos$pos > curr_xo_pos) {
        # If left position isn't smaller than curr_xo_pos, stop and return message
        stop("Left flanking position > current crossover position, this isn't right. 
             Please investigate this more closely before proceeding.")
    }
    # Right position should be greater than current xo position
    if (right_pos$pos < curr_xo_pos) {
        # If right position isn't greater than curr_xo_pos, stop and return message
        stop("Right flanking position < current crossover position, this isn't right.
             Please investigate this more closely before proceeding.")
    }
    
    # Combine closest markers into a single data frame
    closest_markers <- rbind(left_pos, right_pos)
    return(closest_markers)
}

# Categorize closest markers as:
#   LP (left of pericentromere), P (pericentromere), RP (right of pericentromere)
CategorizeMarkers <- function(closest_markers, curr_pcent) {
    # Prepare column names of output data frame
    df <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df) <- c(colnames(closest_markers), "pcent_cat")
    
    for (mrow in 1:nrow(closest_markers)) {
        curr_row <- closest_markers[mrow, ]
        if (curr_row$pos <= curr_pcent$pstartMb) {
            curr_row_cat <- "LP"
        } else if (curr_row$pos > curr_pcent$pstartMb & curr_row$pos < curr_pcent$pendMb) {
            curr_row_cat <- "P"
        } else if (curr_row$pos >= curr_pcent$pendMb) {
            curr_row_cat <- "RP"
        }
        new_row <- data.frame(curr_row, pcent_cat = curr_row_cat)
        df <- rbind(df, new_row)
    }
    return(df)
}

MakePhenoTable <- function(dat, num_chr, lxodf, pcent, samp_name, out_dir) {
    # Reformat physical positions
    marker_pos <- PullMarkerPhysPos(dat)
    
    pcolnames <- c("sampleID")
    # Create column names
    for (c in 1:num_chr) {
        lp <- paste("chr", c, "_LP", sep="")
        p <- paste("chr", c, "_P", sep="")
        rp <- paste("chr", c, "_RP", sep="")
        # Append to existing vector of column names
        pcolnames <- c(pcolnames, lp, p, rp)
    }
    
    # Prepare data frame to store phenotype table
    pheno_df <- as.data.frame(matrix(nrow = 0, ncol = length(pcolnames),
                                     dimnames = list(NULL, c(pcolnames))))
    new_lxodf <- as.data.frame(matrix(nrow = 0, ncol = ncol(lxodf)))
    colnames(new_lxodf) <- colnames(lxodf)
    
    for (i in unique(lxodf$ind)) {
        tmp_row <- c(i)
        for (chrom in 1:num_chr) {
            # Subset xo locations for current individual and current chromosome
            tmp_pchr <- lxodf[lxodf$ind == i & lxodf$chr == chrom, ]$xo_pos
            # Check if xo location is inbetween two markers that span a breakpoint
            # In this case, we can't be sure if the xo falls on one of the chr arms
            #   or in the pericentromere.
            # The counts are based on the number of crossovers stored in tmp_pchr
            #   so, we'll remove xo positions from tmp_pchr if they are set to missing
            new_tmp_pchr <- c()
            curr_pcent <- pcent[pcent$chr == chrom, ] # pericentromere for current chromosome
            for (xopos in tmp_pchr) {
                # Identify the nearest 2 markers to the xo location
                closest_markers <- ClosestFlankingMarkers(marker_pos, curr_xo_pos=xopos, 
                                                            curr_chrom=chrom, curr_indv = i)
                # Categorize closest markers into LP, P, RP
                marker_categories <- CategorizeMarkers(closest_markers, curr_pcent)
                # If the two flanking markers for a xo position has the same pcent_cat, then
                #   we are able to clearly place the xo on LP, P, or RP
                # If the two flanking markers for a xo position don't have the same pcent_cat, then
                #   the xo location is between two markers that span a breakpoint and can't be categorized
                #   as LP, P, or RP with certainty. We'll set this to missing and not count this XO.
                if (length(unique(sort(marker_categories$pcent_cat))) == 1) {
                    # We are able to clearly place the xo on LP, P, or RP
                    new_tmp_pchr <- c(new_tmp_pchr, xopos)
                    # Get current row from lxodf data frame
                    curr_lxodf_row <- lxodf[lxodf$ind == i & lxodf$chr == chrom & lxodf$xo_pos == xopos, ]
                    # This will be used for plotting purposes after setting xo to missing
                    new_lxodf <- rbind(new_lxodf, curr_lxodf_row)
                }
            }
            
            # Get the start position of the pericentromere in Mbp
            pcent_start_Mb <- pcent[pcent$chr == chrom, ]$pstartMb
            # Get the end position of the pericentromere in Mbp
            pcent_end_Mb <- pcent[pcent$chr == chrom, ]$pendMb
            
            # XO Count - Left of pericentromere
            lp <- sum(new_tmp_pchr <= pcent_start_Mb)
            # XO Count - In pericentromere
            p <- sum(new_tmp_pchr > pcent_start_Mb & new_tmp_pchr < pcent_end_Mb)
            # XO Count - Right of pericentromere
            rp <- sum(new_tmp_pchr >= pcent_end_Mb)
            # Combine counts into single row
            tmp_row <- c(tmp_row, lp, p, rp)
        }
        pheno_df[i, ] <- matrix(tmp_row, ncol = length(pcolnames))
    }
    outputs <- list(pheno_df, new_lxodf)
    return(outputs)
}

# After setting ambiguous crossovers to missing, we'll want to re-do our
#   total/sample and total/chromosome counts
UpdateCounts <- function(pheno_df, n_chrom) {
    new_pheno_df <- data.frame(sampleID = pheno_df$sampleID)
    sample_xo_total <- pheno_df %>% dplyr::select(starts_with("chr")) %>%
        mutate_if(is.character, as.numeric) %>%
        mutate(sample_total = rowSums(.)) %>%
        select(sample_total)
    rownames(new_pheno_df) <- rownames(pheno_df)
    new_pheno_df <- cbind(new_pheno_df, sample_xo_total)
    
    for (c in 1:n_chrom) {
        curr_chr <- paste0("chr", c)
        curr_chr_sum <- pheno_df %>% dplyr::select(starts_with(curr_chr)) %>%
            mutate_if(is.character, as.numeric) %>%
            mutate(chr_xo_count = rowSums(.))
        colnames(curr_chr_sum)[4] <- paste0(curr_chr, "_total")
        # Add columns to new phenotype table
        new_pheno_df <- cbind(new_pheno_df, curr_chr_sum)
    }
    
    # Sort by sampleID
    new_pheno_df <- new_pheno_df[order(new_pheno_df$sampleID), ]
    return(new_pheno_df)
}

RunXOAnalysis <- function(dat, pcent, samp_name, userdef_err_prob, out_dir) {
    # Set the total number of chromosomes
    n_chrom <- length(dat$is_x_chr)
    # Plot percent missing
    PlotMissing(dat, samp_name, out_dir)
    # Look for sample duplicates
    cg <- compare_geno(dat, cores=0)
    # Prepare subdirectories
    out_dir_dups <- PrepOutDir(out_dir, "other_summaries")
    capture.output(summary(cg), file = paste0(out_dir_dups, "/", samp_name, "_duplicates_summary.txt"))
    
    # Calculate genotype probabilities first
    # Error probability based on SNP array error rate
    # Reported in: https://www.nature.com/articles/nmeth842 (Steemers et al. 2006 Nature Methods)
    #bpr <- calc_genoprob(dat, error_prob=0.002, map_function="c-f", cores=0)
    bpr <- calc_genoprob(dat, error_prob=as.numeric(userdef_err_prob), map_function="c-f", cores=0)
    # Identify most probable genotype at each position then count exchanges
    bm <- maxmarg(bpr, minprob=0.95, cores=0)
    
    # Crossover counts
    # Returns counts of crossovers on each chromosome (as columns) in each individual
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

    # Plot of # of crossovers
    PlotNumXO(dat, btotxo, "Individuals", samp_name, out_dir)
    
    # Locate crossovers
    blxo <- locate_xo(bm, map = dat$gmap, cores = 0)
    # Try using physical map
    blxo_phys <- locate_xo(bm, map = dat$pmap, cores = 0)
    
    # Make crossover plots
    # Genetic map plot
    out_dir_gmap <- PrepOutDir(out_dir, "genetic_map_plots")
    out_dir_gmap <- paste0(out_dir, "/genetic_map_plots", sep = '')
    GeneticMapPlotting(dat, blxo, samp_name, out_dir_gmap)
    # Physical map plot
    out_dir_pmap <- PrepOutDir(out_dir, "physical_map_plots")
    lxodf <- PhysicalMapPlotting(dat, blxo_phys, pcent, samp_name, out_dir_pmap)
    
    # Create and save phenotype table
    pheno_out_list <- MakePhenoTable(dat, num_chr = n_chrom, lxodf, pcent, samp_name, out_dir)
    pheno_df <- pheno_out_list[[1]]
    new_lxodf <- pheno_out_list[[2]]
    
    # Create %noin%
    `%notin%` <- Negate(`%in%`)
    
    # Remove xo positions set to missing from blxo_phys for plotting purposes
    new_blxo_phys <- blxo_phys
    for (chr in 1:n_chrom) {
        for (indv in unique(new_lxodf$ind)) {
            curr_indv <- new_lxodf[new_lxodf$ind == indv, ]
            if (chr %in% curr_indv$chr) {
                curr_new_lxodf <- new_lxodf[new_lxodf$chr == chr & new_lxodf$ind == indv, ]
                if (nrow(curr_new_lxodf) != 0) {
                    tmp_xo_pos <- c()
                    for (row_num in 1:nrow(curr_new_lxodf)) {
                        curr_xo_pos <- curr_new_lxodf[row_num, ]$xo_pos
                        if (curr_xo_pos %in% new_blxo_phys[[chr]][[indv]]) {
                            # Add curr_xo_pos to vector
                            tmp_xo_pos <- c(tmp_xo_pos, curr_xo_pos)
                        }
                    }
                    # Reassign vector of xo positions to new_blxo_phys
                    new_blxo_phys[[chr]][[indv]] <- tmp_xo_pos
                }
            } else if (chr %notin% curr_indv$chr) {
                # Remove xo position set to missing from blxo_phys
                tmp_xo_pos <- c()
                # Reassign vector of xo positions to new_blxo_phys
                new_blxo_phys[[chr]][[indv]] <- tmp_xo_pos
            }
            
        }
    }
    
    # Physical map plot
    # After setting ambiguous crossovers to missing
    out_dir_miss <- PrepOutDir(out_dir, "physical_map_plots_miss")
    lxodf_miss <- PhysicalMapPlotting(dat, new_blxo_phys, pcent, samp_name, out_dir_miss)
    
    ############### Add feature ############### 
    # Combine total xo counts and xo counts by chr with pheno table
    # Since we want these counts after setting ambiguous crossovers to missing,
    #   we'll sum different cuts of the pheno_df data frame
    new_pheno_df <- UpdateCounts(pheno_df, n_chrom)
    ###########################################
    
    # Prepare subdirectories
    out_dir_pheno <- PrepOutDir(out_dir, "phenotype_tables")
    # Save phenotype table
    write.table(x = new_pheno_df,
                file = paste(out_dir_pheno, "/", samp_name, "_pheno.txt", sep = ""),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
}

Main <- function() {
    # Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    # User provided command line arguments
    yaml_fp <- args[1]
    pcent_fp <- args[2]
    userdef_err_prob <- args[3]
    out_dir <- args[4]
    
    # Read in files
    dat <- ReadFile(yaml_fp)
    pcent <- ReadPericentromerePos(pcent_fp)
    
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
        # Prep subdirectories
        out_dir_too_few <- PrepOutDir(out_dir, "too_few_markers")
        # Send summary output to a file
        capture.output(summary(dat_omit_null), file = paste0(out_dir_too_few, "/", "too_few_markers_per_chr-", samp_name, ".txt"))
    } else {
        # Proceed with analysis
        RunXOAnalysis(dat_omit_null, pcent, samp_name, userdef_err_prob, out_dir)
    }
}

# Run the program
Main()
