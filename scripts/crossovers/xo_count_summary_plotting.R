#!/usr/bin/env Rscript

library(qtl2)
library(data.table)
library(plyr)
library(ggplot2)
library(quantreg)
library(ggbeeswarm)
library(ggpubr) # ggplot correlation visualization

xo_data_dir <- "~/Dropbox/Projects/barley_recombination/Analyses/rqtl2/num_xo_data"
out_dir <- "~/Dropbox/Projects/barley_recombination/Analyses/rqtl2"
pedigree_fp <- "~/Dropbox/Projects/barley_recombination/Data/T3_Full_Pedigree.csv"
yaml_dir <- "~/Dropbox/Projects/barley_recombination/Data/split_by_family_fillIn/qtl2inputs"

#----------------------------------

# Read in full pedigree
pedigree <- read.csv(file = pedigree_fp, header = FALSE)
# Change column names
colnames(pedigree) <- c("progeny", "parent1", "parent2")

# Prepping filepaths
files <- list.files(path = xo_data_dir, pattern = ".txt")
fp <- paste0(xo_data_dir, "/", files)

# Read in SNP data for each family
yaml_files <- list.files(path = yaml_dir, pattern = ".yaml")
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

# Restructure data into long format for beeswarm plot
xo_long <- list()
for (i in 1:length(fp)) {
    tmp_family <- strsplit(xocount[[i]]$ind, split = "-")
    tmp_fam_name <- tmp_family[[1]][1]
    xo_long[[i]] <- data.frame(
        ind = xocount[[i]]$ind,
        total_nxo = xocount[[i]]$xo_count,
        family = rep(tmp_fam_name, times = length(xocount[[i]]$xo_count)),
        all = rep("together", times = length(xocount[[i]]$xo_count))
    )
}
# Concatenate list of dataframes
all_xo_long <- ldply(xo_long, rbind)
# Split by cycles
c1_xo_long <- all_xo_long[grep("^MS10", all_xo_long$family), ]
c2_xo_long <- all_xo_long[grep("^MS11", all_xo_long$family), ]
c3_xo_long <- all_xo_long[grep("^MS12", all_xo_long$family), ]

# Plotting
# Version 1 ggplot
# Cycle 1
c1plot <- ggplot(c1_xo_long, aes(family, total_nxo)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title = element_text(size = 18)) +
    xlab("Cycle 1 Families") +
    ylab("Number of Crossovers") +
    ylim(0, 48)
c1plot
ggsave("c1_nxo_summary_plot.pdf", plot = c1plot, device = "pdf", path = out_dir, width = 10)


# Cycle 2
c2plot <- ggplot(c2_xo_long, aes(family, total_nxo)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title = element_text(size = 18)) +
    xlab("Cycle 2 Families") +
    ylab("Number of Crossovers") +
    ylim(0, 48)
c2plot
ggsave("c2_nxo_summary_plot.pdf", plot = c2plot, device = "pdf", path = out_dir, width = 10)

# Cycle 3
c3plot <- ggplot(c3_xo_long, aes(family, total_nxo)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title = element_text(size = 18)) +
    xlab("Cycle 3 Families") +
    ylab("Number of Crossovers") +
    ylim(0, 48)
c3plot
ggsave("c3_nxo_summary_plot.pdf", plot = c3plot, device = "pdf", path = out_dir, width = 10)


########################################
# Identify parents that contribute to higher mean crossover rate across families
c1_families <- data.frame(family=unique(c1_xo_long$family))
c2_families <- data.frame(family=unique(c2_xo_long$family))
c3_families <- data.frame(family=unique(c3_xo_long$family))

# Cycle 1
c1_fam_mean <- data.frame(family=c1_families$family, xo_mean=NA)
for (i in 1:length(c1_families$family)) {
    c1_fam_mean[i, 2] <- mean(c1_xo_long[c1_xo_long$family %in% c1_families[i, ], 2])
}

c1_fam_mean_wParent <- merge(x = c1_fam_mean, y = pedigree, by.x = "family", by.y = "progeny", all.x = FALSE)
c1_fam_mean_wParent_sorted <- c1_fam_mean_wParent[order(c1_fam_mean_wParent$xo_mean, c1_fam_mean_wParent$parent1, decreasing = TRUE), ]
# Compare to mean of all of Cycle 1
summary(c1_xo_long$total_nxo)
c1_mean <- mean(c1_xo_long$total_nxo)
c1_fam_mean_wParent_sorted[c1_fam_mean_wParent_sorted$xo_mean >= c1_mean, ]
# Identify parents that contribute to higher xo rate
c1p1 <- data.frame(founder=unique(c1_fam_mean_wParent_sorted[c1_fam_mean_wParent_sorted$xo_mean >= c1_mean, 3]))
c1p2 <- data.frame(founder=unique(c1_fam_mean_wParent_sorted[c1_fam_mean_wParent_sorted$xo_mean >= c1_mean, 4]))
c1p_high_mean <- rbind(c1p1, c1p2)
sc1p_high_mean <- data.frame(founder=c1p_high_mean[order(c1p_high_mean$founder), ])
summary(sc1p_high_mean)

# Cycle 2
c2_fam_mean <- data.frame(family=c2_families$family, xo_mean=NA)
for (i in 1:length(c2_families$family)) {
    c2_fam_mean[i, 2] <- mean(c2_xo_long[c2_xo_long$family %in% c2_families[i, ], 2])
}

c2_fam_mean_wParent <- merge(x = c2_fam_mean, y = pedigree, by.x = "family", by.y = "progeny", all.x = FALSE)
c2_fam_mean_wParent_sorted <- c2_fam_mean_wParent[order(c2_fam_mean_wParent$xo_mean, c2_fam_mean_wParent$parent1, decreasing = TRUE), ]
# Compare to mean of all of Cycle 2
summary(c2_xo_long$total_nxo)
c2_mean <- mean(c2_xo_long$total_nxo)
c2_fam_mean_wParent_sorted[c2_fam_mean_wParent_sorted$xo_mean >= c2_mean, ]
# Identify parents that contribute to higher xo rate
c2p1 <- data.frame(founder=unique(c2_fam_mean_wParent_sorted[c2_fam_mean_wParent_sorted$xo_mean >= c2_mean, 3]))
c2p2 <- data.frame(founder=unique(c2_fam_mean_wParent_sorted[c2_fam_mean_wParent_sorted$xo_mean >= c2_mean, 4]))
c2p_high_mean <- rbind(c2p1, c2p2)
sc2p_high_mean <- data.frame(founder=c2p_high_mean[order(c2p_high_mean$founder), ])
summary(sc2p_high_mean)


# Cycle 3
c3_fam_mean <- data.frame(family=c3_families$family, xo_mean=NA)
for (i in 1:length(c3_families$family)) {
    c3_fam_mean[i, 2] <- mean(c3_xo_long[c3_xo_long$family %in% c3_families[i, ], 2])
}
c3_fam_mean_wParent <- merge(x = c3_fam_mean, y = pedigree, by.x = "family", by.y = "progeny", all.x = FALSE)
c3_fam_mean_wParent_sorted <- c3_fam_mean_wParent[order(c3_fam_mean_wParent$xo_mean, c3_fam_mean_wParent$parent1, decreasing = TRUE), ]
# Compare mean to all of Cycle 3
summary(c3_xo_long$total_nxo)
c3_mean <- mean(c3_xo_long$total_nxo)
c3_fam_mean_wParent_sorted[c3_fam_mean_wParent_sorted$xo_mean >= c3_mean, ]
# Figure out pedigree of lowest crossover number individuals
#----------
# MS12_2124
# MS12_2124  2.541667 MS11S2026-019 MS11S2053-023
# Follow parent 1
pedigree[pedigree$progeny == "MS11S2026", ]
# MS11S2026 MS10S3018-001 MS10S3021-013
pedigree[pedigree$progeny == "MS10S3018", ]
# MS10S3018 FEG175-57 ND25986
pedigree[pedigree$progeny == "MS10S3021", ]
# MS10S3021 ND24906 FEG183-52

# Follow parent 2
pedigree[pedigree$progeny == "MS11S2053", ]
# MS11S2053 MS10S3032-016 MS10S3004-004
pedigree[pedigree$progeny == "MS10S3032", ]
# MS10S3032 6B03-4304 FEG175-57
pedigree[pedigree$progeny == "MS10S3004", ]
# MS10S3004 FEG153-25 FEG175-57

#----------
# MS12_2123
# MS12_2123  2.875000 MS11S2026-006 MS11S2075-019
# Follow parent 1
pedigree[pedigree$progeny == "MS11S2026", ]
# MS11S2026 MS10S3018-001 MS10S3021-013
pedigree[pedigree$progeny == "MS10S3018", ]
# MS10S3018 FEG175-57 ND25986
pedigree[pedigree$progeny == "MS10S3021", ]
# MS10S3021 ND24906 FEG183-52

# Follow parent 2
pedigree[pedigree$progeny == "MS11S2075", ]
    # MS11S2075 MS10S3050-024 MS10S3049-003
pedigree[pedigree$progeny == "MS10S3050", ]
# MS10S3050 ND26036 6B05-0922
pedigree[pedigree$progeny == "MS10S3049", ]
# MS10S3049 ND25728 6B01-2218

#----------
# MS12_2125
# MS12_2125  4.583333 MS11S2029-011 MS11S2052-012
# Follow parent 1
pedigree[pedigree$progeny == "MS11S2029", ]
# MS11S2029 MS10S3018-006 MS10S3055-014
pedigree[pedigree$progeny == "MS10S3018", ]
# MS10S3018 FEG175-57 ND25986
pedigree[pedigree$progeny == "MS10S3055", ]
# MS10S3055 6B04-0290 ND26036

# Follow parent 2
pedigree[pedigree$progeny == "MS11S2052", ]
# MS11S2052 MS10S3032-016 MS10S3019-015
pedigree[pedigree$progeny == "MS10S3032", ]
# MS10S3032 6B03-4304 FEG175-57
pedigree[pedigree$progeny == "MS10S3019", ]
# MS10S3019 FEG175-57 ND26104

#----------
# Look at one highest crossover rate in cycle 3
# MS12_2164
# MS12_2164 33.33333 MS11S2076-021 MS11S2058-019
# Follow parent 1
pedigree[pedigree$progeny == "MS11S2076", ]
# MS11S2076 MS10S3053-015 MS10S3019-015
pedigree[pedigree$progeny == "MS10S3053", ]
# MS10S3053 6B03-4478 ND25986
pedigree[pedigree$progeny == "MS10S3019", ]
# MS10S3019 FEG175-57 ND26104

# Follow parent 2
pedigree[pedigree$progeny == "MS11S2058", ]
# MS11S2058 MS10S3034-018 MS10S3029-013
pedigree[pedigree$progeny == "MS10S3034", ]
# MS10S3034 6B03-4478 FEG175-57
pedigree[pedigree$progeny == "MS10S3029", ]
# MS10S3029 FEG183-52 6B01-2218

#---------------------------------

# Check if number of crossovers is correlated with number of informative markers
# Concatenate list of dataframes
tmarkers <- ldply(total_markers, rbind)
shapiro.test(tmarkers$tot_markers)
hist(tmarkers$tot_markers)
tmarkersc1 <- tmarkers[grep("^MS10", tmarkers$family), ]
tmarkersc2 <- tmarkers[grep("^MS11", tmarkers$family), ]
tmarkersc3 <- tmarkers[grep("^MS12", tmarkers$family), ]
shapiro.test(tmarkersc1$tot_markers)
shapiro.test(tmarkersc2$tot_markers)
shapiro.test(tmarkersc3$tot_markers)

# Cycle 1
shapiro.test(c1_fam_mean$xo_mean)
hist(c1_fam_mean$xo_mean)
c1_mark <- merge(x = c1_fam_mean, y = tmarkers, by = "family", all.x = FALSE)
# Cycle 1 Family Means
ggqqplot(c1_mark$xo_mean, ylab = "Crossover Mean")
    
# Cycle 1 Family Markers
ggqqplot(c1_mark$tot_markers, ylab = "Total Number of Markers")
# Cycle 1
c1plot <- ggscatter(c1_mark, x = "xo_mean", y = "tot_markers", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean # of Crossovers", ylab = "Total # of Markers/Family") +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 18))
c1plot
ggsave("cycle1_corr_plot-pearson.pdf", plot = c1plot, path = out_dir)


# Cycle 2
shapiro.test(c2_fam_mean$xo_mean)
c2_mark <- merge(x = c2_fam_mean, y = tmarkers, by = "family", all.x = FALSE)
# Cycle 2 Family Means
ggqqplot(c2_mark$xo_mean, ylab = "Crossover Mean")
# Cycle 2 Family Markers
ggqqplot(c2_mark$tot_markers, ylab = "Total Number of Markers")
# Cycle 2
c2plot <- ggscatter(c2_mark, x = "xo_mean", y = "tot_markers", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mean # of Crossovers", ylab = "Total # of Markers/Family") +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 18))
c2plot
ggsave("cycle2_corr_plot-spearman.pdf", plot = c2plot, path = out_dir)

# Cycle 3
shapiro.test(c3_fam_mean$xo_mean)
c3_mark <- merge(x = c3_fam_mean, y = tmarkers, by = "family", all.x = FALSE)
# Cycle 3 Family Means
ggqqplot(c3_mark$xo_mean, ylab = "Crossover Mean")
# Cycle 3 Family Markers
ggqqplot(c3_mark$tot_markers, ylab = "Total Number of Markers")
# Cycle 3
c3plot <- ggscatter(c3_mark, x = "xo_mean", y = "tot_markers", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mean # of Crossovers", ylab = "Total # of Markers/Family") +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 18))
c3plot
ggsave("cycle3_corr_plot-spearman.pdf", plot = c3plot, path = out_dir)
