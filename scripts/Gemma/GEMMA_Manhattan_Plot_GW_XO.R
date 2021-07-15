# Make a Manhattan plot of -log(P_lrt) from the GEMMA results. These will be
# a little weird because the population is highly structured and the founders
# are pretty closely related.

# Take command-line arguments
args <- commandArgs(trailingOnly = TRUE)
gwxo_filename <- args[1]

# Read the GEMMA association results
gw.xo <- read.table(gwxo_filename, header = TRUE)

# Digest up the plotting a little bit. We want to plot the seven chromosomes
# side-by-side and give them distinct colors
# This is "Dark2" with seven levels from RColorBrewer
color_vector <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")

# We will re-define the data to plot so that the chromosomes will be placed
# end-to-end. Divide by 1M to put it in megabase scale
chr1H_end <- as.numeric(tail(gw.xo$ps[gw.xo$chr == "chr1H"], n=1))
chr2H_end <- as.numeric(tail(gw.xo$ps[gw.xo$chr == "chr2H"], n=1)) + chr1H_end
chr3H_end <- as.numeric(tail(gw.xo$ps[gw.xo$chr == "chr3H"], n=1)) + chr2H_end
chr4H_end <- as.numeric(tail(gw.xo$ps[gw.xo$chr == "chr4H"], n=1)) + chr3H_end
chr5H_end <- as.numeric(tail(gw.xo$ps[gw.xo$chr == "chr5H"], n=1)) + chr4H_end
chr6H_end <- as.numeric(tail(gw.xo$ps[gw.xo$chr == "chr6H"], n=1)) + chr5H_end

gw.xo_plt <- data.frame(
    X=c(
        gw.xo$ps[gw.xo$chr == "chr1H"],
        gw.xo$ps[gw.xo$chr == "chr2H"] + chr1H_end,
        gw.xo$ps[gw.xo$chr == "chr3H"] + chr2H_end,
        gw.xo$ps[gw.xo$chr == "chr4H"] + chr3H_end,
        gw.xo$ps[gw.xo$chr == "chr5H"] + chr4H_end,
        gw.xo$ps[gw.xo$chr == "chr6H"] + chr5H_end,
        gw.xo$ps[gw.xo$chr == "chr7H"] + chr6H_end) / 1000000,
    Y=-log10(p.adjust(gw.xo$p_lrt, method="BH")),
    Chrom=c(
        rep("Chr1H", sum(gw.xo$chr == "chr1H")),
        rep("Chr2H", sum(gw.xo$chr == "chr2H")),
        rep("Chr3H", sum(gw.xo$chr == "chr3H")),
        rep("Chr4H", sum(gw.xo$chr == "chr4H")),
        rep("Chr5H", sum(gw.xo$chr == "chr5H")),
        rep("Chr6H", sum(gw.xo$chr == "chr6H")),
        rep("Chr7H", sum(gw.xo$chr == "chr7H"))
        )
    )

# pdf(file="GEMMA_BOPA_Manhattan.pdf", height=4, width=8)
png(file="GEMMA_Manhattan_bopa_gw_xo_count.png", res=150, height=600, width=1200)
#par(mfrow=c(3, 1), mar=c(1, 3, 0.75, 0.1), mgp=c(2, 1, 0))
par(mfrow=c(1,1))
chroms <- c("Chr1H", "Chr2H", "Chr3H", "Chr4H", "Chr5H", "Chr6H", "Chr7H")
# Make a plot for genome-wide crossover count
plot(0, type="n", axes=FALSE, ylim=c(0, 4), xlim=c(0, 4600), xlab="", ylab="-log10(P)", main="Genome-wide Crossover Count")
sapply(
    seq_along(chroms),
    function(x) {
        cname <- chroms[x]
        color <- color_vector[x]
        points(
            gw.xo_plt$Y[gw.xo_plt$Chrom == cname] ~ gw.xo_plt$X[gw.xo_plt$Chrom == cname],
            pch=19,
            cex=0.25,
            col=color)
    })
abline(h=-log10(0.01), col="red", lwd=0.75, lty=3)
abline(h=-log10(0.1), col="blue", lwd=0.75, lty=3)
abline(
    v=c(chr1H_end, chr2H_end, chr3H_end, chr4H_end, chr5H_end, chr6H_end)/1000000,
    lwd=0.25,
    col="grey",
    lty=2)
axis(side=2)

mtext(
    sub("Chr", "", chroms),
    side=1,
    at=c(
        mean(c(0, chr1H_end)),
        mean(c(chr1H_end, chr2H_end)),
        mean(c(chr2H_end, chr3H_end)),
        mean(c(chr3H_end, chr4H_end)),
        mean(c(chr4H_end, chr5H_end)),
        mean(c(chr5H_end, chr6H_end)),
        mean(c(chr6H_end, 4600000000))
        )/1000000,
    cex=0.75
    )

dev.off()

# Part below is for exploratory purposes
gw.xo$corrected_p = -log10(p.adjust(gw.xo$p_lrt, method="BH"))
head(gw.xo[order(gw.xo$corrected_p, decreasing = T), ])
