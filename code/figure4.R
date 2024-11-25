#!/usr/bin/env Rscript

# Code to generate technology yield figure for NRG SV review
# Ryan Collins, 2024

# Setup
options(scipen=1000, stringsAsFactors=F)
setwd("~/Desktop/Collins/Talkowski/misc/NRG_SV_review_Fall2020/NRG_SV_review_figures/tech_yield_figure/")
colors.omim <- c("lrwgs" = "#A85347",
                "srwgs" = "#BA8C57",
                "ogm" = "#366B90",
                "cma" = "#66568F",
                "wes" = "#42694A")
colors.all <- c("lrwgs" = "#F1CDCB",
                 "srwgs" = "#F6EBC6",
                 "ogm" = "#C4E1F2",
                 "cma" = "#E3CFE4",
                 "wes" = "#D7ECDB")

# Read data
df <- read.table("NRG.tech_yield_figure_data.tsv", header=T, sep="\t")

# Plot
pdf("tech_yield_point_placement.pdf", height=2.55, width=5.75)
par(mar=rep(0, 4))
plot(log10(df$sv_per_genome), log2(df$genes_disrupted_per_genome),
     pch=19, col=colors.all[df$technology], ylim=c(0, 10), cex=0.8)
points(log10(df$sv_per_genome), log2(df$omim_disrupted_per_genome), pch=19, 
       col=colors.omim[df$technology], cex=0.8)
model <- lm(log2(c(df$genes_disrupted_per_genome, df$omim_disrupted_per_genome)) ~ rep(log10(df$sv_per_genome), 2))
abline(model)
summary(model)
# abline(lm(log2() ~ log10(df$sv_per_genome)))
# summary(lm(log2(df$omim_disrupted_per_genome) ~ log10(df$sv_per_genome)))
axis(1, tck=0.03, line=-3)
axis(2, tck=0.03, line=-3)
dev.off()