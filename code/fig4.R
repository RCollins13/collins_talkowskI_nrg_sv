#!/usr/bin/env Rscript
##################################################
## Project/Purpose: NGR Fig 4
## Date: 8/16/2023
## Author: Stephanie Hao
## sphao@bu.edu
##################################################


source("~/Dropbox (Partners HealthCare)/MGH/helperScripts/utils.R")
library(tidyverse)
library(data.table)
library(ggplot2)
WD<-"~/Desktop/Talkowski_NGR/fig4/"
setwd(WD)

cnvs<-fread("./cnvs.txt")
cnvs$NGenes<-as.numeric(cnvs$NGenes)

# Create NGenesCategory column
cnvs <- cnvs %>%
  mutate(NGenesCategory = case_when(
    NGenes < 2 ~ "<2",
    NGenes >= 2 & NGenes < 10 ~ "2-10",
    NGenes >= 10 ~ ">10"
  ))

# Calculate total counts per CNVType, NGenesCategory, and Chrom
total_counts <- cnvs %>%
  group_by(CNVType, NGenesCategory, Chrom) %>%
  summarise(TotalCount = n())

# Calculate total counts per CNVType and NGenesCategory
total_counts_all <- cnvs %>%
  group_by(CNVType, NGenesCategory) %>%
  summarise(TotalCountAll = n())

# Calculate percentage per CNVType, NGenesCategory, and Chrom
summary_data <- cnvs %>%
  group_by(CNVType, NGenesCategory, Chrom) %>%
  summarise(Count = n()) %>%
  left_join(total_counts, by = c("CNVType", "NGenesCategory", "Chrom")) %>%
  left_join(total_counts_all, by = c("CNVType", "NGenesCategory")) %>%
  mutate(Percent = (Count / TotalCount) * 100,
         PercentAll = (Count / TotalCountAll) * 100)

print(summary_data)


# Define a color palette
cnv_color<- c("DEL" = "#dd807a", "DUP"= "#366b90")

sum_percent_all <- summary_data %>%
  group_by(CNVType, NGenesCategory) %>%
  summarise(SumPercentAll = sum(PercentAll))



#Plot Percentage

# Create the combined bar plot with stacked DEL and dodged DUP bars
single_plot <- ggplot(data=subset(summary_data, NGenesCategory=="<2"), aes(x = Chrom, y = PercentAll, fill = CNVType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cnv_color) +
  labs(x = "Chromosome", y = "Percentage (%)", title = "Single Gene CNVs") +
  theme_classic() +
  scale_x_continuous(breaks = unique(cnvs$Chrom), labels = unique(cnvs$Chrom)) +
  scale_y_continuous(limits = c(0,35))+
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "CNV"))

multi_plot <- ggplot(data=subset(summary_data, NGenesCategory=="2-10"), aes(x = Chrom, y = PercentAll, fill = CNVType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cnv_color) +
  labs(x = "Chromosome", y = "Percentage (%)", title = "2-10 Genes CNVs") +
  theme_classic() +
  scale_x_continuous(breaks = unique(cnvs$Chrom), labels = unique(cnvs$Chrom)) +
  scale_y_continuous(limits = c(0,35)) +
  theme(
    legend.position = "right",
        ) +
  guides(fill = guide_legend(title = "CNV"))

max_plot <- ggplot(data=subset(summary_data, NGenesCategory==">10"), aes(x = Chrom, y = PercentAll, fill = CNVType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cnv_color) +
  labs(x = "Chromosome", y = "Percentage (%)", title = ">10 Genes CNVs") +
  theme_classic() +
  #geom_text(aes(label=signif(PercentAll, digits=2 )), hjust = 0.5, vjust = -0.5, position = "stack", size=3.5)+
  scale_x_continuous(breaks = unique(cnvs$Chrom), labels = unique(cnvs$Chrom)) +
  scale_y_continuous(limits = c(0,35))+
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "CNV"))

# # Create the plot
# plot <- ggplot(summary_data, aes(x = Chrom, y = Percent, fill = NGenesCategory)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(x = "Chromosome", y = "Percentage (%)", title = "Percentage of NGenesCategory by Chromosome") +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# 
# # Show the plot
# print(plot)

single_plot
multi_plot
max_plot

pdf("./single_plot.pdf", height=6, width=8)
single_plot
dev.off()

pdf("./multi_plot.pdf", height=6, width=8)
multi_plot
dev.off()

pdf("./max_plot.pdf", height=6, width=8)
max_plot
dev.off()

#########
# 08/23/2023
# 
# Left panel: histogram of CNV size (in Mb) for all genomic disorders. May need to log10-scale this if it's too skewed.
# Middle panel: histogram of # of genes per CNV for all genomic disorders.
# Right panel: histogram of # of OMIM genes per CNV for all genomic disorders.

# Left Panel: histogram of CNV size (in Mb) for all genomic disorders. May need to log10-scale this if it's too skewed.

cnvs<-fread("./cnvs_v2.txt")
cnvs$Ngenes<-as.numeric(cnvs$Ngenes)
cnvs$Size<-as.numeric(cnvs$Size)
cnvs$OMIM_genes<-as.numeric(cnvs$OMIM_genes)
cnvs$lsize<-log10(cnvs$Size)

# Create the bar plot
plot <- ggplot(data = cnvs, aes(x = lsize)) +
  geom_bar(stat = "count") +
  labs(x = "lsize", y = "Count", title = "Bar Plot of lsize vs. Count")

ggplot(data=cnvs, aes(x = lsize)) +
  geom_bar(stat = "bin",  fill="#366b90") + #bins = 40,
  labs(x = "Size (Log10)", y = "Count", title = "CNVs of all GDs") +
  #scale_x_continuous(labels = scales::comma_format(scale = 1e-6)) +
  theme_classic()

ggplot(data=cnvs, aes(x = Size)) +
  geom_bar(stat = "bin", fill="#366b90") + #bins = 40, 
  labs(x = "Size (Mb)", y = "Count", title = "CNVs of all GDs") +
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6)) +
  theme_classic()

size_plot <- ggplot(data=cnvs, aes(x = lsize)) +
  geom_bar(stat = "bin", bins = 30, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  labs(x = "Size (Log10)", y = "Count", title = "CNVs of all GDs") +
  #scale_x_continuous(limits=c(5.3,7)) +
  #scale_y_continuous(limits=c(0,15))+
  theme_classic()


# size_plot <- ggplot(data=cnvs, aes(x = lsize)) +
#   geom_bar(stat = "bin", bins = 30, fill="#9ccdea", color = "#4f6877",size = 0.25) +
#   labs(x = "Size (Log10)", y = "Count", title = "CNVs of all GDs") +
#   #scale_x_continuous(limits=c(5.3,7)) +
#   #scale_y_continuous(limits=c(0,15))+
#   theme_classic()
size_plot 

pdf("./size_plot.pdf", height=4, width=6)
size_plot
dev.off()

library(scales)

median_size <- median(cnvs$Size)  # Calculate median size
mean_size <- mean(cnvs$Size)

nl_size_plot<- ggplot(data=cnvs, aes(x = Size)) +
  geom_bar(stat = "bin", bins = 30, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  # geom_vline(xintercept = median_size, color = "#736f68", linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = mean_size, color = "#736f68", linetype = "dashed", size = 0.5) +
  annotate("text", x = mean_size + 850000, y = 25, label = paste("Mean:", comma(mean_size, scale = 1e-6), "Mb"), color = "#736f68") +
  labs(x = "Size (Mb)", y = "Count", title = "CNVs of all GDs") +
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6)) +
  theme_classic()

nl_size_plot

pdf("./size_nonlog_plot.pdf", height=4, width=6)
nl_size_plot
dev.off()

ave<-mean(cnvs$lsize)

size_plot <- ggplot(data=cnvs, aes(x = lsize)) +
  geom_bar(stat = "bin", bins = 30, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  geom_vline(xintercept = ave, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Size (Log10)", y = "Count", title = "CNVs of all GDs") +
  annotate("text", x = ave - 0.15, y = 15, label = paste("Mean:", comma(ave)), color = "#736f68") +
  #scale_x_continuous(limits=c(5.3,7)) +
  #scale_y_continuous(limits=c(0,15))+
  theme_classic()
size_plot 

pdf("./size_plot_mean.pdf", height=4, width=6)
size_plot
dev.off()

# Calculate the mean of the log-transformed sizes
mean_lsize <- mean(log10(cnvs$Size))

# Create the size plot
size_plot <- ggplot(data = cnvs, aes(x = lsize)) +
  geom_bar(stat = "bin", bins = 30, fill = "#9ccdea") + #, color = "#4f6877", size = 0.25
  geom_vline(xintercept = mean_lsize, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Size (Log10)", y = "Count", title = "CNVs of all GDs") +
  annotate("text", x = mean_lsize - 0.15, y = 15, label = "Mean", color = "#736f68") +
  # scale_x_continuous(limits = c(5.2, 7)) +
  # scale_y_continuous(limits = c(0, 20)) +
  theme_classic()
size_plot

pdf("./size_plot_mean.pdf", height=4, width=6)
size_plot
dev.off()


# Middle panel: histogram of # of genes per CNV for all genomic disorders.
ngene_plot <- ggplot(data=cnvs, aes(x = Ngenes)) +
  geom_bar(stat = "bin", bins = 30, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  # geom_vline(xintercept = ave, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of Genes", y = "Count", title = "CNVs of all GDs") +
  # annotate("text", x = ave + 0.15, y = 15, label = paste("Mean:", comma(ave)), color = "#736f68") +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  # scale_y_continuous(limits=c(0,15))+
  theme_classic()
ngene_plot 

pdf("./ngene_plot.pdf", height=4, width=6)
ngene_plot
dev.off()

# Calculate the mean of number of genes in CNVs across GDs
mean_ngene <- mean(cnvs$Ngenes)

ngene_plot <- ggplot(data=cnvs, aes(x = Ngenes)) +
  geom_bar(stat = "bin", bins = 30, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  geom_vline(xintercept = mean_ngene, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of Genes", y = "Count", title = "CNVs of all GDs") +
  annotate("text", x = ave + 14, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  # scale_y_continuous(limits=c(0,15))+
  theme_classic()
ngene_plot 

pdf("./ngene_plot_mean.pdf", height=4, width=6)
ngene_plot
dev.off()

# Right panel: histogram of # of OMIM genes per CNV for all genomic disorders.

omim_plot <- ggplot(data=cnvs, aes(x = OMIM_genes)) +
  geom_bar(stat = "bin", bins = 25, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  #geom_vline(xintercept = mean_ngene, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of OMIM Genes", y = "Frequency", title = "OMIM Genes of all GDs") +
  #annotate("text", x = ave + 13, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits=c(0,80))+
  theme_classic()
omim_plot 

pdf("./omim_plot_25.pdf", height=4, width=6)
omim_plot
dev.off()

omim_ave<-mean(cnvs$OMIM_genes)
omim_plot <- ggplot(data=cnvs, aes(x = OMIM_genes)) +
  geom_bar(stat = "bin", bins = 25, fill="#9ccdea") + #, color = "#4f6877",size = 0.25
  geom_vline(xintercept = omim_ave, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of OMIM Genes", y = "Frequency", title = "OMIM Genes of all GDs") +
  #annotate("text", x = ave + 13, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits=c(0,80))+
  theme_classic()
omim_plot 

######
# 08/24/2023
#Stacked CNVType

cnv_color <- c("DEL" = "#dd807a", "DUP" = "#366b90")
omim_plot <- ggplot(data=cnvs, aes(x = OMIM_genes, fill = CNVType)) +
  geom_bar(position = "stack", bins = 25) + #, fill="#9ccdea", color = "#4f6877",size = 0.25
  geom_vline(xintercept = omim_ave, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of OMIM Genes in the CNV", y = "Genomice Disorder CNVs", title = "OMIM Genes of all GDs") +
  #annotate("text", x = ave + 13, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_fill_manual(values = cnv_color) + 
  scale_x_continuous(breaks = seq(0, 31, by = 5)) +
  scale_y_continuous(limits=c(0,80))+
  theme_classic()
omim_plot 

pdf("./omim_plot_stacked.pdf", height=4, width=6)
omim_plot
dev.off()




nl_size_plot<- ggplot(data=cnvs, aes(x = Size, fill=CNVType)) +
  geom_bar(stat = "bin", bins = 30) + #, color = "#4f6877",size = 0.25
  # geom_vline(xintercept = median_size, color = "#736f68", linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = mean_size, color = "#736f68", linetype = "dashed", size = 0.5) +
  #annotate("text", x = mean_size + 850000, y = 25, label = paste("Mean:", comma(mean_size, scale = 1e-6), "Mb"), color = "#736f68") +  
  scale_fill_manual(values = cnv_color) +
  labs(x = "Size (Mb)", y = "Genomice Disorder CNVs", title = "CNVs Sizes of all GDs") +
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6), breaks = seq(0, 10e6, by = 2e6), limits=c(0,10e6)) +
  scale_y_continuous(breaks = seq(0, 65, by = 10)) +
  theme_classic()

nl_size_plot

pdf("./size_nonlog_stacked_v3.pdf", height=4, width=6)
nl_size_plot
dev.off()

# Calculate the mean of number of genes in CNVs across GDs
mean_ngene <- mean(cnvs$Ngenes)

ngene_plot <- ggplot(data=cnvs, aes(x = Ngenes, fill = CNVType)) +
  geom_bar(stat = "bin", bins = 30) + #, color = "#4f6877",size = 0.25
  geom_vline(xintercept = mean_ngene, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of Genes in the CNV", y = "Genomice Disorder CNVs", title = "Number of Genes in CNVs") +
  #annotate("text", x = ave + 14, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_fill_manual(values = cnv_color) +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  scale_y_continuous(limits=c(0,32))+
  theme_classic()
ngene_plot 

pdf("./ngene_plot_stacked.pdf", height=4, width=6)
ngene_plot
dev.off()


# 08/26/2023
# Filter cnvs segment IDs to just NAHR CNVS

# continue using cnvs_v2.txt

cnvs<-fread("./cnvs_v2.txt")
cnvs$Ngenes<-as.numeric(cnvs$Ngenes)
cnvs$Size<-as.numeric(cnvs$Size)
cnvs$OMIM_genes<-as.numeric(cnvs$OMIM_genes)
cnvs$lsize<-log10(cnvs$Size)

nahr_segment_ids <- c(
  "merged_DEL_segment_1p36.32-p36.33",
  "merged_DEL_segment_1q21.1",
  "merged_DUP_segment_1q21.1-q21.2",
  "merged_DEL_segment_1q21.1-q21.2",
  "merged_DEL_segment_1q43-q44",
  "LC_GD_4_DUP_2q11.1-q11.2",
  "merged_DUP_segment_2q12.2-q12.3",
  "merged_DEL_segment_2q13_B",
  "merged_DUP_segment_2q13",
  "merged_DEL_segment_2q13_A",
  "MC_GD_1_DUP_2q13",
  "merged_DUP_segment_2q37.1",
  "merged_DUP_segment_3q29",
  "merged_DEL_segment_3q29",
  "merged_DUP_segment_5p15.33",
  "HC_GD_6_DEL_5q35.2-q35.3",
  "merged_DUP_segment_7p22.1",
  "merged_DEL_segment_7q11.23",
  "merged_DUP_segment_7q11.23",
  "merged_DEL_segment_8p23.1-p23.3",
  "merged_DUP_segment_8p23.1",
  "merged_DEL_segment_8p23.1",
  "MC_GD_8_DEL_10q11.22-q11.23",
  "LC_GD_10_DUP_10q22.3-q23.2",
  "merged_DUP_segment_10q26.3",
  "merged_DEL_segment_15q11.2",
  "merged_DUP_segment_15q11.2-q13.3",
  "merged_DEL_segment_15q11.2-q13.1",
  "merged_DEL_segment_15q13.1-q13.2",
  "merged_DEL_segment_15q13.2-q13.3",
  "HC_GD_11_DEL_15q24.1-q24.2",
  "LC_GD_15_DEL_15q24.2-q24.3",
  "merged_DUP_segment_16p13.11",
  "merged_DEL_segment_16p13.11",
  "merged_DUP_segment_16p12.2",
  "merged_DEL_segment_16p12.2",
  "merged_DUP_segment_16p11.2_B",
  "merged_DEL_segment_16p11.2_B",
  "merged_DEL_segment_16p11.2_A",
  "merged_DUP_segment_16p11.2_A",
  "merged_DEL_segment_17p12",
  "merged_DUP_segment_17p12",
  "merged_DEL_segment_17p11.2",
  "merged_DUP_segment_17p11.2",
  "merged_DEL_segment_17q11.2",
  "MC_GD_8_DUP_17q11.2",
  "merged_DUP_segment_17q12",
  "merged_DEL_segment_17q12",
  "merged_DEL_segment_17q21.31",
  "merged_DUP_segment_19p12_B",
  "merged_DUP_segment_22q11.21",
  "merged_DEL_segment_22q11.21",
  "HC_GD_22_DEL_22q11.21-q11.23",
  "merged_DUP_segment_22q11.22-q11.23"
)

nahr_cnvs<-cnvs[cnvs$SegmentID %in% nahr_segment_ids, ]

cnv_color <- c("DEL" = "#dd807a", "DUP" = "#366b90")

omim_ave<-mean(nahr_cnvs$OMIM_genes)

omim_plot <- ggplot(data=nahr_cnvs, aes(x = OMIM_genes, fill = CNVType)) +
  geom_bar(position = "stack", bins = 25) + #, fill="#9ccdea", color = "#4f6877",size = 0.25
  geom_vline(xintercept = omim_ave, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of OMIM Genes in the CNV", y = "Genomice Disorder CNVs", title = "OMIM Genes of all GDs") +
  #annotate("text", x = ave + 13, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_fill_manual(values = cnv_color) + 
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits=c(0,16))+
  theme_classic()
omim_plot 

pdf("./nahr_omim_plot_stacked.pdf", height=4, width=6)
omim_plot
dev.off()

mean_size <- mean(nahr_cnvs$Size)

nl_size_plot <- ggplot(data = nahr_cnvs, aes(x = Size, fill = CNVType)) +
  geom_bar(stat = "bin", bins = 30) +
  geom_vline(xintercept = mean_size, color = "#736f68", linetype = "dashed", size = 0.5) +
  scale_fill_manual(values = cnv_color) +
  labs(x = "Size (Mb)", y = "Genomic Disorder CNVs", title = "CNV Sizes of all GDs") +
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6), breaks = seq(0, 10e6, by = 2e6), limits=c(0,10e6)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  theme_classic()
nl_size_plot

pdf("./nahr_size_nonlog_stacked.pdf", height=4, width=6)
nl_size_plot
dev.off()

# Calculate the mean of number of genes in CNVs across GDs
mean_ngene <- mean(nahr_cnvs$Ngenes)

ngene_plot <- ggplot(data=cnvs, aes(x = Ngenes, fill = CNVType)) +
  geom_bar(stat = "bin", bins = 30) + #, color = "#4f6877",size = 0.25
  geom_vline(xintercept = mean_ngene, color = "#736f68", linetype = "dashed", size = 0.5) +
  labs(x = "Number of Genes in the CNV", y = "Genomice Disorder CNVs", title = "Number of Genes in CNVs") +
  #annotate("text", x = ave + 14, y = 30, label = paste("Mean:", comma(mean_ngene)), color = "#736f68") +
  scale_fill_manual(values = cnv_color) +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  scale_y_continuous(limits=c(0,32))+
  theme_classic()
ngene_plot 

pdf("./nahr_ngene_plot_stacked.pdf", height=4, width=6)
ngene_plot
dev.off()


##########
# Sept 7
# Going babck to all CNVs
# cnv_v2 file from line 129


cnv_color <- c("DEL" = "#dd807a", "DUP" = "#366b90")

nl_size_plot<- ggplot(data=cnvs, aes(x = Size, fill=CNVType)) +
  geom_bar(stat = "bin", bins = 30) + #, color = "#4f6877",size = 0.25
  # geom_vline(xintercept = median_size, color = "#736f68", linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = mean_size, color = "#736f68", linetype = "dashed", size = 0.5) +
  #annotate("text", x = mean_size + 850000, y = 25, label = paste("Mean:", comma(mean_size, scale = 1e-6), "Mb"), color = "#736f68") +  
  scale_fill_manual(values = cnv_color) +
  labs(x = "Size (Mb)", y = "Genomice Disorder CNVs", title = "CNVs Sizes of all GDs") +
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6),breaks = seq(0, 12e6, by = 2e6)) +
  #scale_x_continuous(labels = scales::comma_format(scale = 1e-6), breaks = seq(0, 11e6, by = 2e6), limits=c(0,10e6)) +
  scale_y_continuous(breaks = seq(0, 80, by = 10)) +
  theme_classic()

nl_size_plot

pdf("./size_nonlog_stacked_v3.pdf", height=4, width=6)
nl_size_plot
dev.off()


# Scatter plot of mean genes in CNVs vs OMIM genes in CNV
scatter_plot <- ggplot(data = cnvs, aes(x = Ngenes, y = OMIM_genes, color = CNVType)) +
  geom_point() +
  labs(x = "Number of Genes in the CNV", y = "Number of OMIM Genes in the CNV", title = "Overall Genes in GD CNVs vs OMIM Genes") +
  scale_color_manual(values = cnv_color) +
theme_classic()

# Display the scatter plot
scatter_plot

scatter_plot <- ggplot(data = cnvs, aes(x = Ngenes, y = OMIM_genes, color = CNVType)) +
  geom_point() +
  labs(x = "Number of Genes in the CNV", y = "Number of OMIM Genes in the CNV", title = "Overall Genes in GD CNVs vs OMIM Genes") +
  scale_color_manual(values = cnv_color) +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey", linetype = "dotted", size = 0.5)

scatter_plot

# Display the scatter plot with separate smooth lines for each CNV type
scatter_plot_with_smooth <- scatter_plot +
  geom_smooth(aes(group = CNVType, fill = CNVType), method = "lm", se = TRUE, linetype = "solid", size = 0.5, alpha=0.2)

# Show the scatter plot with separate smooth lines
scatter_plot_with_smooth

pdf("./gene_vs_omim_with_smooth_shade_stacked_v3.pdf", height=4, width=6)
scatter_plot_with_smooth
dev.off()


# Violin Plots

# Create a list of CNV types
cnv_types <- c("DEL", "DUP")

# Create a list to store the CNV data
cnv_data <- list()

# Loop through CNV types
for (cnv_type in cnv_types) {
  # Filter data for CNV type
  cnv_data[[cnv_type]] <- cnvs[cnvs$CNVType == cnv_type, ]
}

# Create a function to generate violin plots
create_violin_plot <- function(data, title) {
  ggplot(data, aes(x = CNVType, y = Ngenes, fill = CNVType)) +
    geom_violin(trim = FALSE, alpha=0.25) +
    geom_jitter(aes(color = CNVType), width = 0.2, height = 0, shape = 15, size = 1.6, alpha = 1.2, stroke=0.7) +  # Add jitter layer
    labs(x = "CNV Type", y = "Number of Genes", title = title) +
    scale_y_continuous(limits=c(0,92),breaks = seq(0, 92, by = 10))+
    scale_fill_manual(values = cnv_color) +
    scale_color_manual(values = cnv_color) +  # Set jitter color to match DEL and DUP colors
    theme_classic() +
    theme(legend.position = "none")
}

# Create the violin plots
cnv_violin_plot <- create_violin_plot(do.call(rbind, cnv_data), "All GD CNVs")

# Filter the cnvs data for NAHR segment IDs
nahr_cnvs <- cnvs[cnvs$SegmentID %in% nahr_segment_ids, ]

# Create a list to store the NAHR CNV data
nahr_cnv_data <- list()

# Loop through CNV types
for (cnv_type in cnv_types) {
  # Filter NAHR CNV data for CNV type
  nahr_cnv_data[[cnv_type]] <- nahr_cnvs[nahr_cnvs$CNVType == cnv_type, ]
}

# Create the violin plots for NAHR CNVs
nahr_cnv_violin_plot <- create_violin_plot(do.call(rbind, nahr_cnv_data), "NAHR CNVs")

# Filter the cnvs data for NAHR segment IDs
not_nahr_cnvs <- cnvs[cnvs$SegmentID %!in% nahr_segment_ids, ]

# Create a list to store the NAHR CNV data
not_nahr_cnv_data <- list()

# Loop through CNV types
for (cnv_type in cnv_types) {
  # Filter NAHR CNV data for CNV type
  not_nahr_cnv_data[[cnv_type]] <- not_nahr_cnvs[not_nahr_cnvs$CNVType == cnv_type, ]
}

# Create the violin plots for NAHR CNVs
not_nahr_cnv_violin_plot <- create_violin_plot(do.call(rbind, not_nahr_cnv_data), "non-NAHR CNVs")

# Arrange the plots side by side
# library(gridExtra)
grid.arrange(cnv_violin_plot, nahr_cnv_violin_plot, ncol = 2)

grid.arrange(not_nahr_cnv_violin_plot, nahr_cnv_violin_plot, ncol = 2)

pdf("./violin_genes_stacked_v3.pdf", height=4, width=6)
grid.arrange(cnv_violin_plot, nahr_cnv_violin_plot, ncol = 2)
dev.off()

pdf("./violin_genes_stacked_v4.pdf", height=4, width=6)
grid.arrange(not_nahr_cnv_violin_plot, nahr_cnv_violin_plot, ncol = 2)
dev.off()

# Perform the Wilcoxon rank-sum test for DEL CNVs
wilcox_DEL <- wilcox.test(not_nahr_cnv_data$DEL$Ngenes, nahr_cnv_data$DEL$Ngenes,
                          alternative = "two.sided", paired = FALSE)

# Perform the Wilcoxon rank-sum test for DUP CNVs
wilcox_DUP <- wilcox.test(not_nahr_cnv_data$DUP$Ngenes, nahr_cnv_data$DUP$Ngenes,
                          alternative = "two.sided", paired = FALSE)

# Extract p-values
p_value_DEL <- wilcox_DEL$p.value
p_value_DUP <- wilcox_DUP$p.value

# Print the p-values
cat("P-Value for DEL CNVs:", p_value_DEL, "\n")
cat("P-Value for DUP CNVs:", p_value_DUP, "\n")

# 
# P-Value for DEL CNVs: 0.004700578 
# 
# P-Value for DUP CNVs: 0.0007873683 
# > scientific(p_value_DUP, digits=3)
# [1] "7.87e-04"
# > scientific(p_value_DEL, digits=3)
# [1] "4.7e-03"
