#!/usr/bin/env Rscript
##################################################
## Project/Purpose: Talkowski NGR Figure 2
## Date: 07/26/2023
## Author: Stephanie Hao
## sphao@mgh.harvard.edu
##################################################

source("~/Dropbox (Partners HealthCare)/MGH/helperScripts/utils.R")
library(tidyverse)
library(data.table)
library(ggplot2)
WD<-"~/Desktop/Talkowski_NGR/"
setwd(WD)

data<-fread("./NRG_SV_review.singleton_counts.tsv")

# A
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
# scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
# scale_colour_manual(values=cbPalette)

##################################################
##  FIGURE 2A
##################################################

singletons <- data %>% filter( svtype == "ALL" & criteria == "all")
colnames(singletons)[1]<-"source"
singletons$perc_singleton<-singletons$N_singleton/singletons$N_all
singletons$perc_singleton<-singletons$perc_singleton*100
s<-singletons

conf1<-print(prop.test(x = singletons$N_singleton[1], n = singletons$N_all[1],
                conf.level = .95,
                correct = FALSE))
conf2<-print(prop.test(x = singletons$N_singleton[2], n = singletons$N_all[2],
                conf.level = .95,
                correct = FALSE))
conf3<-print(prop.test(x = singletons$N_singleton[3], n = singletons$N_all[3],
                conf.level = .95,
                correct = FALSE))
conf4<-print(prop.test(x = singletons$N_singleton[4], n = singletons$N_all[4],
                conf.level = .95,
                correct = FALSE))

# Variance
singletons$variance[1] <- ((conf1$conf.int[2] - conf1$conf.int[1]) / 3.919928) ^ 2
singletons$variance[2] <- ((conf2$conf.int[2] - conf2$conf.int[1]) / 3.919928) ^ 2
singletons$variance[3] <- ((conf3$conf.int[2] - conf3$conf.int[1]) / 3.919928) ^ 2
singletons$variance[4] <- ((conf4$conf.int[2] - conf4$conf.int[1]) / 3.919928) ^ 2

# Weights
singletons$weight<-1/singletons$variance
library(metafor)
library(meta)
# Meta
s<-singletons%>% select(N_singleton, N_all, source)
colnames(s)<-c("events", "n", "study")
s$proportion<-s$events/s$n

meta_result <- metagen(proportion, n, data = s)
summary(meta_result)  
forest(meta_result)

# Weighted Mean
wm<-weighted.mean(singletons$perc_singleton, singletons$weight)
# 0.4995907
# ci = 0.42

# install.packages("metafor")

# singletons$ci[1]<-meta_result$TE[1]
# singletons$ci[2]<-meta_result$TE[2]
# singletons$ci[3]<-meta_result$TE[3]
# singletons$ci[4]<-meta_result$TE[4]

for (i in 1:nrow(singletons)){
  singletons$ci[i]<-meta_result$TE[i]
}


# basic bar plot without meta analysis

p <-ggplot(singletons, aes(x=source, y=perc_singleton, fill=source)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())+
  geom_errorbar(data=singletons, aes(ymin = perc_singleton - ci, ymax = perc_singleton + ci),
                width = 0.5,
                position = position_dodge(), stat="identity")+
  theme_classic()+
  ggtitle("All SVs",subtitle = "Baseline") +
  labs(y="% Singleton") +
  scale_fill_manual(values=cbPalette)+
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5,face = "bold"),
    plot.subtitle = element_text(hjust = 0.5,vjust=1),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank())
p

# install.packages("ggside")
# library(ggside)

si<-singletons
new_row <- list(
  source = "summary",
  svtype = "ALL",
  criteria = "all",
  N_all = 0,
  N_singleton = 0,
  N_polymorphic = 0,
  perc_singleton = wm,
  variance = 0,
  weight = 0,
  ci = meta_result$TE.common
)
si <- rbindlist(list(as.data.table(new_row),si))

# adding diamond 

q <-ggplot(si, aes(x=source, y=perc_singleton, fill=source)) + 
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_errorbar(data=subset(si, source!='summary'), aes(ymin = perc_singleton - ci, ymax = perc_singleton + ci),
                width = 0.23,
                position = position_dodge(), stat="identity")+
  theme_classic()+
  geom_point(data=subset(si, source=='summary'), color='black', shape=18, size=4)+
  ggtitle("All SVs",subtitle = "Baseline") +
  labs(y="% Singleton") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73","#FFFFFF"))+
  theme(
    axis.text = element_text(size = 11, face="bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5,face = "bold"),
    plot.subtitle = element_text(hjust = 0.5,vjust=1),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none")
q

pdf(file="./fig2a.pdf")
q 
dev.off()

fig2a<-ggplot(si, aes(x=source, y=perc_singleton, fill=source)) + 
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_errorbar(data=subset(si, source!='summary'), aes(ymin = perc_singleton - ci, ymax = perc_singleton + ci),
                width = 0.23,
                position = position_dodge(), stat="identity")+
  theme_classic()+
  geom_point(data=subset(si, source=='summary'), color='black', shape=18, size=4)+
  ggtitle("All SVs",subtitle = "Baseline") +
  labs(y="% Singleton") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73","#FFFFFF"))+
  theme(
    axis.text = element_text(size = 11, face="bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5,face = "bold"),
    plot.subtitle = element_text(hjust = 0.5,vjust=1),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank())+
  # Set override.aes to remove the diamonds from the legend
  guides(fill = guide_legend(override.aes = list(shape = NA)))

pdf(file="./fig2a_legend_label.pdf")
plot_grid(fig2a, labels=c("A"))
dev.off()


pdf(file="./fig2a_legend.pdf")
fig2a
dev.off()






##################################################
##  FIGURE 2B
##################################################

# baseline percentages of overall SV % singleton
baseline<-si %>% filter(source!="summary") %>% select(source, perc_singleton, ci) 
colnames(data)[1]<-"source"
svtype<-data %>% filter(svtype!="ALL" & criteria == "all") %>% 
  group_by(source, svtype) %>% 
  summarise(all_n = sum(N_all), all_single = sum(N_singleton))
svtype<-as.data.frame(svtype)

singletons%>%select(source, N_all, N_singleton, N_polymorphic)

# source  N_all N_singleton N_polymorphic
# 1:    1000G30X 153198       69323         83875
# 2:    1000GPh3  63627       23835         39792
# 3:        CCDG 142650       82996         59654
# 4: gnomAD_v2.1 291931      145846        146085

svtype %>% 
  group_by(source) %>% 
  summarise(all_n = sum(all_n), all_single = sum(all_single))

# Define the custom colors for each 'svtype'
svtype_colors <- c("DEL" = "#D43925",
                   "DUP" = "#2376B2",
                   "CNV" = "#7459B2",
                   "INS" = "#D474E0",
                   "INV" = "#FA931E",
                   "CPX" = "#71E38C",
                   "OTH" = "#397246")

# Define the custom shades for each 'study'
study_shades <- c("1000G30X" = "#999999",
                  "1000GPh3" = "#E69F00",
                  "CCDG" = "#56B4E9",
                  "gnomAD_v2.1" = "#009E73",
                  "summary" = "#FFFFFF")

svtype$perc_single<-svtype$all_single/svtype$all_n*100
study<-baseline$source

for (i in 1:nrow(svtype)){
  if (svtype$source[i]==study[1]){
    svtype$delta[i]<-svtype$perc_single[i]-baseline$perc_singleton[1]
  }
  else if (svtype$source[i]==study[2]){
    svtype$delta[i]<-svtype$perc_single[i]-baseline$perc_singleton[2]
  }
  else if (svtype$source[i]==study[3]){
    svtype$delta[i]<-svtype$perc_single[i]-baseline$perc_singleton[3]
    
  }else if (svtype$source[i]==study[4]){
    svtype$delta[i]<-svtype$perc_single[i]-baseline$perc_singleton[4]
  } else {}
}

sv<-svtype %>% arrange(svtype)


#  Generate Meta
sv_categories<-unique(sv$svtype)

# Create an empty list to store the data frames
subtype_data_list <- list()

# Create a for loop to iterate through each subtype
for (subtype in sv_categories) {
  # Create a data frame for each subtype and store it in the list
  subtype_data_list[[subtype]] <- data.frame()
}

for (i in 1:length(sv_categories)){
  type<-sv %>% 
    filter(svtype==sv_categories[i]) %>%
    select(all_single, all_n, source)
  colnames(type)<-c("events", "n", "study")
  type$proportion<-type$events/type$n
  subtype_data_list[[sv_categories[i]]]<-type
}

library(meta)
# Create an empty list to store the metagen results
metagen_results_list <- list()

for (subtype in sv_categories) {
  # Extract the data frame for the current subtype
  current_data <- subtype_data_list[[subtype]]
  
  # Run the metagen function on the current data frame
  meta_result <- metagen(proportion, n, data = current_data)
  
  # Store the metagen result in the 'metagen_results_list' with a specific name
  metagen_results_list[[paste0("meta_result_", subtype)]] <- meta_result
}


# Split the dataframe by 'svtype'
split_svs <- split(sv, sv$svtype)


# Adding CIs
for (subtype in sv_categories) {
  split_svs[[subtype]]$ci<-metagen_results_list[[paste0("meta_result_", subtype)]]$TE
}


# Calculate variance and weight for each sv type
# Function to calculate variance and weight for each row of a data frame
calculate_variance <- function(df) {
  # Apply prop.test to each row and calculate variance and weight
  df$variance <- numeric(nrow(df))
  df$weight <- numeric(nrow(df))
  
  for (i in 1:nrow(df)) {
    prop_results <- prop.test(x = df$all_single[i], n = df$all_n[i], conf.level = 0.95, correct = FALSE)
    if (!is.null(prop_results$p.value)) {
      df$variance[i] <- ((prop_results$conf.int[2] - prop_results$conf.int[1]) / 3.919928) ^ 2
      df$weight[i] <- 1 / df$variance[i]
    } else {
      df$variance[i] <- NA
      df$weight[i] <- NA
    }
  }
  
  # Return the modified data frame
  return(df)
}

# Apply calculate_variance function to each data frame in the split_svs list
split_svs_with_variance_weight <- lapply(split_svs, calculate_variance)


## Calculate Weighted Mean per SV Type

weighted_mean_result_list<-list()

for (subtype in sv_categories) {
  current_data <- split_svs_with_variance_weight[[subtype]]
  weight_mean_result<-weighted.mean(current_data$perc_single, current_data$weight)
  weighted_mean_result_list[[paste0("weight_result_", subtype)]] <- weight_mean_result
}

# Adding WM
for (subtype in sv_categories) {
  split_svs_with_variance_weight[[subtype]]$wm<-weighted_mean_result_list[[paste0("weight_result_", subtype)]]
}

# Create an empty list to store the modified data frames
split_svs_with_summary <- list()

# Loop through each data frame in the split_svs list
for (i in seq_along(split_svs_with_variance_weight)) {
  # Add the corresponding summary_row to the current data frame
  current_df <- split_svs_with_variance_weight[[i]]
  current_summary_row <- list(
    source = "summary",
    svtype = "ALL",
    all_n = 0,
    all_single = 0,
    perc_single = split_svs_with_variance_weight[[i]]$wm[1],
    delta=0, #NEED TO CHANGE
    ci = metagen_results_list[[i]]$TE.common,
    variance = 0,
    weight = 0,
    wm = split_svs_with_variance_weight[[i]]$wm[1]
  )
  current_df <- rbind(current_df, current_summary_row)
  
  # Append the modified data frame to the split_svs_with_summary list
  split_svs_with_summary[[i]] <- current_df
}

# Fill in the empty studies for CPX 
# Get the unique studies present in both data frames
studies_in_svs1 <- unique(split_svs_with_summary[[1]]$source)
studies_in_svs2 <- unique(split_svs_with_summary[[2]]$source)

# Identify the missing studies in split_svs[[1]]
missing_studies <- setdiff(studies_in_svs2, studies_in_svs1)

# Create empty rows for the missing studies in split_svs[[1]]
empty_rows <- data.frame(source = missing_studies,
                         svtype = character(length(missing_studies)),
                         all_n = numeric(length(missing_studies)),
                         all_single = numeric(length(missing_studies)),
                         perc_single = numeric(length(missing_studies)),
                         delta = numeric(length(missing_studies)),
                         ci = numeric(length(missing_studies)),
                         variance = numeric(length(missing_studies)),
                         weight = numeric(length(missing_studies)),
                         wm = numeric(length(missing_studies))
                         )

# Combine split_svs[[1]] and the empty rows using rbind
split_svs_with_summary[[1]] <- rbind(split_svs_with_summary[[1]], empty_rows)

# Function to reorder the data frame with "source == 'summary'" row last
reorder_summary_last <- function(df) {
  # Reorder rows using order() function to put "source == 'summary'" last
  df <- df[order(df$source == "summary"), ]
  
  # Return the modified data frame
  return(df)
}

# Apply reorder_summary_last function to each data frame in the split_svs_with_summary list
split_svs_with_summary <- lapply(split_svs_with_summary, reorder_summary_last)

#Fix Delta for Summary
for (i in seq_along(split_svs_with_summary)) {
  # Add the corresponding summary_row to the current data frame
  split_svs_with_summary[[i]]$delta[5] <- split_svs_with_summary[[i]]$perc_single[5]-49.39729
}


create_bar_plot <- function(df) {
  ggplot(df, aes(x=source, y=delta, fill=source)) + 
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_errorbar(data=subset(df, source!='summary'), aes(ymin = delta - abs(ci), ymax = delta + abs(ci)),
                 width = 0.23,
                 position = position_dodge(), stat="identity")+
  # geom_errorbar(data=df, aes(ymin = delta - abs(ci), ymax = delta + abs(ci)),
  #                width = 0.23,
  #                position = position_dodge(), stat="identity")+
  theme_classic()+
  geom_point(data=subset(df, source=='summary'), color='black', shape=18, size=3)+
  ggtitle(df$svtype[1]) +
  labs(y="%") +
  scale_fill_manual(values=study_shades)+
    scale_y_continuous(limits=c(-60, 60))+
  theme(
    axis.text = element_text(size = 10, face="bold"),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5,vjust = 1.5,face = "bold",color=svtype_colors[[df$svtype[1]]]),
    plot.subtitle = element_text(hjust = 0.5,vjust=1),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none")
    #+ geom_hline(yintercept = 0, color = "black")  # Add the horizontal line at y = 0
}


# Use lapply to create a list of plots
plots_list <- lapply(split_svs_with_summary, create_bar_plot)

#Print the individual plots
for (i in seq_along(plots_list)) {
  print(plots_list[[i]])
}

# install.packages("gridExtra")
library("gridExtra")
# install.packages("cowplot")
library("cowplot")


# x<-plot_grid(q+theme(
#   axis.text = element_text(size = 10, face="bold")), print(plots_list[[1]]), print(plots_list[[2]]), print(plots_list[[3]]), print(plots_list[[4]]), print(plots_list[[5]]), labels=c("A", "B", "", "", "", ""), ncol = 6, nrow = 2)

#individual plot
pdf(file="./fig2b.pdf")
for (i in seq_along(plots_list)) {
   print(plots_list[[i]])
 }
dev.off()

#combo plot
pdf(file="./fig2b_combo.pdf", height=5, width=12)
plot_grid(print(plots_list[[1]]), print(plots_list[[2]]), print(plots_list[[3]]), print(plots_list[[4]]), print(plots_list[[5]]), labels=c("B", "", "", "", ""), ncol = 5, nrow = 1)
dev.off()

grouped_svs<-split_svs_with_summary

# Add the "group" column to each data frame in the list
for (i in seq_along(grouped_svs)) {
  grouped_svs[[i]]$group <- grouped_svs[[i]]$svtype[1]
}

# Combine the data frames into one data frame with the "group" column
combined_df <- bind_rows(grouped_svs)

# Print the combined data frame
print(combined_df)

# Create a named vector for group colors
group_colors <- c("#71E38C", "#D43925", "#2376B2", "#D474E0", "#FA931E")

z<-ggplot(combined_df, aes(x=source, y=delta, fill=source)) + 
  geom_bar(stat="identity", 
           position=position_dodge())+
  facet_wrap(~group, scales="free", nrow=1)+ 
  geom_errorbar(data=subset(combined_df, source!='summary'), aes(ymin = delta - abs(ci), ymax = delta + abs(ci)),
                width = 0.23,
                position = position_dodge(), stat="identity")+
  theme_classic()+
  geom_point(data=subset(combined_df, source=='summary'), color='black', shape=18, size=3)+
  scale_fill_manual(values=study_shades)+
  scale_y_continuous(limits=c(-60, 30))+
  theme(
    axis.text = element_text(size = 10, face="bold"),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 8, colour = c("#71E38C", "#D43925", "#2376B2", "#D474E0", "#FA931E")))

# Remove the y-axis ticks labels from all facets except the first one
z_no_y_axis <- z + theme(axis.text.y = element_blank(), 
                         axis.ticks.y = element_blank(), 
                         axis.title.y = element_blank(),
                         axis.line = element_blank()) + 
              geom_hline(yintercept = 0, color = "black")

# Create an empty plot with only the y-axis labels on the left
y_axis_labels <- ggplot(data.frame(y = c(-60, 30)), aes(x = 0, y = y)) +
  theme_classic()+
  scale_y_continuous(limits = c(-60, 30)) +
  labs(y= "Percent (%)") +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(7.5, 7, 2, 5, "mm"))
  
library(cowplot)
# Combine the plots using plot_grid
plot_grid(y_axis_labels, z_no_y_axis, ncol = 2, align = "b", axis = "t", rel_widths = c(0.1, 1))


pdf(file="./fig2b_facets_V2.pdf", height=5, width=12)
plot_grid(y_axis_labels, z_no_y_axis, ncol = 2, align = "b", axis = "t", rel_widths = c(0.1, 1), labels=c("B"))
dev.off()

pdf(file="./fig2b_facets_no_letter_V2.pdf", height=5, width=12)
plot_grid(y_axis_labels, z_no_y_axis, ncol = 2, align = "b", axis = "t", rel_widths = c(0.1, 1))
dev.off()


##################################################
##  Figure 2C
##################################################
# Finally, for panel C, I think what we want to show is the odds ratio of comparing different subsets of deletions to all deletions per study, again in a similar way as panel B. Instead of subtracting the overall % singletons per study from each value, though, you'll want to run a fisher.test() in R for each study with data formatted in the following 2x2 table:
# +-------------------------------------------+-------------------------+
# | all_polymorphic - polymorphic_in_category | polymorphic_in_category |
# +-------------------------------------------+-------------------------+
# |   all_singleton - singleton_in_category   |   singleton_in_category |
# +-------------------------------------------+-------------------------+

or<-data

stats_all<-data %>% filter(svtype=="DEL", criteria =="all")

or<-or %>% filter(criteria!="all")
or<-or[order(or$source),]

for (i in 1:nrow(or)){
  or$ap[i]<-stats_all$N_polymorphic[stats_all$source==or$source[i]]-or$N_polymorphic[i]
  or$as[i]<-stats_all$N_singleton[stats_all$source==or$source[i]]-or$N_singleton[i]
}

fisher_prep<-or%>%select(source, criteria, N_singleton, N_polymorphic, ap, as)

fisher_df <- as.data.frame(fisher_prep)

# Split the data by both "source" and "criteria"
split_data <- split(fisher_df, list(fisher_df$source, fisher_df$criteria))

#print(split_data)

convert_to_2x2 <- function(df) {
  ap_npolymorphic <- df[c("ap", "N_polymorphic")]
  as_nsingleton <- df[c("as", "N_singleton")]
  ap_npolymorphic <- setNames(ap_npolymorphic, c("delta", "value"))
  as_nsingleton <- setNames(as_nsingleton, c("delta", "value"))
  return(rbind(ap_npolymorphic, as_nsingleton))
}

# Apply the function to each data frame in split_data
split_data_2x2 <- lapply(split_data, convert_to_2x2)

#print(split_data_2x2)

# Define a function to apply fisher.test to each data frame
run_fisher_test <- function(df) {
  fisher_result <- fisher.test(df)
  result <- list(#odds_ratio = fisher_result$estimate,
                 log_or = log2(fisher_result$estimate),
                 conf_int = fisher_result$conf.int[1:2],
                 p_value = fisher_result$p.value)
  return(result)
}

# Apply the fisher.test to each data frame in split_data_2x2
fisher_values <- lapply(split_data_2x2, run_fisher_test)

print(fisher_values)

# Create an empty data frame with column names "LL", "UP", "LogOR", and "source"
result_df <- data.frame(LL = numeric(), UP = numeric(), LogOR = numeric(), source = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Function to convert fisher_values for a specific suffix into a data table
convert_to_data_table <- function(suffix) {
  # Get the list names with the specified suffix
  list_names <- grep(suffix, names(fisher_values), value = TRUE)
  
  # Combine the odds ratios and confidence intervals into a data table
  data_table <- rbindlist(lapply(list_names, function(name) {
    list_data <- fisher_values[[name]]
    data.frame(LL = log2(list_data$conf_int[1]), UP = log2(list_data$conf_int[2]), LogOR = list_data$log_or, source = name, p_value = list_data$p_value) }
    ))
  
  return(data_table)
}

# Suffixes to iterate through
suffixes <- c("exonic", "intergenic", "intronic", "promoter", "utr")

# Convert fisher_values for each suffix into data tables and rbind them together
for (suffix in suffixes) {
  suffix_data <- convert_to_data_table(suffix)
  result_df <- rbind(result_df, suffix_data)
}

# Print the resulting data frame
print(result_df)

result<-result_df

# Split the source column into two separate columns using the '.' as the separator
df <- cbind(result_df, str_split_fixed(result_df$source, "\\.(?=[^.]+$)", 2))
df <- df %>% select(LL, UP, LogOR, V1, V2, p_value)
colnames(df)[4:5]<-c('source', 'criteria')
#df

pl<-ggplot(df, aes(x=source, y=LogOR, fill=source)) + 
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_errorbar(data=subset(df, source!='summary'), aes(ymin = LL, ymax = UP),
                width = 0.23,
                position = position_dodge(), stat="identity")+
  facet_wrap(~criteria, scales="free", nrow=1)+ 
  scale_y_continuous(limits=c(-1, 2))+
  # geom_errorbar(data=df, aes(ymin = delta - abs(ci), ymax = delta + abs(ci)),
  #                width = 0.23,
  #                position = position_dodge(), stat="identity")+
  theme_classic()+
  #geom_point(data=subset(df, source=='summary'), color='black', shape=18, size=3)+
  scale_fill_manual(values=study_shades)+
  theme(
    axis.text = element_text(size = 10, face="bold"),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none")
pl

# Remove the y-axis ticks labels from all facets except the first one
pl_no_axis <- pl + theme(axis.text.y = element_blank(), 
                         axis.ticks.y = element_blank(), 
                         axis.title.y = element_blank(),
                         axis.line = element_blank()) + 
  geom_hline(yintercept = 0, color = "black")

# Create an empty plot with only the y-axis labels on the left
y_axis_labels <- ggplot(data.frame(y = c(-1, 2)), aes(x = 0, y = y)) +
  theme_classic()+
  scale_y_continuous(limits = c(-1, 2)) +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(7.5, 7, 2, 5, "mm"))

pdf(file="./fig2c_v3.pdf", height=5, width=12)
plot_grid(y_axis_labels, pl_no_axis, ncol = 2, align = "b", axis = "t", rel_widths = c(0.1, 1))
dev.off()

## TODO: Add Meta for c (weight & variance), sort by affect size 

#  Generate Meta
del_type<-unique(df$criteria)


#df2<-df
df$LL <- as.numeric(df$LL)
df$UP <- as.numeric(df$UP)
df$LogOR <- as.numeric(df$LogOR)

df$variance <- ((df$UP - df$LL) / 3.919928)^ 2
df$SE <- (df$UP - df$LL) / (2 * 1.96)
df$weight <- 1 / df$variance

# Create an empty list to store the data frames
subtype_list <- list()

# Create a for loop to iterate through each subtype
for (subtype in del_type) {
  print(subtype)
  # Create a data frame for each subtype and store it in the list
  subtype_list[[subtype]] <- data.frame()
}

for (i in 1:length(del_type)){
  type<-df %>% 
    filter(criteria==del_type[i])
  subtype_list[[del_type[i]]]<-type
}

# Create an empty list to store the metagen results
metagen_or_list <- list()

for (subtype in del_type) {
  # Extract the data frame for the current subtype
  current_data <- subtype_list[[subtype]]
  
  # Run the metagen function on the current data frame
  meta_result <- metagen(
    TE = current_data$LogOR,
    seTE =current_data$SE,
    studlab = current_data$source,
    weights = current_data$weight,
    comb.fixed = TRUE
  )
  
  # Store the metagen result in the 'metagen_results_list' with a specific name
  metagen_or_list[[paste0("meta_or_", subtype)]] <- meta_result
}

# Meta, Calculate Weighted Mean per Criteria
metagen_or_list[[i]]$TE.common

# Split the dataframe by 'svtype'
# split_svs <- split(sv, sv$svtype)

del_subtype_list <- subtype_list
# Adding WM
for (subtype in del_type) {
  del_subtype_list[[subtype]]$wm<-metagen_or_list[[paste0("meta_or_", subtype)]]$TE.common
}

del_subtype_list

# Create an empty list to store the modified data frames
split_or_with_summary <- list()

# Loop through each data frame in the split_svs list
for (i in seq_along(del_subtype_list)) {
  # Add the corresponding summary_row to the current data frame
  current_df <- del_subtype_list[[i]]
  current_summary_row <- list(
    source = "summary",
    criteria = del_subtype_list[[i]]$criteria[1],
    LL = metagen_or_list[[i]]$lower.common,
    UP = metagen_or_list[[i]]$upper.common,
    wm = del_subtype_list[[i]]$wm[1],
    LogOR = del_subtype_list[[i]]$wm[1],
    p_value = metagen_or_list[[i]]$pval.common,
    variance = 0,
    SE = 0,
    weight = 0
  )
  current_df <- rbind(current_df, current_summary_row)
  
  # Append the modified data frame to the split_svs_with_summary list
  split_or_with_summary[[i]] <- current_df
}

split_or_with_summary

grouped_or<-split_or_with_summary

# Add the "group" column to each data frame in the list
for (i in seq_along(grouped_or)) {
  grouped_or[[i]]$group <- grouped_or[[i]]$criteria[1]
}

# Combine the data frames into one data frame with the "group" column
combined_or <- bind_rows(grouped_or)

# Print the combined data frame
print(combined_or)

# Define a custom order for criteria
criteria_order <- c("exonic", "utr", "promoter", "intronic", "intergenic")

study_order<-c("1000G30X","1000GPh3","CCDG", "gnomAD_v2.1", "summary")

# Convert 'group' column to a factor with custom levels
combined_or$group <- factor(combined_or$group, levels = criteria_order)

# Convert 'criteria' column to a factor with custom levels
combined_or$criteria <- factor(combined_or$criteria, levels = study_order)

# Sort the data frame based on 'criteria', 'group', and 'LogOR'
sorted_data <- combined_or[order(combined_or$criteria, combined_or$group), ]

# Print the sorted data
print(sorted_data)


por<-ggplot(sorted_data, aes(x=source, y=LogOR, fill=source)) + 
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_errorbar(data=subset(sorted_data, source!='summary'), aes(ymin = LL, ymax = UP),
                width = 0.23,
                position = position_dodge(), stat="identity")+
  facet_wrap(~group, scales="free", nrow=1)+ 
  scale_y_continuous(limits=c(-0.5, 1.5))+
  theme_classic()+
  geom_point(data=subset(sorted_data, source=='summary'), color='black', shape=18, size=3)+
  scale_fill_manual(values=study_shades)+
  theme(
    axis.text = element_text(size = 10, face="bold"),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none")
por

# Remove the y-axis ticks labels from all facets except the first one
por_no_axis <- por + theme(axis.text.y = element_blank(), 
                         axis.ticks.y = element_blank(), 
                         axis.title.y = element_blank(),
                         axis.line = element_blank()) + 
  geom_hline(yintercept = 0, color = "black")

# Create an empty plot with only the y-axis labels on the left
y_axis_labels <- ggplot(data.frame(y = c(-0.5, 1.5)), aes(x = 0, y = y)) +
  theme_classic()+
  scale_y_continuous(limits = c(-0.5, 1.5)) +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(7.5, 7, 2, 5, "mm"))

pdf(file="./fig2c_meta_v4.pdf", height=5, width=12)
plot_grid(y_axis_labels, por_no_axis, ncol = 2, align = "b", axis = "t", rel_widths = c(0.1, 1))
dev.off()

