#Creation of plots for the alpha diversity, beta diversity and PERMANOVA results from Qiime2
library(tidyverse)
library(dplyr)
library(ggplot2)
#loading the shannon diversity file into R
metadata <- read_csv("~/Desktop/Project/Results/Metadata_processing.csv")
View(metadata)
#filtering out NA values/rows
clean_metadata <- metadata %>%
  drop_na()
View(clean_metadata)
#loading in household data:
household <- read_csv("~/Desktop/Project/Results/Household_data.csv")
View(household)
#merging with data
final_metadata <- clean_metadata %>% left_join(household)
View(final_metadata)
#running Mann-Whitney U test (because two separate groups that are unpaired aka each individual is measured once and can only belong to one group)
#and the kruskal-wallis test
wilcox <- wilcox.test(shannon_entropy ~ Malaria, data=final_metadata)
kruskal.test(shannon_entropy ~ Malaria, data=final_metadata)
#also running linear model:
alpha_model <- lm(shannon_entropy ~ Malaria + Age + Sex + Household, data = final_metadata)
summary(alpha_model)
#creation of boxplots for the clean presentation of data
library(ggsignif)

library(ggbeeswarm)

#creation of the boxplots 
final_metadata %>%
  ggplot(aes(x=Malaria, y=shannon_entropy)) +
  geom_boxplot(aes(colour = Malaria), fill = "white", alpha = 0.001) +
  geom_beeswarm(aes(colour = Malaria))+
  scale_color_manual(values = c("Control" = "#4C72B0", "Infected" = "#C44E52", "Control" = "#4C72B057"),
                     labels = c("Control" = "Control (N = 15)", "Infected" = "Infected (N = 14)")) +
  geom_signif(
    comparisons = list(c("Control", "Infected")),
    map_signif_level = FALSE,
    annotations = "Wilcoxon p = 0.0697\nAdjusted model p = 0.119",
    y_position = max(final_metadata$shannon_entropy, na.rm = TRUE) + 0.05,
    tip_length = 0.02,
    colour = "black",
    alpha = 1
  ) +
  ylim(NA, max(final_metadata$shannon_entropy, na.rm = TRUE) + 0.2) +
  xlab("Malaria Status") +
  ylab("Shannon Diversity Index") +
  theme(legend.position = "none")
  
#loading bray-curtis into R 
pcoa_df <- read_csv("~/Desktop/Project/Results/pcoa_df.csv")
#loading the variance explained file containing the percentage of variance explained by each PCoA axis
variance_explained <- read.csv("~/Desktop/Project/Results/pcoa_variance_explained.csv")
head(variance_explained)
#extracting each variance variable from the file 
pc1_var <- variance_explained$PC1
pc2_var <- variance_explained$PC2
#manually storing hte PERMANOVA R2 value and p-value
r2_lab <- 0.036
p_lab_fmt <- "= 0.580"
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Malaria)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(type = "t", linewidth = 0.6) +
  scale_color_manual(
    values = c("Control" = "#4C72B0", "Infected" = "#C44E52"),
    labels = c("Control" = "Control (N = 15)", "Infected" = "Infected (N = 14)")) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    caption = paste0("PERMANOVA: R\u00B2 = ", round(r2_lab, 3), ", p ", p_lab_fmt)
  ) +
  theme_bw()
#the R2 value for this shows how much variance is explained by infection status Malaria/infection status explains only 3.6% of the variance in Bray-Curtis community composition, and this is nowhere near significant (p = 0.580).
