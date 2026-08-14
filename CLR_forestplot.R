#CLR calculation:
install.packages("compositions")
library(tidyverse)
library(compositions)
#loading up data files - metadata and feature/OTU file created in qiime2
metadata <- read.delim("~/Desktop/Project/Results/results_metadata.tsv")
head(metadata)
household <- read.csv("~/Desktop/Project/Results/Household_data.csv")
View(household)
meta1 <- metadata %>% rename(SampleID = sample.id)
meta <- meta1 %>% left_join(household)
View(meta)
colnames(household)
colnames(meta)
"Household" %in% colnames(meta)
otu <- read.delim("~/Desktop/Project/Results/exported-table/feature-table.tsv", skip = 1, row.names = 1, check.names = FALSE)
#creating a dataframe for clr - transposing the table:
otu_t <- as.data.frame(t(otu))
#adding the pseudocount (pseudocount to prevent logging 0 values)
otu_pseudo <- otu_t + 1
#transforming abundance values onto a scale (clr = centered log-ratio) (https://microbiome.github.io/OMA/docs/devel/pages/transformation.html, https://www.reddit.com/r/bioinformatics/comments/16yloem/which_transformation_method_to_use_with/)
otu_clr <- clr(otu_pseudo)
#converting back to plain data format
otu_clr <- as.data.frame(unclass(otu_clr))
rownames(otu_clr) <- rownames(otu_pseudo)
colnames(otu_clr)
otu_clr$SampleID <- rownames(otu_clr)
#merging with metadata
merged <- left_join(otu_clr, meta, by = c("SampleID"))
View(merged)
#identifying species of interest:
species_of_interest <- colnames(otu_clr)[colnames(otu_clr) %in% c(
  "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__",
  "d__Bacteria;p__Bacillota;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__",
  "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella;s__",
  "d__Bacteria;p__Actinomycetota;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium;s__",
  "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__",
  "d__Bacteria;p__Bacillota;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__",
  "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__",
  "d__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__",
  "d__Bacteria;p__Bacillota;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Peptostreptococcaceae;g__Intestinibacter;s__"
)]
length(species_of_interest)
species_of_interest
#running the linear model per species with an lapply loop:
results_list <- lapply(species_of_interest, function(sp) {
  f <- as.formula(paste0("`", sp, "` ~ Malaria + Age + Sex + (1|Household)")) #defining the model and setting household as a random effect
  model <- lm(f, data = merged)
  coefs <- summary(model)$coefficients
  row <- grep("^Malaria", rownames(coefs))
  
  data.frame(
    species = sp,
    estimate = coefs[row, "Estimate"],
    se = coefs[row, "Std. Error"],
    p_value = coefs[row, "Pr(>|t|)"]
  )
})

#row bind/structuring a dataframe with one row per genus
results_df <- do.call(rbind, results_list)
#applying the benjamini-hochberg correction to adjust p-values (accounting for the fact we are testing multiple genera)
results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
#extracting only the genus name to make the graphs easier to read
results_df$genus <- sub(".*g__([^;]+);s__$", "\\1", results_df$species)
#calculating the 95% confidence interval for each genus 
results_df$lower <- results_df$estimate - 1.96 * results_df$se
results_df$upper <- results_df$estimate + 1.96 * results_df$se
#labelling genera as significant or not based on p value (i.e. more or less than 0.05)
results_df$significant <- ifelse(results_df$p_adj < 0.05, "Significant (q < 0.05)", "Not significant")

#plot creation
#ordering the data from smallest to largest
results_df <- results_df %>% arrange(estimate)
#converting genus column into an ordered factor:
results_df$genus <- factor(results_df$genus, levels = results_df$genus)

results_df$significant <- factor(
  results_df$significant,
  levels = c("Significant (q < 0.05)", "Not significant")
)
p_forest <- ggplot(results_df, aes(x = estimate, y = genus, color = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 0.7, show.legend = TRUE) +
  geom_point(size = 3, show.legend = TRUE) +
  scale_color_manual(values = c("Significant (q < 0.05)" = "#C44E52", "Not significant" = "grey40"), drop = FALSE) +
  labs(
    x = "CLR-scale effect size (95% CI)\n(adjusted for age, sex, household)",
    y = NULL,
    color = NULL
  ) +
  theme_bw()
#show the plot:
p_forest



                    