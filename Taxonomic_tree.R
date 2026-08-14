#Creation of taxonomic tree (References for this section; Lin, H., 2023, "ANCOM-BC2 Tutorial", available:https://bioconductor.uib.no/packages/3.18/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html ; Lin, H., 2025, "ANCOM-BC2 Tutorial", available: https://bioconductor.posit.co/packages/3.22/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html)
devtools::install_github('grunwaldlab/heattree')
packageVersion(heattree)
library(jsonlite)
library(dplyr)
library(qiime2R)
library(tidyverse)
library(metacoder)
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

#extracting the log-fold change and q values (adjusted p-values) from the ancom-bc files to enable visualisation of the changes and significance in species 
#turning the json fils from ancom-bc into a data table to enable the creation of a taxonomic tree (https://rdrr.io/cran/jsonlite/man/stream_in.html):
lfc <- stream_in(file("~/Desktop/Project/Results/ancombc_exported/lfc.jsonl"))
qval <- stream_in(file("~/Desktop/Project/Results/ancombc_exported/q.jsonl"))
#filtering the datatable to remove any NA values:
head(lfc)
lfc_clean <- lfc %>%
  filter(!is.na(taxon)) %>%
  select(taxon, lfc_malaria = `Malaria.status::Infected`)

head(lfc_clean)
head(qval)

qval_clean <- qval %>%
  filter(!is.na(taxon)) %>%
  select(taxon, q_malaria = `Malaria.status::Infected`)
head(qval_clean)
summary(qval_clean$q_malaria)

#joining files together using the taxon column
ancombc_results <- left_join(lfc_clean, qval_clean, by = "taxon")

#extracting the taxonomic data from the taxonomic classification step in the QIIME2 workflow
taxonomy <- read_qza("~/Desktop/Project/Results/results_taxonomy.qza")$data %>%
  rename(OTU_ID = Feature.ID, lineage = Taxon)
head(taxonomy)
#joining the taxonomy file with the existing file 
ancombc_results <- ancombc_results %>%
  left_join(taxonomy, by = c("taxon" = "OTU_ID"))
head(ancombc_results)
#importing the OTU data and joining with existing data-frame
otu_table <- read_qza("~/Desktop/Project/Results/table-filtered.qza")$data %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") #making sure that the format and column names match so that the two dataframes can be joined 
obj_input <- ancombc_results %>%
  rename(OTU_ID = taxon) %>%
  left_join(otu_table, by = "OTU_ID")
#looking at the dimensions of the object
dim(obj_input)
head(obj_input)
#filtering out taxa that only appear a few times in the dataset so that the taxonomic tree is not too big 
obj_input <- obj_input %>% filter(str_count(lineage, ";") == 6)
nrow(obj_input)
head(obj_input)
View(obj_input)
#importing the metadata
meta <- read.delim("~/Desktop/Project/Results/results_metadata.tsv") %>% rename(SampleID = sample.id)
#parsing the taxonomic classification data
obj <- parse_tax_data(
  obj_input,
  class_cols = "lineage",  #naming the column that contains the classification
  class_sep = ";",       #defining the separator
  class_key = c(tax_rank = "info", tax_name = "taxon_name"),
  class_regex = "^([a-z])__(.*)$"     #shortening the taxonomic names to make them more easily readable 
 )
dim(obj)
head(obj)
names(obj$data)[names(obj$data) == "tax_data"] <- "otu_data"
meta_filtered <- meta %>% 
  filter(SampleID %in% colnames(obj$data$otu_data))
head(meta_filtered)
View(meta_filtered)
#calculating the taxon abundance - conversion of the per observation counts (OTU) to per-taxon counts (https://www.rdocumentation.org/packages/metacoder/versions/0.3.8/topics/calc_taxon_abund)
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_data", cols = meta_filtered$SampleID)
#calculating and storing one lfc value for each taxon
obj$data$lfc_summary <- obj %>%
  obs("otu_data", subset = TRUE) %>%
  map_dbl(function(rows) {             #looping through the OTU rows and applying function to each
    if (length(rows) == 0) return(0) #handling taxa with no OTU rows
    current_otus <- obj$data$otu_data$OTU_ID[rows]
    sub_data <- obj_input %>% filter(OTU_ID %in% current_otus) #finding the OTUs in the differential abundance results
    mean(sub_data$lfc_malaria, na.rm = TRUE) #calculation of the mean log-fold change
  })
#applying the same sort of function but to find the significant taxa
obj$data$sig_summary <- obj %>%
  obs("otu_data", subset = TRUE) %>%
  map_lgl(function(rows) {
    if (length(rows) == 0) return(FALSE)
    current_otus <- obj$data$otu_data$OTU_ID[rows]
    sub_data <- obj_input %>% filter(OTU_ID %in% current_otus)
    any(sub_data$q_malaria < 0.05, na.rm = TRUE)
  })
true_ranks <- obj$taxon_ranks()

obj$data$sig_summary[true_ranks != "g" & true_ranks != "s"] <- FALSE
#again applying the same sort of function but for abundance values 
n_samples <- nrow(meta_filtered)
obj$data$prevalence <- obj %>%
  obs("otu_data") %>%
  map_dbl(function(rows) {
    if (length(rows) == 0) return(0)
    abund_sub <- obj$data$otu_data[rows, meta_filtered$SampleID, drop = FALSE]
    mean(colSums(abund_sub > 0) > 0)
  })
set.seed(1)
#replace missing lfc values i.e. if the function was unable to calculate an LFC value it should be left as 0
obj$data$lfc_summary[is.na(obj$data$lfc_summary)] <- 0
obj$data$sig_summary[is.na(obj$data$sig_summary)] <- FALSE
taxon_labels <- taxon_names(obj)
#highlighting the significant genera with a star
taxon_labels_starred <- ifelse(obj$data$sig_summary, paste0(taxon_labels, " *"), taxon_labels)
#checking which (if any genera have been starred)
head(taxon_labels_starred)
#creation of the taxonomic tree (Reference: https://github.com/grunwaldlab/heat-tree)
heat_tree(
  obj,
  node_label = taxon_labels_starred,
  node_size = obj$data$prevalence,
  node_color = obj$data$lfc_summary,
  node_color_range = c("#08519c", "#f0f0f0", "#a50f15"),
  node_color_trans = "linear",
  node_color_interval = c(-2, 2),
  node_size_axis_label = "Prevalence",
  node_color_axis_label = "ANCOM-BC2 log fold change (Infected vs Control)",
  layout = "davidson-harel",
  initial_layout = "reingold-tilford"
)
#saving the tree to the local device 
ggsave("~/Desktop/Project/heat_tree_infected_vs_control_22.pdf", width = 12, height = 12)

