#Maaslin3 analysis
#the maaslin3 analysis was run closely following the tutorial produced by the creators of the tool (https://github.com/biobakery/biobakery/wiki/maaslin3)
#Nickols, W. A., Kuntz, T., Shen, J. et al. (2026), "MaAsLin3: Refining and extending generalized multivariable linear models for meta-omic association discovery", Nature Methods, 23, p.554-564, doi:10.1038/s41592-025-02923-9  
library(maaslin3)
packageVersion("maaslin3")
taxa_table <- read.delim("table-genus-relative/feature-table.tsv", row.names = 1, check.names = FALSE)
#extract genus names and put them into their own columns 
rownames(taxa_table) <- sub(".*g__([^;]+)$", "\\1", taxa_table[, 1])
#aggregate any duplicate genera together
taxa_table_summed <- aggregate(. ~ genus, data = taxa_table[, -1], FUN = sum)
#set unique row names and remove null values
rownames(taxa_table_summed) <- taxa_table_summed$genus
taxa_table_summed$genus <- NULL
#loading the metadata into the environment and then exporting it as .tsv and re-importing it 
meta <- read.csv("results_metadata.csv")
write.table(meta, file = "metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
metadata <- read.delim("metadata.tsv", row.names = 1, check.names = FALSE)

fit_out <- maaslin3(
  input_data = taxa_table,
  input_metadata = metadata,
  output = "maaslin3_output_household_final",
  formula = "~ malaria_status + Age + Sex + (1|Household)",
  max_significance = 0.05
  )
fit_out <- maaslin3(
  input_data = taxa_table,
  input_metadata = metadata,
  output = "maaslin3_output_household",
  formula = "~ Malaria.status + Age + Sex + (1|Household)"
  )

holdemanella <- read_tsv("~/Desktop/Project/Results/significant_results.tsv")
View(holdemanella)

