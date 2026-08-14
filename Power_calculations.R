#mPower calculations
#the mPower power analysis was done following the guide created by the creators of the tool (https://github.com/chloelulu/mPower)
#Yang, L, & Chen, J., (2026), "mPower: a real data-based power analysis tool for microbiome study design", Microbiome, 14:196, doi:10.1186/s40168-026-02427-4  
library(mPower)
packageVersion("mPower")
install.package"mPower"install.packages("mPower")
#loading the file 
feature.data <- read.delim("~/Desktop/Project/Results/exported-table/feature-table.tsv", skip = 1, row.names = 1, check.names = FALSE)
feature.data <- as.matrix(feature.data)
dim(feature.data)
head(rownames(feature.data))
head(colnames(feature.data))
#preprocessing, excluding features present in fewer than two samples
feature.data <- feature.data[rowSums(feature.data != 0) > 2, ]
#estimating parameters
model.paras <- EstPara(ref.otu.tab = feature.data)
dim(feature.data)
sample_depths <- colSums(feature.data)
mean(sample_depths)
sd(sample_depths)
ls("package:mPower")
#calculation of community level power 
res1 <- mPower(feature.dat = feature.data,
               model.paras = model.paras, 
               test = 'Community', 
               design = 'CaseControl', 
               nSams = c(10, 29), 
               grp.ratio = 0.5, 
               iters = 500, 
               alpha = 0.05, 
               distance = 'BC', 
               diff.otu.pct = 0.1, 
               covariate.eff.min = 0, 
               covariate.eff.maxs = 2, 
               diff.otu.direct = 'balanced', 
               diff.otu.mode = 1, 
               confounder = 'no', 
               depth.mu = 18398.86, 
               depth.sd = 3841.693, 
               verbose = F)
#outputting the value for the community level power
knitr::kable(res1$power, format = "markdown")
res1$plot
feature.data <- as.matrix(feature.data)
feature.data <- feature.data[rowSums(feature.data != 0) > 2, ]
model.paras <- EstPara(ref.otu.tab = feature.data)
sample_depths <- colSums(feature.data)
#calculating the taxa level power
res2 <- mPower(
  feature.dat = feature.data,
  model.paras = model.paras,
  test = 'Taxa',
  design = 'CaseControl',
  nSams = c(15, 29),
  grp.ratio = 0.5,
  iters = 500, 
  alpha = 0.05,
  diff.otu.direct = 'balanced',
  diff.otu.mode = 1,
  covariate.eff.min = 0,
  covariate.eff.maxs = 2,
  prev.filter = 0.1,
  max.abund.filter = 0.002,
  confounder ='no',
  depth.mu = mean(sample_depths),
  depth.sd = sd(sample_depths),
  verbose = F
)
#outputting the values for the taxa level power 
knitr::kable(res2$aTPR, format = "markdown")
knitr::kable(res2$pOCR, format = "markdown")

