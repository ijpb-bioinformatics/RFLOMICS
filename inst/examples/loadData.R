# load ecoseed data
data(ecoseed)

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = list(ecoseed$RNAtest, ecoseed$metatest, ecoseed$protetest),
  omicsNames  = c("RNAtest", "metatest", "protetest"),
  omicsTypes  = c("RNAseq","metabolomics","proteomics"),
  ExpDesign   = ecoseed$design,
  factorRef   = ecoseed$factorRef)

# generate upset plot
MultiAssayExperiment::upsetSamples(MAE)

# generate data overview plot
plotDataOverview(MAE)

# generate plot of coverage of condition by data
plotConditionsOverview(MAE)