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

# generate all statistical model formulae
formulae <- generateModelFormulae(MAE)

# chose and set model formula to rflomicsMAE object
MAE <- setModelFormula(MAE, formulae[[1]])

# Generate expression of contrasts from chosen model
contrastList <- generateExpressionContrast(MAE, "averaged")

# Set the contrasts List and choose the first 3 contrasts of type averaged
MAE <- setSelectedContrasts(MAE, contrastList = contrastList[c(1,2,3),])
