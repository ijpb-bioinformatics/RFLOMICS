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

# Set the statistical model
formulae <- generateModelFormulae(MAE)
MAE <- setModelFormula(MAE, formulae[[1]])  

# set the contrast list
contrastList <- generateExpressionContrast(MAE, "averaged")
MAE <- setSelectedContrasts(MAE, contrastList = contrastList[c(1, 2, 3),])

# Data processing of RNAseq dataset : RNAtest
## using data processing functions for RNAseq data
### filter low RNAseq count 
# MAE <- filterLowAbundance(MAE, SE.name = "RNAtest",
#                           filterStrategy = "NbReplicates", 
#                           cpmCutoff = 1)
# ### filter outlier samples
# MAE <- runSampleFiltering(MAE, SE.name = "RNAtest", 
#                           samples = colnames(MAE[["RNAtest"]])[-1])
# ### data normalisation outlier samples
# MAE <- runNormalization(MAE, SE.name = "RNAtest", 
#                         normMethod = "TMM")

## use runDataProcessing function that combines the previous three functions
MAE <- runDataProcessing(MAE, SE.name = "RNAtest", 
                         samples = colnames(MAE[["RNAtest"]])[-1], 
                         lowCountFiltering_strategy = "NbReplicates", 
                         lowCountFiltering_CPM_Cutoff = 1, 
                         normMethod = "TMM") 

## check completness of RNAtest data
checkExpDesignCompleteness(MAE, omicName = "RNAtest")$messages

# Data processing of proteimics dataset : protetest
# ## transform data
# MAE <- runTransformData(MAE, SE.name = "protetest",  transformMethod = "log2")
# ## normalise data
# MAE <- runNormalization(MAE, SE.name = "protetest", normMethod = "median") 

## use runDataProcessing function
MAE <- runDataProcessing(MAE, SE.name = "protetest", 
                         normMethod = "median", 
                         transformMethod = "log2")

plotExpDesignCompleteness(MAE[["RNAtest"]])

# plot Library Size
plotDataDistribution(MAE[["RNAtest"]], raw=TRUE)
plotDataDistribution(MAE[["RNAtest"]], raw=FALSE)

# plot gene expression distribution
plotDataDistribution(MAE[["RNAtest"]], raw=TRUE, plot = "boxplot")
plotDataDistribution(MAE[["RNAtest"]], raw=FALSE, plot = "boxplot")

# plot PCA 
plotOmicsPCA(MAE[["RNAtest"]], raw="raw", groupColor = "imbibition")
plotOmicsPCA(MAE[["RNAtest"]], raw="norm", groupColor = "imbibition")

