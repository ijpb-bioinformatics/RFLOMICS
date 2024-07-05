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

# Set the statistical model and contrasts to test
formulae <- generateModelFormulae(MAE)
MAE <- setModelFormula(MAE, formulae[[1]])  

# Get the contrasts List and choose the first 3 contrasts of type averaged
contrastList <- generateExpressionContrast(MAE, "averaged")

MAE <- setSelectedContrasts(MAE, contrastList = contrastList[c(1, 2, 3),])

# Run the data preprocessing and perform the differential analysis 
MAE <- runDataProcessing(MAE, SE.name = "protetest",  
                         transformMethod = "log2",
                         normMethod = "median")
MAE <- runDiffAnalysis(MAE, SE.name = "protetest")

# Run co-expression analysis
MAE <- runCoExpression(MAE, SE.name = "protetest", 
                       K = 2:5, replicates = 5, 
                       merge = "union")

# get parametres used to run co-expression analysis
coExp.set.list <- getCoexpSettings(MAE[["protetest"]])
coExp.set.list$method

# get results
clusters <- getCoseqClusters(MAE[["protetest"]])
proteins_c1 <- getClusterEntities(MAE[["protetest"]], 
                                      clusterName = "cluster_1")

# plots
#plotCoExpression(MAE[["protetest"]])

#plotCoExpressionProfile(MAE[["protetest"]], cluster = 2)

#plotCoseqContrasts(MAE[["protetest"]])
