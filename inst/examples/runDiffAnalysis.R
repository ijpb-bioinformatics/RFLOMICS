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

MAE <- runDiffAnalysis(MAE, SE.name = "protetest", 
                       method = "limmalmFit", 
                       contrastList = contrastList)
# or
# MAE[["protetest"]] <- runDiffAnalysis(MAE[["protetest"]], 
#                                       method = "limmalmFit", 
#                                       contrastList = contrastList)

# Filter the results of the differential analysis with new cut-off values 
# for p-value and fold change.
MAE <- filterDiffAnalysis(MAE, SE.name = "protetest",
                          p.adj.cutoff = 0.01,
                          logFC.cutoff = 1)
# or
# MAE[["protetest"]] <- filterDiffAnalysis(MAE[["protetest"]], 
#                                          p.adj.cutoff = 0.01, 
#                                          logFC.cutoff = 1)

# Access to the diff analysis settings
## Get DE matrix from DiffExpAnalysis
head(getDEMatrix(MAE[["protetest"]]))

## Get union or intersection from list of contrasts
getDEList(MAE[["protetest"]], contrasts = ) 

## Get diff setting
getDiffSettings(MAE[["protetest"]])


# generate plot results of a differential analysis
thiscontrast <- "(temperatureMedium - temperatureLow) in mean"

## generate MAplot from diff analysis
plotDiffAnalysis(MAE[["protetest"]], 
                 contrastName = thiscontrast, 
                 typeofplots = "MA.plot")

## plot the heatmap
plotHeatmapDesign(MAE[["protetest"]], 
                  contrastName = thiscontrast)

## plot boxplot with feature expression
plotBoxplotDE(MAE[["protetest"]], 
              features = "AT1G47128", 
              groupColor = "temperature")
plotBoxplotDE(MAE[["protetest"]], 
              features = "AT1G79550", 
              groupColor = "imbibition")

