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
                       method = "limmalmFit")

# Run GO annotation (enrichGO)
MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
                               list_args = list(OrgDb = "org.At.tair.db",
                                                keyType = "TAIR",
                                                pvalueCutoff = 0.05),
                               from = "DiffExp", database = "GO",
                               domain = "CC")

# Run KEGG annotation (enrichKEGG)
MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
                               list_args = list(organism = "ath",
                                                keyType = "kegg",
                                                pvalueCutoff = 0.05),
                               from = "DiffExp", database = "KEGG")


# Search for the pvalue cutoff:
sumORA(MAE[["protetest"]], from = "DiffExp", database = "KEGG")

# need internet connection
#plotKEGG(MAE[["protetest"]], pathway_id = "ath00710", species = "ath",
#         contrastName = "cluster.4", from = "Coexp")

# From differential analysis proteins lists:
plotClusterProfiler(MAE[["protetest"]],
                    contrastName = "(temperatureElevated - temperatureMedium) in mean",
                    database = "KEGG", from = "DiffExp",
                    plotType = "heatplot", p.adj.cutoff = 0.05,
                    domain = "no-domain")

plotEnrichComp(MAEtest[["protetest"]], from = "DiffExp",
               database = "KEGG", matrixType = "FC")

# Get all results from KEGG on differential expression lists:
results <- getEnrichRes(MAE[["protetest"]],
                        from = "diffexp", database = "KEGG")

# Search for the pvalue cutoff:
usedPvalue <- 
  getEnrichPvalue(MAE[["protetest"]], from = "diffexp", database = "KEGG")
settings <- 
  getEnrichSettings(MAE[["protetest"]], from = "diffexp", database = "KEGG")

