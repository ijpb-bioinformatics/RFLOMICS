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

formulae <- RFLOMICS::generateModelFormulae(MAE)
MAE <- RFLOMICS::setModelFormula(MAE, modelFormula = formulae[[1]])

selectedContrasts <- 
  RFLOMICS::generateExpressionContrast(MAE, contrastType="simple")

MAE <- setSelectedContrasts(MAE, contrastList = selectedContrasts)

## data processing
MAE <- RFLOMICS::runDataProcessing(
  object = MAE, 
  SE.name = "protetest", 
  samples=NULL, 
  normMethod="none", 
  transformMethod="none")

## diff analysis
MAE <- RFLOMICS::runDiffAnalysis(
  object = MAE, 
  SE.name = "protetest", 
  contrastList = 
    selectedContrasts, 
  p.adj.method="BH", 
  method = "limmalmFit",  
  p.adj.cutoff = 0.05, 
  logFC.cutoff = 0)

## Enrichment
MAE <- RFLOMICS::runAnnotationEnrichment(
  object = MAE, 
  SE.name = "protetest", 
  database = "GO", 
  domain = c("MF"), 
  list_args = list(OrgDb = "org.At.tair.db", 
                   keyType = "TAIR", 
                   pvalueCutoff = 0.05))

# get name of performed analysis
getAnalyzedDatasetNames(MAE)

# generate report
#generateReport(object = MAE, fileName = ecoseed_report.html)