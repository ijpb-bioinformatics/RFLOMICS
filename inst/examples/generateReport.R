ExpDesign <- RFLOMICS::readExpDesign(
  file = paste0(system.file(package = "RFLOMICS"), 
                "/ExamplesFiles/ecoseed/condition.txt"))
factorRef <- data.frame(
  factorName  = c("Repeat", "temperature" , "imbibition"),
  factorRef   = c("rep1",   "Low",          "DS"),
  factorType  = c("batch",  "Bio",          "Bio"),
  factorLevels= c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI"))

omicsData <- list(
  RFLOMICS::readOmicsData(
    file = paste0(system.file(package = "RFLOMICS"), 
                  "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")))

omicsNames  = c("RNAtest", "metatest", "protetest")
omicsTypes  = c("RNAseq","metabolomics","proteomics")
MAE <- RFLOMICS::createRflomicsMAE(projectName = "Tests",
                                   omicsData   = omicsData,
                                   omicsNames  = omicsNames,
                                   omicsTypes  = omicsTypes,
                                   ExpDesign   = ExpDesign,
                                   factorRef   = factorRef)

formulae <- RFLOMICS::generateModelFormulae(MAE)
MAE <- RFLOMICS::setModelFormula(MAE, modelFormula = formulae[[1]])

Contrasts.List <- RFLOMICS::generateExpressionContrast(MAE)
selectedContrasts <- Contrasts.List$simple[1:4,]

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

#generateReport(object = MAE, fileName = ecoseed_report.html)