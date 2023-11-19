library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----

RNAdat <- RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))
corresp <- read.table(file = paste0(system.file(package = "RFLOMICS"),"/ExamplesFiles/ecoseed/transcript_genes.txt"), sep = "\t", header = TRUE)

# Use gene id to ease use of clusterprofiler
RNAdat <- RNAdat[corresp$ensembl_transcript_id,]
rownames(RNAdat) <- corresp$ensembl_gene_id[match(corresp$ensembl_transcript_id, rownames(RNAdat))]

## ---- Construction MAE RFLOMICS ready for CPR analysis : ----
ExpDesign <- RFLOMICS::read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))
factorRef <- data.frame(factorName  = c("Repeat", "temperature" , "imbibition"),
                        factorRef   = c("rep1",   "Low",          "DS"),
                        factorType  = c("batch",  "Bio",          "Bio"),
                        factorLevels= c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI"))

omicsData <- list(
  RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
  RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
  RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")))

MAE <- RFLOMICS::FlomicsMultiAssay.constructor(projectName = "Tests", 
                                               omicsData   = omicsData,
                                               omicsNames  = c("RNAtest", "metatest", "protetest"),
                                               omicsTypes  = c("RNAseq","metabolomics","proteomics"),
                                               ExpDesign   = ExpDesign,
                                               factorRef   = factorRef)
names(MAE) <- c("RNAtest", "metatest", "protetest")


formulae <- RFLOMICS::GetModelFormulae(MAE = MAE) 
MAE <- MAE |>
  RFLOMICS::getExpressionContrast(model.formula = formulae[[1]]) 
MAE <- MAE  |> RFLOMICS::getContrastMatrix(contrastList = c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                                            "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3",
                                                            "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" )) 

MAE2 <- MAE

MAE <- MAE |>
  TransformData(     SE.name = "metatest",  transformMethod = "log2")          |>
  RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")            |>
  RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")                 |>
  RunNormalization(  SE.name = "protetest", NormMethod = "median")              |>
  FilterLowAbundance(SE.name = "RNAtest") 

MAE[["metatest"]]@metadata$DataProcessing$done  <- TRUE
MAE[["protetest"]]@metadata$DataProcessing$done <- TRUE
MAE[["RNAtest"]]@metadata$DataProcessing$done   <- TRUE

MAE <- MAE |>
  RunDiffAnalysis(   SE.name = "metatest",  DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "protetest", DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "RNAtest",   DiffAnalysisMethod = "edgeRglmfit") |>
  FilterDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 1.5) |>
  runCoExpression(   SE.name = "RNAtest",   K = 2:12, replicates = 2)   |>
  runCoExpression(   SE.name = "protetest", K = 2:12, replicates = 2)   |>
  runCoExpression(   SE.name = "metatest",  K = 2:12, replicates = 2)


# ---- Annotation test function - DiffExpEnrichment ----
# 

test_that("it's running from diffExpAnal - GO - RNASeq", {
  
  # Selecting only one contrast
  expect_no_error({
   MAE <- runAnnotationEnrichment(MAE, nameList = "H1", SE.name = "RNAtest", dom.select = "GO",
                                     list_args = list(OrgDb = "org.At.tair.db", 
                                                      keyType = "TAIR", 
                                                      pvalueCutoff = 0.05),
                                     Domain = c("BP", "MF", "CC"))
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H1", ont = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  # All contrasts
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", dom.select = "GO",
                                     list_args = list(OrgDb = "org.At.tair.db", 
                                                      keyType = "TAIR", 
                                                      pvalueCutoff = 0.05),
                                     Domain = c("BP", "MF", "CC"))
                  
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H2", ont = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  

  
})

test_that("it's running from diffExpAnal - Custom - RNASeq", {
  
  df_custom <- vroom::vroom(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/GO_annotations/Arabidopsis_thaliana_Ensembl_55.txt"))
  
  MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", dom.select = "custom",
                                     list_args = list(pvalueCutoff = 0.05),
                                     col_term = "GO term accession", 
                                     col_gene = "Gene stable ID",
                                     col_name = "GO term name",
                                     col_domain = "GO domain",
                                     annot = df_custom)
  
  expect(!is.null(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$custom$summary), failure_message = "Custom RNAseq didn't work (summary not present)")
  
  # Selecting only one contrast, custom file
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, nameList = "H1", SE.name = "RNAtest", dom.select = "custom",
                                       list_args = list(pvalueCutoff = 0.05),
                                       col_term = "GO term accession", 
                                       col_gene = "Gene stable ID",
                                       col_name = "GO term name",
                                       col_domain = "GO domain",
                                       annot = df_custom)
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H1", ont = "custom", domain = "biological_process")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  
  # All contrasts, GO database
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", dom.select = "GO",
                                       list_args = list(OrgDb = "org.At.tair.db", 
                                                        keyType = "TAIR", 
                                                        pvalueCutoff = 0.05),
                                       Domain = c("BP", "MF", "CC"))
    
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H2", ont = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  
  
})


# ---- Annotation test function - CoExpression enrichment ----

test_that("it's running from CoExpAnal - GO - RNASeq", {
  
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", from = "CoExpAnal", dom.select = "GO",
                                       list_args = list(OrgDb = "org.At.tair.db", 
                                                        keyType = "TAIR", 
                                                        pvalueCutoff = 0.05),
                                       Domain = c("BP", "MF", "CC"))
    
  })
  
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "cluster.1", ont = "GO", domain = "BP", from = "CoExpAnal")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from CoExp) - There is no result in the enrichment metadata part.")

  
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", nameList = c("cluster.1", "cluster.2") ,
                                       from = "CoExpAnal", dom.select = "GO",
                                       list_args = list(OrgDb = "org.At.tair.db", 
                                                        keyType = "TAIR", 
                                                        pvalueCutoff = 0.05),
                                       Domain = c("BP", "MF", "CC"))
    
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "cluster.1", ont = "GO", domain = "BP", from = "CoExpAnal")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from CoExp) - There is no result in the enrichment metadata part.")
  
})



plotCPR(MAE[["RNAtest"]], 
        contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
        ont = "GO",  
        Domain = "BP", 
        type = "cnetplot", searchExpr = "bou")

.doNotPlot(plotCPR(MAE[["RNAtest"]], 
                   contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
                   ont = "GO",  
                   Domain = "BP", 
                   type = "cnetplot", searchExpr = "bou"))

plotCPR(MAE[["RNAtest"]], 
        contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
        ont = "GO",  
        Domain = "BP", 
        type = "cnetplot", searchExpr = "bou")


outheat <- tryCatch(plotCPR(MAE[["RNAtest"]], 
                            contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
                            ont = "GO",  
                            Domain = "BP", 
                            type = "heatplot"),
                    error = function(e) e,
                    warnings = function(w) w)

outheat
class(outheat)

outdot <- tryCatch(plotCPR(MAE[["RNAtest"]], 
                            contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
                            ont = "GO",  
                            Domain = "BP", 
                            type = "dotplot"),
                    error = function(e) e,
                    warnings = function(w) w)

outdot
class(outdot)

is(outdot, "gg")

outcnet <- tryCatch(plotCPR(MAE[["RNAtest"]], 
                 contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
                 ont = "GO",  
                 Domain = "BP", 
                 type = "cnetplot"),
         error = function(e) e,
         warnings = function(w) w)

outcnet

outcnete <- tryCatch(plotCPR(MAE[["RNAtest"]], 
                            contrast = "(temperatureElevated - temperatureLow) in imbibitionDS",
                            ont = "GO",  
                            Domain = "BP", searchExpr = "ERR" ,
                            type = "cnetplot"),
                    error = function(e) e,
                    warnings = function(w) w)
outcnete
