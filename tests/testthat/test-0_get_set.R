library(testthat)
library(RFLOMICS)

# ---- Construction MAE RFLOMICS ready for differential analysis : ----
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

contrastList <- RFLOMICS::getExpressionContrast(object = MAE, model.formula = formulae[[1]]) |> purrr::reduce(rbind) |>
  dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
                                "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))


MAE2 <- MAE

MAE <- MAE |>
  TransformData(     SE.name = "metatest",  transformMethod = "log2")          |>
  RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")            |>
  RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")                 |>
  RunNormalization(  SE.name = "protetest", NormMethod = "median")              |>
  FilterLowAbundance(SE.name = "RNAtest")                                       |>
  RunDiffAnalysis(   SE.name = "metatest",  DiffAnalysisMethod = "limmalmFit", contrastList = contrastList, modelFormula = formulae[[1]])  |>
  RunDiffAnalysis(   SE.name = "protetest", DiffAnalysisMethod = "limmalmFit", contrastList = contrastList, modelFormula = formulae[[1]])  |>
  RunDiffAnalysis(   SE.name = "RNAtest",   DiffAnalysisMethod = "edgeRglmfit", contrastList = contrastList, modelFormula = formulae[[1]]) |>
  FilterDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 1.5) |>
  integrationWrapper(omicsToIntegrate = c("protetest", "metatest"), silent = TRUE, cmd = FALSE) |>
  integrationWrapper(omicsToIntegrate = c("RNAtest", "metatest"),   silent = TRUE, cmd = FALSE, method = "mixOmics")

# ----- TESTS -----

test_that("Design accessors", {
  expect_identical(getFactorTypes(MAE), MAE@metadata$design@Factors.Type)
  
})

# ---- factor types ----

test_that("Factors types", {
  
  expect_error(getFactorTypes(MAE[["RNAtest"]]))
  
  expect_identical(bioFactors(MAE), c("temperature", "imbibition"))
  expect_identical(batchFactors(MAE), c("Repeat"))
  
})

# ---- omics dictionnary ----

# Test of internal function
test_that("Omics dictionnary", {
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "RNAtest"), list(variableName = "transcripts", valueType = "counts"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "metatest"), list(variableName = "metabolites", valueType = "XIC"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "protetest"), list(variableName = "proteins", valueType = "XIC"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "metatest"), RFLOMICS:::omicsDic(MAE[["metatest"]]))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "RNAtest"), RFLOMICS:::omicsDic(MAE[["RNAtest"]]))
  
  expect_error(RFLOMICS:::omicsDic(MAE))

})

# ---- opDEList - operation on DE lists ----
 
test_that("Operation DE lists", {
  
  expect(!identical(opDEList(MAE, SE.name = "RNAtest", contrasts = c("H1", "H2"), operation = "union"),
                         opDEList(MAE, SE.name = "RNAtest", contrasts = c("H1"), operation = "union")), 
         failure_message = "union of (H1, H2) should be different from H1 alone")
  
  expect(!identical(opDEList(MAE, SE.name = "RNAtest", contrasts = c("H1", "H2"), operation = "intersection"),
                    opDEList(MAE, SE.name = "RNAtest", contrasts = c("H1"), operation = "intersection")), 
         failure_message = "Intersection of (H1, H2) should be different from H1 alone")

  H1H2DE <- getDE(object = MAE[["RNAtest"]], contrast = c("H1", "H2"))
  expect_identical(opDEList(MAE, SE.name = "RNAtest", contrasts = c("H1", "H2"), operation = "union"),
                   H1H2DE$DEF)
  expect_identical(opDEList(MAE, SE.name = "RNAtest", contrasts = c("H1", "H2"), operation = "intersection"),
                   H1H2DE$DEF[rowSums(H1H2DE[,-1]) == 2])
  
  H1H2DE <- getDE(object = MAE[["protetest"]], contrast = c("H1", "H2"))
  expect_identical(opDEList(MAE, SE.name = "protetest", contrasts = c("H1", "H2"), operation = "union"),
                   H1H2DE$DEF)
  expect_identical(opDEList(MAE, SE.name = "protetest", contrasts = c("H1", "H2"), operation = "intersection"),
                   H1H2DE$DEF[rowSums(H1H2DE[,-1]) == 2])
  

})
 


