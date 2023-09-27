library(testthat)
library(RFLOMICS)

## ---- Construction MAE RFLOMICS ready for integration analysis : ----
MAE <- RFLOMICS::FlomicsMultiAssay.constructor(
  list("RNAtest"     = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
                            "omicType" = "RNAseq"),
       "metatest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
                         "omicType" = "metabolomics"),
       "protetest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")), 
                          "omicType" = "proteomics")
  ),
  projectName = "Tests", 
  ExpDesign = RFLOMICS::read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt")),
  refList = c("Repeat" = "rep1", "temperature" = "Low", "imbibition" = "DS"),
  typeList = c("Repeat" = "batch", "temperature" = "Bio", "imbibition" = "Bio"))

formulae <- RFLOMICS::GetModelFormulae(MAE = MAE) 
MAE <- MAE |>
  RFLOMICS::getExpressionContrast(model.formula = formulae[[1]]) 
MAE <- MAE  |> RFLOMICS::getContrastMatrix(contrastList = c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                                            "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3",
                                                            "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" )) 
contrastsDF <- RFLOMICS::getSelectedContrasts(MAE)

MAE2 <- MAE

MAE <- MAE |>
  TransformData(     SE.name = "metatest",  transform_method = "log2")          |>
  RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")            |>
  RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")                 |>
  RunNormalization(  SE.name = "protetest", NormMethod = "median")              |>
  FilterLowAbundance(SE.name = "RNAtest")                                       |>
  RunDiffAnalysis(   SE.name = "metatest",  DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "protetest", DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "RNAtest",   DiffAnalysisMethod = "edgeRglmfit") |>
  FilterDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 1.5) |>
  integrationWrapper(omicsToIntegrate = c("protetest", "metatest"), silent = TRUE, cmd = FALSE) |>
  integrationWrapper(omicsToIntegrate = c("RNAtest", "metatest"),   silent = TRUE, cmd = FALSE, method = "mixOmics")

# ----- TESTS -----

test_that("Design accessors", {
  expect_identical(getFactorTypes(MAE), MAE@metadata$design@Factors.Type)
  
})


# ---- omics dictionnary ----

test_that("Omics dictionnary", {
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "RNAtest"), list(variableName = "genes", valueType = "counts"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "metatest"), list(variableName = "metabolites", valueType = "XIC"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "protetest"), list(variableName = "proteins", valueType = "XIC"))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "metatest"), RFLOMICS:::omicsDic(MAE[["metatest"]]))
  expect_identical(RFLOMICS:::omicsDic(MAE, SE.name = "RNAtest"), RFLOMICS:::omicsDic(MAE[["RNAtest"]]))
  
  expect_error(RFLOMICS:::omicsDic(MAE))

})

# ---- factor types ----

test_that("Factors types", {
  
  expect_error(getFactorTypes(MAE[["RNAtest"]]))
  
  expect_identical(bioFactors(MAE), c("temperature", "imbibition"))
  expect_identical(batchFactors(MAE), c("Repeat"))
  
})
