library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----

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
  FilterDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 1.5)

# ----- TESTS -----
 
test_that("Working?", code = {
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest"))
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest"))
  
  getMOFA(MAE)
  
  
  expect(!is.null(MAE@metadata$MOFA$MOFA_results), failure_message = "MOFA didn't run as it should have run")
  
  expect_error(integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "proteotest"), contrasts_names = c("H1", "H2")))
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "protetest"), contrasts_names = c("H1", "H2"))
  

  
  # TODO tests: with feature selection, without, with only one response variable, with several, etc. 
  # TODO compare: with the same methods outside of RFLOMICS. (normalisation, data treatment, find the equivalent pipeline)
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest"))
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest"), silent = TRUE, cmd = TRUE)  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest"), silent = FALSE, cmd = TRUE)  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("protetest", "metatest"), silent = TRUE, cmd = TRUE)  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("protetest", "metatest"), silent = FALSE, cmd = TRUE)  
  
  MAE@metadata$MOFA$MOFA_results
  
  MOFA2::plot_data_overview(MAE@metadata$MOFA$MOFA_results)
  
  resMOFA <- MAE@metadata$MOFA$MOFA_results
  MOFA2::plot_factor_cor(resMOFA)
  MOFA2::plot_variance_explained(resMOFA, plot_total = TRUE)[[2]]
  MOFA2::views_names(resMOFA)
  resMOFA@dimensions$K
  
}) 


