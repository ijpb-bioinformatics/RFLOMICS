library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for integration analysis : ----
MAE <- generateExample(
  annotation = FALSE,
  integration = FALSE
) 

# ----- TESTS -----
 
test_that("Working?", code = {
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest"))
  
  
  expect(!is.null(MAE@metadata$MOFA$MOFA_results), failure_message = "MOFA didn't run as it should have run")
  expect_error(integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "proteotest"), contrasts_names = c("H1", "H2")))
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "protetest"), contrasts_names = c("H1", "H2"))

# 
#   MAE2 <- MAE |>
#     FilterDiffAnalysis(SE.name = "protetest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 0)
  
  # MAE2 <- MAE
  # MAE2 <- integrationWrapper(MAE2, omicsToIntegrate = c("metatest", "protetest"), contrasts_names = c("H1", "H2"), type = "intersection")
  
  MAE2 <- integrationWrapper(MAE2, omicsToIntegrate = c("metatest", "protetest"), contrasts_names = c("H1", "H2"), type = "intersection")
  getPossibleContrasts(MAE[["RNAtest"]], returnTable = TRUE)
  getValidContrasts(MAE[["RNAtest"]])
  MAE[["metatest"]] <- setValidContrasts(MAE[["metatest"]], c("H1", "H2"))
  MAE[["protetest"]] <- setValidContrasts(MAE[["protetest"]], c("H1", "H2"))
  
  getValidContrasts(MAE[["protetest"]])

  RFLOMICS:::convertTagToContrast(MAE, c("H1", "H2"))
  RFLOMICS:::convertContrastToTag(MAE, c("(temperatureElevated - temperatureLow) in mean", "(imbibitionEI - imbibitionDS) in mean" ))
  
  
  allrownames <- lapply( c("metatest", "protetest"), FUN = function(nam){
    try_rflomics(opDEList(MAE, SE.name = nam, contrasts = c("H1", "H2"), operation = 'intersection'))
  })
  names(allrownames) <- c("metatest", "protetest")
  allrownames
  
  allrownames <- lapply( c("metatest", "protetest"), FUN = function(nam){
    tryCatch(opDEList(MAE, SE.name = nam, contrasts = c("H1", "H2"), operation = 'intersection'), 
             error = function(e) e, 
             warning = function(w) w)
  })
  names(allrownames) <- c("metatest", "protetest")
  any(sapply(allrownames, class) != "character")
  
  if (any(sapply(allrownames, class) != "character")) {
    probOmics <- names(allrownames)[sapply(allrownames, class) != "character"]
    stop("It seems there  is a problem with : ", sapply(probOmics, FUN = function(namesError) paste("\n", probOmics, " -- error message:\n", allrownames[[probOmics]])))
  }
  
  MAE@metadata$MOFA$MOFA_results

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

