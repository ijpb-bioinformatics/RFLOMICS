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
  FilterDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 2)

# ----- TESTS -----

test_that("Working?", code = {
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest"), method = "mixomics") 
  expect({!is.null( getMixOmics(MAE, response = "temperature"))}, failure_message = "There is no temperature results here")
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("metatest"), method = "mixomics")
  expect({!is.null( getMixOmics(MAE, response = "temperature"))}, failure_message = "There is no temperature results here")

  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest"), method = "mixomics")
  expect({!is.null( getMixOmics(MAE, response = "temperature"))}, failure_message = "There is no temperature results here")
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "protetest"), 
                            contrasts_names = c("H1", "H2"), method = "mixOmics")
  
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "protetest"), 
                            contrasts_names = c("H1", "H2"), method = "mixOmics", 
                            selectedResponse = c("temperature", "imbibition"))
  
  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("RNAtest", "metatest", "protetest"), 
                            contrasts_names = c("H1", "H2"), method = "mixOmics", 
                            selectedResponse = c("temperature"))
  
  expect({is.null( getMixOmics(MAE, response = "imbibition"))}, failure_message = "There is imbibition results here")
  
  expect(identical(names(MAE@metadata$mixOmics), c("temperature", "imbibition"), attrib.as.set = FALSE), 
         failure_message = "Taking only two responses for mixOmics does not work")
  

  MAE <- integrationWrapper(MAE, omicsToIntegrate = c("metatest", "protetest"), 
                            contrasts_names = c("H1", "H2"), method = "mixOmics", 
                            selectedResponse = c("temperature", "imbibition"), sparsity = TRUE, 
                            cmd = TRUE, ncomp = 2,
                            cases_to_try = 10)
  
  mixOmics::plotDiablo(MAE@metadata$mixOmics$temperature$MixOmics_results)
  
})

mixOmics::circosPlot(MAE@metadata$mixOmics$temperature$MixOmics_results, cutoff = 0.95)

cor_mat <- mixOmics::network(MAE@metadata$mixOmics$temperature$MixOmics_results,
                             plot.graph = FALSE)
try_rflomics(mixOmics::circosPlot(MAE@metadata$mixOmics$temperature$MixOmics_results, cutoff = 1))



MAE <- integrationWrapper(MAE, omicsToIntegrate = c("protetest", "metatest"), method = "mixomics")

outw <- tryCatch({mixOmics::circosPlot(MAE@metadata$mixOmics$temperature$MixOmics_results, cutoff = 1)},
                error = function(e) e,
                warning = function(w) w
                )
class(outw)

out <- tryCatch({mixOmics::circosPlot(MAE@metadata$mixOmics$temperature$MixOmics_results, cutoff = 0.9)},
                error = function(e) e,
                warning = function(w) w
)
out

class(out)
out[1:5,1:5]
?mixOmics::circosPlot

.doNotPlot <- function(expr){
  pdf(file = NULL)
  out <- tryCatch({capture.output(eval(expr))},
    error = function(e) e,
    warning = function(w) w
    )
  dev.off()
  return(out)
}

outMat <- .doNotPlot(mixOmics::circosPlot(MAE@metadata$mixOmics$temperature$MixOmics_results, cutoff = 0.75))
outMat 

outW <- .doNotPlot(mixOmics::circosPlot(MAE@metadata$mixOmics$temperature$MixOmics_results, cutoff = 1))
outW 

out <- mixOmics::network(mat = MAE@metadata$mixOmics$temperature$MixOmics_results, 
                         blocks = 1:2,
                         cutoff = 0.5, 
                         shape.node = "rectangle")

outElse <- .doNotPlot({
  mixOmics::network(mat = MAE@metadata$mixOmics$temperature$MixOmics_results, 
                    blocks = 1:2,
                    cutoff = 0.5, 
                    shape.node = "rectangle")
})
class(outElse)

outE <- .doNotPlot({
  mixOmics::network(mat = MAE@metadata$mixOmics$temperature$MixOmics_results, 
                    blocks = 1:2,
                    cutoff = 1, 
                    shape.node = "rectangle")
})

class(outE)


outE <- RFLOMICS:::.doNotSpeak({
  mixOmics::network(mat = MAE@metadata$mixOmics$temperature$MixOmics_results, 
                    blocks = 1:2,
                    cutoff = 1, 
                    shape.node = "rectangle")
})

class(outE)
