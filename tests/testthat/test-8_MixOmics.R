library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for integration analysis : ----
MAE <- generateExample(integration = FALSE, annotation = FALSE)

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
