library(testthat)
library(RFLOMICS)

# ---- Construction MAE RFLOMICS ready for differential analysis : ----

MAE <- generateExample(annotation = FALSE, coexp = FALSE, integration = FALSE)

# ----- TESTS -----

test_that("Design accessors", {
  expect_identical(getFactorTypes(MAE), metadata(MAE)$design$Factors.Type)
  
})

# ---- factor types ----

test_that("Factors types", {
  
  expect_identical(getBioFactors(MAE), c("temperature", "imbibition"))
  expect_identical(getBatchFactors(MAE), c("Repeat"))
  
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



