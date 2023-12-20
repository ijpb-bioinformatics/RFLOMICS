library(testthat)
library(RFLOMICS)

# ---- Construction MAE RFLOMICS ready for differential analysis : ----

MAE <- generateExample(annotation = FALSE, coexp = FALSE, integration = FALSE)

# ----- TESTS -----

test_that("Design accessors", {
  expect_identical(getFactorTypes(MAE), MAE@metadata$design@Factors.Type)
  
})

# ---- factor types ----

test_that("Factors types", {
  
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



