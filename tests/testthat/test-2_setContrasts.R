library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in the code or package versions affect the expected results. 

# ---------------- RUN RFLOMICS ---------------

# load ecoseed data
data(ecoseed)

# create rflomicsMAE object with ecoseed data
MAE <- createRflomicsMAE(
    projectName = "Tests",
    omicsData   = list(ecoseed$RNAtest, ecoseed$metatest, ecoseed$protetest),
    omicsNames  = c("RNAtest", "metatest", "protetest"),
    omicsTypes  = c("RNAseq","metabolomics","proteomics"),
    ExpDesign   = ecoseed$design,
    factorRef   = ecoseed$factorRef)


formulae <- generateModelFormulae(object = MAE) 
MAE <- setModelFormula(MAE, modelFormula = formulae[[1]])

# ----- 
test_that("generateExpressionContrast", {
    
    Contrasts.List.m <- generateExpressionContrast(MAE)
    Contrasts.List.f <- RFLOMICS:::.getExpressionContrastF(ExpDesign = colData(MAE), 
                                                           factorBio = c("temperature", "imbibition"), 
                                                           modelFormula = formulae[[1]])
    
    expect_equal(Contrasts.List.m, Contrasts.List.f)
})


