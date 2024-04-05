library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in the code or package versions affect the expected results. 

# ---------------- RUN RFLOMICS ---------------

## construct rflomics + PCA raw + check completness
ExpDesign <- RFLOMICS::readExpDesign(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))
factorRef <- data.frame(factorName  = c("Repeat", "temperature" , "imbibition"),
                        factorRef   = c("rep1",   "Low",          "DS"),
                        factorType  = c("batch",  "Bio",          "Bio"),
                        factorLevels= c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI"))

omicsData <- list(
    RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
    RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
    RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")))

MAE <- RFLOMICS::createRflomicsMAE(projectName = "Tests", 
                                   omicsData   = omicsData,
                                   omicsNames  = c("RNAtest.raw", "metatest.raw", "protetest.raw"),
                                   omicsTypes  = c("RNAseq","metabolomics","proteomics"),
                                   ExpDesign   = ExpDesign,
                                   factorRef   = factorRef)

formulae <- RFLOMICS::generateModelFormulae(object = MAE) 
MAE <- RFLOMICS::setModelFormula(MAE, modelFormula = formulae[[1]])

# ----- 
test_that("generateExpressionContrast", {
  
  Contrasts.List.m <- RFLOMICS::generateExpressionContrast(MAE)
  Contrasts.List.f <- RFLOMICS:::.getExpressionContrastF(ExpDesign = MAE@colData, factorBio=c("temperature", "imbibition"), modelFormula = formulae[[1]])
  
  expect_equal(Contrasts.List.m, Contrasts.List.f)
})


