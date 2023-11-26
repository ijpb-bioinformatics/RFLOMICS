library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in the code or package versions affect the expected results. 

# ---------------- RUN RFLOMICS ---------------

## construct rflomics + PCA raw + check completness
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

## choice of model formulae
formulae <- RFLOMICS::GetModelFormulae(MAE = MAE)

## list of all hypothesis grouped by contarst type (run on MAE or SE)
Contrasts.List <- RFLOMICS::getExpressionContrast(MAE, modelFormula = formulae[[1]])

## choice of hypothesis
selcetedContrasts <- rbind(Contrasts.List$simple[1:3,],
                           Contrasts.List$averaged[1:3,],
                           Contrasts.List$interaction[1:3,])

## data processing
sampleToKeep <- dplyr::filter(as.data.frame(MAE@colData), temperature != "Elevated") %>% row.names()
MAE <- MAE |> RFLOMICS::runDataProcessing(SE.name = "RNAtest"  , samples=sampleToKeep, lowCountFiltering_strategy="NbReplicates", lowCountFiltering_CPM_Cutoff=1, normalisation_method="TMM") |>
              RFLOMICS::runDataProcessing(SE.name = "protetest", samples=NULL, normalisation_method="none", transformation_method="none") |>
              RFLOMICS::runDataProcessing(SE.name = "metatest" , samples=NULL, normalisation_method=NULL, transformation_method="log2")


selcetedContrasts <- getExpressionContrast(MAE[["RNAtest"]], modelFormula=formulae[[1]]) %>% purrr::reduce(rbind) %>% dplyr::filter(contrast %in% selcetedContrasts$contrast)

## get contrast coef
MAE <- MAE |> RFLOMICS::getContrastMatrix(SE.name = "RNAtest",   modelFormula = formulae[[1]], contrastList = selcetedContrasts) |>
              RFLOMICS::getContrastMatrix(SE.name = "protetest", modelFormula = formulae[[1]], contrastList = selcetedContrasts) |>
              RFLOMICS::getContrastMatrix(SE.name = "metatest",  modelFormula = formulae[[1]], contrastList = selcetedContrasts)
  




RFLOMICS::getExpressionContrast(MAE[["RNAtest"]], modelFormula = formulae[[1]])$simple
RFLOMICS::getExpressionContrastF(ExpDesign = MAE[["RNAtest"]]@colData, factorBio=c("temperature", "imbibition"), modelFormula = formulae[[1]])$simple



test_that("getExpressionContrast", {
  
  Contrasts.List.m <- RFLOMICS::getExpressionContrast(MAE, modelFormula = formulae[[1]])
  Contrasts.List.f <- RFLOMICS::getExpressionContrastF(ExpDesign = MAE@colData, factorBio=c("temperature", "imbibition"), modelFormula = formulae[[1]])
  
  expect_equal(Contrasts.List.m, Contrasts.List.f)
})


