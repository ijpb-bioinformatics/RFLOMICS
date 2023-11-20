library(RFLOMICS)
library(testthat)

# ---- Construction of objects for the tests ----

# ---- Construction MAE RFLOMICS ready for differential analysis : ----
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
names(MAE) <- c("RNAtest", "metatest", "protetest")

formulae <- RFLOMICS::GetModelFormulae(MAE = MAE) 

contrastList <- RFLOMICS::getExpressionContrast(object = MAE, model.formula = formulae[[1]]) |> purrr::reduce(rbind) |>
  dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
                                "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))


MAE <- MAE |>
  TransformData(     SE.name = "metatest",  transformMethod = "log2")          |>
  RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")            |>
  RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")                 |>
  RunNormalization(  SE.name = "protetest", NormMethod = "median")              |>
  FilterLowAbundance(SE.name = "RNAtest")                                       |>
  RunDiffAnalysis(   SE.name = "metatest",  DiffAnalysisMethod = "limmalmFit", contrastList = contrastList, modelFormula = formulae[[1]])  |>
  RunDiffAnalysis(   SE.name = "protetest", DiffAnalysisMethod = "limmalmFit", contrastList = contrastList, modelFormula = formulae[[1]])  |>
  RunDiffAnalysis(   SE.name = "RNAtest",   DiffAnalysisMethod = "edgeRglmfit", contrastList = contrastList, modelFormula = formulae[[1]]) |>
  FilterDiffAnalysis(SE.name = "RNAtest",  Adj.pvalue.cutoff = 0.05, logFC.cutoff = 2)

# ----- Test equivalence geneList selection : ----


test_that("runCoExpression - Former code and new code equivalence", {
  # Former code for genelist inside runCoExpression
  # Test if equivalent to new code
  
  object <- MAE[["RNAtest"]]
  nameList <- colnames(object@metadata$DiffExpAnal[["mergeDEF"]])[-1]
  merge <- "union"
  
  geneList <- dplyr::select(object@metadata$DiffExpAnal[["mergeDEF"]], DEF, tidyselect::all_of(nameList)) %>% 
    dplyr::mutate(intersection = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) == length(nameList), "YES", "NO"), 
                  union = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) != 0 , "YES", "NO")) %>% 
    dplyr::filter(union != "NO", get(merge) == "YES") 
  geneList <- geneList$DEF
  
  geneList2 <- opDEList(object = MAE[["RNAtest"]], contrasts = nameList, operation = merge)
  
  expect_equal(geneList, geneList2)
  
  # Change merge type
  merge <- "intersection"
  
  geneList <- dplyr::select(object@metadata$DiffExpAnal[["mergeDEF"]], DEF, tidyselect::all_of(nameList)) %>% 
    dplyr::mutate(intersection = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) == length(nameList), "YES", "NO"), 
                  union = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) != 0 , "YES", "NO")) %>% 
    dplyr::filter(union != "NO", get(merge) == "YES") 
  geneList <- geneList$DEF
  
  geneList2 <- opDEList(object = MAE[["RNAtest"]], contrasts = nameList, operation = merge)
  
  expect_equal(geneList, geneList2)
  
  # change contrasts
  nameList <- c("H1")
  
  geneList <- dplyr::select(object@metadata$DiffExpAnal[["mergeDEF"]], DEF, tidyselect::all_of(nameList)) %>% 
    dplyr::mutate(intersection = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) == length(nameList), "YES", "NO"), 
                  union = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) != 0 , "YES", "NO")) %>% 
    dplyr::filter(union != "NO", get(merge) == "YES") 
  geneList <- geneList$DEF
  
  geneList2 <- opDEList(object = MAE[["RNAtest"]], contrasts = nameList, operation = merge)
  
  expect_equal(geneList, geneList2)

})





