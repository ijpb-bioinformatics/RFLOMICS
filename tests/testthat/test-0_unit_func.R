library(RFLOMICS)
library(testthat)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for coseq analysis : ----
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

MAE <- MAE |>
  TransformData(     SE.name = "metatest",  transform_method = "log2")          |>
  RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")            |>
  RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")                 |>
  RunNormalization(  SE.name = "protetest", NormMethod = "median")              |>
  FilterLowAbundance(SE.name = "RNAtest")                                       |>
  RunDiffAnalysis(   SE.name = "metatest",  DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "protetest", DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "RNAtest",   DiffAnalysisMethod = "edgeRglmfit") |>
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





