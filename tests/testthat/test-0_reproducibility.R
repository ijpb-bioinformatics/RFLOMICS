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

MAE <- RFLOMICS::getExpressionContrast(object = MAE, model.formula = formulae[[1]])

## choice of hypothesis
selcetedContrasts <- c(MAE@metadata$design@Contrasts.List$simple$contrast[1:3],
                       MAE@metadata$design@Contrasts.List$averaged$contrast[1:3],
                       MAE@metadata$design@Contrasts.List$interaction$contrast[1:3])

MAE <- RFLOMICS::getContrastMatrix(object = MAE, contrastList = selcetedContrasts)
contrastsDF <- RFLOMICS::getSelectedContrasts(MAE)


## Process data RNAseq
MAE <- MAE |> RFLOMICS::runDataProcessing(SE.name = "RNAtest"  , samples=colnames(MAE[["RNAtest.raw"]])[-1], lowCountFiltering_strategy="NbReplicates", lowCountFiltering_CPM_Cutoff=1, normalisation_method="TMM") |>
              RFLOMICS::runDataProcessing(SE.name = "protetest", samples=NULL, normalisation_method="none", transformation_method="none") |>
              RFLOMICS::runDataProcessing(SE.name = "metatest" , samples=NULL, normalisation_method=NULL, transformation_method="log2")   

## diff analysis
MAE <- MAE |> RFLOMICS::RunDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.method="BH", DiffAnalysisMethod = "edgeRglmfit", Adj.pvalue.cutoff = 0.05, logFC.cutoff = 0) |>
              RFLOMICS::RunDiffAnalysis(SE.name = "protetest", Adj.pvalue.method="BH", DiffAnalysisMethod = "limmalmFit",  Adj.pvalue.cutoff = 0.05, logFC.cutoff = 0) |>
              RFLOMICS::RunDiffAnalysis(SE.name = "metatest",  Adj.pvalue.method="BH", DiffAnalysisMethod = "limmalmFit",  Adj.pvalue.cutoff = 0.05, logFC.cutoff = 0)

## co expression
MAE <- MAE |> RFLOMICS::runCoExpression(SE.name = "RNAtest",   nameList = "(temperatureMedium - temperatureLow) in imbibitionDS", K = 2:10, replicates = 5, merge = "union", model = "normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin", normFactors = "TMM") |>
              RFLOMICS::runCoExpression(SE.name = "protetest", nameList = "(temperatureMedium - temperatureLow) in imbibitionDS", K = 2:10, replicates = 5, merge = "union", model = "normal") |>
              RFLOMICS::runCoExpression(SE.name = "metatest",  nameList = "(temperatureMedium - temperatureLow) in imbibitionDS", K = 2:10, replicates = 5, merge = "union", model = "normal")

## Enrichment
MAE <- MAE |> RFLOMICS::runAnnotationEnrichment(nameList = "H1", SE.name = "RNAtest", dom.select = "GO", Domain = c("BP", "MF", "CC"), list_args = list(OrgDb = "org.At.tair.db", keyType = "TAIR", pvalueCutoff = 0.05))
              # RFLOMICS::runAnnotationEnrichment(nameList = "H1", SE.name = "RNAtest", dom.select = "KEGG", list_args = list(..., pvalueCutoff = 0.05))

# ---------------- TESTS ---------------

# comparaison avec un autre objet
test_that("colData", {
  
  colData <- data.frame(Repeat      = rep(c("rep1", "rep2", "rep3"), 9),
                        temperature = c(rep("Elevated", 9), rep("Low", 9), rep("Medium", 9)),
                        imbibition  = rep(c(rep("DS",3), rep("EI", 3), rep("LI", 3)), 3))
  
  rownames(colData) <- c("Elevated_DS_1", "Elevated_DS_2", "Elevated_DS_3", 
                         "Elevated_EI_1", "Elevated_EI_2", "Elevated_EI_3", 
                         "Elevated_LI_1", "Elevated_LI_2", "Elevated_LI_3", 
                         "Low_DS_1", "Low_DS_2", "Low_DS_3", 
                         "Low_EI_1", "Low_EI_2", "Low_EI_3", 
                         "Low_LI_1", "Low_LI_2", "Low_LI_3", 
                         "Medium_DS_1", "Medium_DS_2", "Medium_DS_3", 
                         "Medium_EI_1", "Medium_EI_2", "Medium_EI_3",
                         "Medium_LI_1", "Medium_LI_2", "Medium_LI_3")
  
  colData$Repeat      <- factor(colData$Repeat,      levels =c("rep1", "rep2", "rep3"))
  colData$temperature <- factor(colData$temperature, levels =c("Low", "Medium", "Elevated"))
  colData$imbibition  <- factor(colData$imbibition,  levels =c("DS", "EI", "LI"))
  
  colData$Repeat      <- relevel(as.factor(colData$Repeat),      ref="rep1")
  colData$temperature <- relevel(as.factor(colData$temperature), ref="Low")
  colData$imbibition  <- relevel(as.factor(colData$imbibition),  ref="DS")
  
  expect_equal(as.data.frame(MAE@colData), colData)
  
})

test_that("contrast", {
  
  Contrasts.names <- c("(temperatureMedium_imbibitionDS - temperatureLow_imbibitionDS)",                                                                                                                                                     
                       "(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",                                                                                                                                                   
                       "(temperatureElevated_imbibitionDS - temperatureMedium_imbibitionDS)",                                                                                                                                              
                       "((temperatureMedium_imbibitionDS - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureLow_imbibitionEI) + (temperatureMedium_imbibitionLI - temperatureLow_imbibitionLI))/3",          
                       "((temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) + (temperatureElevated_imbibitionLI - temperatureLow_imbibitionLI))/3",
                       "((temperatureElevated_imbibitionDS - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureMedium_imbibitionEI) + (temperatureElevated_imbibitionLI - temperatureMedium_imbibitionLI))/3",
                       "((temperatureMedium_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureMedium_imbibitionDS - temperatureLow_imbibitionDS))",
                       "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))",
                       "((temperatureElevated_imbibitionEI - temperatureMedium_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureMedium_imbibitionDS))")
  
  Contrasts.Coeff <- rbind(
    c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0, 0, 0.333333333333333, 0, 0.333333333333333, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0.333333333333333, 0, 0.333333333333333),
    c(0, 0, 0, -1, 1, 0, 0, -0.333333333333333, 0.333333333333333, -0.333333333333333, 0.333333333333333),
    c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
  ) %>% as.data.frame()
  
  names(Contrasts.Coeff) <- c("Intercept", "Repeatrep2", "Repeatrep3", 
                                "temperatureMedium", "temperatureElevated", 
                                "imbibitionEI", "imbibitionLI", 
                                "temperatureMedium:imbibitionEI", "temperatureElevated:imbibitionEI",
                                "temperatureMedium:imbibitionLI", "temperatureElevated:imbibitionLI")
  
  row.names(Contrasts.Coeff) <- Contrasts.names
  
  expect_equal(MAE@metadata$design@Contrasts.Coeff, Contrasts.Coeff)
  
})


test_that("data processing", {
  
  ## low count Filtering
  # remplacer par le bon getter
  nb_lowGene <- length(MAE[["RNAtest"]]@metadata$DataProcessing$Filtering$results$filteredFeatures)
  expect_equal(nb_lowGene, 6725)
  expect_equal(MAE[["protetest"]]@metadata$DataProcessing$Filtering, NULL)
  expect_equal(MAE[["metatest"]]@metadata$DataProcessing$Filtering, NULL)
  
  ## sample filtering
  expect_equal(MAE[["RNAtest"]]@metadata$DataProcessing$filteredSamples, "Elevated_DS_1")
  
  ## transformation
  
  ## Normalisation
  
})


test_that("diff analysis", {
  
  ## RNAtest
  stats.RNAseq <- data.frame(
    All = c(2589, 5925, 1579, 5535, 9735, 4901, 1473, 5110,  817),
    Up  = c(1323, 3017,  790, 2729, 4778, 2324,  568, 2453,  344),
    Down= c(1266, 2908,  789, 2806, 4957, 2577,  905, 2657,  473)
  ) %>% as.matrix()
  rownames(stats.RNAseq) <- c(
    "(temperatureMedium - temperatureLow) in imbibitionDS",                                                              
    "(temperatureElevated - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureMedium) in imbibitionDS",
    "(temperatureMedium - temperatureLow) in mean",
    "(temperatureElevated - temperatureLow) in mean",
    "(temperatureElevated - temperatureMedium) in mean",
    "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS"
  )
  
  expect_equal(MAE[["RNAtest"]]@metadata$DiffExpAnal$stats, stats.RNAseq)
  
  ## protetest
  stats.prot <- data.frame(
    All = c(36, 132,  53,  99, 238, 130,   0,   1,   0),
    Up  = c(24,  67,  21,  59, 124,  47,   0,   1,   0),
    Down= c(12,  65,  32,  40, 114,  83,   0,   0,   0)
  ) %>% as.matrix()
  rownames(stats.prot) <- c(
    "(temperatureMedium - temperatureLow) in imbibitionDS",                                                              
    "(temperatureElevated - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureMedium) in imbibitionDS",
    "(temperatureMedium - temperatureLow) in mean",
    "(temperatureElevated - temperatureLow) in mean",
    "(temperatureElevated - temperatureMedium) in mean",
    "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS"
  )
  
  expect_equal(MAE[["protetest"]]@metadata$DiffExpAnal$stats, stats.prot)
  
  ## metatest
  stats.met <- data.frame(
    All = c(54, 67, 45, 73, 87, 69, 12, 18,  5),
    Up  = c(41, 44, 23, 62, 57, 33,  6,  8,  1),
    Down= c(13, 23, 22, 11, 30, 36,  6, 10,  4)
  ) %>% as.matrix()
  rownames(stats.met) <- c(
    "(temperatureMedium - temperatureLow) in imbibitionDS",                                                              
    "(temperatureElevated - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureMedium) in imbibitionDS",
    "(temperatureMedium - temperatureLow) in mean",
    "(temperatureElevated - temperatureLow) in mean",
    "(temperatureElevated - temperatureMedium) in mean",
    "(temperatureMedium - temperatureLow) in imbibitionEI - (temperatureMedium - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureLow) in imbibitionEI - (temperatureElevated - temperatureLow) in imbibitionDS",
    "(temperatureElevated - temperatureMedium) in imbibitionEI - (temperatureElevated - temperatureMedium) in imbibitionDS"
  )
  
  expect_equal(MAE[["metatest"]]@metadata$DiffExpAnal$stats, stats.met)
})

test_that("diff analysis", {
  
  expect_equal(MAE[["RNAtest"]]@metadata$CoExpAnal$cluster.nb[[1]],   4)
  expect_equal(MAE[["protetest"]]@metadata$CoExpAnal$cluster.nb[[1]], 3)
  expect_equal(MAE[["metatest"]]@metadata$CoExpAnal$cluster.nb[[1]],  6)
})



