library(testthat)
library(RFLOMICS)

# These tests will allow us to see if changes in the code or package versions affect the expected results. 

# ---- Create RflomicsMAE object ----

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
                                   omicsNames  = c("RNAtest", "metatest", "protetest"),
                                   omicsTypes  = c("RNAseq","metabolomics","proteomics"),
                                   ExpDesign   = ExpDesign,
                                   factorRef   = factorRef)

SE <- getRflomicsSE(MAE, "RNAtest.raw")

# ---- check design ----
plotConditionsOverview(MAE)

test_that("test if RflomicsMAE / RflomicsSE", {
  
  check <- checkExpDesignCompleteness(SE)
  
  expect_false(check$error)
  expect_equal(check$message, "The experimental design is complete and balanced.")
})

## choice of model formulae
formulae <- RFLOMICS::generateModelFormulae(MAE)

MAE <- RFLOMICS::setModelFormula(MAE, modelFormula = formulae[[1]])

## list of all hypothesis grouped by contarst type (run on MAE or SE)
Contrasts.List <- RFLOMICS::generateExpressionContrast(MAE)

# utilisation en dehors de MAE ou SE

## choice of hypothesis
selectedContrasts <- rbind(Contrasts.List$simple[1:3,],
                           Contrasts.List$averaged[1:3,],
                           Contrasts.List$interaction[1:3,])

selectedContrasts <- Contrasts.List$simple[1:4,]

MAE <- setSelectedContrasts(MAE, contrastList = selectedContrasts)

## data processing
sampleToKeep <- colnames(MAE[["RNAtest.raw"]])[-1]
MAE <- MAE |> RFLOMICS::runDataProcessing(SE.name = "RNAtest"  , samples=sampleToKeep, lowCountFiltering_strategy="NbReplicates", lowCountFiltering_CPM_Cutoff=1, normMethod="TMM",transformMethod = "none") |>
  RFLOMICS::runDataProcessing(SE.name = "protetest", samples=NULL, normMethod="none", transformMethod="none") |>
  RFLOMICS::runDataProcessing(SE.name = "metatest" , samples=NULL, normMethod=NULL, transformMethod="log2")

## diff analysis

MAE <- MAE |> RFLOMICS::runDiffAnalysis(SE.name = "RNAtest", contrastList = selectedContrasts, p.adj.method="BH", method = "edgeRglmfit", p.adj.cutoff = 0.05, logFC.cutoff = 0) |>
  RFLOMICS::runDiffAnalysis(SE.name = "protetest", contrastList = selectedContrasts, p.adj.method="BH", method = "limmalmFit",  p.adj.cutoff = 0.05, logFC.cutoff = 0) |>
  RFLOMICS::runDiffAnalysis(SE.name = "metatest", contrastList = selectedContrasts, p.adj.method="BH", method = "limmalmFit",  p.adj.cutoff = 0.05, logFC.cutoff = 0)

## co expression
MAE <- MAE |> RFLOMICS::runCoExpression(SE.name = "RNAtest",   contrastNames = "H1" , K = 2:10, replicates = 5, merge = "union", model = "normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin", normFactors = "TMM") |>
  RFLOMICS::runCoExpression(SE.name = "protetest", contrastNames = "H1" , K = 2:10, replicates = 5, merge = "union", model = "normal") |>
  RFLOMICS::runCoExpression(SE.name = "metatest",  contrastNames = "H1" , K = 2:10, replicates = 5, merge = "union", model = "normal")


## Enrichment
MAE <- MAE |> RFLOMICS::runAnnotationEnrichment(SE.name = "RNAtest", database = "GO", domain = c("BP", "MF", "CC"), list_args = list(OrgDb = "org.At.tair.db", keyType = "TAIR", pvalueCutoff = 0.05)) |>
  RFLOMICS::runAnnotationEnrichment(SE.name = "protetest", database = "GO", domain = c("BP", "MF", "CC"), list_args = list(OrgDb = "org.At.tair.db", keyType = "TAIR", pvalueCutoff = 0.05))
 #RFLOMICS::runAnnotationEnrichment(nameList = "H1", SE.name = "RNAtest", database = "KEGG", list_args = list(..., pvalueCutoff = 0.05))

## runReport
# generateReport(object=MAE, tmpDir ="/Users/dcharif/Desktop/", fileName = "projet",export=TRUE)

# ---- test accessors ----
## ---- Rflomics class ----
test_that("test if RflomicsMAE / RflomicsSE", {
  
  expect_true(is(MAE, "RflomicsMAE"))
})

## ---- getDatasetNames ----
test_that("test getDatasetNames", {
  
  expect_true(all(getDatasetNames(MAE) %in% c("RNAtest", "metatest", "protetest")))
})

## ---- getOmicsTypes ----
test_that("test getOmicsTypes", {
  
  exp.res <- c("RNAseq", "metabolomics", "proteomics")
  names(exp.res) <- c("RNAtest", "metatest", "protetest")
  expect_equal(getOmicsTypes(MAE), exp.res)
})

## ---- colData ----
test_that("colData", {
  
  Repeat      <- factor(rep(c("rep1", "rep2", "rep3"), 9), 
                        levels =c("rep1", "rep2", "rep3"))
  temperature <- factor(c(rep("Elevated", 9), rep("Low", 9), rep("Medium", 9)), 
                        levels =c("Low", "Medium", "Elevated"))
  imbibition  <- factor(rep(c(rep("DS",3), rep("EI", 3), rep("LI", 3)), 3),  
                        levels =c("DS", "EI", "LI"))
  Repeat      <- relevel(as.factor(Repeat),      ref="rep1")
  temperature <- relevel(as.factor(temperature), ref="Low")
  imbibition  <- relevel(as.factor(imbibition),  ref="DS")
  
  # groups
  groups.level <- character(0)
  # Boucles imbriquées pour créer les combinaisons
  for (v1 in c("Low", "Medium", "Elevated")) {
    for (v2 in c("DS", "EI", "LI")) {
      groups.level <- c(groups.level, paste(v1, v2, sep = "_"))
    }
  }
  groups <- factor(paste(temperature, imbibition, sep = "_"),
                   levels =groups.level)
  
  # samples
  samples.level <- character(0)
  # Boucles imbriquées pour créer les combinaisons
  for (v1 in groups.level) {
    for (v2 in 1:3) {
      samples.level <- c(samples.level, paste(v1, v2, sep = "_"))
    }
  }
  
  samples <- factor(paste(temperature, imbibition, sub("rep", "", Repeat), sep = "_"),
                    levels = samples.level)
  
  
  colData <- data.frame(Repeat      = Repeat, 
                        groups      = groups,
                        temperature = temperature,
                        imbibition  = imbibition,
                        samples     = samples)
  
  rownames(colData) <- samples
  
  expect_equal(as.data.frame(getDesignMat(MAE)), colData)
})


## ---- getRflomicsSE ----
test_that("test getRflomicsSE", {
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_true(is(SE, "RflomicsSE"))
  
  expect_null(getRflomicsSE(MAE))
  expect_null(getRflomicsSE(MAE, "toto"))
  
  expect_equal(getDatasetNames(SE), "RNAtest")
})


## ---- get design factors ----
test_that("test getFactorNames", {
  
  expect_equal(getFactorNames(MAE), c("Repeat", "temperature", "imbibition" ))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getFactorNames(SE), c("Repeat", "temperature", "imbibition" ))
})

test_that("test getFactorModalities", {
  
  expect_equal(getFactorModalities(MAE, "imbibition"), c("DS", "EI", "LI"))
  expect_equal(getFactorModalities(SE, "imbibition"),  c("DS", "EI", "LI"))
})



test_that("test getFactorTypes", {
  
  vec <- c("batch", "Bio", "Bio")
  names(vec) <- c("Repeat", "temperature", "imbibition" )
  expect_equal(getFactorTypes(MAE), vec)
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getFactorTypes(SE), vec)
})

test_that("test getBioFactors", {
  
  expect_equal(getBioFactors(MAE), c("temperature", "imbibition"))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getBioFactors(SE), c("temperature", "imbibition"))
})

test_that("test getBatchFactors", {
  
  expect_equal(getBatchFactors(MAE), c("Repeat"))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_equal(getBatchFactors(SE), c("Repeat"))
})

test_that("test getMetaFactors", {
  
  expect_null(getMetaFactors(MAE))
  
  SE <- getRflomicsSE(MAE, "RNAtest")
  expect_null(getMetaFactors(SE))
})

## ---- sub set of RflomicsMAE ----
test_that("test subRflomicsMAE", {
  
  miniMAE <- MAE[,, c("RNAtest", "metatest")]
  miniMAE.rf <- subRflomicsMAE(MAE, c("RNAtest", "metatest"))
  expect_equal(miniMAE.rf, miniMAE)
  
})

# ---- contrasts ----

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
  
  expect_equal(MAE[["RNAtest"]]@metadata$design$Contrasts.Coeff, Contrasts.Coeff)
  expect_equal(MAE[["protetest"]]@metadata$design$Contrasts.Coeff, Contrasts.Coeff)
  expect_equal(MAE[["metatest"]]@metadata$design$Contrasts.Coeff, Contrasts.Coeff)
  
  expect_equal(CheckExpDesignCompleteness(MAE[["RNAtest"]])$messages, "The experimental design is complete but not balanced.")
  
  
})

# ---- processing ----

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

# ---- diff analysis ----
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

# ---- co-expression ----
test_that("coseq analysis", {
  
  expect_equal(MAE[["RNAtest"]]@metadata$CoExpAnal$cluster.nb[[1]],   4)
  expect_equal(MAE[["protetest"]]@metadata$CoExpAnal$cluster.nb[[1]], 3)
  expect_equal(MAE[["metatest"]]@metadata$CoExpAnal$cluster.nb[[1]],  6)
})

# ---- enrichment ----
test_that("Enrichment analysis", {
  
  expect_true(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$GO$summary["BP"] == 1)
  expect_true(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$GO$summary["MF"] == 6)
  expect_true(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$GO$summary["CC"] == 11)
})



