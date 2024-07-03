library(RFLOMICS)
library(testthat)

# ---- Tests read_exp_design ----
### Comment when not run locally, need the txt files.
# 
# expect_no_error(readExpDesign(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt")))
expect_error(readExpDesign())

# ---- Tests read_omics_data ----

# expect_no_error(RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), 
#                                                       "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")))
# expect_error(RFLOMICS::readOmicsData(file = "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))

# ---- Tests RFLOMICS constructor ----

# load data
# load ecoseed data
data(ecoseed)

# create rflomicsMAE object with ecoseed data
MAE <- RFLOMICS::createRflomicsMAE(
  projectName = "Tests",
  omicsData   = list(ecoseed$RNAtest, ecoseed$metatest, ecoseed$protetest),
  omicsNames  = c("RNAtest", "metatest", "protetest"),
  omicsTypes  = c("RNAseq","metabolomics","proteomics"),
  ExpDesign   = ecoseed$design,
  factorRef   = ecoseed$factorRef)
names(MAE) <- c("RNAtest", "metatest", "protetest")

test_that("FlomicsMultiAssay.constructor fonction return MultiAssayExperiment object", {

  # test type if MAE class
  expect_true("MultiAssayExperiment" %in% is(MAE))
  
  # test type if element of MAE class
  for (SE in names(MAE)) {
    expect_true("SummarizedExperiment" %in% is(MAE[[SE]]))
  }
})

test_that("test if RflomicsMAE / RflomicsSE", {
    expect_true(is(MAE, "RflomicsMAE"))
})

## ---- getDatasetNames ----
test_that("test getDatasetNames", {
    expect_true(all(getDatasetNames(MAE) %in% c("RNAtest", "metatest", "protetest")))
})


# ---- Tests order of samples when all samples are the same ----

test_that("All omics data are ordred in same way", {
  
  expect_equal(colnames(MAE[[1]]), colnames(MAE[[2]]))
  expect_equal(colnames(MAE[[3]]), colnames(MAE[[2]]))
  expect_equal(colnames(MAE[[3]]), colnames(MAE[[2]]))
})

test_that("Test if samples in data matrix and rownames in design are ordered in same way", {
  
  expect_equal(colnames(MAE[[1]]), as.character(MAE[[1]]$samples))
  expect_equal(colnames(MAE[[2]]), as.character(MAE[[2]]$samples))
  expect_equal(colnames(MAE[[3]]), as.character(MAE[[3]]$samples))
})

test_that("Test if samples in data matrix and rownames in design are orderd in same way (with levels)", {
    
    expect_equal(colnames(MAE[[1]]), as.vector(MAE[[1]]$samples))
    expect_equal(colnames(MAE[[2]]), as.vector(MAE[[2]]$samples))
    expect_equal(colnames(MAE[[3]]), as.vector(MAE[[3]]$samples))
})

# ---- Tests order of samples when samples are not all the same ----

# load data
ExpDesign <- ecoseed$design
factorRef <- data.frame(factorName  = c("Repeat", "temperature" , "imbibition"),
                        factorRef   = c("rep1",   "Low",          "DS"),
                        factorType  = c("batch",  "Bio",          "Bio"),
                        factorLevels= c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI"))

omicsData <- list(
    ecoseed$RNAtest,
    ecoseed$metatest,
    ecoseed$protetest)

omicsData[[1]] <- omicsData[[1]][,-5]
omicsData[[2]] <- omicsData[[2]][,-10]
ExpDesign      <- ExpDesign[-20, ]

MAE <- RFLOMICS::createRflomicsMAE(projectName = "Tests", 
                                   omicsData   = omicsData,
                                   omicsNames  = c("RNAtest", "metatest", "protetest"),
                                   omicsTypes  = c("RNAseq","metabolomics","proteomics"),
                                   ExpDesign   = ExpDesign,
                                   factorRef   = factorRef)

test_that("Test if samples in data matrix and rownames in design are orderd in same way (with levels)", {
  
  expect_equal(colnames(MAE[[1]]), as.vector(MAE[[1]]$samples))
  expect_equal(colnames(MAE[[2]]), as.vector(MAE[[2]]$samples))
  expect_equal(colnames(MAE[[3]]), as.vector(MAE[[3]]$samples))
})

