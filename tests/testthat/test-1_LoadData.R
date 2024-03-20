library(RFLOMICS)
library(testthat)

# ---- Creating files to test ----

# # Create files with missing values, special char, etc.
# condEcoseed <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))
# set.seed(2000)
# modDelete1 <- sample(x = 1:nrow(condEcoseed), size = 3)
# set.seed(2001)
# modDelete2 <- sample(x = 1:nrow(condEcoseed), size = 3)
# # set.seed(2002)
# # modDelete3 <- sample(x = 1:nrow(condEcoseed), size = 3)
# 
# condEcoseedNA <- condEcoseed
# condEcoseedNA$Repeat[modDelete1] <- NA
# condEcoseedNA <- condEcoseedNA %>% tibble::rownames_to_column(var = "sample") %>% dplyr::relocate(sample)
# write.table(condEcoseedNA, file = "inst/testFiles/EcoseedCondNA.txt", sep = "\t", row.names = FALSE)
# 
# condEcoseedspecChar <- condEcoseed
# levels(condEcoseedspecChar$temperature) <- c("Ele_vated", "L&w", "Médium")
# condEcoseedspecChar <- condEcoseedspecChar %>% tibble::rownames_to_column(var = "sample") %>% dplyr::relocate(sample)
# write.table(condEcoseedspecChar, file = "inst/testFiles/EcoseedCondSpecChar.txt", sep = "\t", row.names = FALSE)

## NEED 
# TODO
# Create file with duplicated colnames
# Create file with duplicated modalities
# Same work with data

# ---- Tests read_exp_design ----

expect_no_error(read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt")))
expect_error(read_exp_design())

# Missing values
#condEcoseedTest <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/testFiles/EcoseedCondNA.txt"))
# TODO verifier le constructeur ET lancer une analyse diff avec des NA

# Special characters
#condEcoseedTest <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/testFiles/EcoseedCondSpecChar.txt"))
# TODO verifier le constructeur ET lancer une analyse diff avec des characteres speciaux

# ---- Tests read_omics_data ----

expect_no_error(RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")))
# expect_error(RFLOMICS::readOmicsData(file = "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))

# ---- Tests RFLOMICS constructor ----

# load data
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
    
    expect_equal(as.factor(colnames(MAE[[1]])), MAE[[1]]$samples)
    expect_equal(as.factor(colnames(MAE[[2]])), MAE[[2]]$samples)
    expect_equal(as.factor(colnames(MAE[[3]])), MAE[[3]]$samples)
})

# ---- Tests order of samples when samples are not all the same ----

# load data
ExpDesign <- RFLOMICS::readExpDesign(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))
factorRef <- data.frame(factorName  = c("Repeat", "temperature" , "imbibition"),
                        factorRef   = c("rep1",   "Low",          "DS"),
                        factorType  = c("batch",  "Bio",          "Bio"),
                        factorLevels= c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI"))

omicsData <- list(
    RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
    RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
    RFLOMICS::readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")))

omicsData[[1]] <- omicsData[[1]][,-5]
omicsData[[2]] <- omicsData[[2]][,-10]
ExpDesign      <- ExpDesign[-20, ]

MAE <- RFLOMICS::createRflomicsMAE(projectName = "Tests", 
                                   omicsData   = omicsData,
                                   omicsNames  = c("RNAtest", "metatest", "protetest"),
                                   omicsTypes  = c("RNAseq","metabolomics","proteomics"),
                                   ExpDesign   = ExpDesign,
                                   factorRef   = factorRef)

test_that("Test if samples in data matrix and rownames in design are orderd in same way", {
  
  expect_equal(as.factor(colnames(MAE[[1]])), MAE[[1]]$samples)
  expect_equal(as.factor(colnames(MAE[[2]])), MAE[[2]]$samples)
  expect_equal(as.factor(colnames(MAE[[3]])), MAE[[3]]$samples)
})

