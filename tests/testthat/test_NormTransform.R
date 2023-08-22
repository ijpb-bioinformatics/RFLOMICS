library(testthat)

##### Checks on transformation and normalization of data through command line functions

# ---- Construction of objects for the tests ----

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

protMat <- RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt"))
rnaMat <- RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))

######################################-
########### PROTEO/METABO ############
######################################-

# ---- PROTEO - Check if yielding the expected result ----
# Test all possibles combination (so far) of data transformation and data normalization for proteomics data
# also true for metabolomics data.
# Also test results from PCA (raw and norm data) 

test_that("Transformation and normalisation combination - proteomics", {
  
  casesMat <- expand.grid(c("none", "log2", "log10", "log1p", "squareroot"), c("none", "median", "totalSum"))
  colnames(casesMat) <- c("Trans", "Norm")
  
  res_equal <- lapply(1:nrow(casesMat), FUN = function(i){
    
    case_vect <- casesMat[i,]
    
    # matrix version
    protMattransnorm <- protMat
    
    pca.raw <- FactoMineR::PCA(t(protMattransnorm), ncp = 5, graph = FALSE)
    
    expect_identical(pca.raw$eig, MAE[["protetest"]]@metadata$PCAlist$raw$eig)
    expect_identical(pca.raw$svd, MAE[["protetest"]]@metadata$PCAlist$raw$svd)
    expect_identical(pca.raw$ind, MAE[["protetest"]]@metadata$PCAlist$raw$ind)
    expect_identical(pca.raw$var, MAE[["protetest"]]@metadata$PCAlist$raw$var)
    # the call is obligatory different between the two
    
    protMattransnorm <- switch(as.character(case_vect[[1]]),
                               "none"       = protMattransnorm,
                               "log2"       = log2(protMattransnorm + 1),
                               "log10"      = log10(protMattransnorm + 1),
                               "log1p"      = log1p(protMattransnorm),
                               "squareroot" = sqrt(protMattransnorm)
    )
    
    protMattransnorm <- switch(as.character(case_vect[[2]]),
                               "none"     = protMattransnorm,
                               "median"   = apply(protMattransnorm, 2, FUN = function(vect) vect - median(vect)),
                               "totalSum" = apply(protMattransnorm, 2, FUN = function(vect) vect/sum(vect^2))
    )
    
    pca.norm <- FactoMineR::PCA(t(protMattransnorm), ncp = 5, graph = FALSE)
    
    # RFLOMICS version
    MAE2 <- MAE
    MAE2 <- RFLOMICS::TransformData(MAE2, SE = "protetest", transform_method = as.character(case_vect[[1]]))
    MAE2 <- RFLOMICS::RunNormalization(MAE2, SE.name = "protetest", NormMethod = as.character(case_vect[[2]]))
    
    MAE2 <- RFLOMICS::RunPCA(MAE2, SE = "protetest")
    expect_identical(pca.norm$eig, MAE2[["protetest"]]@metadata$PCAlist$norm$eig)
    expect_identical(pca.norm$svd, MAE2[["protetest"]]@metadata$PCAlist$norm$svd)
    expect_identical(pca.norm$ind, MAE2[["protetest"]]@metadata$PCAlist$norm$ind)
    expect_identical(pca.norm$var, MAE2[["protetest"]]@metadata$PCAlist$norm$var)
    
    MAE2[["protetest"]] <- RFLOMICS:::checkTransNorm(MAE2[["protetest"]]) 
    
    expect_identical(SummarizedExperiment::assay(MAE2[["protetest"]]), as.matrix(protMattransnorm))
    
  })
  
})

# ---- PROTEO - Check if transformation working without arguments ----
test_that("Transformation - no method - no modification", {
  
  MAE2 <- RFLOMICS::TransformData(MAE, SE = "protetest")
  
  expect_message(RFLOMICS::TransformData(MAE, SE = "protetest"))
  expect(MAE2[["protetest"]]@metadata$transform$transform_method == "log2", failure_message = "The transformation is not log2 by default.")
  expect_equal(SummarizedExperiment::assay(MAE[["protetest"]]), SummarizedExperiment::assay(MAE2[["protetest"]]))
})

test_that("Transformation - no method - modification", {
  
  MAE2 <- RFLOMICS::TransformData(MAE, SE = "protetest", modify_assay = TRUE)
  
  # Message if argument method is forgotten
  expect_message(RFLOMICS::TransformData(MAE, SE = "protetest"))
  
  # Default transformation must be log2 + 1 for proteomics data
  expect_equal(log2(SummarizedExperiment::assay(MAE[["protetest"]]) + 1), SummarizedExperiment::assay(MAE2[["protetest"]]))
  expect(MAE2[["protetest"]]@metadata$transform$transform_method == "log2", failure_message = "The transformation is not log2 by default.")
  
  # Modification of the assay: expect transformed metadata to be TRUE 
  expect(MAE2[["protetest"]]@metadata$transform$transformed, failure_message = "The assay was not transformed.")
})

######################################-
########### RNASEQ ####################
######################################-

# ---- RNAseq - Check if yielding the expected result ----

test_that("RNAseq - none + TMM + log2", {
  
  
  # SumarizedExperiment::assay(MAE[["RNAtest"]]) 
  
  # matrix version (filtering is done in the FlomicsMultiAssay constructor)
  rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))

  pca.raw <- FactoMineR::PCA(t(log2(rnaSeqMat + 1)), ncp = 5, graph = FALSE)
  
  expect_identical(pca.raw$eig, MAE[["RNAtest"]]@metadata$PCAlist$raw$eig)
  expect_identical(pca.raw$svd, MAE[["RNAtest"]]@metadata$PCAlist$raw$svd)
  expect_identical(pca.raw$ind, MAE[["RNAtest"]]@metadata$PCAlist$raw$ind)
  expect_identical(pca.raw$var, MAE[["RNAtest"]]@metadata$PCAlist$raw$var)
  # the call is obligatory different between the two
  
  # Manually transforming and normalizing rnaSeqMat
  normFactors <- edgeR::calcNormFactors(rnaSeqMat, method = "TMM")
  libSize <-  colSums(rnaSeqMat)
  
  tnDat <- scale(rnaSeqMat + 1, center = FALSE, scale = normFactors*libSize)
  
  pca.norm <- FactoMineR::PCA(t(log2(tnDat)), ncp = 5, graph = FALSE)
  
  # RFLOMICS version
  MAE2 <- MAE
  MAE2 <- RFLOMICS::TransformData(MAE2, SE = "RNAtest", transform_method = "none")
  MAE2 <- RFLOMICS::RunNormalization(MAE2, SE.name = "RNAtest", NormMethod = "TMM")
  
  MAE2 <- RFLOMICS::RunPCA(MAE2, SE = "RNAtest")
  expect_identical(pca.norm$eig, MAE2[["RNAtest"]]@metadata$PCAlist$norm$eig)
  expect_identical(pca.norm$svd, MAE2[["RNAtest"]]@metadata$PCAlist$norm$svd)
  expect_identical(pca.norm$ind, MAE2[["RNAtest"]]@metadata$PCAlist$norm$ind)
  expect_identical(pca.norm$var, MAE2[["RNAtest"]]@metadata$PCAlist$norm$var)
  
})

# ---- RNAseq - behaviour of transforData and RunNormalization ---

test_that("RNAseq - correct behaviour of normalization and transformation",  {
  
  MAE2 <- MAE
  MAE2 <- RFLOMICS::TransformData(MAE2, SE = "RNAtest")
  expect(MAE2[["RNAtest"]]@metadata$transform$transform_method == "none", failure_message = "TransformData does not put 'none' as default for RNAseq data")
  expect(!MAE2[["RNAtest"]]@metadata$transform$transformed, failure_message = "TransformData transformed the data when not asked to")

  MAE2 <- RFLOMICS::RunNormalization(MAE2, SE = "RNAtest")
  expect(MAE2[["RNAtest"]]@metadata$Normalization$methode == "TMM", failure_message = "Normalization is not defaulted to TMM for RNAseq data.")
  expect(!MAE2[["RNAtest"]]@metadata$Normalization$normalized, failure_message = "RunNormalization transformed the data when not asked to")
  
  MAE2 <- RFLOMICS::TransformData(MAE2, SE = "RNAtest", transform_method = "log2") # shouldn't be at all
  MAE2 <- RFLOMICS::RunNormalization(MAE2, SE = "RNAtest", NormMethod = "median") # shouldn't be at all.
  
  expect(MAE2[["RNAtest"]]@metadata$transform$transform_method == "none", failure_message = "TransformData did not replace by 'none' for RNAseq data ")
  expect(!MAE2[["RNAtest"]]@metadata$transform$transformed, failure_message = "TransformData transformed the data when not asked to")
  
  expect(MAE2[["RNAtest"]]@metadata$Normalization$methode == "TMM", failure_message = "Normalization did not replace by TMM for RNAseq data.")
  expect(!MAE2[["RNAtest"]]@metadata$Normalization$normalized, failure_message = "RunNormalization transformed the data when not asked to")

}) 


test_that("RNAseq - correct behaviour of checkTransNorm",  {
  
  # Compare with the MAE built buy Flomics constructor 
  MAE2 <- MAE # none, none at the building
  
  MAE2[["RNAtest"]] <- RFLOMICS:::checkTransNorm(MAE2[["RNAtest"]])
  
  MAE3 <- MAE
  MAE3 <- RFLOMICS::TransformData(MAE3, SE = "RNAtest", transform_method = "none")  
  MAE3 <- RFLOMICS::RunNormalization(MAE3, SE = "RNAtest", NormMethod = "TMM")

  MAE3[["RNAtest"]] <- RFLOMICS:::checkTransNorm(MAE3[["RNAtest"]])

  expect_identical(SummarizedExperiment::assay(MAE2[["RNAtest"]]), SummarizedExperiment::assay(MAE3[["RNAtest"]]))  
  
  # Compare with the transformation directely on raw data
  rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))
  
  normFactors <- edgeR::calcNormFactors(rnaSeqMat, method = "TMM")
  libSize <-  colSums(rnaSeqMat)
  
  tnDat <- log2(scale(rnaSeqMat + 1, center = FALSE, scale = normFactors*libSize))
  
  expect_identical(as.data.frame(SummarizedExperiment::assay(MAE2[["RNAtest"]])), data.frame(tnDat))  

})

# ---- RNAseq - Check if transformation working without arguments ----
# 
test_that("RNAseq - no arguments - methods", {
  
  MAE2 <- RFLOMICS::TransformData(MAE, SE = "RNAtest", modify_assay = FALSE)
  
  # Message if argument method is forgotten
  expect_message(RFLOMICS::TransformData(MAE, SE = "RNAtest"))

  # Default transformation must be "none" for RNAseq data when no transformation is asked.
  expect(MAE2[["RNAtest"]]@metadata$transform$transform_method == "none", failure_message = "The transformation is not 'none' by default.")
  
}) 



