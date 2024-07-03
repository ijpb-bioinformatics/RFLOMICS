library(testthat)
library(RFLOMICS)

##### Checks on transformation and normalization of data through command line functions
##### Essentially testing the .checkTransNorm pipeline. 

# ---- Construction of objects for the tests ---- 

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
names(MAE) <- c("RNAtest", "metatest", "protetest")

# Comparison data:
protMat <- ecoseed$protetest
rnaMat <- ecoseed$RNAtest

rnaMat <- rnaMat[, match(colnames(assay(MAE[["RNAtest"]])), colnames(rnaMat))]
protMat <- protMat[, match(colnames(assay(MAE[["protetest"]])), colnames(protMat))]

# ---- Some functions ----

isTransformed <- function(object, SE.name) RFLOMICS:::.isTransformed(object[[SE.name]])
isNorm        <- function(object, SE.name) RFLOMICS:::.isNorm(object[[SE.name]])

######################################-
########### FUNCTIONS TESTS ###########
######################################-

# ---- transformData, apply_transform ----

test_that("transformData and apply_transform yield expected results", {
    
    MAE2 <- MAE3 <- MAE4 <- MAE5 <- MAE5b <- MAE6 <- MAE 
    rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))
    
    ####
    # --- no transformation, no modification asked, nothing is supposed to happen.
    
    
    MAE2 <- runTransformData(MAE, SE.name = "protetest", modifyAssay = FALSE, transformMethod = "none")
    expect_equal( assay(MAE2[["protetest"]]), as.matrix(protMat))
    expect(!RFLOMICS:::.isTransformed(MAE2[["protetest"]]), failure_message = "It was transformed, it shouldn't be.")
    
    MAE2 <- runTransformData(MAE, SE.name = "RNAtest", modifyAssay = FALSE, transformMethod = "none")
    expect_equal( assay(MAE2[["RNAtest"]]), as.matrix(rnaSeqMat))
    expect(!RFLOMICS:::.isTransformed(MAE2[["RNAtest"]]), failure_message = "It was transformed, it shouldn't be.")
    
    ####
    # --- modifyAssay = FALSE, transform method is wrong. 
    
    MAE3 <- runTransformData(MAE, SE.name = "protetest", modifyAssay = FALSE, transformMethod = "nothing")
    expect_equal( assay(MAE3[["protetest"]]), as.matrix(protMat))
    expect(!RFLOMICS:::.isTransformed(MAE3[["protetest"]]), failure_message = "It was transformed, it shouldn't be.")
    # TODO It should have a warning... or something
    
    MAE3 <- runTransformData(MAE, SE.name = "RNAtest", modifyAssay = FALSE, transformMethod = "nothing")
    expect_message(runTransformData(MAE, SE.name = "RNAtest", modifyAssay = FALSE, transformMethod = "nothing")) # forced to "none".
    expect_equal( assay(MAE3[["RNAtest"]]), as.matrix(rnaSeqMat))
    expect(!RFLOMICS:::.isTransformed(MAE3[["RNAtest"]]), failure_message = "It was transformed, it shouldn't be.")
    
    ####
    # --- modifyAssay = TRUE, wrong transformation
    
    MAE4 <- runTransformData(MAE, SE.name = "protetest", modifyAssay = TRUE, transformMethod = "nothing")
    expect_message(runTransformData(MAE, SE.name = "protetest", modifyAssay = TRUE, transformMethod = "nothing"))
    expect_equal( assay(MAE4[["protetest"]]), as.matrix(protMat)) # no transformation applied
    expect(!RFLOMICS:::.isTransformed(MAE4[["protetest"]]), failure_message = "It was transformed, it shouldn't be.")
    
    MAE4 <- runTransformData(MAE, SE.name = "RNAtest", modifyAssay = TRUE, transformMethod = "nothing")
    expect_message(runTransformData(MAE, SE.name = "RNAtest", modifyAssay = TRUE, transformMethod = "nothing"))
    expect_equal( assay(MAE4[["RNAtest"]]), as.matrix(rnaSeqMat)) # no transformation applied
    expect(!RFLOMICS:::.isTransformed(MAE4[["RNAtest"]]), failure_message = "It was transformed, it shouldn't be.")
    
    ####
    # --- right transformation, no assay modification.
    
    MAE5 <- MAE5b <- runTransformData(MAE, SE.name = "protetest", modifyAssay = FALSE, transformMethod = "log2")
    expect_equal( assay(MAE5[["protetest"]]), as.matrix(protMat)) 
    expect(!RFLOMICS:::.isTransformed(MAE5[["protetest"]]), failure_message = "It was transformed, it shouldn't be.")
    
    ####
    # --- apply transformation:
    
    MAE5b[["protetest"]] <- RFLOMICS:::.applyTransformation(MAE5[["protetest"]])
    expect_equal( assay(MAE5b[["protetest"]]), as.matrix(log2(protMat + 10^-10))) 
    expect(RFLOMICS:::.isTransformed(MAE5b[["protetest"]]), failure_message = "It wasn't transformed, it should be.")
    
    ####
    # --- apply transformation directly: 
    
    MAE6 <- runTransformData(MAE, SE.name = "protetest", modifyAssay = TRUE, transformMethod = "log2")
    expect_equal( assay(MAE6[["protetest"]]), as.matrix(log2(protMat + 10^-10))) 
    expect(RFLOMICS:::.isTransformed(MAE6[["protetest"]]), failure_message = "It wasn't transformed, it should be.")
    
    
})

# ---- RunNormalization, apply_norm ----

test_that("RunNormalization and apply_norm yield expected results", {
    
    MAE2 <- MAE3 <- MAE4 <- MAE5 <- MAE5b <- MAE6 <- MAE 
    rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))
    
    ####
    # --- Missing norm argument
    
    MAE2 <- runNormalization(MAE, SE.name = "protetest", modifyAssay = FALSE)
    expect_equal( assay(MAE2[["protetest"]]), as.matrix(protMat))
    expect( !isNorm(MAE2,"protetest"), failure_message = "It was normalized, it shouldn't be.")
    
    MAE2 <- runNormalization(MAE, SE.name = "RNAtest", modifyAssay = FALSE)
    expect_equal( assay(MAE2[["RNAtest"]]), as.matrix(rnaSeqMat))
    expect(getNormSettings(MAE2[["RNAtest"]])$method == "TMM", failure_message = "TMM was not forced on RNAseq data")
    expect(!isNorm(MAE2,"RNAtest"), failure_message = "It was normalized, it shouldn't be.")
    
    ####
    # --- no norm, no modification asked, nothing is supposed to happen.
    
    MAE2 <- runNormalization(MAE, SE.name = "protetest",normMethod = "none",  modifyAssay = FALSE)
    expect_equal( assay(MAE2[["protetest"]]), as.matrix(protMat))
    expect(!isNorm(MAE2,"protetest"), failure_message = "It was normalized, it shouldn't be.")
    
    MAE2 <- runNormalization(MAE, SE.name = "RNAtest", normMethod = "none", modifyAssay = FALSE)
    expect_equal( assay(MAE2[["RNAtest"]]), as.matrix(rnaSeqMat))
    expect(getNormSettings(MAE2[["RNAtest"]])$method == "TMM", failure_message = "TMM was not forced on RNAseq data")
    expect(!isNorm(MAE2,"RNAtest"), failure_message = "It was normalized, it shouldn't be.")
    
    ####
    # --- modifyAssay = FALSE
    
    MAE3 <- runNormalization(MAE, SE.name = "protetest", modifyAssay = FALSE, normMethod = "nothing")
    expect_equal( assay(MAE3[["protetest"]]), as.matrix(protMat))
    expect(!isNorm(MAE3,"protetest"), failure_message = "It was normalized, it shouldn't be.")
    # TODO It should have a warning... or something
    
    MAE3 <- runNormalization(MAE, SE.name = "RNAtest", modifyAssay = FALSE, normMethod = "nothing")
    expect_message(runNormalization(MAE, SE.name = "RNAtest", modifyAssay = FALSE, normMethod = "nothing")) # forced to "none".
    expect_equal( assay(MAE3[["RNAtest"]]), as.matrix(rnaSeqMat))
    expect(!isNorm(MAE3, "RNAtest"), failure_message = "It was normalized, it shouldn't be.")
    
    ####
    # --- modifyAssay = TRUE, wrong transformation
    
    MAE4 <- runNormalization(MAE, SE.name = "protetest", modifyAssay = TRUE, normMethod = "nothing")
    expect_message(runNormalization(MAE, SE.name = "protetest", modifyAssay = TRUE, normMethod = "nothing"))
    expect_equal( assay(MAE4[["protetest"]]), as.matrix(protMat)) # no transformation applied
    expect(!isNorm(MAE4, "protetest"), failure_message = "It was normalized, it shouldn't be.")
    
    MAE4 <- runNormalization(MAE, SE.name = "RNAtest", modifyAssay = TRUE, normMethod = "nothing")
    expect_message(runNormalization(MAE, SE.name = "RNAtest", modifyAssay = TRUE, normMethod = "nothing"))
    
    scales_factors <- getCoeffNorm(MAE4[["RNAtest"]])$lib.size*getCoeffNorm(MAE4[["RNAtest"]])$norm.factors
    assayTransform <-  assay(MAE[["RNAtest"]])
    RNAnorm <- scale(assayTransform + 1, center = FALSE, scale = scales_factors)
    
    expect_equal(as.data.frame( assay(MAE4[["RNAtest"]])), as.data.frame(as.matrix(RNAnorm))) # transformation applied: "TMM"
    expect(isNorm(MAE4, "RNAtest"), failure_message = "It wasn't normalized.")
    
    ####
    # --- right transformation.
    
    MAE5 <- MAE5b <- runNormalization(MAE, SE.name = "protetest", modifyAssay = FALSE, normMethod = "median")
    expect_equal( assay(MAE5[["protetest"]]), as.matrix(protMat)) 
    expect(!isNorm(MAE5, "protetest"), failure_message = "It was normalized, it shouldn't be.")
    
    ####
    # --- apply transformation:
    
    MAE5b[["protetest"]] <- RFLOMICS:::.applyNorm(MAE5[["protetest"]])
    # getNorm(MAE5b, "protetest") # median
    protMed <- apply(protMat, 2, FUN = function(vect) vect - median(vect))
    expect_equal( assay(MAE5b[["protetest"]]), as.matrix(protMed)) 
    expect(isNorm(MAE5b, "protetest"), failure_message = "It wasn't normalized.")
    
    ####
    # --- apply transformation directly: 
    
    MAE6 <- runNormalization(MAE, SE.name = "protetest", modifyAssay = TRUE, normMethod = "median")
    protMed <- apply(protMat, 2, FUN = function(vect) vect - median(vect))
    expect_equal( assay(MAE6[["protetest"]]), as.matrix(protMed)) 
    expect(isNorm(MAE6, "protetest"), failure_message = "It wasn't normalized.")
    
    
})



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
        # print(case_vect)
        # matrix version
        protMattransnorm <- protMat
        
        pca.raw <- FactoMineR::PCA(t(protMattransnorm), ncp = 5, graph = FALSE)
        
        expect_equal(pca.raw$eig, MAE[["protetest"]]@metadata$PCAlist$raw$eig)
        expect_equal(pca.raw$svd, MAE[["protetest"]]@metadata$PCAlist$raw$svd)
        expect_equal(pca.raw$ind, MAE[["protetest"]]@metadata$PCAlist$raw$ind)
        expect_equal(pca.raw$var, MAE[["protetest"]]@metadata$PCAlist$raw$var)
        # the call is obligatory different between the two
        
        protMattransnorm <- switch(as.character(case_vect[[1]]),
                                   "none"       = protMattransnorm,
                                   "log2"       = log2(protMattransnorm + 10^-10),
                                   "log10"      = log10(protMattransnorm + 10^-10),
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
        MAE2 <- runTransformData(MAE2, SE = "protetest", transformMethod = as.character(case_vect[[1]]))
        MAE2 <- runNormalization(MAE2, SE.name = "protetest", normMethod = as.character(case_vect[[2]]))
        
        MAE2 <- RFLOMICS::runOmicsPCA(MAE2, SE = "protetest")
        expect_equal(pca.norm$eig, MAE2[["protetest"]]@metadata$PCAlist$norm$eig)
        expect_equal(pca.norm$svd, MAE2[["protetest"]]@metadata$PCAlist$norm$svd)
        expect_equal(pca.norm$ind, MAE2[["protetest"]]@metadata$PCAlist$norm$ind)
        expect_equal(pca.norm$var, MAE2[["protetest"]]@metadata$PCAlist$norm$var)
        
        MAE2[["protetest"]] <- RFLOMICS:::.checkTransNorm(MAE2[["protetest"]]) 
        
        expect_equal( assay(MAE2[["protetest"]]), as.matrix(protMattransnorm))
        
    })
    
})

# ---- PROTEO - Check if transformation working without arguments ----
test_that("Transformation - no method - no modification", {
    
    MAE2 <- runTransformData(MAE, SE = "protetest")
    
    expect_message(runTransformData(MAE, SE = "protetest"))
    expect(getTransSettings(MAE2[["protetest"]])$method == "log2", failure_message = "The transformation is not log2 by default.")
    expect_equal( assay(MAE[["protetest"]]),  assay(MAE2[["protetest"]]))
})

test_that("Transformation - no method - modification", { # 
    
    MAE2 <- runTransformData(MAE, SE = "protetest", transformMethod = NULL, modifyAssay = TRUE)
    
    # Message if argument method is forgotten
    expect_message(runTransformData(MAE, SE = "protetest"))
    
    # Default transformation must be log2 + 1 for proteomics data
    expect_equal(log2( assay(MAE[["protetest"]]) + 10^-10),  assay(MAE2[["protetest"]]))
    expect(getTransSettings(MAE2[["protetest"]])$method == "log2", failure_message = "The transformation is not log2 by default.")
    
    # Modification of the assay: expect transformed metadata to be TRUE 
    expect(RFLOMICS:::.isTransformed(MAE2[["protetest"]]), failure_message = "The assay was not transformed.")
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
    
    expect_equal(pca.raw$eig, MAE[["RNAtest"]]@metadata$PCAlist$raw$eig)
    expect_equal(pca.raw$svd, MAE[["RNAtest"]]@metadata$PCAlist$raw$svd)
    expect_equal(pca.raw$ind, MAE[["RNAtest"]]@metadata$PCAlist$raw$ind)
    expect_equal(pca.raw$var, MAE[["RNAtest"]]@metadata$PCAlist$raw$var)
    # the call is obligatory different between the two
    
    # Manually transforming and normalizing rnaSeqMat
    normFactors <- edgeR::calcNormFactors(rnaSeqMat, method = "TMM")
    libSize <-  colSums(rnaSeqMat)
    
    tnDat <- scale(rnaSeqMat + 1, center = FALSE, scale = normFactors*libSize)
    
    pca.norm <- FactoMineR::PCA(t(log2(tnDat)), ncp = 5, graph = FALSE)
    
    # RFLOMICS version
    MAE2 <- MAE
    MAE2 <- runTransformData(MAE2, SE = "RNAtest", transformMethod = "none")
    MAE2 <- runNormalization(MAE2, SE.name = "RNAtest", normMethod = "TMM")
    
    MAE2 <- RFLOMICS::runOmicsPCA(MAE2, SE = "RNAtest")
    expect_equal(pca.norm$eig, MAE2[["RNAtest"]]@metadata$PCAlist$norm$eig)
    expect_equal(pca.norm$svd, MAE2[["RNAtest"]]@metadata$PCAlist$norm$svd)
    expect_equal(pca.norm$ind, MAE2[["RNAtest"]]@metadata$PCAlist$norm$ind)
    expect_equal(pca.norm$var, MAE2[["RNAtest"]]@metadata$PCAlist$norm$var)
    
})

# ---- RNAseq - behaviour of transforData and RunNormalization ---

test_that("RNAseq - correct behaviour of normalization and transformation",  {
    
    MAE2 <- MAE
    MAE2 <- runTransformData(MAE2, SE = "RNAtest")
    expect(getTransSettings(MAE2[["RNAtest"]])$method == "none", 
           failure_message = "TransformData does not put 'none' as default for RNAseq data")
    expect(!RFLOMICS:::.isTransformed(MAE2[["RNAtest"]]), 
           failure_message = "TransformData transformed the data when not asked to")
    
    MAE2 <- runNormalization(MAE2, SE = "RNAtest")
    expect(getNormSettings(MAE2[["RNAtest"]])$method == "TMM", 
           failure_message = "Normalization is not defaulted to TMM for RNAseq data.")
    expect(!RFLOMICS:::.isNorm(MAE2[["RNAtest"]]), 
           failure_message = "RunNormalization transformed the data when not asked to")
    
    MAE2 <- runTransformData(MAE2, SE = "RNAtest", transformMethod = "log2") # shouldn't be at all
    MAE2 <- runNormalization(MAE2, SE = "RNAtest", normMethod = "median") # shouldn't be at all.
    
    expect(getTransSettings(MAE2[["RNAtest"]])$method == "none",
           failure_message = "TransformData did not replace by 'none' for RNAseq data ")
    expect(!RFLOMICS:::.isTransformed(MAE2[["RNAtest"]]), 
           failure_message = "TransformData transformed the data when not asked to")
    
    expect(getNormSettings(MAE2[["RNAtest"]])$method == "TMM", 
           failure_message = "Normalization did not replace by TMM for RNAseq data.")
    expect(!RFLOMICS:::.isNorm(MAE2[["RNAtest"]]),
           failure_message = "RunNormalization transformed the data when not asked to")
    
}) 


test_that("RNAseq - correct behaviour of checkTransNorm",  {
    
    # Compare with the MAE built buy Flomics constructor 
    MAE2 <- MAE # none, none at the building
    
    MAE2[["RNAtest"]] <- RFLOMICS:::.checkTransNorm(MAE2[["RNAtest"]])
    
    MAE3 <- MAE
    MAE3 <- runTransformData(MAE3, SE = "RNAtest", transformMethod = "none")  
    MAE3 <- runNormalization(MAE3, SE = "RNAtest", normMethod = "TMM")
    
    MAE3[["RNAtest"]] <- RFLOMICS:::.checkTransNorm(MAE3[["RNAtest"]])
    
    expect_equal( assay(MAE2[["RNAtest"]]), 
                  assay(MAE3[["RNAtest"]]))  
    
    # Compare with the transformation directely on raw data
    rnaSeqMat <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))
    
    normFactors <- edgeR::calcNormFactors(rnaSeqMat, method = "TMM")
    libSize <-  colSums(rnaSeqMat)
    
    tnDat <- log2(scale(rnaSeqMat + 1, center = FALSE, scale = normFactors*libSize))
    
    expect_equal(as.data.frame( assay(MAE2[["RNAtest"]])), 
                 data.frame(tnDat))  
    
})

# ---- RNAseq - Check if transformation working without arguments ----
# 
test_that("RNAseq - no arguments - methods", {
    
    MAE2 <- runTransformData(MAE, SE = "RNAtest", modifyAssay = FALSE)
    
    # Message if argument method is forgotten
    expect_message(runTransformData(MAE, SE = "RNAtest"))
    
    # Default transformation must be "none" for RNAseq data when no transformation is asked.
    expect(getTransSettings(MAE2[["RNAtest"]])$method == "none", 
           failure_message = "The transformation is not 'none' by default.")
    
}) 



