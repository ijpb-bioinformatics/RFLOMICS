library(testthat)
library(RFLOMICS)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for coseq analysis : ----
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


formulae <- RFLOMICS::GetModelFormulae(MAE = MAE) 
MAE <- MAE |>
  RFLOMICS::getExpressionContrast(model.formula = formulae[[1]]) 
MAE <- MAE  |> RFLOMICS::getContrastMatrix(contrastList = c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                                            "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS))/3",
                                                            "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" )) 
contrastsDF <- RFLOMICS::getSelectedContrasts(MAE)

MAE2 <- MAE

MAE <- MAE |>
  TransformData(     SE.name = "metatest",  transformMethod = "log2")          |>
  RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")            |>
  RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")                 |>
  RunNormalization(  SE.name = "protetest", NormMethod = "median")              |>
  FilterLowAbundance(SE.name = "RNAtest")                                       |>
  RunDiffAnalysis(   SE.name = "metatest",  DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "protetest", DiffAnalysisMethod = "limmalmFit")  |>
  RunDiffAnalysis(   SE.name = "RNAtest",   DiffAnalysisMethod = "edgeRglmfit") |>
  FilterDiffAnalysis(SE.name = "RNAtest",   Adj.pvalue.cutoff = 0.05, logFC.cutoff = 1.5)

## ---- Construction of data tables differential analysis : ----

protMat <- RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt"))
rnaMat  <- RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))
metMat  <- RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt"))
condMat <- RFLOMICS::read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))

condMat$Repeat      <- relevel(condMat$Repeat, ref = "rep1")
condMat$imbibition  <- relevel(condMat$imbibition, ref = "DS")
condMat$temperature <- relevel(condMat$temperature, ref = "Low")

orderNames <- rownames(condMat)

protMat <- protMat[match(orderNames, colnames(protMat))]
rnaMat  <- rnaMat[match(orderNames, colnames(rnaMat))]
metMat  <- metMat[match(orderNames, colnames(metMat))]

# Contrasts
design <- model.matrix(~Repeat + temperature + imbibition + temperature:imbibition, data = condMat)

# Not checking if the coefficients are ok in there.
# taking the ones computed by RFLOMICS functions.
contrastsCoeff <- MAE@metadata$design@Contrasts.Coeff

## ---- Data Transformation ----
rnaMat2 <- rnaMat %>% dplyr::filter(rownames(.) %in% rownames(MAE[["RNAtest"]]))
protMat2 <- apply(protMat, 2, FUN = function(vect) vect - median(vect))
metMat2  <- apply(log2(metMat + 1), 2, FUN = function(vect) vect/sum(vect^2))


# ---- Tests unitaires -----

test_that("Everything works as expected", {
  
  # MAE2 don't have any differential analysis results. 
  expect_error(runCoExpression(object = MAE2, SE.name = "RNAtest", K = 2:10, replicates = 5, merge = "union", 
                               model = "normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin", 
                               normFactors = "TMM")) 
  
  # Wrong SE.name
  expect_error(runCoExpression(object = MAE, SE.name = "NA", K = 2:10, replicates = 5, merge = "union",
                               model = "normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin",
                               normFactors = "TMM"))
})

# ---- Tests seed -----

test_that("Two runs, same results - seed is working - RNAseq", {
  
  res1 <- runCoExpression (object = MAE, SE.name = "RNAtest", K = 2:10, replicates = 5, merge = "union", 
                               model = "normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin", 
                               normFactors = "TMM", nameList = c("H1"))
  
  # Wrong SE.name
  res2 <- runCoExpression(object = MAE, SE.name = "RNAtest", K = 2:10, replicates = 5, merge = "union",
                               model = "normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation = "arcsin",
                               normFactors = "TMM", nameList = c("H1"))
  
  expect_identical(coseq::clusters(res1[["RNAtest"]]@metadata$CoExpAnal$coseqResults),
                   coseq::clusters(res2[["RNAtest"]]@metadata$CoExpAnal$coseqResults))
  
  
})

test_that("Two runs, same results - seed is working - proteomics", {
  
  res1 <- runCoExpression(object = MAE, SE.name = "protetest", K = 2:20, replicates = 5, merge = "union", 
                          model = "normal")
  
  res2 <- runCoExpression(object = MAE, SE.name = "protetest", K = 2:20, replicates = 5, merge = "union",
                          model = "normal")
  
  expect_identical(coseq::clusters(res1[["protetest"]]@metadata$CoExpAnal$coseqResults),
                   coseq::clusters(res2[["protetest"]]@metadata$CoExpAnal$coseqResults))
  
  # Check GaussianModel
  expect_equal(res1[["protetest"]]@metadata$CoExpAnal$coseqResults@metadata$GaussianModel, 
               "Gaussian_pk_Lk_Bk")
  expect_equal(res2[["protetest"]]@metadata$CoExpAnal$coseqResults@metadata$GaussianModel, 
               "Gaussian_pk_Lk_Bk")
  
  # Check transformation
  expect_equal(res1[["protetest"]]@metadata$CoExpAnal$coseqResults@transformation, "none")
  expect_equal(res2[["protetest"]]@metadata$CoExpAnal$coseqResults@transformation, "none")
  
  # Check normFactors
  expect_equal(res1[["protetest"]]@metadata$CoExpAnal$coseqResults@normFactors, 
               rep(1, ncol(SummarizedExperiment::assay(res1[["protetest"]]))))
  expect_equal(res2[["protetest"]]@metadata$CoExpAnal$coseqResults@normFactors, 
               rep(1, ncol(SummarizedExperiment::assay(res2[["protetest"]]))))
  
  
})

# ------ TESTS --------###########

## ----- RNASeq -----

test_that("Coseq on RNAseq equivalence", {
  
  # Parameters for the three analyses
  merge = "union"
  
  K = 2:10
  replicates = 2
  iter <-  rep(K, each = replicates)
  geneList <- opDEList(MAE[["RNAtest"]], operation = merge)
  
  param.list = list(model = "normal",
                    GaussianModel = "Gaussian_pk_Lk_Ck",
                    transformation = "arcsin",
                    normFactors = "TMM",
                    meanFilterCutoff = 50)
  
  # RFLOMICS function:
  MAE <- runCoExpression(object = MAE, SE.name = "RNAtest", 
                         K = K, replicates = replicates, merge = merge, 
                         model = param.list$model, 
                         GaussianModel = param.list$GaussianModel, 
                         transformation = param.list$transformation, 
                         normFactors = param.list$normFactors, 
                         meanFilterCutoff = param.list$meanFilterCutoff, silent = FALSE)
  
  clustersMAE <- coseq::clusters(MAE[["RNAtest"]]@metadata$CoExpAnal$coseqResults)
  resultsMAE  <- MAE[["RNAtest"]]@metadata$CoExpAnal$coseqResults
  tcountsMAE <- resultsMAE@tcounts
  yprofMAE <- resultsMAE@y_profiles
  
  # Test run coseqLocal only:
  countMat <- SummarizedExperiment::assay(MAE[["RNAtest"]])[geneList,]
  countMat <- countMat[, match(rownames(MAE[["RNAtest"]]@metadata$Groups), colnames(countMat))]
  identical(rownames(MAE[["RNAtest"]]@metadata$Groups), colnames(countMat), attrib.as.set = FALSE)
  
  test_rcl <- runCoseq_local(countMat, 
                             conds = MAE[["RNAtest"]]@metadata$Groups$groups,
                             K = K, replicates = replicates, 
                             param.list = param.list)
  
  
  clustersRCL <- coseq::clusters(test_rcl$coseqResults)
  tcountsRCL <- test_rcl$coseqResults@tcounts
  yprofRCL <- test_rcl$coseqResults@y_profiles
  
  # table(clustersRCL, clustersMAE)
  
  expect_equal(tcountsMAE, tcountsRCL)
  expect_equal(yprofRCL, yprofMAE)
  expect_identical(clustersRCL, clustersMAE)
  
  # Equivalent pipeline outside (supposed to be equivalent)
  
  countMat <- rnaMat[geneList,]
  countMat <- countMat[, match(rownames(MAE[["RNAtest"]]@metadata$Groups), colnames(countMat))]
  identical(rownames(MAE[["RNAtest"]]@metadata$Groups), colnames(countMat), attrib.as.set = FALSE)
  
  # set.seed(12345)
  co <- capture.output(suppressMessages(
    coseq_res <-  lapply(1:replicates, function(x){
      coseq::coseq(countMat, K = K, 
                   model = param.list$model, 
                   transformation = param.list$transformation, 
                   GaussianModel = param.list$GaussianModel, 
                   normFactors = param.list$normFactors, 
                   meanFilterCutoff = param.list$meanFilterCutoff, 
                   parallel = TRUE, seed = x)
    })
  ))
  names(coseq_res) <- c(1:replicates)

  ICL.vec <- lapply(1:length(coseq_res), function(x){ coseq::ICL(coseq_res[[x]]) }) %>% unlist()
  ICL.tab <- data.frame(K = stringr::str_replace(names(ICL.vec), "K=", ""), ICL = ICL.vec) %>%
    dplyr::mutate(K = as.numeric(K))
  
  ICL.n <- ICL.tab  %>% 
    dplyr::group_by(., K) %>% dplyr::filter(!is.na(ICL)) %>%
    dplyr::summarise(median = median(ICL, na.rm = TRUE), n = dplyr::n()) %>%
    dplyr::mutate(K = as.numeric(K))
  
  K.ICL.median.min <- ICL.n[which.min(ICL.n$median),]$K
  K.ICL.min <- min(ICL.vec[names(ICL.vec) == paste0("K=", K.ICL.median.min)], na.rm = TRUE)
  
  index <- sapply(names(coseq_res), function(x){(TRUE %in% (coseq::ICL(coseq_res[[x]]) == K.ICL.min))})
  
  coseq.res <- coseq_res[index][[1]]
  
  clustersRES <- coseq::clusters(coseq.res)
  yprofRES <- coseq.res@y_profiles
  tcountsRES <- coseq.res@tcounts

  expect_equal(tcountsMAE, tcountsRES)
  expect_equal(yprofRES, yprofMAE)
  expect_identical(clustersRES, clustersMAE)

})

## ---- RNAseq - ClusterMQ ----
# 
# test_that("Coseq on RNAseq equivalence - clusterMQ", {
# 
#   # Parameters for the three analyses
#   merge = "union"
# 
#   K = 2:10
#   replicates = 2
#   iter <-  rep(K, each = replicates)
# 
#   rep(1:replicates, replicates)
# 
#   param.list = list(model = "normal",
#                     GaussianModel = "Gaussian_pk_Lk_Ck",
#                     transformation = "arcsin",
#                     normFactors = "TMM",
#                     meanFilterCutoff = 50)
# 
#   # RFLOMICS function:
#   MAE <- runCoExpression(object = MAE, SE.name = "RNAtest",
#                          K = K, replicates = replicates, merge = merge,
#                          model = param.list$model,
#                          GaussianModel = param.list$GaussianModel,
#                          transformation = param.list$transformation,
#                          normFactors = param.list$normFactors,
#                          meanFilterCutoff = param.list$meanFilterCutoff,
#                          clustermq = TRUE,
#                          silent = FALSE, cmd = TRUE)
# 
#   clustersMAE <- coseq::clusters(MAE[["RNAtest"]]@metadata$CoExpAnal$coseqResults)
#   resultsMAE  <- MAE[["RNAtest"]]@metadata$CoExpAnal$coseqResults
#   tcountsMAE <- resultsMAE@tcounts
#   yprofMAE <- resultsMAE@y_profiles
# 
#   # Equivalent pipeline outside (supposed to be equivalent)
# 
#   countMat <- data.frame(SummarizedExperiment::assay(MAE[["RNAtest"]])) %>%
#     dplyr::filter(rownames(.) %in% opDEList(MAE, SE.name = "RNAtest", operation = merge))
# 
#   fx <- function(x, seed_arg){
#     co <- suppressMessages(capture.output(
#       res <- coseq::coseq(object = param.list[["object"]],
#                           K = x,
#                           model = param.list$model,
#                           transformation = param.list$transformation,
#                           GaussianModel = param.list$GaussianModel,
#                           normFactors = param.list$normFactors,
#                           meanFilterCutoff = param.list$meanFilterCutoff,
#                           seed = seed_arg)  ))
# 
#     return(res)
#   }
# 
#   iter <-  rep(K, each = replicates)
#   seed_arg = rep(1:replicates, max(K) - 1)
# 
#   df_args <- cbind(iter, seed_arg)
#   colnames(df_args) = c("x", "seed_arg")
#   nbr_iter <- length(iter)
# 
#   coseq.res.list <- list()
#   coseq.res.list <- clustermq::Q_rows(fun = fx,
#                                  df = df_args,
#                                  export = param.list, n_jobs = nbr_iter, pkgs = "coseq")
# 
#   names(coseq.res.list) <- c(1:nbr_iter)
# 
# })
# 
# 
## ----- Proteomics -----

test_that("Coseq on Proteomics equivalence", {
  
  # Parameters for the three analyses
  merge = "union"
  
  geneList <- opDEList(MAE[["protetest"]], operation = merge)
  
  K = 2:20
  replicates = 2
  iter <-  rep(K, each = replicates)
  
  param.list = list(model = "normal",
                    GaussianModel = "Gaussian_pk_Lk_Bk",
                    transformation = "none",
                    normFactors = "none",
                    meanFilterCutoff = NULL)
  
  # RFLOMICS function:
  MAE <- runCoExpression(object = MAE, SE.name = "protetest", 
                         K = K, replicates = replicates, merge = merge, 
                         model = param.list$model, 
                         GaussianModel = param.list$GaussianModel, 
                         transformation = param.list$transformation, 
                         normFactors = param.list$normFactors, 
                         meanFilterCutoff = param.list$meanFilterCutoff)
  
  clustersMAE <- coseq::clusters(MAE[["protetest"]]@metadata$CoExpAnal$coseqResults)
  resultsMAE  <- MAE[["protetest"]]@metadata$CoExpAnal$coseqResults
  tcountsMAE <- resultsMAE@tcounts
  yprofMAE <- resultsMAE@y_profiles
  
  # Test run coseqLocal only:
  
  countMat <- data.frame(SummarizedExperiment::assay(MAE[["protetest"]])) 
  countMat <- apply(countMat, 2, FUN = function(x) x - median(x))
  countMat <- countMat[geneList,]
  countMat <- countMat[, match(rownames(MAE[["protetest"]]@metadata$Groups), colnames(countMat))]
  identical(rownames(MAE[["protetest"]]@metadata$Groups), colnames(countMat), attrib.as.set = FALSE)
  
  countMat2 <- t(apply(countMat , 1, function(x){scale(x, center = TRUE, scale = TRUE) }))
  colnames(countMat2) <- colnames(countMat)
  
  test_rcl <- runCoseq_local(countMat2, 
                             conds = MAE[["protetest"]]@metadata$Groups$groups,
                             K = K, replicates = replicates, 
                             param.list = param.list)
  
  
  clustersRCL <- coseq::clusters(test_rcl$coseqResults)
  tcountsRCL <- test_rcl$coseqResults@tcounts
  yprofRCL <- test_rcl$coseqResults@y_profiles

  # table(clustersRCL, clustersMAE)
  
  expect_equal(tcountsMAE, tcountsRCL)
  expect_equal(data.frame(tcountsRCL), data.frame(countMat2)) # no transformation, no filtering
  expect_equal(yprofRCL, yprofMAE)
  expect_identical(clustersRCL, clustersMAE)
  
  # Equivalent pipeline outside (supposed to be equivalent)
  
  countMat <- protMat2[geneList,]
  countMat <- countMat[, match(rownames(MAE[["protetest"]]@metadata$Groups), colnames(countMat))]
  identical(rownames(MAE[["protetest"]]@metadata$Groups), colnames(countMat), attrib.as.set = FALSE)
  
  countMat2 <- t(apply(countMat , 1, function(x){scale(x, center = TRUE, scale = TRUE) }))
  colnames(countMat2) <- colnames(countMat)
  
  # set.seed(12345)
  co <- capture.output(suppressMessages(
    coseq_res <-  lapply(1:replicates, function(x){
      coseq::coseq(countMat2, K = K, 
                   model = param.list$model, 
                   transformation = param.list$transformation, 
                   GaussianModel = param.list$GaussianModel, 
                   normFactors = param.list$normFactors, 
                   meanFilterCutoff = param.list$meanFilterCutoff, 
                   parallel = TRUE, seed = x)
    })
  ))
  names(coseq_res) <- c(1:replicates)
  
  ICL.vec <- lapply(1:length(coseq_res), function(x){ coseq::ICL(coseq_res[[x]]) }) %>% unlist()
  ICL.tab <- data.frame(K = stringr::str_replace(names(ICL.vec), "K=", ""), ICL = ICL.vec) %>%
    dplyr::mutate(K = as.numeric(K))
  
  ICL.n <- ICL.tab  %>% 
    dplyr::group_by(., K) %>% dplyr::filter(!is.na(ICL)) %>%
    dplyr::summarise(median = median(ICL, na.rm = TRUE), n = dplyr::n()) %>%
    dplyr::mutate(K = as.numeric(K))
  
  K.ICL.median.min <- ICL.n[which.min(ICL.n$median),]$K
  K.ICL.min <- min(ICL.vec[names(ICL.vec) == paste0("K=", K.ICL.median.min)], na.rm = TRUE)
  
  index <- sapply(names(coseq_res), function(x){(TRUE %in% (coseq::ICL(coseq_res[[x]]) == K.ICL.min))})
  
  coseq.res <- coseq_res[index][[1]]
  
  clustersRES <- coseq::clusters(coseq.res)
  yprofRES <- coseq.res@y_profiles
  tcountsRES <- coseq.res@tcounts
  
  # table(clustersRES, clustersMAE)
  expect_equal(tcountsMAE, tcountsRES) 
  expect_equal(data.frame(tcountsRES), data.frame(countMat2)) # no transformation, no filtering
  expect_equal(yprofRES, yprofMAE)
  expect_identical(clustersRES, clustersMAE)
  
})