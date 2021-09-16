context("test-flomics-generics")

# data(count)
# data(QCmat)

# test_that("FlomicsExperiment Object can be instantiate without QCmatrix", {
#   FE <- FlomicsExperiment(count)
#   expect_true( class(FE) == "FlomicsExperiment")
# })


################################################
## tests for "getExpressionContrast" method ####

#######################
## case : 1 bio factors

# create ExpDesign class object 
ExpDesign.tbl <- data.frame(Treatment = c(rep("NoSi",3), rep("Si",3)),
                            Replicate = paste("R",rep(1:3, 2), sep=""), 
                            row.names = c("NoSi_R1", "NoSi_R2", "NoSi_R3", "Si_R1", "Si_R2", "Si_R3"))
dF.List.ref  <- c("NoSi", "R1");   names(dF.List.ref)  <- names(ExpDesign.tbl)
dF.Type.dFac <- c("Bio", "batch"); names(dF.Type.dFac) <- names(ExpDesign.tbl)
Design <- ExpDesign.constructor(ExpDesign = ExpDesign.tbl, projectName = "Brassica", refList = dF.List.ref, typeList = dF.Type.dFac)


########
test_that("getExpressionContrast method return object from ExpDesign class (in 1 bio factor cases)", {
  
  GLM_model <- "~Replicate + Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_true( class(Design) == "ExpDesign")
})

########
test_that("getExpressionContrast generate 1 types of contrasts (in 1 bio factor cases)", {
 
  GLM_model <- "~Replicate + Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_equal(names(Design@Contrasts.List), c("simple"))
})

########
test_that("getExpressionContrast generate correct expression contrast (in 1 bio factor cases)", {
  
  GLM_model <- "~Replicate + Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_equal(Design@Contrasts.List[["simple"]]$contrast, c("(TreatmentSi - TreatmentNoSi)"))
})

#######################
## case : 2 bio factors

# create ExpDesign class object
ExpDesign.tbl <- data.frame(tissue   =c(rep("MatureLeaf", 6), rep("Root", 6)), 
                            Treatment=c(rep("NoSi",3), rep("Si",3), rep("NoSi",3), rep("Si",3)),
                            Replicate=paste("R",rep(1:3, 4), sep=""), 
                            row.names = c("MatureLeaf_NoSi_R1", "MatureLeaf_NoSi_R2", "MatureLeaf_NoSi_R3", "MatureLeaf_Si_R1", "MatureLeaf_Si_R2", "MatureLeaf_Si_R3", 
                                          "Root_NoSi_R1", "Root_NoSi_R2", "Root_NoSi_R3", "Root_Si_R1", "Root_Si_R2", "Root_Si_R3"))
dF.List.ref <- c("Root","NoSi", "R1");  names(dF.List.ref) <- names(ExpDesign.tbl)
dF.Type.dFac <- c("Bio", "Bio", "batch"); names(dF.Type.dFac) <- names(ExpDesign.tbl)
Design <- ExpDesign.constructor(ExpDesign = ExpDesign.tbl, projectName = "Brassica", refList = dF.List.ref, typeList = dF.Type.dFac)

########
test_that("getExpressionContrast method return object from ExpDesign class (in 2 bio factor case with interaction)", {
  
  GLM_model <- "~Replicate + tissue + Treatment + tissue:Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_true( class(Design) == "ExpDesign")
})

########
test_that("getExpressionContrast method return object from ExpDesign class (in 2 bio factor case without interaction)", {
  
  GLM_model <- "~Replicate + tissue + Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_true( class(Design) == "ExpDesign")
})

########
test_that("getExpressionContrast generate 3 types of contrasts (in 2 bio factor case with interaction)", {
  
  GLM_model <- "~Replicate + tissue + Treatment + tissue:Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_equal(names(Design@Contrasts.List), c("simple", "averaged", "interaction"))
})

########
test_that("getExpressionContrast generate correct contrast expressions (in 2 bio factor case with interaction)", {
  
  GLM_model <- "~Replicate + tissue + Treatment + tissue:Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  
  exp.results <- c("(tissueRoot - tissueMatureLeaf) in TreatmentNoSi", "(tissueRoot - tissueMatureLeaf) in TreatmentSi", 
                   "(TreatmentSi - TreatmentNoSi) in tissueMatureLeaf", "(TreatmentSi - TreatmentNoSi) in tissueRoot",
                   "(tissueRoot - tissueMatureLeaf) in mean", "(TreatmentSi - TreatmentNoSi) in mean",
                   "(tissueRoot - tissueMatureLeaf) in TreatmentSi - (tissueRoot - tissueMatureLeaf) in TreatmentNoSi") 
  obs.results <- lapply(c("simple", "averaged", "interaction"), function(i){ x <- Design@Contrasts.List[[i]]$contrastName }) %>% purrr::reduce(., c)
  
  expect_equal(exp.results, obs.results)
})

########
test_that("getExpressionContrast generate 1 types of contrasts (in 2 bio factor case without interaction)", {
  
  GLM_model <- "~Replicate + tissue + Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  expect_equal(names(Design@Contrasts.List), c("simple", "averaged")) #########################################
})

########
test_that("getExpressionContrast generate correct contrast expressions (in 2 bio factor case without interaction)", {
  
  GLM_model <- "~Replicate + tissue + Treatment"
  Design <- getExpressionContrast(Design, GLM_model)
  
  exp.results <- c("(tissueRoot - tissueMatureLeaf) in TreatmentNoSi", "(tissueRoot - tissueMatureLeaf) in TreatmentSi", 
                   "(TreatmentSi - TreatmentNoSi) in tissueMatureLeaf", "(TreatmentSi - TreatmentNoSi) in tissueRoot",
                   "(tissueRoot - tissueMatureLeaf) in mean", "(TreatmentSi - TreatmentNoSi) in mean") 
  obs.results <- lapply(c("simple", "averaged"), function(i){ x <- Design@Contrasts.List[[i]]$contrastName }) %>% purrr::reduce(., c)
  
  expect_equal(exp.results, obs.results)
})


