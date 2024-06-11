library(testthat)
library(RFLOMICS)
library(MOFA2)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for integration analysis : ----
MAE <- generateExample(
    annotation = FALSE,
    coexp = FALSE,
    integration = FALSE
) 

MAE0 <- MAE

protMat <- readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt"))
metMat  <- readOmicsData(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt"))
condMat <- readExpDesign(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))

condMat$Repeat      <- factor(condMat$Repeat, levels = c("rep1", "rep2", "rep3"))
condMat$imbibition  <- factor(condMat$imbibition, levels = c("DS", "EI", "LI"))
condMat$temperature <- factor(condMat$temperature, levels = c("Low", "Medium", "Elevated"))

condMat$Repeat      <- relevel(condMat$Repeat, ref = "rep1")
condMat$imbibition  <- relevel(condMat$imbibition, ref = "DS")
condMat$temperature <- relevel(condMat$temperature, ref = "Low")

orderNames <- rownames(colData(MAE))
condMat <- condMat[match(orderNames, rownames(condMat)),]
protMat <- protMat[match(orderNames, colnames(protMat))]
metMat  <- metMat[match(orderNames, colnames(metMat))]

identical(orderNames, colnames(protMat), attrib.as.set = FALSE)
identical(orderNames, colnames(metMat), attrib.as.set = FALSE)

metMat2  <- apply(log2(metMat + 10^-10), 2, FUN = function(vect) vect/sum(vect^2))
protMat2 <- apply(protMat, 2, FUN = function(vect) vect - median(vect))

# from .rbeFunction
colBatch <- getBatchFactors(MAE)
newFormula <-
    gsub(pattern = paste(paste(colBatch, "[+]"),  collapse = "|"),
         "", getModelFormula(MAE))
colData <- getDesignMat(MAE)
designToPreserve <-
    model.matrix(as.formula(newFormula), data = colData)


metMat3 <- limma::removeBatchEffect(metMat2, batch = colData[, colBatch],
                                    design = designToPreserve)
protMat3 <- limma::removeBatchEffect(protMat2, batch = colData[, colBatch],
                                     design = designToPreserve)

# ----- TESTS -----

test_that("Equivalence", {
    selectedData <- c("metatest", "protetest")
    selectMethod <- c("metatest" = "diff", "protetest" = "none")
    operation <- c("metatest" = "union", "protetest" = "union")
    
    variableList <- 
        lapply(selectedData, function(set) {
            switch(
                selectMethod[[set]],
                "diff"  = getDEList(
                    object = MAE[[set]],
                    contrasts = getSelectedContrasts(MAE),
                    operation = operation[[set]]
                ),
                "none"  = names(MAE[[set]])
            )
        })
    names(variableList) <- selectedData
    
    MAE2 <- prepareForIntegration(object           = MAE,
                                  omicsNames       = selectedData,
                                  variableLists    = variableList,
                                  method           = "MOFA", 
                                  transformData    = TRUE
    )
    
    # ---- Equivalence after preparation : ----
    expect(is(MAE2, "MOFA"), failure_message = "Prepared MAE is not a MOFA object")
    expect_equal(get_dimensions(MAE2)$D, lengths(variableList))
    
    protMat4 <- protMat3[variableList[["protetest"]],]
    protRes <- MAE2@data$protetest$group1
    
    expect_equal(dim(protMat4), dim(protRes))
    expect(identical(rownames(protMat4), rownames(protRes), attrib.as.set = FALSE), 
           failure_message = "proteins rownames are not identical")
    expect(identical(colnames(protMat4), colnames(protRes), attrib.as.set = FALSE), 
           failure_message = "proteins colnames are not identical")
    expect_identical(protMat4, protRes)
    
    metMat4 <- metMat3[variableList[["metatest"]],]
    metaRes <- MAE2@data$metatest$group1
    
    expect_equal(dim(metMat4), dim(metaRes))
    expect(identical(rownames(metMat4), rownames(metaRes), attrib.as.set = FALSE), 
           failure_message = "metabolites rownames are not identical")
    expect(identical(colnames(metMat4), colnames(metaRes), attrib.as.set = FALSE), 
           failure_message = "metabolites colnames are not identical")
    expect_identical(metMat4, metaRes)
    
    # ---- Equivalence on results: ----
    
    MAE3 <- runOmicsIntegration(MAE, preparedObject = MAE2, 
                                method = "MOFA", scale_views = TRUE, maxiter = 1000)
    
    # equivalence
    mofaobject <- create_mofa(data = list("protetest" = protMat4, "metatest" = metMat4), 
                              extract_metadata = TRUE)
    data_opts  <- get_default_data_options(mofaobject)
    model_opts <- get_default_model_options(mofaobject)
    train_opts <- get_default_training_options(mofaobject)
    
    data_opts$scale_views  <- TRUE
    train_opts$maxiter     <- 1000
    train_opts$verbose     <- FALSE
    model_opts$num_factors <- 10
    MOFAObject.untrained <- prepare_mofa(
        object           = mofaobject,
        data_options     = data_opts,
        model_options    = model_opts,
        training_options = train_opts
    )
    
    MOFAObject.trained <-
        run_mofa(MOFAObject.untrained,
                 use_basilisk = FALSE,
                 save_data = TRUE)
    
    resRFLOMICS <- get_factors(getMOFA(MAE3))$group1
    resEquivalence <- get_factors(MOFAObject.trained)$group1
    
    resRFLOMICSW <- get_weights(getMOFA(MAE3))$group1
    resEquivalenceW <- get_weights(MOFAObject.trained)$group1
    
    expect_equal(resRFLOMICSW, resEquivalenceW)
    expect_equal(resRFLOMICS, resEquivalence, tolerance = 10^-5)
})

