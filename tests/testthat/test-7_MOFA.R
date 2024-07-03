library(testthat)
library(RFLOMICS)
library(MOFA2)

# ---- Construction of objects for the tests ----

## ---- Construction MAE RFLOMICS ready for integration analysis : ----
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

formulae <- generateModelFormulae( MAE) 
MAE <- setModelFormula(MAE, formulae[[1]])

contrastList <- generateExpressionContrast(object = MAE) |> 
    purrr::reduce(rbind) |>
    dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                  "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
                                  "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))
MAE <- MAE |>
    setSelectedContrasts(contrastList) |>
    runTransformData(SE.name = "metatest", transformMethod = "log2") |>
    runNormalization(SE.name = "metatest", normMethod = "median")  |>
    runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")     |>
    runNormalization(SE.name = "protetest", normMethod = "median")  |>
    filterLowAbundance(SE.name = "RNAtest")                          |>
    runNormalization(SE.name = "RNAtest", normMethod = "TMM")        |>
    runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit")

MAE0 <- MAE

# ---- Comparison data: ----

protMat <- ecoseed$protetest
metMat <- ecoseed$metatest
condMat <- ecoseed$design

condMat$Repeat      <- factor(condMat$Repeat, 
                              levels = c("rep1", "rep2", "rep3"))
condMat$imbibition  <- factor(condMat$imbibition, 
                              levels = c("DS", "EI", "LI"))
condMat$temperature <- factor(condMat$temperature, 
                              levels = c("Low", "Medium", "Elevated"))

condMat$Repeat      <- relevel(condMat$Repeat, ref = "rep1")
condMat$imbibition  <- relevel(condMat$imbibition, ref = "DS")
condMat$temperature <- relevel(condMat$temperature, ref = "Low")

orderNames <- rownames(colData(MAE))
condMat <- condMat[match(orderNames, rownames(condMat)),]
protMat <- protMat[match(orderNames, colnames(protMat))]
metMat  <- metMat[match(orderNames, colnames(metMat))]

identical(orderNames, colnames(protMat), attrib.as.set = FALSE)
identical(orderNames, colnames(metMat), attrib.as.set = FALSE)

metMat2  <- apply(log2(metMat + 10^(-10)), 2, 
                  FUN = function(vect) vect - median(vect))
protMat2 <- apply(protMat, 2, FUN = function(vect) vect - median(vect))

# Select data
selectedData <- c("metatest", "protetest")
selectMethod <- c("metatest" = "diff", "protetest" = "none")
operation <- c("metatest" = "union", "protetest" = "union")

variableList <- 
    lapply(selectedData, function(set) {
        switch(
            selectMethod[[set]],
            "diff"  = getDEList(
                object = MAE[[set]],
                contrasts = getSelectedContrasts(MAE)$contrastName,
                operation = operation[[set]]
            ),
            "none"  = names(MAE[[set]])
        )
    })
names(variableList) <- selectedData


protMat2 <- protMat2[variableList[["protetest"]],]
metMat2 <- metMat2[variableList[["metatest"]],]

# Transform
metMat3 <- t(scale(t(metMat2), center = TRUE, scale = TRUE))
protMat3 <- t(scale(t(protMat2), center = TRUE, scale = TRUE))
# protMat3 <- protMat2
# metMat3 <- metMat2

# from .rbeFunction
colBatch <- getBatchFactors(MAE)
newFormula <-
    gsub(pattern = paste(paste(colBatch, "[+]"),  collapse = "|"),
         "", getModelFormula(MAE))
colData <- getDesignMat(MAE)
designToPreserve <-
    model.matrix(as.formula(newFormula), data = colData)

metMat4 <- limma::removeBatchEffect(metMat3, batch = colData[, colBatch],
                                    design = designToPreserve)
protMat4 <- limma::removeBatchEffect(protMat3, batch = colData[, colBatch],
                                     design = designToPreserve)

# ----- TESTS -----

test_that("Equivalence", {
    
    MAE2 <- prepareForIntegration(object           = MAE,
                                  omicsNames       = selectedData,
                                  variableLists    = variableList,
                                  method           = "MOFA", 
                                  transformData    = TRUE
    )
    
    # ---- Equivalence after preparation : ----
    expect(is(MAE2, "MOFA"), failure_message = "Prepared MAE is not a MOFA object")
    expect_equal(get_dimensions(MAE2)$D, lengths(variableList))
    
    protRes <- MAE2@data$protetest$group1
    
    expect_equal(dim(protMat4), dim(protRes))
    expect(identical(rownames(protMat4), rownames(protRes),
                     attrib.as.set = FALSE), 
           failure_message = "proteins rownames are not identical")
    expect(identical(colnames(protMat4), colnames(protRes),
                     attrib.as.set = FALSE), 
           failure_message = "proteins colnames are not identical")
    expect_identical(as.data.frame(protMat4), as.data.frame(protRes))
    
    metaRes <- MAE2@data$metatest$group1
    
    expect_equal(dim(metMat4), dim(metaRes))
    expect(identical(rownames(metMat4), rownames(metaRes), 
                     attrib.as.set = FALSE), 
           failure_message = "metabolites rownames are not identical")
    expect(identical(colnames(metMat4), colnames(metaRes),
                     attrib.as.set = FALSE), 
           failure_message = "metabolites colnames are not identical")
    expect_identical(as.data.frame(metMat4), as.data.frame(metaRes))
    
    # ---- Equivalence on results: ----
    
    MAE3 <- runOmicsIntegration(MAE, preparedObject = MAE2, 
                                method = "MOFA", scale_views = TRUE, 
                                maxiter = 1000, num_factors = 5)
    
    # equivalence
    mofaobject <- create_mofa(data = list("protetest" = protMat4, 
                                          "metatest" = metMat4), 
                              extract_metadata = TRUE)
    data_opts  <- get_default_data_options(mofaobject)
    model_opts <- get_default_model_options(mofaobject)
    train_opts <- get_default_training_options(mofaobject)
    
    data_opts$scale_views  <- TRUE
    train_opts$maxiter     <- 1000
    train_opts$verbose     <- FALSE
    model_opts$num_factors <- 5
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

