# load ecoseed data
library(RFLOMICS)
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
contrastList <- Reduce(rbind, generateExpressionContrast(MAE)) 

MAE <- MAE |>
    setSelectedContrasts(contrastList[c(3,6,25)]) |>
    runTransformData(SE.name = "metatest", transformMethod = "log2") |>
    runNormalization(SE.name = "metatest", normMethod = "median")    |>
    runNormalization(SE.name = "protetest", normMethod = "median")   |>
    runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")     |>
    runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")    

# Integration using MOFA
# Prepare mofa object:
mofaObj <- prepareForIntegration(MAE,
                                 omicsNames = c("protetest", "metatest"),
                                 variableLists = rownames(MAE), 
                                 method = "MOFA")

# Perform integration:
MAEtest <- runOmicsIntegration(MAE, 
                               preparedObject = mofaObj, 
                               method = "MOFA", num_factors = 5)

# Integration using MixOmics
mixObj <- prepareForIntegration(MAE,
                                omicsNames = c("protetest", "metatest"),
                                variableLists = rownames(MAE),
                                method = "mixOmics")
MAEtest <- runOmicsIntegration(MAEtest, preparedObject = mixObj, 
                               method = "mixOmics")


# Access mixOmics results:
getMixOmics(MAEtest, response = "temperature")
getMixOmicsSettings(MAEtest)
mixOmics::plotIndiv(getMixOmics(MAEtest, response = "imbibition"))

# Access MOFA2 results:
getMOFA(MAEtest)
getMOFASettings(MAEtest)
MOFA2::plot_variance_explained(getMOFA(MAEtest))

