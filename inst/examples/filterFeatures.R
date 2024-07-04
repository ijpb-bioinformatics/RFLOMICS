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

contrastList <- generateExpressionContrast(object = MAE) |> 
    purrr::reduce(rbind) |>
    dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
                                  "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
                                  "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))
MAE <- MAE |>
    setSelectedContrasts(contrastList) |>
    runTransformData(SE.name = "metatest", transformMethod = "log2") |>
    runNormalization(SE.name = "metatest", normMethod = "median")    |>
    runNormalization(SE.name = "protetest", normMethod = "median")   |>
    runDiffAnalysis(SE.name = "metatest", method = "limmalmFit")     |>
    runDiffAnalysis(SE.name = "protetest", method = "limmalmFit")    |>
    runCoExpression(SE.name = "protetest", K = 2:10, replicates = 5, 
                    merge = "union")

# Define the selections options
selOpt = list("protetest" = c("cluster.1", "H3"), "metatest" = c("DE"))
MAE1 <- filterFeatures(MAE, selOpt)
MAE1

selOpt2 = list("protetest" = c("cluster.2", "H2", "H1"), metatest = c("DE"))
MAE2 <- filterFeatures(MAE, selOpt2,
                       type = c("metatest" = "intersection",
                                "protetest" = 'union'))
MAE2
