# ---- Create Example MAE ----
#' @title Initializing an Example MultiAssayExperiment compatible with RFLOMICS
#'  functions.
#'
#' @return a multiAssayExperiment ready for use in RFLOMICS
#' @export
#' @rdname Example-functions
#' @examples
#' initExampleMAE()
#'
initExampleMAE <- function() {
    datPath <- paste0(system.file(package = "RFLOMICS"),
                      "/ExamplesFiles/ecoseed/")
    
    ExpDesign <-
        readExpDesign(file = paste0(datPath, "condition.txt"))
    facRef <-
        data.frame(
            factorName   = c("Repeat", "temperature" , "imbibition"),
            factorRef    = c("rep1",   "Low",          "DS"),
            factorType   = c("batch",  "Bio",          "Bio"),
            factorLevels = c("rep1,rep2,rep3",
                             "Low,Medium,Elevated",
                             "DS,EI,LI")
        )
    
    omicsData <- list(
        readOmicsData(file = paste0(
            datPath, "transcriptome_ecoseed.txt"
        )),
        readOmicsData(file = paste0(datPath, "metabolome_ecoseed.txt")),
        readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt"))
    )
    
    RMAE <- createRflomicsMAE(
        projectName = "Tests",
        omicsData   = omicsData,
        omicsNames  = c("RNAtest", "metatest",
                        "protetest"),
        omicsTypes  = c("RNAseq", "metabolomics",
                        "proteomics"),
        ExpDesign   = ExpDesign,
        factorRef   = facRef
    )
    names(RMAE) <- c("RNAtest", "metatest", "protetest")
    
    return(RMAE)
}

#' @title singleOmicsExample
#'
#' @return a Rflomics SummarizedExperiment, result of some RFLOMICS analyses
#' @importFrom vroom vroom
#' @export
#' @rdname Example-functions
#' @examples
#' singleOmicsExample(omicType = "RNAseq")
#'
singleOmicsExample <- function(omicType = "RNAseq") {
    if (missing(omicType))
        stop("Please provide an omic type.")
    
    datPath <- paste0(system.file(package = "RFLOMICS"),
                      "/ExamplesFiles/ecoseed/")
    
    ExpDesign <-
        readExpDesign(file = paste0(datPath, "condition.txt"))
    
    design <- c(
        "Repeat" = "batch",
        "temperature" = "bio",
        "imbibition" = "bio"
    )
    
    omicsData <-
        switch(
            toupper(omicType),
            "RNASEQ" = {
                readOmicsData(file = paste0(datPath, "transcriptome_ecoseed.txt"))
            },
            "PROTEOMICS" = {
                readOmicsData(file = paste0(datPath, "metabolome_ecoseed.txt"))
            },
            "METABOLOMICS" = {
                readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt"))
            }
        )
    
    RSE <- createRflomicsSE(
        omicData  = omicsData,
        omicType  = omicType,
        ExpDesign = ExpDesign,
        design    = design
    )
    
    return(RSE)
}


#' @title generateExample
#'
#' @return a multiAssayExperiment, results of some RFLOMICS analyses
#' @importFrom vroom vroom
#' @export
#' @rdname Example-functions
#' @importFrom purrr reduce
#' @importFrom dplyr filter
#' @examples
#' MAE <- generateExample(coexp = FALSE, annotation = FALSE, integration = FALSE)
#'


generateExample <- function(processing   = TRUE,
                            diffanalysis = TRUE,
                            coexp        = TRUE,
                            annotation   = TRUE,
                            integration  = TRUE) {
    MAE <- initExampleMAE()
    
    formulae <- generateModelFormulae(MAE)
    MAE <- setModelFormula(MAE, formulae[[1]])
    
    contrastList <- rbind(
        getPossibleContrasts(
            MAE,
            formula = formulae[[1]],
            typeContrast = "simple",
            returnTable = TRUE
        )[c(1, 2, 3),],
        getPossibleContrasts(
            MAE,
            formula = formulae[[1]],
            typeContrast = "averaged",
            returnTable = TRUE
        )[c(1, 2, 3),],
        getPossibleContrasts(
            MAE,
            formula = formulae[[1]],
            typeContrast = "interaction",
            returnTable = TRUE
        )[c(1, 2, 3),]
    )
    
    if (processing) {
        MAE <- MAE |>
            runTransformData(SE.name = "metatest",  transformMethod = "log2") |>
            runNormalization(SE.name = "metatest",  normMethod = "totalSum")  |>
            runNormalization(SE.name = "RNAtest",   normMethod = "TMM")       |>
            runNormalization(SE.name = "protetest", normMethod = "median")    |>
            filterLowAbundance(SE.name = "RNAtest")
    }
    
    if (diffanalysis) {
        MAE <- MAE |>
            runDiffAnalysis(SE.name = "metatest",
                            method = "limmalmFit",
                            contrastList = contrastList)  |>
            runDiffAnalysis(SE.name = "protetest",
                            method = "limmalmFit",
                            contrastList = contrastList)  |>
            runDiffAnalysis(SE.name = "RNAtest",
                            method = "edgeRglmfit",
                            contrastList = contrastList)  |>
            filterDiffAnalysis(
                SE.name = "RNAtest",
                p.adj.cutoff = 0.05,
                logFC.cutoff = 1.5
            )
    }
    
    if (coexp) {
        MAE <- MAE |>
            runCoExpression(SE.name = "RNAtest",
                            K = 2:12,
                            replicates = 2) |>
            runCoExpression(SE.name = "protetest",
                            K = 2:12,
                            replicates = 2) |>
            runCoExpression(SE.name = "metatest",
                            K = 2:12,
                            replicates = 2)
    }
    
    if (annotation) {
        annotPath <- paste0(
            system.file(package = "RFLOMICS"),
            "/ExamplesFiles/ecoseed/AT_GOterm_EnsemblPlants.txt"
        )
        
        df_custom <- vroom(file = annotPath)
        MAE <- MAE |>
            runAnnotationEnrichment(
                SE.name = "RNAtest",
                nameList = getSelectedContrasts(MAE[["RNAtest"]])$tag,
                database = "custom",
                list_args = list(pvalueCutoff = 0.05),
                col_term = "GO term accession",
                col_gene = "Gene stable ID",
                col_name = "GO term name",
                col_domain = "GO domain",
                annot = df_custom
            ) |>
            runAnnotationEnrichment(
                SE.name = "protetest",
                nameList = getSelectedContrasts(MAE[["protetest"]])$tag,
                database = "custom",
                list_args = list(pvalueCutoff = 0.05),
                col_term = "GO term accession",
                col_gene = "Gene stable ID",
                col_name = "GO term name",
                col_domain = "GO domain",
                annot = df_custom
            )
        
        
    }
    
    
    if (annotation && coexp) {
        df_custom <- vroom(file = annotPath)
        MAE <- MAE |>
            runAnnotationEnrichment(
                SE.name = "RNAtest",
                database = "custom",
                from = "coexp",
                list_args = list(pvalueCutoff = 0.05),
                col_term = "GO term accession",
                col_gene = "Gene stable ID",
                col_name = "GO term name",
                col_domain = "GO domain",
                annot = df_custom
            ) |>
            runAnnotationEnrichment(
                SE.name = "protetest",
                database = "custom",
                from = "coexp",
                list_args = list(pvalueCutoff = 0.05),
                col_term = "GO term accession",
                col_gene = "Gene stable ID",
                col_name = "GO term name",
                col_domain = "GO domain",
                annot = df_custom
            )
        
    }
    
    if (integration) {
        MAE <- integrationWrapper(MAE,
                                  omicsNames = c("RNAtest", "metatest",
                                                 "protetest"))
        MAE <- integrationWrapper(
            MAE,
            omicsNames = c("RNAtest", "metatest",
                           "protetest"),
            method = "mixOmics"
        )
    }
    
    return(MAE)
    
}