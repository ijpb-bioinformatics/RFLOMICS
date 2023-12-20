# ---- Create Example MAE ----
#' @title Initializing an Example MultiAssayExperiment compatible with RFLOMICS functions.
#'
#' @return a multiAssayExperiment ready for use in RFLOMICS
#' @export
#' @rdname Example-functions
initExampleMAE <- function(){
  
  ExpDesign <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), 
                                             "/ExamplesFiles/ecoseed/condition.txt"))
  factorRef <- data.frame(factorName  = c("Repeat", "temperature" , "imbibition"),
                          factorRef   = c("rep1",   "Low",          "DS"),
                          factorType  = c("batch",  "Bio",          "Bio"),
                          factorLevels= c("rep1,rep2,rep3", 
                                          "Low,Medium,Elevated", 
                                          "DS,EI,LI"))
  
  omicsData <- list(
    read_omics_data(file = paste0(system.file(package = "RFLOMICS"), 
                                  "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
    read_omics_data(file = paste0(system.file(package = "RFLOMICS"), 
                                  "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
    read_omics_data(file = paste0(system.file(package = "RFLOMICS"), 
                                  "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")))
  
  MAE <- FlomicsMultiAssay.constructor(projectName = "Tests", 
                                       omicsData   = omicsData,
                                       omicsNames  = c("RNAtest", "metatest", "protetest"),
                                       omicsTypes  = c("RNAseq","metabolomics","proteomics"),
                                       ExpDesign   = ExpDesign,
                                       factorRef   = factorRef)
  names(MAE) <- c("RNAtest", "metatest", "protetest")
  
  return(MAE)  
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
#' initExampleMAE()
#' MAE <- generateExample(coexp = FALSE, integration = FALSE)
#' 


generateExample <- function(processing   = TRUE,
                            diffanalysis = TRUE,
                            coexp        = TRUE,
                            annotation   = TRUE,
                            integration  = TRUE){
  
  MAE <- initExampleMAE()
  
  formulae <- GetModelFormulae(MAE = MAE) 
  contrastList <- rbind(getPossibleContrasts(MAE, 
                                       formula = formulae[[1]], 
                                       typeContrast = "simple",
                                       returnTable = TRUE)[1:3,],
                        getPossibleContrasts(MAE, 
                                             formula = formulae[[1]], 
                                             typeContrast = "averaged",
                                             returnTable = TRUE)[1:3,],
                        getPossibleContrasts(MAE, 
                                             formula = formulae[[1]], 
                                             typeContrast = "interaction",
                                             returnTable = TRUE)[1:3,])
  
  if (processing) {
    MAE <- MAE |>
      TransformData(     SE.name = "metatest",  transformMethod = "log2") |>
      RunNormalization(  SE.name = "metatest",  NormMethod = "totalSum")  |>
      RunNormalization(  SE.name = "RNAtest",   NormMethod = "TMM")       |>
      RunNormalization(  SE.name = "protetest", NormMethod = "median")    |>
      FilterLowAbundance(SE.name = "RNAtest")                                    
  }
  
  if (diffanalysis) {
    MAE <- MAE |>
      RunDiffAnalysis(SE.name = "metatest",  
                      DiffAnalysisMethod = "limmalmFit", 
                      contrastList = contrastList, 
                      modelFormula = formulae[[1]])  |>
      RunDiffAnalysis(SE.name = "protetest", 
                      DiffAnalysisMethod = "limmalmFit", 
                      contrastList = contrastList, 
                      modelFormula = formulae[[1]])  |>
      RunDiffAnalysis(SE.name = "RNAtest",   
                      DiffAnalysisMethod = "edgeRglmfit", 
                      contrastList = contrastList, 
                      modelFormula = formulae[[1]])  |>
      FilterDiffAnalysis(SE.name = "RNAtest",   
                         Adj.pvalue.cutoff = 0.05, 
                         logFC.cutoff = 1.5)
  }
  
  if (coexp) {
    MAE <- MAE |> 
      runCoExpression(SE.name = "RNAtest",   K = 2:12, replicates = 2) |>
      runCoExpression(SE.name = "protetest", K = 2:12, replicates = 2) |>
      runCoExpression(SE.name = "metatest",  K = 2:12, replicates = 2)
  }
  
  if (annotation) {
    df_custom <- vroom(file = paste0(system.file(package = "RFLOMICS"), 
                                            "/ExamplesFiles/GO_annotations/Arabidopsis_thaliana_Ensembl_55.txt"))
    MAE <- MAE |> 
      runAnnotationEnrichment(SE.name = "RNAtest", 
                              ontology = "custom",
                              list_args = list(pvalueCutoff = 0.05),
                              col_term = "GO term accession", 
                              col_gene = "Gene stable ID",
                              col_name = "GO term name",
                              col_domain = "GO domain",
                              annot = df_custom) |>
      runAnnotationEnrichment(SE.name = "protetest", 
                              ontology = "custom",
                              list_args = list(pvalueCutoff = 0.05),
                              col_term = "GO term accession", 
                              col_gene = "Gene stable ID",
                              col_name = "GO term name",
                              col_domain = "GO domain",
                              annot = df_custom) 
    
    
  }
  
  if (annotation && coexp) {
    df_custom <- vroom(file = paste0(system.file(package = "RFLOMICS"), 
                                     "/ExamplesFiles/GO_annotations/Arabidopsis_thaliana_Ensembl_55.txt"))
    MAE <- MAE |> 
      runAnnotationEnrichment(SE.name = "RNAtest", 
                              ontology = "custom",
                              from = "coexp",
                              list_args = list(pvalueCutoff = 0.05),
                              col_term = "GO term accession", 
                              col_gene = "Gene stable ID",
                              col_name = "GO term name",
                              col_domain = "GO domain",
                              annot = df_custom) |>
      runAnnotationEnrichment(SE.name = "protetest", 
                              ontology = "custom",
                              from = "coexp",
                              list_args = list(pvalueCutoff = 0.05),
                              col_term = "GO term accession", 
                              col_gene = "Gene stable ID",
                              col_name = "GO term name",
                              col_domain = "GO domain",
                              annot = df_custom) 
    
  }
  
  if (integration) {
    
    MAE <- integrationWrapper(MAE, 
                              omicsToIntegrate = c("RNAtest", "metatest", "protetest"))
    MAE <- integrationWrapper(MAE, 
                              omicsToIntegrate = c("RNAtest", "metatest", "protetest"),
                              method = "mixOmics")
  }
  
  return(MAE)
  
}