
################################### DIFF-ANALYSIS #############################


###### Statistical METHOD

## METHOD to perform differential analysis

#' @title runDiffAnalysis
#' @description This is an interface method which run a differential 
#' analysis on omics datasets stored in an object of class 
#' \link{RflomicsSE} or \link{RflomicsMAE}. According to the type 
#' of omics and to a list of contrasts, a differential analysis  
#' is performef for each contrasts.
#' Two methods are available according to the type of object:
#' \itemize{
#' \item{For RNAseq data: }{the \code{glmFit} function/model of the 
#' \code{edgeR} package} is applied.
#' \item{For proteomic and metabolomic data: }{the \code{\link{lmFit}} 
#' function/model of the \code{limma} package} is applied.
#' }
#' Functions and parameters used for RNAseq are those recommended in DiCoExpress 
#' workflow (see the paper in reference).
#' Functions and parameters used for proteomic and metabolomic data are those 
#' recommended in the 
#' @return
#' All the results are stored as a named list \code{DiffExpAnal} 
#' in the metadata slot of a given \link{RflomicsSE} object.
#' 
#' Objects are:
#' \itemize{
#' \item{stats:}{ data.frame giving a summary of the differential statistical analysis results by contrast:
#'  number of DE features, number of up and down regulated features}
#' 
#' \item{setting: }{Parameters used for the differential analysis}
#' \itemize{
#' \item{method: }{ The method used for the differential analysis}
#' \item{p.adj.method:}{ The applied p-value correction method}
#' \item{p.adj.cutoff:}{ The cut-off applied for the adjusted p-value}
#' \item{logFC.cutoff:}{ The absolute log FC cut-off}
#' }
#' 
#' \item{RawDEFres: }{a list giving for each contrast the raw results of 
#' the differential analysis method}
#' \item{DEF: }{a list giving for each contrast a data.frame of non filtered 
#' differential expressed features with their statistics}
#' \item{TopDEF: }{a list giving for each contrast a data.frame of 
#' differential expressed features ordered and filtered by p.adj.cutoff 
#' with their statistics}
#' \item{mergeDEF: }{a data frame of 0/1 indicating for each features in row, 
#' if it is DE in a given contrasts in column}
#' \item{contrasts: }{a data.table of the contrasts used for the differential 
#' analysis}
#' }
#' @param object an object of class \link{RflomicsSE} or \link{RflomicsMAE} 
#' @param SE.name the name of the dataset to fetch if the object is 
#' of class \link{RflomicsMAE}
#' @param method A character vector giving the name of the differential 
#' analysis method to run. Either "edgeRglmfit" or "limmalmFit".
#' @param contrastList data.frame of contrast from getExpressionContrastF()
#' @param p.adj.method The method choosen to adjust pvalue. Takes the same 
#' values as the ones of adj.p.adjust method.
#' @param p.adj.cutoff The adjusted pvalue cut-off
#' @param logFC.cutoff The lof2FC cutoff
#' @param clustermq A boolean indicating whether the constrasts have to be 
#' computed in local or in a distant machine
#' @param parallel boolean. If TRUE, uses mclapply function from the parallel
#' package. 
#' @param nworkers integers. The number of CPU to use.
#' @param cmd boolean. Used in the interface. If TRUE, print cmd for the user.
#' @return An object of class \link{RflomicsSE}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool 
#' to process multifactorial RNAseq experiments from quality controls to 
#' co-expression analysis through differential analysis based on contrasts 
#' inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runDiffAnalysis
#' @importFrom dplyr filter
#' @rdname runDiffAnalysis
#' @aliases runDiffAnalysis
#' @seealso \code{\link{filterDiffAnalysis}}
#' @seealso \code{\link{getDEMatrix}}
#' @seealso \code{\link{getDEList}}
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")), 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest", "protetest"),
#'                          omicsTypes  = c("RNAseq","proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest", "protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  filterLowAbundance(SE.name = "RNAtest")                           |>    
#'  runNormalization(SE.name = "RNAtest",   normMethod = "TMM")       |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |>
#'  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit", 
#'                  contrastList = contrastList)  |>
#'  runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit", 
#'                  contrastList = contrastList) 
#'  
#'  # Access to the diff analysis settings
#'  getDiffSettings(experiments(MAE)[["RNAtest"]])
#'  getDiffSettings(experiments(MAE)[["protetest"]])
#' 
#'  getDEList(experiments(MAE)[["RNAtest"]]) 

setMethod(f         = "runDiffAnalysis",
          signature = "RflomicsSE",
          definition = function(object, p.adj.method="BH",
                                contrastList = NULL, method = NULL, 
                                p.adj.cutoff=0.05, logFC.cutoff=0, 
                                clustermq=FALSE, parallel = FALSE, 
                                nworkers = 1,
                                cmd = FALSE){
              
              modelFormula <- getModelFormula(object)
              
              # check args
              if (is.null(contrastList) || nrow(contrastList) == 0) {
                  stop("contrastList arg is mandatory.")
              }
              
              if (is.null(method) || isFALSE(method %in% c("edgeRglmfit", "limmalmFit"))) {
                  switch(getOmicsTypes(object),
                         "RNAseq" = { method <- "edgeRglmfit"},
                         { method <- "limmalmFit" }
                  )
                  warning("DiffAnalyseMethod was missing. Detected omic type is ", 
                          getOmicsTypes(object)," using ", 
                          method, " for differential analysis.")
              }
              
              ## check completness
              Completeness <- checkExpDesignCompleteness(object)
              if (isTRUE(Completeness[["error"]])) {
                  stop(Completeness[["messages"]])
              }
              
              ## getcontrast
              contrastList <- generateExpressionContrast(object) %>% reduce(rbind) %>% 
                  filter(contrast %in% contrastList$contrast)
              
              object <- generateContrastMatrix(object, contrastList = contrastList)
              Contrasts.Sel <- mutate(contrastList, 
                                      tag = paste0("H", seq_len(nrow(contrastList))))
              object <- setSelectedContrasts(object, Contrasts.Sel)
              
              object@metadata$DiffExpAnal <- list()
              object@metadata$DiffExpAnal[["contrasts"]] <- Contrasts.Sel
              
              # remplacera à terme les lignes ci-dessus
              object@metadata$DiffExpAnal[["setting"]][["method"]] <- method
              object@metadata$DiffExpAnal[["setting"]][["p.adj.method"]] <- p.adj.method
              object@metadata$DiffExpAnal[["setting"]][["p.adj.cutoff"]] <- p.adj.cutoff
              object@metadata$DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
              
              # transform and norm if needed
              if (method == "limmalmFit") {
                  object2 <- object
                  
                  if (!.isTransformed(object2)) {
                      object2 <- .applyTransformation(object2)
                  }
                  if (!.isNorm(object2)) { 
                      object2 <- .applyNorm(object2)
                  }
              }
              
              # move in ExpDesign Constructor
              model_matrix <- model.matrix(as.formula(paste(modelFormula, collapse = " ")), data = getDesignMat(object))
              
              ListRes <- switch(
                  method,
                  "edgeRglmfit" = .tryRflomics(.edgeRAnaDiff(count_matrix  = assay(object),
                                                             model_matrix    = model_matrix[colnames(object),],
                                                             group           = getCoeffNorm(object)$group,
                                                             lib.size        = getCoeffNorm(object)$lib.size,
                                                             norm.factors    = getCoeffNorm(object)$norm.factors,
                                                             Contrasts.Sel   = Contrasts.Sel,
                                                             Contrasts.Coeff = object@metadata$design$Contrasts.Coeff,
                                                             FDR             = 1,
                                                             clustermq       = clustermq,
                                                             parallel        = parallel,
                                                             nworkers        = nworkers,
                                                             cmd             = cmd)),
                  "limmalmFit" = .tryRflomics(.limmaAnaDiff(count_matrix    = assay(object2),
                                                            model_matrix      = model_matrix[colnames(object2),],
                                                            Contrasts.Sel     = Contrasts.Sel,
                                                            Contrasts.Coeff   = object@metadata$design$Contrasts.Coeff,
                                                            p.adj.cutoff = 1,
                                                            p.adj.method = p.adj.method,
                                                            clustermq         = clustermq,
                                                            cmd               = cmd)))
              
              if (!is.null(ListRes$value)) {
                  if (!is.null(ListRes$value[["RawDEFres"]])) {
                      object@metadata$DiffExpAnal[["results"]] <- TRUE
                      object@metadata$DiffExpAnal[["RawDEFres"]] <- ListRes$value[["RawDEFres"]]
                      object@metadata$DiffExpAnal[["DEF"]] <- ListRes$value[["TopDEF"]]
                  }else{
                      object@metadata$DiffExpAnal[["results"]]    <- FALSE
                      object@metadata$DiffExpAnal[["ErrorStats"]] <- ListRes$value[["ErrorTab"]]
                  }
              }else{
                  object@metadata$DiffExpAnal[["results"]]    <- FALSE
                  object@metadata$DiffExpAnal[["Error"]]      <- ListRes$error
                  object@metadata$DiffExpAnal[["ErrorStats"]] <- NULL
                  
                  
                  #return(object)
              }
              
              ## filtering
              object <-  filterDiffAnalysis(object = object, 
                                            p.adj.cutoff = p.adj.cutoff, 
                                            logFC.cutoff = logFC.cutoff)
              
              return(object)
          })


#' @rdname runDiffAnalysis
#' @title runDiffAnalysis
#' @aliases runDiffAnalysis
#' @exportMethod runDiffAnalysis
#' 
setMethod(f          = "runDiffAnalysis",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, p.adj.method="BH",
                                contrastList = NULL, method = NULL,
                                p.adj.cutoff=0.05, 
                                logFC.cutoff=0, clustermq=FALSE, 
                                parallel = FALSE, nworkers = 1, 
                                cmd = FALSE){
              
              # all verifications are done in this method
              object[[SE.name]] <-  runDiffAnalysis(object = object[[SE.name]],
                                                    p.adj.method = p.adj.method,
                                                    contrastList = contrastList,
                                                    method = method,
                                                    p.adj.cutoff = p.adj.cutoff,
                                                    logFC.cutoff = logFC.cutoff,
                                                    clustermq = clustermq,
                                                    parallel = parallel,
                                                    nworkers = nworkers,
                                                    cmd = cmd
              )
              return(object)
          })

## METHOD to filter differential analysis

#' @title filterDiffAnalysis
#' Filter differential analysis
#'
#' @param object A RflomicsSE object
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a \link{RflomicsMAE}
#' @param p.adj.cutoff adjusted pvalue cutoff. Default is the parameter from 
#' the differential analysis.
#' @param logFC.cutoff cutoff for absolute value of log2FC. Default is the 
#' parameter from the differential analysis. 
#' @return A \link{RflomicsSE} object or a \link{RflomicsMAE}, depending on the object type, 
#' where the differential analysis results have been actualized with the 
#' new parameters.
#' @exportMethod filterDiffAnalysis
#' @rdname filterDiffAnalysis
#' @aliases filterDiffAnalysis
#' @importFrom dplyr full_join filter if_else mutate_at
#' @importFrom data.table data.table
#' @importFrom purrr reduce
#' @examples
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")), 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest", "protetest"),
#'                          omicsTypes  = c("RNAseq","proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest", "protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  filterLowAbundance(SE.name = "RNAtest")                           |>    
#'  runNormalization(SE.name = "RNAtest",   normMethod = "TMM")       |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |>
#'  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit", 
#'                  contrastList = contrastList)  |>
#'  runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit", 
#'                  contrastList = contrastList) |> 
#'  filterDiffAnalysis(SE.name = "protetest", p.adj.cutoff = 0.01, logFC.cutoff = 1) |> 
#'  filterDiffAnalysis(SE.name = "RNAtest", p.adj.cutoff = 0.01, logFC.cutoff = 1)
setMethod(f          = "filterDiffAnalysis",
          signature  = "RflomicsSE",
          definition <- function(object, 
                                 p.adj.cutoff = NULL, 
                                 logFC.cutoff = NULL){
              
              if (is.null(object@metadata$DiffExpAnal[["RawDEFres"]])) {
                  stop("can't filter the DiffExpAnal object because it doesn't exist")
              }
              
              if (is.null(p.adj.cutoff)) 
                  p.adj.cutoff <- getDiffSettings(object)[["p.adj.cutoff"]]
              
              if (is.null(logFC.cutoff))
                 logFC.cutoff <- getDiffSettings(object)[["abs.logFC.cutoff"]]
              
              # remplacera à terme les lignes ci-dessus
              metadata(object)$DiffExpAnal[["setting"]][["p.adj.cutoff"]] <- p.adj.cutoff
              metadata(object)$DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
              
              ## TopDEF: Top differential expressed features
              DEF_filtred <- lapply(seq_len(length(object@metadata$DiffExpAnal[["DEF"]])), function(x){
                  res <- object@metadata$DiffExpAnal[["DEF"]][[x]]
                  keep <- (res$Adj.pvalue < p.adj.cutoff) & (abs(res$logFC) > logFC.cutoff)
                  res <- res[keep,]
                  return(res)
              })
              names(DEF_filtred) <- names(object@metadata$DiffExpAnal[["RawDEFres"]])
              object@metadata$DiffExpAnal[["TopDEF"]] <- DEF_filtred
              
              ## stats
              object@metadata$DiffExpAnal[["stats"]] <- sumDiffExp(object)
              
              ## merge results in bin matrix
              DEF_list <- list()
              for (x in names(object@metadata$DiffExpAnal[["TopDEF"]])){
                  res <- object@metadata$DiffExpAnal[["TopDEF"]][[x]]
                  tmp <- data.frame(DEF = rownames(res), bin = rep(1,length(rownames(res))))
                  colnames(tmp) <- c("DEF", filter(object@metadata$DiffExpAnal$contrasts, contrastName == x)$tag)
                  
                  if(dim(tmp)[1] != 0){ DEF_list[[x]] <- tmp }
              }
              
              object@metadata$DiffExpAnal[["mergeDEF"]] <- NULL
              
              if (length(DEF_list) != 0) {
                  object@metadata$DiffExpAnal[["mergeDEF"]] <- DEF_list %>% 
                      reduce(full_join, by="DEF") %>%
                      mutate_at(.vars = 2:(length(DEF_list)+1),
                                .funs = function(x){
                                    if_else(is.na(x), 0, 1)}) %>%
                      data.table()
              }
              
              return(object)
          })

#' @rdname filterDiffAnalysis
#' @aliases filterDiffAnalysis
#' @title filterDiffAnalysis
#' @exportMethod filterDiffAnalysis
setMethod(f          = "filterDiffAnalysis",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, 
                                p.adj.cutoff = NULL, 
                                logFC.cutoff = NULL){
              
              if (!SE.name %in% names(object))
                  stop(SE.name, " isn't the name of an experiment in ", object)
              
              if (is.null(p.adj.cutoff) && is.null(logFC.cutoff)) {
                  
                  message("Parameter p.adj.cutoff and |logFC.cutoff| are both NULL. 
                          Not changing anything")
                  return(object)
                  
              }else{
                  
                  if (is.null(p.adj.cutoff)) 
                      p.adj.cutoff <- getDiffSettings(object[[SE.name]])$p.adj.cutoff
                  
                  if (is.null(logFC.cutoff))
                    logFC.cutoff <- getDiffSettings(object[[SE.name]])$abs.logFC.cutoff
                  
                  object[[SE.name]] <-  filterDiffAnalysis(object = object[[SE.name]],
                                                           p.adj.cutoff = p.adj.cutoff,
                                                           logFC.cutoff = logFC.cutoff)
                  
                  return(object)
              }
              
          })

###### Graphical METHOD

## Method to plot results of a differential analysis

#' @title plotDiffAnalysis
#' @description
#' This is an interface method which draw a MAplot, a volcano plot and the 
#' pvalues distribution from the results of a differential analysis
#' performed on omic datasets stored in an object of class \link{RflomicsSE}
#' @param object An object of class \link{RflomicsSE}
#' @param contrastName The contrastName for which the plots has to be drawn
#' @param typeofplots The plots you want to return. Default is all possible 
#' plots: MA plot, Volcano plot and non adjusted pvalues histogram.
#' @return plot
#' @exportMethod plotDiffAnalysis
#' @rdname plotDiffAnalysis
#' @aliases plotDiffAnalysis
#' @export
#' 
setMethod(f="plotDiffAnalysis",
          signature="RflomicsSE",
          definition <- function(object, 
                                 contrastName,
                                 typeofplots = c("MA.plot", "volcano", "histogram")){

              if (.isTagName(object, contrastName)) {
                  contrastName <- .convertTagToContrast(object, contrastName)
              }
              
              plots <- list()
              
              resTable <- object@metadata$DiffExpAnal[["DEF"]][[contrastName]]
              
              logFC.cutoff <- getDiffSettings(object)[["abs.logFC.cutoff"]]
              p.adj.cutoff <- getDiffSettings(object)[["p.adj.cutoff"]]
              
              if ("MA.plot" %in% typeofplots) plots[["MA.plot"]]        <-  .plotMA(data = resTable, p.adj.cutoff = p.adj.cutoff, logFC.cutoff = logFC.cutoff, contrastName=contrastName)
              if ("volcano" %in% typeofplots) plots[["Volcano.plot"]]   <-  .plotVolcanoPlot(data = resTable, p.adj.cutoff = p.adj.cutoff, logFC.cutoff = logFC.cutoff, contrastName=contrastName)
              if ("histogram" %in% typeofplots) plots[["Pvalue.hist"]]  <-  .plotPValue(data =resTable, contrastName=contrastName)
              return(plots)
          })

#' @rdname plotDiffAnalysis
#' @title plotDiffAnalysis
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @exportMethod plotDiffAnalysis
#' @aliases plotDiffAnalysis
setMethod(f          = "plotDiffAnalysis",
          signature  = "RflomicsMAE",
          definition = function(object,
                                SE.name, 
                                contrastName, 
                                typeofplots = c("MA.plot", "volcano", "histogram")){
              
              if (.isTagName(object, contrastName)) { 
                  contrastName <-  .convertTagToContrast(object, contrastName)
              }
              
              return(plotDiffAnalysis(object      = object[[SE.name]],
                                      contrastName  = contrastName,
                                      typeofplots = typeofplots))
              
          })


#' @title plotHeatmapDesign
#' @description
#' This mehod draws a heatmap from the results of a 
#' differential analysis performed on omic datasets stored in an object of 
#' class \link{RflomicsSE}
#' @param object An object of class \link{RflomicsSE}
#' @param contrastName The contrastName for which the MAplot has to be drawn
#' @param splitFactor characters. Default to none. Name of a feature in the 
#' design matrix, splits the samples on the heatmap according to its modalities.  
#' @param title characters. Title of the heatmap. 
#' @param annotNames vector. Names of the annotations to keep in the Heatmap. 
#' Default takes all available information.
#' @param modalities named list of vectors of modalities to subset and print 
#' on the heatmap. 
#' @param drawArgs,heatmapArgs  named lists. Any additional parameter passed 
#' to ComplexHeatmap::Heatmap or ComplexHeatmap::draw
#' @return plot
#' @exportMethod plotHeatmapDesign
#' @export
#' @importFrom dplyr arrange select
#' @importFrom tidyselect any_of
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importClassesFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importMethodsFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @rdname plotHeatmapDesign
#' @aliases plotHeatmapDesign
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"),
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"),omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |>
#'  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit", 
#'                   contrastList = contrastList)
#'  
#'  # plot the heatmap
#'  plotHeatmapDesign(experiments(MAE)[[1]], 
#'                    getSelectedContrasts(experiments(MAE)[[1]])$contrastName[1])

 
setMethod(f          = "plotHeatmapDesign",
          signature  = "RflomicsSE",
          definition = function(object, 
                                contrastName, 
                                splitFactor="none", 
                                title = "", 
                                annotNames = NULL, 
                                modalities = NULL, 
                                drawArgs = list(), 
                                heatmapArgs = list()){
              
              Groups     <- getDesignMat(object)
              
              if (is.null(object@metadata$DiffExpAnal[["TopDEF"]][[contrastName]])) {
                  stop("no DE variables")
              }
              
              resTable <- arrange(object@metadata$DiffExpAnal[["TopDEF"]][[contrastName]], Adj.pvalue)
              
              if (dim(resTable)[1] == 0) {
                  stop("no differentially expressed variables...")
              }
              
              if (dim(resTable)[1] > 2000) {
                  message("differentially expressed variables exceeding 2000 variables, only the first 2000 will be displayed")
                  resTable <- resTable[seq_len(2000),]
                  title <- ifelse(title == "", 
                                  paste0(title, "plot only 2000 TOP DE variables"),
                                  paste0(title, "\nplot only 2000 TOP DE variables"))
              }
              
              object2 <- .checkTransNorm(object, raw = FALSE)
              m.def  <- assay(object2)
              
              m.def <- as.data.frame(m.def) %>%
                  select(any_of(Groups$samples))
              
              # filter by DE
              m.def.filter <- subset(m.def, 
                                     rownames(m.def) %in% row.names(resTable))
              
              # normalize count
              
              # Center
              m.def.filter.center <- t(scale(t(m.def.filter), center = TRUE, scale = FALSE))
              
              # Annotations datatable
              df_annotation <- Groups %>% select(!samples & !groups)  
              df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),] 
              
              # Subset the dataset to print only interesting modalities
              if (!is.null(modalities)) {
                  if (is.null(names(modalities))) {
                      message("In heatmapPlot, modalities argument needs a named list. Not subsetting")
                  }else{ 
                      samplesToKeep <- Reduce("intersect", lapply(
                          seq_len(length(modalities)),
                          FUN = function(i){
                              col_nam <- names(modalities)[i]
                              rownames(df_annotation[which(df_annotation[[col_nam]] %in% modalities[[i]]),])
                          }
                      ))
                      
                      df_annotation <- df_annotation[which(rownames(df_annotation) %in% samplesToKeep),]
                      m.def.filter.center <- m.def.filter.center[, which(colnames(m.def.filter.center) %in% samplesToKeep)]
                      
                      df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),]
                  }
              }
              
              # Split management
              column_split.value <- if (splitFactor != "none") { df_annotation[, splitFactor] } else { NULL }
              
              # Select the right columns
              if (!is.null(annotNames)) {
                  df_annotation <- df_annotation %>% 
                      select(any_of(annotNames))
              }
              
              # Color annotations
              nAnnot <- ncol(df_annotation)
              selectPal <- rownames(brewer.pal.info)[seq_len(nAnnot)]
              
              color_list <- lapply(seq_len(nAnnot), 
                                   FUN = function(i){
                                       annot_vect <- unique(df_annotation[,i])
                                       
                                       col_vect <-  colorRampPalette(
                                           brewer.pal(n = min(length(annot_vect), 8), 
                                                      name = selectPal[i])
                                       )(length(annot_vect)) 
                                       names(col_vect) <- annot_vect 
                                       col_vect[!is.na(names(col_vect))] 
                                   })
              names(color_list) <- colnames(df_annotation)
              
              column_ha <- HeatmapAnnotation(df = df_annotation, 
                                             col = color_list)
              
              namArg <- ifelse(getOmicsTypes(object) == "RNAseq", 
                               "normalized counts", "XIC")
              
              # Arguments for Heatmap
              heatmapArgs <- c(
                  list(matrix = m.def.filter.center,
                       name = namArg,
                       show_row_names = ifelse(dim(m.def.filter.center)[1] > 50, FALSE, TRUE),
                       row_names_gp = gpar(fontsize = 8),
                       column_names_gp = gpar(fontsize = 12),
                       row_title_rot = 0 ,
                       clustering_method_columns = "ward.D2",
                       cluster_column_slice = FALSE,
                       column_split = column_split.value,
                       top_annotation = column_ha,
                       column_title = title),
                  heatmapArgs)
              
              # Arguments for drawing the heatmap
              drawArgs <- c(list(merge_legend = TRUE),
                            drawArgs)            
              
              # Drawing heatmap in a null file to not plot it
              pdf(file = NULL)
              ha <- do.call(Heatmap, heatmapArgs)
              
              drawArgs$object <- ha
              ha <- do.call(draw, drawArgs)
              
              dev.off()
              
              return(ha)
          })


#' @rdname plotHeatmapDesign
#' @title plotHeatmapDesign
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @exportMethod plotHeatmapDesign
#' @aliases plotHeatmapDesign
setMethod(f          = "plotHeatmapDesign",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, 
                                contrastName, 
                                splitFactor="none", 
                                title = "", annotNames = NULL, 
                                modalities = NULL, 
                                drawArgs = list(), 
                                heatmapArgs = list()){
              
              
              if (.isTagName(object, contrastName)) {
                  contrastName <- .convertTagToContrast(object, contrastName)
              }
              
              return(plotHeatmapDesign(object        = object[[SE.name]],
                                       contrastName    = contrastName,
                                       splitFactor     = splitFactor,
                                       title         = title,
                                       annotNames = annotNames,
                                       modalities   = modalities,
                                       drawArgs     = drawArgs,
                                       heatmapArgs  = heatmapArgs))
              
          })




#' @title plotBoxplotDE
#'
#' @param object An object of class \link{RflomicsSE}
#' @param features variable name (gene/protein/metabolite name)
#' @param groupColor default to groups, indicate a variable in the design to
#' color the boxplots accordingly. 
#' @param raw Boolean. Plot the raw data or the transformed ones (TRUE)
#' @return a ggplot graph
#' @exportMethod plotBoxplotDE
#' @importFrom ggplot2  ggplot aes geom_boxplot  element_text theme guides 
#' guide_legend xlab ylab theme_void ggtitle
#' @importFrom dplyr full_join  arrange
#' @rdname plotBoxplotDE
#' @aliases plotBoxplotDE
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3",
#'                                       "Low,Medium,Elevated",
#'                                       "DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", 
#'                          omicsData   = omicsData,
#'                          omicsNames  = c("protetest"),
#'                          omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |>
#'  runDiffAnalysis(SE.name = "protetest", method = "limmalmFit", contrastList = contrastList)
#'  
#'  # plot the boxplot
#'  plotBoxplotDE(experiments(MAE)[[1]], features = "AT1G15670", groupColor = "temperature")
#'  plotBoxplotDE(experiments(MAE)[[1]], features = "AT1G15670", groupColor = "imbibition")

setMethod(f          = "plotBoxplotDE",
          signature  = "RflomicsSE",
          definition = function(object, 
                                features = NULL, 
                                groupColor="groups", 
                                raw = FALSE){
              
              # check variable name
              if (is.null(features) || features == "" || length(features) != 1) {
                  message("set variable name")
                  
                  p <- ggplot() + theme_void() + ggtitle("set variable name")
                  
                  return(p)
              }
              
              Groups <- getDesignMat(object)
              object <- .checkTransNorm(object, raw = raw)
              
              # check presence of variable in SE
              object.DE <- tryCatch(object[features], error = function(e) e)
              if (!is.null(object.DE$message)) {
                  message(object.DE$message)
                  
                  p <- ggplot() + theme_void() + ggtitle(object.DE$message) 
                  
                  return(p)
              }
              
              if (raw) {
                  if (object.DE@metadata$omicType != "RNAseq") {
                      
                      pseudo <- assay(object.DE)
                      x_lab  <- features
                      title  <- features
                      
                  } else {
                      pseudo <- log2(assay(object.DE) + 1)
                      
                      x_lab  <- paste0("log2(", DE, " data)")
                      title  <- features
                  }
              } else{ 
                  if (object.DE@metadata$omicType != "RNAseq") {
                      
                      title <- features
                      pseudo <- assay(object.DE)
                      x_lab  <- paste0(features, " data")
                      
                      if (.isTransformed(object.DE) && getTransSettings(object.DE)$method != "none") {
                          title  <- paste0("Transformed (", 
                                           getTransSettings(object.DE)$method,
                                           ") ", title)
                      }
                      if (.isNorm(object.DE) && getNormSettings(object.DE)$method != "none") {
                          title <- paste0(title, " - normalization: ", 
                                          getNormSettings(object.DE)$method)
                      }  
                  } else {
                      
                      pseudo <- assay(object.DE) 
                      title  <- features
                      x_lab  <- paste0("log2(", features, " data)") 
                      
                  }
              }
              
              pseudo.gg <- melt(pseudo)
              colnames(pseudo.gg) <- c("features", "samples", "value")
              
              pseudo.gg <- pseudo.gg %>% 
                  full_join(Groups, by = "samples") %>%
                  arrange(groups)
              
              pseudo.gg <- arrange(pseudo.gg, get(groupColor))
              
              pseudo.gg$groups <- factor(pseudo.gg$groups, 
                                         levels = unique(pseudo.gg$groups))
              
              p <-  ggplot(pseudo.gg, 
                           aes(x = groups, y = value, label = features)) +
                  geom_boxplot( aes(fill = get(groupColor))) +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  guides(fill = guide_legend(title = "condition")) + 
                  xlab("") + 
                  ylab(x_lab) + 
                  ggtitle(title) #+
              #geom_point(alpha = 1/100,size=0)
              
              return(p)
              
          }
)

#' @rdname plotBoxplotDE
#' @title plotBoxplotDE
#' @param SE.name the name of the \link{RflomicsSE} object to fetch if the object is of class \link{RflomicsMAE}
#' @exportMethod plotBoxplotDE
#' @aliases plotBoxplotDE
setMethod(f          = "plotBoxplotDE",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name, 
                                features = NULL, 
                                groupColor = "groups",
                                raw = FALSE){
              
              plotBoxplotDE(object = object[[SE.name]], 
                            features = features,
                            groupColor = groupColor, raw = raw)
              
          })




###### ACCESSEUR: getteur and setteur


# ---- Get DE matrix from DiffExpAnalysis ----


#' @title Get DE matrix
#'
#' @param object of class \link{RflomicsSE} or \link{RflomicsMAE}
#' @param SE.name the name of the data to fetch in the object if the object is  \link{RflomicsMAE}
#' @return a matrix of results from the differential analyses.
#' @exportMethod getDEMatrix
#' @aliases getDEMatrix
#' @rdname getDEMatrix

setMethod(f          = "getDEMatrix",
          signature  = "RflomicsSE",
          definition = function(object){
              if (!is.null(metadata(object)$DiffExpAnal$mergeDEF)) {
                  metadata(object)$DiffExpAnal$mergeDEF
              } else {
                  stop("There is no DE matrix in this object.")
              }
          })

#' @title getDEMatrix
#' @rdname getDEMatrix
#' @exportMethod getDEMatrix
#' @aliases getDEMatrix

setMethod(f          = "getDEMatrix",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getDEMatrix(object = object[[SE.name]])
          })

# ---- Get union or intersection from list of contrasts ----

# very similar to filter_DE_from_SE but returns a vector instead of a SE.

#' @title Operation on differential analyses lists. Get union vector of DE 
#' entities from list of contrasts
#'
#' @param object an object of class \link{RflomicsSE}. Expects to find a slot with 
#' differential analyses results.
#' @param contrasts Vector of characters, expect to be contrast names. 
#' Default is null, the operation (union) is performed
#' on every contrasts found.
#' @param operation character. Either union or intersection.
#' Defines the operation to perform on the DE lists from the contrasts.
#' @return vector of unique DE entities
#' @rdname getDEList
#' @seealso [runDiffAnalysis()], [filterDiffAnalysis()]
#' @exportMethod getDEList
#' @importFrom tidyselect starts_with any_of
#' @importFrom dplyr select mutate filter
#' @aliases getDEList

setMethod(f          = "getDEList",
          signature  = "RflomicsSE",
          definition = function(object, contrasts = NULL, operation = "union"){
              
              if (is.null(contrasts) || length(contrasts) == 0) 
                  contrasts <- getSelectedContrasts(object)[["tag"]]
              if (.isContrastName(object, contrasts)) 
                  contrasts <- .convertContrastToTag(object, contrasts)
              
              if (!is.null(object@metadata$DiffExpAnal$Validcontrasts)) {
                  validTags <- .convertContrastToTag(object, getValidContrasts(object)$contrastName)
              } else {
                  validTags <- contrasts
              }
              
              tagsConcerned <- intersect(contrasts, validTags)
              
              if (length(tagsConcerned) == 0) 
                  stop("It seems there is no contrasts to select DE entities from.")
              
              df_DE <- getDEMatrix(object) |>
                  select(c("DEF", any_of(tagsConcerned)))
              
              if (operation == "intersection") {
                  DETab <- df_DE %>%
                      mutate(SUMCOL = select(., starts_with("H")) %>%
                                 rowSums(na.rm = TRUE)) %>%
                      filter(SUMCOL == length(tagsConcerned))
                  
              } else {
                  DETab <- df_DE %>%
                      mutate(SUMCOL = select(., starts_with("H")) %>%
                                 rowSums(na.rm = TRUE)) %>%
                      filter(SUMCOL >= 1)
              }
              
              return(unique(DETab$DEF))
          })

#' @rdname getDEList
#' @title getDEList
#' @aliases getDEList
#' @param SE.name the name of the data to fetch in the object if the object is 
#' a RflomicsMAE
#' @exportMethod getDEList

setMethod(f          = "getDEList",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, contrasts = NULL, 
                                operation = "union"){
              
              getDEList(object = object[[SE.name]], 
                        contrasts = contrasts,
                        operation = operation)
          })

# ---- Get diff setting ----

#' @title Get differential analysis setting parameters
#'
#' @param object of class RflomicsSE
#' @return List of differential analysis setting parametres.
#' @exportMethod getDiffSettings
#' @rdname runDiffAnalysis
#' @aliases getDiffSettings


setMethod(f          = "getDiffSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
              return(metadata(object)$DiffExpAnal$setting)   
          })

#' @title getDiffSettings
#' @rdname runDiffAnalysis
#' @aliases getDiffSettings
#' @exportMethod getDiffSettings

setMethod(f          = "getDiffSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getDiffSettings(object = object[[SE.name]])
          })

