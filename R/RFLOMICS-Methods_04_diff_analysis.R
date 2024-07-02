#==== DIFF-ANALYSIS METHODS ====
# D. Charif, 
# N. Bessoltane, 
# A. Hulot

##==== STAT METHOD ====

###==== METHOD to perform differential analysis ====

#' @title Run Differential Expression Analysis
#' @name runDiffAnalysis
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
#' @param object An object of class \link{RflomicsSE} or 
#' class \link{RflomicsMAE}
#' @param method A character vector giving the name of the differential 
#' analysis method to run. Either "edgeRglmfit" or "limmalmFit".
#' @param contrastList data.frame of contrast from generateExpressionContrast().
#' if NULL, it takes all selected contrasts.
#' @param p.adj.method The method choosen to adjust pvalue. Takes the same 
#' values as the ones of adj.p.adjust method.
#' @param p.adj.cutoff The adjusted pvalue cut-off
#' @param logFC.cutoff The lof2FC cutoff
#' @param clustermq A boolean indicating whether the constrasts have to be 
#' computed in local or in a distant machine
#' @param cmd boolean. Used in the interface. If TRUE, print cmd for the user.
#' @param ... Additional arguments.
#' @details
#' Functions and parameters used for RNAseq are those recommended in DiCoExpress 
#' workflow (see the paper in reference).
#' Functions and parameters used for proteomic and metabolomic data are those 
#' recommended in the (Efstathiou *et al.*, 2017)
#' @return A \link{RflomicsSE} or a \link{RflomicsMAE} object.
#' All the results are stored as a named list \code{DiffExpAnal} 
#' in the metadata slot of a given \link{RflomicsSE} object.
#' Objects are:
#' \itemize{
#' \item{stats:}{ data.frame giving a summary of the differential 
#'  statistical analysis results by contrast:
#'  number of DE features, number of up and down regulated features}
#' \item{setting: }{Parameters used for the differential analysis}
#' \itemize{
#' \item{method: }{ The method used for the differential analysis}
#' \item{p.adj.method:}{ The applied p-value correction method}
#' \item{p.adj.cutoff:}{ The cut-off applied for the adjusted p-value}
#' \item{logFC.cutoff:}{ The absolute log FC cut-off}
#' }
#' \item{RawDEFres: }{a list giving for each contrast the raw results of 
#'  the differential analysis method}
#' \item{DEF: }{a list giving for each contrast a data.frame of non filtered 
#'  differential expressed features with their statistics}
#' \item{TopDEF: }{a list giving for each contrast a data.frame of 
#'  differential expressed features ordered and filtered by p.adj.cutoff 
#'  with their statistics}
#' \item{mergeDEF: }{a data frame of 0/1 indicating for each features in row, 
#'  if it is DE in a given contrasts in column}
#' \item{contrasts: }{a data.table of the contrasts used for the differential 
#'  analysis}
#' }
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool 
#' to process multifactorial RNAseq experiments from quality controls to 
#' co-expression analysis through differential analysis based on contrasts 
#' inside GLM models. Plant Methods 16, 68 (2020).
#' 
#' Efstathiou G, Antonakis AN, Pavlopoulos GA, et al. ProteoSign: 
#' an end-user online differential proteomics statistical analysis platform. 
#' Nucleic Acids Res. 2017;45(W1):W300-W306.
#' @exportMethod runDiffAnalysis
#' @importFrom dplyr filter
#' @rdname runDiffAnalysis-methods
#' @seealso \code{\link{getDiffSettings}}, \code{\link{getDEList}},
#'  \code{\link{getDEMatrix}}
#' @seealso \code{\link{plotDiffAnalysis}}, \code{\link{plotHeatmapDesign}},
#' \code{\link{plotBoxplotDE}}
#' @example inst/examples/runDiffAnalysis.R
setMethod(f         = "runDiffAnalysis",
          signature = "RflomicsSE",
          definition = function(object, 
                                contrastList = NULL,
                                method = NULL,
                                p.adj.method="BH",
                                p.adj.cutoff=0.05, 
                                logFC.cutoff=0, 
                                clustermq=FALSE,
                                cmd = FALSE,
                                ...){
            
            DiffExpAnal <- list()
            
            modelFormula <- getModelFormula(object)
            if(length(modelFormula) == 0)
              stop("No model defined in the ", getDatasetNames(object), " object.")
            
            contrast.sel <- getSelectedContrasts(object)
            if(nrow(contrast.sel) == 0 || is.null(contrast.sel))
              stop("No contrasts defined in the ", getDatasetNames(object), " object.")
            
            if(is.null(contrastList)){
              contrastList <- getSelectedContrasts(object)
            }
            else{
              contrastList <- intersect(contrastList, contrast.sel)
              if(length(contrastList) == 0)
                stop("The specified contrasts do not match the selected contrasts")
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
            object <- generateContrastMatrix(object, contrastList = contrastList)
            
            
            DiffExpAnal[["contrasts"]] <- contrastList
            
            # remplacera à terme les lignes ci-dessus
            DiffExpAnal[["setting"]][["method"]] <- method
            DiffExpAnal[["setting"]][["p.adj.method"]] <- p.adj.method
            DiffExpAnal[["setting"]][["p.adj.cutoff"]] <- p.adj.cutoff
            DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
            
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
            model_matrix <- 
              model.matrix(
                as.formula(paste(modelFormula, collapse = " ")), 
                data = getDesignMat(object))
            
            ListRes <- 
              switch(
                method,
                "edgeRglmfit" = 
                  .tryRflomics(
                    .edgeRAnaDiff(count_matrix    = assay(object),
                                  model_matrix    = model_matrix[colnames(object),],
                                  group           = getCoeffNorm(object)$group,
                                  lib.size        = getCoeffNorm(object)$lib.size,
                                  norm.factors    = getCoeffNorm(object)$norm.factors,
                                  Contrasts.Sel   = contrastList,
                                  Contrasts.Coeff = object@metadata$design$Contrasts.Coeff,
                                  FDR             = 1,
                                  clustermq       = clustermq,
                                  cmd             = cmd)),
                "limmalmFit" = 
                  .tryRflomics(
                    .limmaAnaDiff(count_matrix    = assay(object2),
                                  model_matrix    = model_matrix[colnames(object2),],
                                  Contrasts.Sel   = contrastList,
                                  Contrasts.Coeff = object@metadata$design$Contrasts.Coeff,
                                  p.adj.cutoff    = 1,
                                  p.adj.method    = p.adj.method,
                                  clustermq       = clustermq,
                                  cmd             = cmd)))
            
            if (!is.null(ListRes$value)) {
              if (!is.null(ListRes$value[["RawDEFres"]])) {
                DiffExpAnal[["results"]] <- TRUE
                DiffExpAnal[["RawDEFres"]] <- ListRes$value[["RawDEFres"]]
                DiffExpAnal[["DEF"]] <- ListRes$value[["TopDEF"]]
              }else{
                DiffExpAnal[["results"]]    <- FALSE
                DiffExpAnal[["ErrorStats"]] <- ListRes$value[["ErrorTab"]]
              }
            }else{
              DiffExpAnal[["results"]]    <- FALSE
              DiffExpAnal[["Error"]]      <- ListRes$error
              DiffExpAnal[["ErrorStats"]] <- NULL
              
              
              #return(object)
            }
            
            object <- 
              setElementToMetadata(object, 
                                   name = "DiffExpAnal", 
                                   content = DiffExpAnal)
            
            ## filtering
            object <-  filterDiffAnalysis(object = object, 
                                          p.adj.cutoff = p.adj.cutoff, 
                                          logFC.cutoff = logFC.cutoff)
            
            return(object)
          })

#' @name runDiffAnalysis
#' @param SE.name SE.name the name of the dataset if the input object 
#' is a \link{RflomicsMAE}
#' @rdname runDiffAnalysis-methods
#' @exportMethod runDiffAnalysis
#' 
setMethod(f          = "runDiffAnalysis",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, 
                                contrastList = NULL,
                                method = NULL,
                                p.adj.method="BH",
                                p.adj.cutoff=0.05,
                                logFC.cutoff=0,
                                clustermq=FALSE, 
                                cmd = FALSE,
                                ...){
            
            # all verifications are done in this method
            object[[SE.name]] <-  
              runDiffAnalysis(object = object[[SE.name]],
                              contrastList = contrastList,
                              p.adj.method = p.adj.method,
                              method = method,
                              p.adj.cutoff = p.adj.cutoff,
                              logFC.cutoff = logFC.cutoff,
                              clustermq = clustermq,
                              cmd = cmd
              )
            return(object)
          })

###==== METHOD to generateContrastMatrix ====

#' @rdname runDiffAnalysis-methods
#' @name generateContrastMatrix
#' @description
#' \itemize{
#'    \item generateContrastMatrix: 
#'  Defines contrast matrix or contrast list with contrast 
#'  name and contrast coefficients}
#' @param contrastList a data.frame of contrasts generated by 
#' \link{generateExpressionContrast}
#' @exportMethod generateContrastMatrix
#' @importFrom stats formula terms.formula
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
setMethod(f          = "generateContrastMatrix",
          signature  = "RflomicsSE",
          definition = function(object, contrastList=NULL){
            
            if(is.null(contrastList))
              contrastList <- getSelectedContrasts(object)
            
            if(is.null(contrastList)) 
              stop("You need to select the contrasts (see ?generateContrastMatrix)")
            
            ExpDesign <- getDesignMat(object)
            
            factorBio <- getBioFactors(object)
            
            modelFormula <- getModelFormula(object)
            object@metadata$design$Contrasts.Coeff <- 
              .getContrastMatrixF(ExpDesign = ExpDesign, 
                                  factorBio = factorBio, 
                                  contrastList = contrastList$contrast, 
                                  modelFormula)
            object@metadata$design$Contrasts.Sel   <- contrastList
            
            return(object)
          })

#' @rdname runDiffAnalysis-methods
#' @exportMethod generateContrastMatrix
setMethod(f          = "generateContrastMatrix",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, 
                                contrastList=NULL){
            
            if (is.null(object[[SE.name]])) {
              stop("no Experiment named ", SE.name, " in MAE object")
            }
            if (is.null(modelFormula)) {
              stop("Model.formula arg is mandatory.")
            }
            if (is.null(contrastList)) stop("contrastList is mandatory.")
            if (any(!c("contrast", "contrastName", "groupComparison", "type") %in% names(contrastList))) {
              stop("contrastList data.frame must contain at least these colomn : contrast, contrastName, groupComparison, type")
            }
            
            object <- setModelFormula(object, modelFormula)
            
            object[[SE.name]] <- 
              generateContrastMatrix(object = object[[SE.name]], 
                                     contrastList = contrastList)
            
            return(object)
          })


###==== METHOD to filter differential analysis ====

#' @rdname runDiffAnalysis-methods
#' @name filterDiffAnalysis
#' @description
#' \itemize{
#'    \item filterDiffAnalysis: The filterDiffAnalysis method allows 
#'    filtering the results of the differential analysis based on a 
#'    new cutoff for p-value and fold change.}
#' @param p.adj.cutoff adjusted pvalue cutoff. Default is the parameter from 
#' the differential analysis.
#' @param logFC.cutoff cutoff for absolute value of log2FC. Default is the 
#' parameter from the differential analysis. 
#' @exportMethod filterDiffAnalysis
#' @importFrom dplyr full_join filter if_else mutate_at
#' @importFrom data.table data.table
#' @importFrom purrr reduce
setMethod(f          = "filterDiffAnalysis",
          signature  = "RflomicsSE",
          definition = function(object, 
                                p.adj.cutoff = 0.05, 
                                logFC.cutoff = 0){
            
            DiffExpAnal <- getAnalysis(object, name = "DiffExpAnal")
            
            if (is.null(DiffExpAnal[["RawDEFres"]])) {
              stop("can't filter the DiffExpAnal object because it doesn't exist")
            }
            
            # remplacera à terme les lignes ci-dessus
            DiffExpAnal[["setting"]][["p.adj.cutoff"]] <- p.adj.cutoff
            DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
            
            ## TopDEF: Top differential expressed features
            DEF_filtred <- lapply(seq_len(length(DiffExpAnal[["DEF"]])), function(x){
              res <- DiffExpAnal[["DEF"]][[x]]
              keep <- (res$Adj.pvalue < p.adj.cutoff) & (abs(res$logFC) > logFC.cutoff)
              res <- res[keep,]
              return(res)
            })
            names(DEF_filtred) <- names(DiffExpAnal[["RawDEFres"]])
            DiffExpAnal[["TopDEF"]] <- DEF_filtred
            
            ## stats
            DiffExpAnal[["stats"]] <- sumDiffExp(object)
            
            ## merge results in bin matrix
            DEF_list <- list()
            for (x in names(DiffExpAnal[["TopDEF"]])){
              res <- DiffExpAnal[["TopDEF"]][[x]]
              tmp <- data.frame(DEF = rownames(res), bin = rep(1,length(rownames(res))))
              colnames(tmp) <- c("DEF", x)
              
              if(dim(tmp)[1] != 0){ DEF_list[[x]] <- tmp }
            }
            
            DiffExpAnal[["mergeDEF"]] <- NULL
            
            if (length(DEF_list) != 0) {
              DiffExpAnal[["mergeDEF"]] <- DEF_list %>% 
                reduce(full_join, by="DEF") %>%
                mutate_at(.vars = 2:(length(DEF_list)+1),
                          .funs = function(x){
                            if_else(is.na(x), 0, 1)}) %>%
                as.data.frame()
            }
            
            object <- 
              setElementToMetadata(object,
                                   name = "DiffExpAnal", 
                                   content = DiffExpAnal)
            
            return(object)
          })

#' @rdname runDiffAnalysis-methods
#' @name filterDiffAnalysis
#' is a \link{RflomicsMAE}
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

###==== Set Valid Contrasts : (after differential analysis) ====

#' @rdname runDiffAnalysis-methods
#' @name setValidContrasts
#' @description
#' \itemize{
#'    \item setValidContrasts: Set the valid contrasts stored in \code{metadata} slot.}
#' @param contrastList A data.frame of contrast
#' @exportMethod setValidContrasts
setMethod(f          = "setValidContrasts",
          signature  = "RflomicsSE",
          definition = function(object, contrastList=NULL){
            
            object@metadata$DiffExpAnal[["Validcontrasts"]] <- contrastList
            
            return(object)
          })

#' @rdname runDiffAnalysis-methods
#' @param omicName a dataset name
#' @exportMethod setValidContrasts
setMethod(f          = "setValidContrasts",
          signature  = "RflomicsMAE",
          definition <- function(object, omicName=NULL, contrastList=NULL){
            
            object[[omicName]] <- 
              setValidContrasts(object[[omicName]], contrastList = contrastList)
            
            return(object)
          })

##==== GRAPHICAL METHOD ====

###==== Method to plot results of a differential analysis ====

#' @name plotDiffAnalysis
#' @rdname Rflomics-plots
#' @description
#' Graphical functions for differential analysis results:
#' \itemize{
#'    \item plotDiffAnalysis method draws a MAplot, a volcano plot and the 
#' p-values distribution from the results of a differential analysis.}
#' @param contrastName The contrastName for which the plots has to be drawn
#' @param typeofplots The plots you want to return. Default is all possible 
#' plots: MA plot, Volcano plot and non adjusted pvalues histogram.
#' @exportMethod plotDiffAnalysis
#' @export
#' @examples
#' # See runDiffAnalysis for an example that includes plotDiffAnalysis
setMethod(f="plotDiffAnalysis",
          signature="RflomicsSE",
          definition <- function(object, 
                                 contrastName,
                                 typeofplots = c("MA.plot", "volcano", "histogram")){
            
            plots <- list()
            DiffExpAnal <- 
              getAnalysis(object, name = "DiffExpAnal")
            
            resTable <- DiffExpAnal[["DEF"]][[contrastName]]
            
            logFC.cutoff <- getDiffSettings(object)[["abs.logFC.cutoff"]]
            p.adj.cutoff <- getDiffSettings(object)[["p.adj.cutoff"]]
            
            if ("MA.plot" %in% typeofplots) 
              plots[["MA.plot"]] <-  
              .plotMA(data = resTable, 
                      p.adj.cutoff = p.adj.cutoff, 
                      logFC.cutoff = logFC.cutoff, 
                      contrastName=contrastName)
            if ("volcano" %in% typeofplots) 
              plots[["Volcano.plot"]]   <-  
              .plotVolcanoPlot(data = resTable, 
                               p.adj.cutoff = p.adj.cutoff, 
                               logFC.cutoff = logFC.cutoff, 
                               contrastName=contrastName)
            if ("histogram" %in% typeofplots) 
              plots[["Pvalue.hist"]]  <-  
              .plotPValue(data =resTable, contrastName=contrastName)
            
            return(plots)
          })

#' @rdname Rflomics-plots
#' @name plotDiffAnalysis
#' @param SE.name SE.name the name of the dataset if the input object 
#' is a \link{RflomicsMAE}
#' @exportMethod plotDiffAnalysis
setMethod(f          = "plotDiffAnalysis",
          signature  = "RflomicsMAE",
          definition = function(object,
                                SE.name, 
                                contrastName, 
                                typeofplots = c("MA.plot", "volcano", "histogram")){

            return(plotDiffAnalysis(object      = object[[SE.name]],
                                    contrastName  = contrastName,
                                    typeofplots = typeofplots))
          })

###==== Method to plot plotHeatmapDesign ====

#' @name plotHeatmapDesign
#' @rdname Rflomics-plots
#' @description
#' \itemize{
#'    \item plotHeatmapDesign method draws a heatmap from the results of a 
#' differential analysis.}
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
#' @examples
#' # See runDiffAnalysis for an example that includes plotHeatmapDesign
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
            DiffExpAnal <- 
              getAnalysis(object, name = "DiffExpAnal")
            
            if (is.null(DiffExpAnal[["TopDEF"]][[contrastName]])) {
              stop("no DE variables")
            }
            
            resTable <- arrange(DiffExpAnal[["TopDEF"]][[contrastName]], Adj.pvalue)
            
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

#' @name plotHeatmapDesign
#' @rdname Rflomics-plots
#' @exportMethod plotHeatmapDesign
setMethod(f          = "plotHeatmapDesign",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, 
                                contrastName, 
                                splitFactor="none", 
                                title = "", annotNames = NULL, 
                                modalities = NULL, 
                                drawArgs = list(), 
                                heatmapArgs = list()){
            
            return(plotHeatmapDesign(object        = object[[SE.name]],
                                     contrastName    = contrastName,
                                     splitFactor     = splitFactor,
                                     title         = title,
                                     annotNames = annotNames,
                                     modalities   = modalities,
                                     drawArgs     = drawArgs,
                                     heatmapArgs  = heatmapArgs))
            
          })


###==== Method to plot BoxplotDE ====

#' @name plotBoxplotDE
#' @rdname Rflomics-plots
#' @description
#' \itemize{
#'    \item plotBoxplotDE method draws a boxplot showing the expression of 
#'    given differentially expressed feature.}
#' @param features variable name (gene/protein/metabolite name)
#' @param groupColor default to groups, indicate a variable in the design to
#' color the boxplots accordingly. 
#' @param raw Boolean. Plot the raw data or the transformed ones (TRUE)
#' @exportMethod plotBoxplotDE
#' @importFrom ggplot2  ggplot aes geom_boxplot  element_text theme guides 
#' guide_legend xlab ylab theme_void ggtitle
#' @importFrom dplyr full_join  arrange
#' @examples
#' # See runDiffAnalysis for an example that includes plotBoxplotDE
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

#' @rdname Rflomics-plots
#' @name plotBoxplotDE
#' @exportMethod plotBoxplotDE
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




##==== ACCESSEUR: getteur and setteur ====

### ---- Get DE matrix from DiffExpAnalysis ----

#' @rdname Rflomics-accessors
#' @name getDEMatrix
#' @description
#' Getters and Setters for Differential Analysis results:
#' \itemize{
#'    \item getDEMatrix: return a matrix of experimental design.}
#' @exportMethod getDEMatrix
#' @examples
#' # See runDiffAnalysis for an example that includes getDEMatrix
setMethod(f          = "getDEMatrix",
          signature  = "RflomicsSE",
          definition = function(object){
            
            if (!is.null(metadata(object)$DiffExpAnal$mergeDEF)){
              return(metadata(object)$DiffExpAnal$mergeDEF)
            } 
            else{
              return(NULL)
            }
          })

#' @rdname Rflomics-accessors
#' @name getDEMatrix
#' @exportMethod getDEMatrix
setMethod(f          = "getDEMatrix",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getDEMatrix(object = object[[SE.name]])
          })

### ---- Get union or intersection from list of contrasts ----

# very similar to filter_DE_from_SE but returns a vector instead of a SE.

#' @rdname Rflomics-accessors
#' @name getDEList
#' @description
#' \itemize{
#'    \item getDEList: return a vector of union or intersection of differential
#'    expressed features from list of contrasts.}
#' @param contrasts Vector of characters, expect to be contrast names. 
#' Default is null, the operation (union) is performed
#' on every contrasts found.
#' @param operation character. Either union or intersection.
#' Defines the operation to perform on the DE lists from the contrasts.
#' @exportMethod getDEList
#' @importFrom tidyselect starts_with any_of
#' @importFrom dplyr select mutate filter
#' @examples
#' # See runDiffAnalysis for an example that includes getDEList
setMethod(f          = "getDEList",
          signature  = "RflomicsSE",
          definition = function(object, contrasts = NULL, operation = "union"){
            
            validContrasts <- getValidContrasts(object)[["contrastName"]]
            if (is.null(validContrasts) || length(validContrasts) == 0){
              validContrasts <- names(getDEMatrix(object))[-1]
              if (is.null(validContrasts) || length(validContrasts) == 0)
                return(NULL)
            }
            
            if (is.null(contrasts) || length(contrasts) == 0){
              contrasts <- validContrasts 
            }
            else{
              contrasts <- intersect(contrasts, validContrasts)
              if (is.null(contrasts) || length(contrasts) == 0)
                return(NULL)
            }
            
            df_DE <- getDEMatrix(object) |>
              select(c("DEF", any_of(contrasts)))
            
            # if (is.null(df_DE) || nrow(df_DE) == 0 || ncol(df_DE) < 2)
            #   stop("")
            
            if (operation == "intersection") {
              DETab <- df_DE %>%
                mutate(SUMCOL = select(., -DEF) %>% rowSums(na.rm = TRUE)) %>%
                filter(SUMCOL == length(contrasts))
              
            } else {
              DETab <- df_DE %>%
                mutate(SUMCOL = select(., -DEF) %>% rowSums(na.rm = TRUE)) %>%
                filter(SUMCOL >= 1)
            }
            
            return(unique(DETab$DEF))
          })

#' @rdname Rflomics-accessors
#' @name getDEList
#' @exportMethod getDEList
setMethod(f          = "getDEList",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, contrasts = NULL, 
                                operation = "union"){
            
            getDEList(object = object[[SE.name]], 
                      contrasts = contrasts,
                      operation = operation)
          })

### ---- Get diff setting ----
#' @rdname Rflomics-accessors
#' @name getDiffSettings
#' @description
#' \itemize{
#'    \item getDiffSettings: retrun a list of differential expression analysis 
#'    settings of a given omic dataset}
#' @exportMethod getDiffSettings
#' @examples
#' # See runDiffAnalysis for an example that includes getDiffSettings
setMethod(f          = "getDiffSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(metadata(object)$DiffExpAnal$setting)   
          })

#' @rdname Rflomics-accessors
#' @name getDiffSettings
#' @exportMethod getDiffSettings
setMethod(f          = "getDiffSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getDiffSettings(object = object[[SE.name]])
          })


### ---- get Valid Contrasts : (after differential analysis) ----

#' @rdname Rflomics-accessors
#' @name getValidContrasts
#' @description
#' \itemize{
#'    \item getValidContrasts: return a data.frame of valided contrasts}
#' @exportMethod getValidContrasts
setMethod(f          = "getValidContrasts",
          signature  = "RflomicsSE",
          definition = function(object){
            
            return(object@metadata$DiffExpAnal[["Validcontrasts"]])
          })

#' @rdname Rflomics-accessors
#' @param omicName a dataset name
#' @exportMethod getValidContrasts
setMethod(f          = "getValidContrasts",
          signature  = "RflomicsMAE",
          definition = function(object, omicName){
            
            res <- getValidContrasts(object[[omicName]])
            
            return(res)
          })

### ---- getDiffAnalysesSummary ----
#' @name getDiffAnalysesSummary
#' @description
#' \itemize{
#'    \item getDiffAnalysesSummary: ...}
#' @param plot FALSE or TRUE
#' @importFrom dplyr mutate
#' @importFrom purrr reduce
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_col geom_text
#' @importFrom ggplot2 facet_grid scale_x_continuous labs position_stack
#' @exportMethod getDiffAnalysesSummary
#' @return a data.frame with differential analyses summary
#' @rdname Rflomics-accessors
setMethod(
  f          = "getDiffAnalysesSummary",
  signature  = "RflomicsMAE",
  definition = function(object, plot = FALSE, 
                        ylabelLength = 30,
                        nbMaxLabel = 20,
                        interface = FALSE) {
    # DataProcessing
    
    omicNames <- getAnalyzedDatasetNames(object, "DiffExpAnal")
    
    df.list <- list()
    for (dataset in omicNames) {
      
      Validcontrasts <- getValidContrasts(object[[dataset]])$contrastName
      if (is.null(Validcontrasts) || length(Validcontrasts) == 0)
        next
      
      p.adj <- getDiffSettings(object[[dataset]])$p.adj.cutoff
      logFC <- getDiffSettings(object[[dataset]])$abs.logFC.cutoff
      
      df.list[[dataset]] <-
        as.data.frame(
          object[[dataset]]@metadata$DiffExpAnal$stats)[Validcontrasts,] %>%
        mutate(dataset = dataset, contrasts = rownames(.), 
               settings = paste0("(p.adj: ", p.adj, " & logFC: ", logFC, ")"))
    }
    
    if (length(df.list) == 0)
      return(NULL)
    
    
    df <- reduce(df.list, rbind) |>
      melt(id = c("dataset", "settings", "contrasts", "All"), value.name = "Up_Down") |>
      filter(All != 0, !is.na(All)) |>
      mutate(percent = Up_Down / All * 100)
    
    if (isFALSE(plot))
      return(df)
    
    # For the interface, the labels (contrastNames) were replaced 
    # by tags if they exceed nbMaxLabel to simplify the figures.
    if(isTRUE(interface) && length(unique(df$contrasts)) > nbMaxLabel){
      
      contrast.df <- getSelectedContrasts(object)
      df$tabel <- lapply(df$contrasts, function(x){
        contrast.df[contrastName == x,]$tag
      }) %>% unlist
      
    }else{
      
      df$tabel <- 
        lapply(df$contrasts, function(x){
          vec <- seq(1, stringr::str_length(x), ylabelLength)
          stringr::str_sub(x, vec, vec+ylabelLength-1) %>% 
            paste(., collapse = "\n")
        }) |> unlist()
    }
    
    df$tabel <- factor(df$tabel, levels = unique(df$tabel))
    
    p <- ggplot(data = df, aes(y = tabel, x = percent, fill = variable)) +
      geom_col() +
      geom_text(aes(label = Up_Down), position = position_stack(vjust = 0.5)) +
      facet_wrap(.~ dataset+settings, ncol = 4) +
      scale_x_continuous(breaks = seq(0, 100, 25),
                         labels = paste0(seq(0, 100, 25), "%")) +
      labs(fill = NULL, x = "", y="")
    
    return(p)
  }
)
