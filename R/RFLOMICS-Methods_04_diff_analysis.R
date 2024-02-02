
################################### DIFF-ANALYSIS #############################


###### Statistical METHOD

## METHOD to perform differential analysis

#' @title RunDiffAnalysis
#' @description This is an interface method which run a differential analysis method on
#' omic datasets stored in an object of class \link{RflomicsSE}.
#' According to the type of omic and to a list of contrasts,
#' a differential analysis method is applied to each contrasts (or hypothesis).
#' Three methods are available according to the type of object:
#' \itemize{
#' \item{For RNAseq data: }{the \code{glmFit} function of the \code{edgeR} package}
#' \item{For proteomic and metabolomic data: }{the \code{lmFit} function of the \code{limma} package}
#' }
#' Parameters used for RNAseq are those recommended in DiCoExpress workflow (see the paper in reference)
#' @return
#' All the results are stored as a named list \code{DiffExpAnal} in the metadata slot of a
#' given \code{RflomicsSE} object.
#' Objects are:
#' \itemize{
#' \item{contrasts: }{The selected contrasts for which the differential analysis has been conducted}
#' \item{method: }{The method used for the differential analysis. }
#' \item{Adj.pvalue.method: The method applied for the pvalue adjustment.}
#' \item{Adj.pvalue.cutoff: The threshold applied for the pvalue adjustment}
#' \item{FDR: }{The false discovery rate given in input}
#' \item{RawDEFres: }{a list giving for each contrast the raw results of the differential analysis method}
#' \item{DEF: }{a list giving for each contrast a data.frame of non filtered differential expressed features}
#' \item{TopDEF: }{a list giving for each contrast a data.frame of differential expressed features by Adj.pvalue.cutoff}
#' \item{mergeDEF: }{A data frame indicating for each features in row, if it is DE in a given contrasts in column}
#' }
#' @param object an object of class \link{RflomicsSE} or \link{RflomicsMAE} 
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE 
#' @param design an object of class \link{ExpDesign-class}
#' @param DiffAnalysisMethod A character vector giving the name of the differential analysis method
#' to run. Either "edgeRglmfit" or "limmalmFit".
#' @param contrastList data.frame of contrast from getExpressionContrastF()
#' @param Adj.pvalue.method The method choosen to adjust pvalue. Takes the same values as the ones of adj.p.adjust method.
#' @param Adj.pvalue.cutoff The adjusted pvalue cut-off
#' @param clustermq A boolean indicating whether the constrasts have to be computed in local or in a distant machine
#' @return An object of class \link{RflomicsSE}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod RunDiffAnalysis
#' @importFrom dplyr filter
#' @rdname RunDiffAnalysis
#' 
methods::setMethod(f         = "RunDiffAnalysis",
                   signature = "RflomicsSE",
                   definition <- function(object, design, Adj.pvalue.method="BH", contrastList = NULL, DiffAnalysisMethod = NULL, 
                                          Adj.pvalue.cutoff=0.05, logFC.cutoff=0, clustermq=FALSE, parallel = FALSE, nworkers = 1,
                                          cmd = FALSE, modelFormula=NULL){
                     
                     # check args
                     if(is.null(contrastList) || nrow(contrastList) == 0) stop("contrastList arg is mandatory.")
                     
                     if(is.null(DiffAnalysisMethod) || isFALSE(DiffAnalysisMethod %in% c("edgeRglmfit", "limmalmFit"))) {
                       switch (object@metadata$omicType,
                               "RNAseq" = { DiffAnalysisMethod <- "edgeRglmfit"},
                               { DiffAnalysisMethod <- "limmalmFit" }
                       )
                       warning("DiffAnalyseMethod was missing. Detected omic type is ", object@metadata$omicType ," using ", DiffAnalysisMethod, " for differential analysis.")
                     }
                     
                     ## check completness
                     Completeness <- CheckExpDesignCompleteness(object)
                     if(isTRUE(Completeness[["error"]])) stop(Completeness[["messages"]])
                     
                     ## getcontrast
                     contrastList <- RFLOMICS::getExpressionContrast(object, modelFormula = modelFormula) %>% purrr::reduce(rbind) %>% 
                       dplyr::filter(contrast %in% contrastList$contrast)
                     
                     object <- getContrastMatrix(object, modelFormula = modelFormula, contrastList = contrastList)
                     object@metadata$design$Contrasts.Sel <- dplyr::mutate(contrastList, tag=paste0("H", 1:nrow(contrastList)))
                     
                     #Contrasts.Sel <- dplyr::filter(design@Contrasts.Sel, contrastName %in% contrastList)
                     
                     object@metadata$DiffExpAnal <- list()
                     object@metadata$DiffExpAnal[["contrasts"]] <- object@metadata$design$Contrasts.Sel
                     
                     # remplacera à terme les lignes ci-dessus
                     object@metadata$DiffExpAnal[["setting"]][["method"]]            <- DiffAnalysisMethod
                     object@metadata$DiffExpAnal[["setting"]][["Adj.pvalue.method"]] <- Adj.pvalue.method
                     object@metadata$DiffExpAnal[["setting"]][["Adj.pvalue.cutoff"]] <- Adj.pvalue.cutoff
                     object@metadata$DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
                     
                     # transform and norm if needed
                     if (DiffAnalysisMethod == "limmalmFit") {
                       
                       object2 <- object
                       
                       if (!isTransformed(object2)) object2 <-  apply_transformation(object2)
                       if (!isNorm(object2))        object2 <-  apply_norm(object2)
                     }
                     
                     # move in ExpDesign Constructor
                     model_matrix <- model.matrix(as.formula(paste(modelFormula, collapse = " ")), data = getDesignMat(object))
                     # model_matrix <- model.matrix(as.formula(paste(design@Model.formula, collapse = " ")), data = as.data.frame(design@ExpDesign))
                     # rownames(model_matrix) <- rownames(design@ExpDesign)
                     
                     ListRes <- switch(DiffAnalysisMethod,
                                       "edgeRglmfit" = try_rflomics(edgeR.AnaDiff(count_matrix  = SummarizedExperiment::assay(object),
                                                                                  model_matrix    = model_matrix[colnames(object),],
                                                                                  group           = getCoeffNorm(object)$group,
                                                                                  lib.size        = getCoeffNorm(object)$lib.size,
                                                                                  norm.factors    = getCoeffNorm(object)$norm.factors,
                                                                                  Contrasts.Sel   = object@metadata$design$Contrasts.Sel,
                                                                                  Contrasts.Coeff = object@metadata$design$Contrasts.Coeff,
                                                                                  FDR             = 1,
                                                                                  clustermq       = clustermq,
                                                                                  parallel        = parallel,
                                                                                  nworkers        = nworkers,
                                                                                  cmd             = cmd)),
                                       "limmalmFit" = try_rflomics(limma.AnaDiff(count_matrix    = SummarizedExperiment::assay(object2),
                                                                                 model_matrix      = model_matrix[colnames(object2),],
                                                                                 Contrasts.Sel     = object@metadata$design$Contrasts.Sel,
                                                                                 Contrasts.Coeff   = object@metadata$design$Contrasts.Coeff,
                                                                                 Adj.pvalue.cutoff = 1,
                                                                                 Adj.pvalue.method = Adj.pvalue.method,
                                                                                 clustermq         = clustermq,
                                                                                 cmd               = cmd)))
                     
                     if(! is.null(ListRes$value)){
                       if(! is.null(ListRes$value[["RawDEFres"]])){
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
                     object <-  FilterDiffAnalysis(object = object, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff)
                     
                     return(object)
                   })


#' @rdname RunDiffAnalysis
#' @title RunDiffAnalysis
#' @exportMethod RunDiffAnalysis
methods::setMethod(f          = "RunDiffAnalysis",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, Adj.pvalue.method="BH",
                                         contrastList = NULL, DiffAnalysisMethod = NULL,
                                         Adj.pvalue.cutoff=0.05, logFC.cutoff=0, clustermq=FALSE, 
                                         parallel = FALSE, nworkers = 1, cmd = FALSE, modelFormula=NULL){
                     
                     # all verifications are done in this method
                     object[[SE.name]] <-  RunDiffAnalysis(object = object[[SE.name]],
                                                           design = object@metadata$design,
                                                           Adj.pvalue.method = Adj.pvalue.method,
                                                           modelFormula = modelFormula,
                                                           contrastList = contrastList,
                                                           DiffAnalysisMethod = DiffAnalysisMethod,
                                                           Adj.pvalue.cutoff = Adj.pvalue.cutoff,
                                                           logFC.cutoff = logFC.cutoff,
                                                           clustermq = clustermq,
                                                           parallel = parallel,
                                                           nworkers = nworkers,
                                                           cmd = cmd
                     )
                     return(object)
                   })

## METHOD to filter differential analysis

#' Filter differential analysis
#'
#' @param object A RflomicsSE object
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @param Adj.pvalue.cutoff adjusted pvalue cutoff. Default is the parameter from the differential analysis.
#' @param logFC.cutoff cutoff for absolute value of log2FC. Default is the parameter from the differential analysis. 
#'
#' @return A RflomicsSE object or a RflomicsMAE, depending on the object type, 
#' where the differential analysis results have been actualized with the new parameters.
#' @exportMethod FilterDiffAnalysis
#' @rdname FilterDiffAnalysis
#' @importFrom dplyr filter if_else mutate_at
#' @importFrom data.table data.table
#' @importFrom purrr reduce
#'
methods::setMethod(f          = "FilterDiffAnalysis",
                   signature  = "RflomicsSE",
                   definition <- function(object, Adj.pvalue.cutoff = NULL, logFC.cutoff = NULL){
                     
                     if(is.null(object@metadata$DiffExpAnal[["RawDEFres"]])){
                       stop("can't filter the DiffExpAnal object because it doesn't exist")
                     }
                     
                     if (is.null(Adj.pvalue.cutoff)) 
                       Adj.pvalue.cutoff <- getDiffSetting(object)$Adj.pvalue.cutoff
                     
                     if (is.null(logFC.cutoff))
                       logFC.cutoff <- getDiffSetting(object)$abs.logFC.cutoff
                     
                     # remplacera à terme les lignes ci-dessus
                     object@metadata$DiffExpAnal[["setting"]][["Adj.pvalue.cutoff"]] <- Adj.pvalue.cutoff
                     object@metadata$DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
                     
                     ## TopDEF: Top differential expressed features
                     DEF_filtred <- lapply(1:length(object@metadata$DiffExpAnal[["DEF"]]), function(x){
                       res <- object@metadata$DiffExpAnal[["DEF"]][[x]]
                       keep <- (res$Adj.pvalue < Adj.pvalue.cutoff) & (abs(res$logFC) > logFC.cutoff)
                       res <- res[keep,]
                       return(res)
                     })
                     names(DEF_filtred) <- names(object@metadata$DiffExpAnal[["RawDEFres"]])
                     object@metadata$DiffExpAnal[["TopDEF"]] <- DEF_filtred
                     
                     ## stats
                     object@metadata$DiffExpAnal[["stats"]] <- sumDiffExp(object)
                     
                     ## merge results in bin matrix
                     DEF_list <- list()
                     for(x in names(object@metadata$DiffExpAnal[["TopDEF"]])){
                       res <- object@metadata$DiffExpAnal[["TopDEF"]][[x]]
                       tmp <- data.frame(DEF = rownames(res), bin = rep(1,length(rownames(res))))
                       colnames(tmp) <- c("DEF", dplyr::filter(object@metadata$DiffExpAnal$contrasts, contrastName == x)$tag)
                       
                       if(dim(tmp)[1] != 0){ DEF_list[[x]] <- tmp }
                     }
                     
                     object@metadata$DiffExpAnal[["mergeDEF"]] <- NULL
                     
                     if(length(DEF_list) != 0){
                       
                       object@metadata$DiffExpAnal[["mergeDEF"]] <- DEF_list %>% purrr::reduce(dplyr::full_join, by="DEF") %>%
                         dplyr::mutate_at(.vars = 2:(length(DEF_list)+1),
                                          .funs = function(x){
                                            dplyr::if_else(is.na(x), 0, 1)}) %>%
                         data.table::data.table()
                     }
                     
                     return(object)
                   })

#' @rdname FilterDiffAnalysis
#' @title FilterDiffAnalysis
#' @exportMethod FilterDiffAnalysis
methods::setMethod(f          = "FilterDiffAnalysis",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, 
                                         Adj.pvalue.cutoff = NULL, logFC.cutoff = NULL){
                     
                     if(!SE.name %in% names(object))
                       stop(SE.name, " isn't the name of an experiment in ", object)
                     
                     if (is.null(Adj.pvalue.cutoff) && is.null(logFC.cutoff)) {
                       
                       message("Parameter Adj.pvalue.cutoff and logFC.cutoff are both NULL. Not changing anything")
                       return(object)
                       
                     }else{
                       
                       if (is.null(Adj.pvalue.cutoff)) 
                         Adj.pvalue.cutoff <- getDiffSetting(object[[SE.name]])$Adj.pvalue.cutoff
                       
                       if (is.null(logFC.cutoff))
                         logFC.cutoff <- getDiffSetting(object[[SE.name]])$abs.logFC.cutoff
                       
                       object[[SE.name]] <-  FilterDiffAnalysis(object = object[[SE.name]],
                                                                Adj.pvalue.cutoff = Adj.pvalue.cutoff,
                                                                logFC.cutoff = logFC.cutoff)
                       
                       return(object)
                     }
                     
                   })

###### Graphical METHOD

## Method to plot results of a differential analysis

#' @title DiffAnal.plot
#' @description
#' This is an interface method which draw a MAplot, a volcano plot and the pvalues distribution from the results of a differential analysis
#' performed on omic datasets stored in an object of class \link{RflomicsSE}
#' @param object An object of class \link{RflomicsSE}
#' @param hypothesis The hypothesis for which the plots has to be drawn
#' @param typeofplots The plots you want to return. Default is all possible plots: MA plot, Volcano plot and non adjusted pvalues histogram.
#' @return plot
#' @exportMethod DiffAnal.plot
#' @rdname DiffAnal.plot
#' @export
#' 
methods::setMethod(f="DiffAnal.plot",
                   signature="RflomicsSE",
                   
                   definition <- function(object, hypothesis, typeofplots = c("MA.plot", "volcano", "histogram")){
                     
                     if ( isTagName(object, hypothesis)) hypothesis <-  convertTagToContrast(object, hypothesis)
                     
                     plots <- list()
                     
                     res      <- object@metadata$DiffExpAnal[["RawDEFres"]][[hypothesis]]
                     resTable <- object@metadata$DiffExpAnal[["DEF"]][[hypothesis]]
                     
                     logFC.cutoff      <- getDiffSetting(object)[["abs.logFC.cutoff"]]
                     Adj.pvalue.cutoff <- getDiffSetting(object)[["Adj.pvalue.cutoff"]]
                     
                     if ("MA.plot" %in% typeofplots) plots[["MA.plot"]]        <-  MA.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff, hypothesis=hypothesis)
                     if ("volcano" %in% typeofplots) plots[["Volcano.plot"]]   <-  Volcano.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff, hypothesis=hypothesis)
                     if ("histogram" %in% typeofplots) plots[["Pvalue.hist"]]  <-  pvalue.plot(data =resTable, hypothesis=hypothesis)
                     return(plots)
                   })

#' @rdname DiffAnal.plot
#' @title DiffAnal.plot
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod DiffAnal.plot
methods::setMethod(f          = "DiffAnal.plot",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, hypothesis, typeofplots = c("MA.plot", "volcano", "histogram")){
                     
                     if (isTagName(object, hypothesis)) hypothesis <-  convertTagToContrast(object, hypothesis)
                     
                     return(DiffAnal.plot(object      = object[[SE.name]],
                                          hypothesis  = hypothesis,
                                          typeofplots = typeofplots))
                     
                   })


#' @title heatmapPlot
#' @description
#' This is an interface method which draw a heatmap from the results of a differential analysis
#' performed on omic datasets stored in an object of class \link{RflomicsSE}
#' @param object An object of class \link{RflomicsSE}
#' @param hypothesis The hypothesis for which the MAplot has to be drawn
#' @param condition characters. Default to none. Name of a feature in the design matrix, splits the samples on the heatmap according to its modalities.  
#' @param title characters. Title of the heatmap. 
#' @param annot_to_show vector. Names of the annotations to keep in the Heatmap. Default takes all available information.
#' @param subset_list named list of vectors of modalities to subset and print on the heatmap. 
#' @param draw_args,heatmap_args  named lists. Any additional parameter passed to ComplexHeatmap::Heatmap or ComplexHeatmap::draw
#' @return plot
#' @exportMethod heatmapPlot
#' @export
#' @importFrom dplyr arrange select
#' @importFrom tidyselect any_of
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importClassesFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importMethodsFrom ComplexHeatmap draw
#' @importFrom grid gpar
#' @rdname heatmapPlot
#' 
methods::setMethod(f          = "heatmapPlot",
                   signature  = "RflomicsSE",
                   definition = function(object, 
                                         hypothesis, 
                                         condition="none", 
                                         title = "", 
                                         annot_to_show = NULL, 
                                         subset_list = NULL, 
                                         draw_args = list(), 
                                         heatmap_args = list()){
                     
                     Groups     <- getDesignMat(object)
                     
                     if (is.null(object@metadata$DiffExpAnal[["TopDEF"]][[hypothesis]])) {
                       stop("no DE variables")
                     }
                     
                     resTable <- dplyr::arrange(object@metadata$DiffExpAnal[["TopDEF"]][[hypothesis]], Adj.pvalue)
                     
                     if (dim(resTable)[1] == 0) {
                       stop("no differentially expressed variables...")
                     }
                     
                     if (dim(resTable)[1] > 2000) {
                       message("differentially expressed variables exceeding 2000 variables, only the first 2000 will be displayed")
                       resTable <- resTable[1:2000,]
                       title <- ifelse(title == "", paste0(title, "plot only 2000 TOP DE variables"),
                                       paste0(title, "\nplot only 2000 TOP DE variables"))
                     }
                     
                     object2 <- checkTransNorm(object, raw = FALSE)
                     m.def  <- assay(object2)
                     
                     m.def <- as.data.frame(m.def) %>%
                       dplyr::select(tidyselect::any_of(Groups$samples))
                     
                     # filter by DE
                     m.def.filter <- subset(m.def, rownames(m.def) %in% row.names(resTable))
                     
                     # normalize count
                     
                     # Center
                     m.def.filter.center <- t(scale(t(m.def.filter), center = TRUE, scale = FALSE))
                     
                     # Annotations datatable
                     df_annotation <- Groups %>% dplyr::select(!samples & !groups)  
                     df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),] 
                     
                     # Subset the dataset to print only interesting modalities
                     if (!is.null(subset_list)) {
                       if (is.null(names(subset_list))) {
                         message("In heatmapPlot, subset_list argument needs a named list. Not subsetting")
                       }else{ 
                         samplesToKeep <- Reduce("intersect", lapply(
                           1:length(subset_list),
                           FUN = function(i){
                             col_nam <- names(subset_list)[i]
                             rownames(df_annotation[which(df_annotation[[col_nam]] %in% subset_list[[i]]),])
                           }
                         ))
                         
                         df_annotation <- df_annotation[which(rownames(df_annotation) %in% samplesToKeep),]
                         m.def.filter.center <- m.def.filter.center[, which(colnames(m.def.filter.center) %in% samplesToKeep)]
                         
                         df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),]
                       }
                     }
                     
                     # Split management
                     column_split.value <- if (condition != "none") { df_annotation[, condition] } else { NULL }
                     
                     # Select the right columns
                     if (!is.null(annot_to_show)) {
                       df_annotation <- df_annotation %>% dplyr::select(tidyselect::any_of(annot_to_show))
                     }
                     
                     # Color annotations
                     set.seed(10000) # seed for chosing palette. Ensure running twice the function get the same annotations colors.
                     selectPal <- sample(rownames(RColorBrewer::brewer.pal.info),  size = ncol(df_annotation), replace = FALSE)
                     
                     color_list <- lapply(1:ncol(df_annotation), FUN = function(i){
                       annot_vect <- unique(df_annotation[,i])
                       
                       col_vect <-  grDevices::colorRampPalette(
                         suppressWarnings({ RColorBrewer::brewer.pal(n = min(length(annot_vect), 8), name = selectPal[i])}))(length(annot_vect)) 
                       names(col_vect) <- annot_vect 
                       col_vect[!is.na(names(col_vect))] # RcolorBrewer::brewer.pal n is minimum 3, remove NA names if only 2 levels
                     })
                     names(color_list) <- colnames(df_annotation)
                     
                     column_ha <- ComplexHeatmap::HeatmapAnnotation(df = df_annotation, col = color_list)
                     
                     namArg <- ifelse(getOmicsTypes(object) == "RNAseq", "normalized counts", "XIC")
                     
                     # Arguments for Heatmap
                     heatmap_args <- c(
                       list(matrix = m.def.filter.center,
                            name = namArg,
                            show_row_names = ifelse( dim(m.def.filter.center)[1] > 50, FALSE, TRUE),
                            row_names_gp = grid::gpar(fontsize = 8),
                            column_names_gp = grid::gpar(fontsize = 12),
                            row_title_rot = 0 ,
                            clustering_method_columns = "ward.D2",
                            cluster_column_slice = FALSE,
                            column_split = column_split.value,
                            top_annotation = column_ha,
                            column_title = title),
                       heatmap_args)
                     
                     # Arguments for drawing the heatmap
                     draw_args <- c(list(merge_legend = TRUE),
                                    draw_args)            
                     
                     # Drawing heatmap in a null file to not plot it
                     pdf(file = NULL)
                     ha <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
                     
                     draw_args$object <- ha
                     ha <- do.call(ComplexHeatmap::draw, draw_args)
                     
                     dev.off()
                     
                     return(ha)
                   })


#' @rdname heatmapPlot
#' @title heatmapPlot
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod heatmapPlot
methods::setMethod(f          = "heatmapPlot",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, hypothesis, condition="none", title = "", annot_to_show = NULL, subset_list = NULL, draw_args = list(), heatmap_args = list()){
                     
                     
                     if (isTagName(object, hypothesis)) hypothesis <- convertTagToContrast(object, hypothesis)
                     
                     return(heatmapPlot(object        = object[[SE.name]],
                                        hypothesis    = hypothesis,
                                        condition     = condition,
                                        title         = title,
                                        annot_to_show = annot_to_show,
                                        subset_list   = subset_list,
                                        draw_args     = draw_args,
                                        heatmap_args  = heatmap_args))
                     
                   })




#' @title boxplot.DE.plot
#'
#' @param object An object of class \link{RflomicsSE}
#' @param DE variable name (gene/protein/metabolite name)
#' @exportMethod boxplot.DE.plot
#' @importFrom ggplot2 geom_density xlab theme_void ggtitle ggplot geom_boxplot guide_legend guides
#' @importFrom dplyr full_join  arrange
#' @noRd

methods::setMethod(f          = "boxplot.DE.plot",
                   signature  = "RflomicsSE",
                   definition = function(object, DE = NULL, condition="groups", raw = FALSE){
                     
                     # check variable name
                     if (is.null(DE) || DE == "" || length(DE) != 1) {
                       message("set variable name")
                       
                       p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("set variable name")
                       
                       return(p)
                     }
                     
                     Groups <- getDesignMat(object)
                     object <- checkTransNorm(object, raw = raw)
                     
                     # check presence of variable in SE
                     object.DE <- tryCatch(object[DE], error = function(e) e)
                     if (!is.null(object.DE$message)) {
                       message(object.DE$message)
                       
                       p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(object.DE$message) 
                       
                       return(p)
                     }
                     
                     if (raw) {
                       if (object.DE@metadata$omicType != "RNAseq") {
                         
                         pseudo <- SummarizedExperiment::assay(object.DE)
                         x_lab  <- DE
                         title  <- DE
                         
                       } else {
                         pseudo <- log2(SummarizedExperiment::assay(object.DE) + 1)
                         
                         x_lab  <- paste0("log2(", DE, " data)")
                         title  <- DE
                       }
                     } else{ 
                       if (object.DE@metadata$omicType != "RNAseq") {
                         
                         title <- DE
                         pseudo <- SummarizedExperiment::assay(object.DE)
                         x_lab  <- paste0(DE, " data")
                         
                         if (isTransformed(object.DE) && getTransSetting(object.DE)$method != "none") {
                           title  <- paste0("Transformed (", getTransSetting(object.DE)$method, ") ", title)
                         }
                         if (isNorm(object.DE) && getNormSetting(object.DE)$method != "none") {
                           title <- paste0(title, " - normalization: ", getNormSetting(object.DE)$method)
                         }  
                       } else {
                         
                         pseudo <- SummarizedExperiment::assay(object.DE) 
                         title  <- DE
                         x_lab  <- paste0("log2(", DE, " data)") 
                         
                       }
                     }
                     
                     pseudo.gg <- pseudo %>% reshape2::melt()
                     colnames(pseudo.gg) <- c("features", "samples", "value")
                     
                     pseudo.gg <- pseudo.gg %>% dplyr::full_join(Groups, by="samples") %>%
                       dplyr::arrange(groups)
                     
                     pseudo.gg <- dplyr::arrange(pseudo.gg, get(condition))
                     
                     pseudo.gg$groups <- factor(pseudo.gg$groups, levels = unique(pseudo.gg$groups))
                     
                     p <- ggplot2::ggplot(pseudo.gg, ggplot2::aes(x=groups, y=value, label = features)) +
                       ggplot2::geom_boxplot(ggplot2::aes(fill=get(condition))) +
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
                       ggplot2::guides(fill=guide_legend(title="condition")) + 
                       ggplot2::xlab("") + 
                       ggplot2::ylab(x_lab) + 
                       ggtitle(title) #+
                     #geom_point(alpha = 1/100,size=0)
                     
                     return(p)
                     
                   }
)

#' @rdname boxplot.DE.plot
#' @title boxplot.DE.plot
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod boxplot.DE.plot
methods::setMethod(f          = "boxplot.DE.plot",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, DE = NULL, condition="groups"){
                     
                     boxplot.DE.plot(object = object[[SE.name]], 
                                     DE = DE,
                                     condition = condition)
                     
                   })




###### ACCESSEUR: getteur and setteur


# ---- Get DE matrix from DiffExpAnalysis ----


#' @title Get DE matrix
#'
#' @param object a RflomicsSE object
#' @return a matrix of results from the differential analyses.
#' @exportMethod getDEMatrix
#' @rdname getDEMatrix

methods::setMethod(f          = "getDEMatrix",
                   signature  = "RflomicsSE",
                   definition = function(object){
                       if (!is.null(object@metadata$DiffExpAnal$mergeDEF)) {
                         object@metadata$DiffExpAnal$mergeDEF
                       } else {
                         stop("There is no DE matrix in this object.")
                       }
                   })

#' @rdname getDEMatrix
#' @title getDEMatrix
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod getDEMatrix

methods::setMethod(f          = "getDEMatrix",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name){
                     getDEMatrix(object = object[[SE.name]])
                   })

#' 
#' #' @param contrast character name (can be a vector of name) for the contrast to select.
#' #' @param union Boolean value. TRUE : union; FALSE : intersection
#' #' @export
#' #' @importFrom tidyselect any_of
#' #' @importFrom dplyr select
#' #' @rdname getDE
#' getDE <- function(object, contrast, union = TRUE) {
#'   
#'   if (isContrastName(object, contrast)) 
#'     contrast <- convertContrastToTag(object, contrast)
#'   
#'   DEmat <- getDEMatrix(object)
#'   DEmat <- DEmat %>% select(any_of(c("DEF", contrast)))
#'   
#'   if (union) return(DEmat[rowSums(DEmat[,-1]) >= 1,])
#'   
#'   return(DEmat[rowSums(DEmat[,-1]) >= length(contrast),])
#' }


# ---- Get union or intersection from list of contrasts ----

# very similar to filter_DE_from_SE but returns a vector instead of a SE.

#' @title Operation on differential analyses lists. Get union vector of DE entities from list of contrasts
#'
#' @param object an object of class RflomicsSE. Expects to find a slot with differential analyses results.
#' @param contrasts Vector of characters, expect to be contrast names. Default is null, the operation (union) is performed
#' on every contrasts found.
#' @param operation character. Either union or intersection.
#' Defines the operation to perform on the DE lists from the contrasts.
#' @return vector of unique DE entities
#' @rdname getDEList
#' @exportMethod getDEList
#' @importFrom tidyselect starts_with any_of
#' @importFrom dplyr select mutate filter

methods::setMethod(f          = "getDEList",
                   signature  = "RflomicsSE",
                   definition = function(object, contrasts = NULL, operation = "union"){
                   
  if (is.null(contrasts) || length(contrasts) == 0) 
    contrasts <- getSelectedContrasts(object)[["tag"]]
  if (isContrastName(object, contrasts)) 
    contrasts <- convertContrastToTag(object, contrasts)
  
  if (!is.null(object@metadata$DiffExpAnal$Validcontrasts)) {
    validTags <- convertContrastToTag(object, getValidContrasts(object)$contrastName)
  } else {
    validTags <- contrasts
  }
  
  tagsConcerned <- intersect(contrasts, validTags)
  
  if (length(tagsConcerned) == 0) 
    stop("It seems there is no contrasts to select DE entities from.")
  
  df_DE <- getDEMatrix(object) %>%
    select(c("DEF", any_of(tagsConcerned)))
  
  if (operation == "intersection") {
    DETab <- df_DE %>%
      mutate(SUMCOL = select(., starts_with("H")) %>%
               rowSums(na.rm = TRUE)) %>%
      filter(SUMCOL == length(validTags))
    
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
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod getDEList

methods::setMethod(f          = "getDEList",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, contrasts = NULL, operation = "union"){
    
                     getDEList(object = object[[SE.name]], 
                                      contrasts = contrasts,
                                      operation = operation)
                   })


#' #' return gene list
#' #'
#' #' @param matrix matrix of DE results
#' #' @param colnames colnames
#' #' @param mergeType either union or intersection.
#' #' @return list of genes
#' #' @export
#' #' @noRd
#' #'
#' getDEGlist_for_coseqAnalysis <- function(matrix, colnames = colnames(matrix)[-1], mergeType="union"){
#'   
#'   if (length(colnames) == 0 ){ return(NULL) }
#'   
#'   matrix_sum <- matrix %>% dplyr::mutate(sum = dplyr::select(., tidyselect::all_of(colnames)) %>% rowSums(.))
#'   
#'   DEG_list <- switch(mergeType,
#'                      
#'                      "union"={        dplyr::filter(matrix_sum, sum != 0) },
#'                      "intersection"={ dplyr::filter(matrix_sum, sum == length(colnames)) }
#'   )
#'   
#'   if (length(DEG_list$DEF) == 0 ){ return(NULL) }
#'   
#'   return(DEG_list$DEF)
#' }



