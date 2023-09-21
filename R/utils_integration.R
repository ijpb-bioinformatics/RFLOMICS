#---- Filter DE entities from SE (returns vector) ----

# TODO A voir pour remplacer par un getter peut-être
# TODO Est-ce que cette fonction est toujorus utilisée ? Remplacée par opDEList logiquement.

#' @title filter_DE_from_SE
#'
#' @param SEobject An object of class \link{SummarizedExperiment}
#' @param contrasts_arg the name (or names) of the contrasts to consider.
#' @param type type of union for the entities, either union or intersection.
#' @return An object of class \link{SummarizedExperiment}
#' @export
#' @noRd
#'

filter_DE_from_SE <- function(SEobject, contrasts_arg = NULL, type = "union"){
  
  tabCorresp <- SEobject@metadata$DiffExpAnal$Validcontrasts %>% 
    dplyr::select(contrastName, tag)
  if (is.null(contrasts_arg))   contrasts_arg <- SEobject@metadata$DiffExpAnal$Validcontrasts$contrastName
  
  tabCorresp <- tabCorresp %>% dplyr::filter(contrastName %in% contrasts_arg)
  contrasts_select <- tabCorresp$tag
  
  tab1 <- SEobject@metadata$DiffExpAnal$mergeDEF %>%
    dplyr::select(tidyselect::any_of(c("DEF", contrasts_select)))
  
  if(type == "intersection"){
    
    DETab <- tab1 %>%
      dplyr::mutate(SUMCOL = dplyr::select(., tidyselect::starts_with("H")) %>% rowSums(na.rm = TRUE))  %>%
      dplyr::filter(SUMCOL==length(contrasts_select))
    
  }else{
    
    DETab <- tab1 %>%
      dplyr::mutate(SUMCOL = dplyr::select(., tidyselect::starts_with("H")) %>% rowSums(na.rm = TRUE))  %>%
      dplyr::filter(SUMCOL >= 1)
  }
  
  SEobject <- SEobject[DETab$DEF,]
  
  return(SEobject)
}


# ---- remove batch effects from omics : ----

# TODO A MODIFIER ON DOIT POUVOIR PRENDRE EN COMPTE TOUS LES BATCH EFFECTS
#' @title rbe_function
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEobject An object of class \link{SummarizedExperiment}
#' @return An object of class \link{SummarizedExperiment}
#' @export
#' @noRd
#'
rbe_function = function(object, SEobject, cmd = FALSE){
  
  assayTransform <- SummarizedExperiment::assay(SEobject)
  ftypes <- RFLOMICS::getFactorTypes(object)
  colBatch <- names(ftypes)[ftypes == "batch"]
  
  if (cmd) {
    print(paste0("#     =>Correction for Batch: ", 
                 paste(colBatch, collapse = " "), 
                 " in ", SEobject@metadata$omicType))
  }
  
  newFormula <- gsub(pattern = paste(paste(colBatch, "[+]"),  collapse = "|"), "", RFLOMICS::getModelFormula(object))
  newFormula <- gsub(pattern = "~ [+] ", "~ ", newFormula) # use ?
  designToPreserve <- model.matrix(as.formula(newFormula), data = SEobject@metadata$Groups)
  
  if (length(colBatch) == 1) {
    rbeRes <- limma::removeBatchEffect(assayTransform, batch = SEobject@metadata$Groups[,colBatch], design = designToPreserve)
  }else if (length(colBatch) >= 2) {
    
    rbeRes <- limma::removeBatchEffect(assayTransform,
                                       batch  = SEobject@metadata$Groups[,colBatch[1]],
                                       batch2 = SEobject@metadata$Groups[,colBatch[2]],
                                       design = designToPreserve)
  }
  if (length(colBatch) > 2) print("sorry, only 2 batches effect for now!") 
  # TODO : find a way to have more than 2 batches (for commandline only)
  
  SEobject@metadata[["correction_batch_method"]] <- "limma (removeBatchEffect)"
  
  SummarizedExperiment::assay(SEobject) <- rbeRes
  
  return(SEobject)
}

# ----- Transform rnaseq assay from SE ----

#' @title Remove batch effect and transform rnaseq data
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEname the name of the rnaseq data to transform. Supposed to be a SummarizedExperiment.
#' @param correctBatch if TRUE, correction of batch effects. 
#' @param transformation the name of the transformation to be applied on counts. Default is limma voom. 
#' No other option for now. 
#' @return An object of class \link{MultiAssayExperiment}
#' @export
#' @noRd
#'

rnaseqRBETransform <- function(object, 
                               SEname, 
                               correctBatch = FALSE, 
                               transformation = "limma (voom)",
                               contrasts_names = NULL,
                               type = "union",
                               choice = "DE",
                               cmd = FALSE
){
  
  # TODO : check if differential analysis has been performed.
  # TODO : in case of removeBatchEffect not working, what do we do?
  # TODO : if no DE (ex: intersection is 0 genes), what do we do?
  
  if (class(object) != "MultiAssayExperiment") stop("object is not a MultiAssyExperiment")
  
  rnaDat <- object[[SEname]] 
  assayTransform <- SummarizedExperiment::assay(rnaDat)
  
  if (!is.integer(assayTransform) && !identical(assayTransform, floor(assayTransform))) {
    message(paste0("You indicated RNASeq data for ", SEname, "but it is not recognized as count data")) 
    print(assayTransform[1:3, 1:3])
  }
  
  DMat      <- getDesignMat(object)
  coefNorm  <- getNormCoeff(rnaDat)
  designMat <- model.matrix(formula(getModelFormula(object)), data = DMat)
  
  DGEObject <- edgeR::DGEList(
    counts       = assayTransform,
    norm.factors = coefNorm$norm.factors,
    lib.size     = coefNorm$lib.size,
    samples      = DMat %>% 
      dplyr::filter(row.names(DMat) %in% colnames(assayTransform))
  )
  
  limmaRes <- limma::voom(DGEObject,
                          design = designMat[which(rownames(designMat) %in% colnames(assayTransform)),]) 
  
  SummarizedExperiment::assay(rnaDat) <- limmaRes$E
  
  if (correctBatch)    rnaDat <- rbe_function(object, rnaDat)
  
  if (choice == "DE")  rnaDat <- rnaDat[opDEList(rnaDat, contrasts = contrasts_names, operation = type),]
  
  rnaDat@metadata[["correction_batch"]]             <- correctBatch
  rnaDat@metadata[["transform_results_all"]]        <- limmaRes 
  rnaDat@metadata[["transform_method_integration"]] <- transformation
  
  object[[SEname]] <- rnaDat
  
  return(object)
} 


# ----- Transform rnaseq assay from SE ----

#' @title RBETransform
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEname the name of the omics data to transform. No counts data.  
#' @param correctBatch if TRUE, correction of batch effects. 
#' @return An object of class \link{MultiAssayExperiment}
#' @export
#' @noRd
#'

RBETransform <- function(object,
                         SEname,
                         correctBatch = TRUE,
                         contrasts_names = NULL,
                         type = "union",
                         choice = "DE",
                         cmd = FALSE){
  
  omicsDat <- object[[SEname]]
  omicsDat@metadata[["correction_batch"]]             <- correctBatch
  omicsDat@metadata[["transform_method_integration"]] <- omicsDat@metadata$transform_method
  
  if (correctBatch)    omicsDat <- rbe_function(object, omicsDat)
  if (choice == "DE")  omicsDat <- omicsDat[opDEList(omicsDat, contrasts = contrasts_names, operation = type),]
  
  object[[SEname]] <- omicsDat
  
  return(object)
}

# ----- MixOmics: plot variance explained ----

#' @title plot_MO_varExp
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param selectedResponse a character string of the response variable to consider
#' @param mode Can be NULL (default), "cumulative" or "comp". Defines the type of graph to return
#' @return An object of class \link{MultiAssayExperiment}
#' @export
#'

plot_MO_varExp <- function(object, selectedResponse, mode = NULL){
  
  if (class(object) != "MultiAssayExperiment") stop("Object is not a MultiAssayExperiment.")
  if (is.null(object@metadata$mixOmics)) stop("It seems this object has no mixOmics results.")
  if (is.null(object@metadata$mixOmics[[selectedResponse]])) stop("It seems you didn't run MixOmics on this particular variable.")
  
  Data_res <- object@metadata$mixOmics[[selectedResponse]]$MixOmics_results
  gg_return <- NULL
  
  if (is.null(mode)) {
    gg_return <- ggpubr::ggarrange(plot_MO_1(Data_res),
                                   plot_MO_2(Data_res), 
                                   ncol = 2)
  }
  else if (tolower(mode) == "cumulative") {
    gg_return <- plot_MO_1(Data_res)
  }
  else if (tolower(mode) == "comp") {
    gg_return <- plot_MO_2(Data_res)
  }
  
  return(gg_return)
}

#' @title Plot cumulative Explained variance for mixOmics results
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string. 
#' @return A ggplot2 graph
#' @keywords internal
#' @noRd

plot_MO_1 <- function(Data_res){
  dat_explained <- reshape2::melt(do.call("rbind", Data_res$prop_expl_var))
  colnames(dat_explained) <- c("Dataset", "Component", "% of explained variance")
  dat_explained$`% of explained variance` <- dat_explained$`% of explained variance`*100
  
  dat_comb <- dat_explained %>% 
    dplyr::group_by(Dataset) %>% 
    dplyr::summarise("Cumulative Explained Variance" = sum(`% of explained variance`))
  
  if (is(Data_res, "block.splsda") || is(Data_res, "block.plsda")) {
    dat_comb <- dat_comb %>% dplyr::filter(Dataset != "Y")
  }
  
  gg1 <- ggplot2::ggplot(dat_comb, aes(x = Dataset, y = `Cumulative Explained Variance`)) +
    ggplot2::geom_col(fill = "darkblue") + 
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 12),
      axis.line  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12)) + 
    ggplot2::ylab("") + 
    ggplot2::ggtitle("Cumulative explained variance")  
  
  return(gg1)
}

plot_MO_2 <- function(Data_res){
  dat_explained <- reshape2::melt(do.call("rbind", Data_res$prop_expl_var))
  colnames(dat_explained) <- c("Dataset", "Component", "% of explained variance")
  dat_explained$`% of explained variance` <- dat_explained$`% of explained variance`*100
  
  if (is(Data_res, "block.splsda") || is(Data_res, "block.plsda")) {
    dat_explained <- dat_explained %>% dplyr::filter(Dataset != "Y")
  }
  
  # Chunk of code to be cohesive with MOFA2::plot_explained_variance
  gg2 <- ggplot2::ggplot(dat_explained, aes(x = Dataset, y = Component)) +
    ggplot2::geom_tile(aes(fill = `% of explained variance`)) + 
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 12),
      axis.line  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12),
    ) + ggplot2::ylab("") +  
    ggplot2::scale_fill_gradientn(colors = c("gray97","darkblue"), guide = "colorbar", 
                                  limits = c(min(dat_explained$`% of explained variance`),
                                             max(dat_explained$`% of explained variance`))) + 
    ggplot2::ggtitle("Percentage of explained variance \n per component per block")
  
  return(gg2)
}


#---- MOFA2 - plot network ----

#' @title Plot correlation network between omics from MOFA2 results
#'
#' @param object a MOFA result
#' @param factor_choice chose weight from this factor
#' @param abs_weight_network threshold of weight to select entities to print
#' @param abs_min_cor_network correlation threshold 
#' @param network_layout one of spring, circle , circle + omics
#' @param omics_colors named list of colors palettes, one for each omics (can be NULL) 
#' @param posCol colors of positive edges
#' @param negCol colors of negative edges
#' 
#' @return A ggplot2 graph
#' @export
#' 
MOFA_cor_network <- function(resMOFA,
                             factor_choice = 1,
                             abs_weight_network = 0.5,
                             abs_min_cor_network = 0.5,
                             network_layout = "Circle + omics",
                             omics_colors = NULL,
                             posCol = "red",
                             negCol = "blue"
){
  
  # Correlation matrix is done on all ZW, not on the selected factor. 
  data_reconst_list <- lapply(MOFA2::get_weights(resMOFA), FUN = function(mat){
    MOFA2::get_factors(resMOFA)$group1 %*% t(mat)})
  data_reconst <- do.call(cbind, data_reconst_list)
  cor_mat <- stats::cor(data_reconst)
  
  
  features_metadata <- do.call(rbind, lapply(1:length(MOFA2::get_weights(resMOFA)), FUN = function(i){
    mat_weights <- data.frame(MOFA2::get_weights(resMOFA, scale = TRUE)[[i]])
    mat_weights$Table <- names(MOFA2::get_weights(resMOFA))[i]
    return(mat_weights)
  }))
  
  factor_selected <- paste0("Factor", factor_choice)
  
  feature_filtered <- features_metadata %>% 
    tibble::rownames_to_column("EntityName") %>%
    dplyr::mutate(F_selected = abs(get(factor_selected))) %>% 
    dplyr::arrange(desc(abs(F_selected))) %>% 
    dplyr::group_by(Table) %>% 
    dplyr::filter(abs(F_selected) > abs_weight_network)
  
  if (nrow(feature_filtered) > 0) {
    
    if(is.null(omics_colors)){
      omics_colors <- as.list(rownames(brewer.pal.info)[1:length(unique(feature_filtered$Table))]) 
      names(omics_colors) <- unique(feature_filtered$Table)
    }
    
    omics_colors <- lapply(unique(feature_filtered$Table), FUN = function(omicTable){
      RColorBrewer::brewer.pal(name = omics_colors[[omicTable]], n = 5)
    })
    names(omics_colors) <- unique(feature_filtered$Table)
    
    feature_filtered <- feature_filtered %>% dplyr::group_by(Table) %>%
      dplyr::mutate(Color = cut(F_selected, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
    feature_filtered$Color2 <- sapply(1:nrow(feature_filtered), 
                                      FUN = function(i) omics_colors[[feature_filtered$Table[i]]][as.numeric(feature_filtered$Color[i])])
    
    # Layout
    layout_arg <- tolower(network_layout)
    if (tolower(layout_arg) == tolower("Circle + omics")) {
      layout_arg <- "groups"
    }
    
    # Network main graph
    cor_display <- cor_mat[rownames(cor_mat) %in% feature_filtered$EntityName, colnames(cor_mat) %in% feature_filtered$EntityName]
    
    if (any(abs(cor_display[upper.tri(cor_display)]) >= abs_min_cor_network)) {
      qgraph_plot <- qgraph::qgraph(cor_display, minimum = abs_min_cor_network, 
                                    cut = 0,
                                    shape = "rectangle", labels = rownames(cor_display), vsize2 = 2, 
                                    vsize = sapply(rownames(cor_display), nchar)*1.1,  layout = layout_arg,
                                    esize = 2,
                                    groups = gsub("[.]filtred", "", features_metadata$Table[match(rownames(cor_display), rownames(features_metadata))]),
                                    posCol = posCol, negCol = negCol, 
                                    details = FALSE,  legend = FALSE,
                                    color = feature_filtered$Color2[match(rownames(cor_display), feature_filtered$EntityName)])
      qgraph_plot <- recordPlot()
      
      # Legend
      legend_matrix <- do.call("rbind", omics_colors)
      colnames(legend_matrix) <- c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8, 1]")
      rownames(legend_matrix) <- gsub("[.]filtred", "", rownames(legend_matrix))
      legend.reshape <- reshape2::melt(legend_matrix)
      
      gg.legend <-  ggplot2::ggplot(legend.reshape, ggplot2::aes(x = Var2, y = Var1)) + 
        ggplot2::geom_tile(fill = legend.reshape$value) + ggplot2::xlab("") + ggplot2::ylab("") + 
        ggplot2::theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(), 
                                    axis.ticks.y = ggplot2::element_blank(),
                                    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
      
      gg.legend <- ggpubr::ggarrange(gg.legend, nrow = 3, ncol = 1) 
      
      # Actual plotting
      ggpubr::ggarrange(plotlist = list(qgraph_plot, gg.legend), nrow = 1, ncol = 2, widths = c(3, 1))
      
      
      
    }else{
      print("There is nothing to plot. Please try to lower the absolute weight, the correlation threshold or change the factor.")
    }
  }else{
    print("There is nothing to plot. Please try to lower the absolute weight, the correlation threshold or change the factor.")
  }
}


