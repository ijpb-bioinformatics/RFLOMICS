#---- Filter DE entities from SE (returns vector) ----

# TODO A voir pour remplacer par un getter peut-être
# TODO Est-ce que cette fonction est toujorus utilisée ? Remplacée par opDEList logiquement.

#' @title filter_DE_from_SE
#'
#' @param SEobject An object of class \link{SummarizedExperiment}
#' @param contrasts_arg
#' @param type
#' @return An object of class \link{SummarizedExperiment}
#' @export
#' @examples
#' @noRd
#'

filter_DE_from_SE <- function(SEobject, contrasts_arg = NULL, type = "union"){
  
  tabCorresp <- SEobject@metadata$DiffExpAnal$Validcontrasts %>% 
    dplyr::select(contrastName, tag)
  if(is.null(contrasts_arg))   contrasts_arg <- SEobject@metadata$DiffExpAnal$Validcontrasts$contrastName
  
  tabCorresp <- tabCorresp %>% dplyr::filter(contrastName %in% contrasts_arg)
  contrasts_select <- tabCorresp$tag
  
  tab1 <- SEobject@metadata$DiffExpAnal$mergeDEF %>%
    dplyr::select(any_of(c("DEF", contrasts_select)))
  
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
#' @examples
#' @noRd
#'
rbe_function = function(object, SEobject){
  
  assayTransform <- SummarizedExperiment::assay(SEobject)
  ftypes <- RFLOMICS::getFactorTypes(object)
  colBatch <- names(ftypes)[ftypes == "batch"]
  
  print(paste0("#     =>Correction for Batch: ", paste(colBatch, collapse = " "), " in ", SEobject@metadata$omicType))
  
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
#' @examples
#' @noRd
#'

rnaseqRBETransform <- function(object, 
                               SEname, 
                               correctBatch = FALSE, 
                               transformation = "limma (voom)",
                               contrasts_names = NULL,
                               type = "union",
                               choice = "DE"
){
  
  # TODO : check if differential analysis has been performed.
  # TODO : in case of removeBatchEffect not working, what do we do?
  # TODO : if no DE (ex: intersection is 0 genes), what do we do?
  
  if (class(object) != "MultiAssayExperiment") stop("object is not a MultiAssyExperiment")
  
  rnaDat <- object[[SEname]] 
  assayTransform <- SummarizedExperiment::assay(rnaDat)
  if (!is.integer(assayTransform)) {
    message(paste0("You indicated RNASeq data for ", SEname, "but it is not recognized as count data")) 
    print(assayTransform[1:3, 1:3])
  }
  
  DMat      <- RFLOMICS::getDesignMat(object)
  coefNorm  <- RFLOMICS::getNormCoeff(rnaDat)
  designMat <- model.matrix(as.formula(RFLOMICS::getModelFormula(object)), data = DMat)
  
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
  
  if (correctBatch)    rnaDat <- RFLOMICS::rbe_function(object, rnaDat)
  if (choice == "DE")  rnaDat <- rnaDat[RFLOMICS::opDEList(rnaDat, contrasts = contrasts_names, operation = type),]
  
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
#' @examples
#' @noRd
#'

RBETransform <- function(object,
                         SEname,
                         correctBatch = TRUE,
                         contrasts_names = NULL,
                         type = "union",
                         choice = "DE"){
  
  omicsDat <- object[[SEname]]
  omicsDat@metadata[["correction_batch"]]             <- correctBatch
  omicsDat@metadata[["transform_method_integration"]] <- omicsDat@metadata$transform_method
  
  if (correctBatch)    omicsDat <- RFLOMICS::rbe_function(object, omicsDat)
  if (choice == "DE")  omicsDat <- omicsDat[RFLOMICS::opDEList(omicsDat, contrasts = contrasts_names, operation = type),]
  
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
#' @examples
#'

plot_MO_varExp <- function(object, selectedResponse, mode = NULL){
  
  if (class(object) != "MultiAssayExperiment") stop("Object is not a MultiAssayExperiment.")
  if (is.null(object@metadata$mixOmics)) stop("It seems this object has no mixOmics results.")
  if (is.null(object@metadata$mixOmics[[selectedResponse]])) stop("It seems you didn't run MixOmics on this particular variable.")
  
  Data_res <- object@metadata$mixOmics[[selectedResponse]]$MixOmics_results
  
  if (is.null(mode)) {
    return(ggpubr::ggarrange(RFLOMICS:::plot_MO_1(Data_res),
                             RFLOMICS:::plot_MO_2(Data_res), 
                             ncol = 2))
  }
  else if (toupper(mode) == "cumulative") {
    return(RFLOMICS:::plot_MO_1(Data_res)) 
  }
  else if (toupper(mode) == "comp") {
    return(RFLOMICS:::plot_MO_2(Data_res))
  }
}

#' @title Plot cumulative Explained variance for mixOmics results
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string. 
#' @return A ggplot2 graph
#' @keywords internal
#' @NoRd

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
}
