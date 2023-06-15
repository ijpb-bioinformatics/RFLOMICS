######## INTEGRATION WITH MOFA #############

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

filter_DE_from_SE <- function(SEobject, contrasts_arg, type = "union"){
  
  tabCorresp <- SEobject@metadata$DiffExpAnal$Validcontrasts %>% dplyr::select(contrastName, tag)
  if("all" %in% contrasts_arg)   contrasts_arg <- SEobject@metadata$DiffExpAnal$Validcontrasts$contrastName
  
  tabCorresp <- tabCorresp %>% dplyr::filter(contrastName %in% contrasts_arg)
  contrasts_select <- tabCorresp$tag
  
  tab1 <- SEobject@metadata$DiffExpAnal$mergeDEF %>%
    dplyr::select(any_of(c("DEF", contrasts_select)))
  
  if(type == "intersection"){
    
    DETab <- tab1 %>%
      dplyr::mutate(SUMCOL = dplyr::select(., starts_with("H")) %>% rowSums(na.rm = TRUE))  %>%
      dplyr::filter(SUMCOL==length(contrasts_select))
    
  }else{
    
    DETab <- tab1 %>%
      dplyr::mutate(SUMCOL = dplyr::select(., starts_with("H")) %>% rowSums(na.rm = TRUE))  %>%
      dplyr::filter(SUMCOL>=1)
  }
  
  SEobject <- SEobject[DETab$DEF,]
  
  return(SEobject)
}


# rbe_function : pour corriger sur le batch effect.
# object: l'objet flomics
# SEtransform : le SE qui contient le tableau transforme (la case metadata[["transform_results"]])
# TODO A MODIFIER ON DOIT POUVOIR PRENDRE EN COMPTE TOUS LES BATCH EFFECTS
#' @title rbe_function
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEobject An object of class \link{SummarizedExperiment}
#' @return An object of class \link{MultiAssayExperiment}
#' @export
#' @examples
#' @noRd
#'
rbe_function = function(object, SEobject){
  
  # object = object
  # SEobject = omicsDat
  
  assayTransform <- SummarizedExperiment::assay(SEobject)
  
  colBatch <- names(object@metadata$design@Factors.Type)[object@metadata$design@Factors.Type=="batch"]
  
  print(paste0("#     =>Correction for Batch: ", paste(colBatch, collapse = " "), " in ", SEobject@metadata$omicType))
  
  newFormula <- gsub(pattern = paste(colBatch, collapse = "[+]|"), "", object@metadata$design@Model.formula)
  newFormula <- gsub(pattern = "~ [+] ", "~ ", newFormula)
  # designToPreserve <- model.matrix(as.formula(newFormula), data = object@metadata$design@ExpDesign)
  designToPreserve <- model.matrix(as.formula(newFormula), data = SEobject@metadata$Groups)
  
  if(length(colBatch)==1){
    rbeRes <- limma::removeBatchEffect(assayTransform, batch = SEobject@metadata$Groups[,colBatch], design = designToPreserve)
  }else if(length(colBatch) >= 2){
    
    rbeRes <- limma::removeBatchEffect(assayTransform,
                                       batch = SEobject@metadata$Groups[,colBatch[1]],
                                       batch2 = SEobject@metadata$Groups[,colBatch[2]],
                                       design = designToPreserve)
  }
  # else{
  if(length(colBatch) > 2) print("sorry, only 2 batches effect for now!!!") # trouver un moyen de prendre en compte automatiquement plusieurs batch factors. C'est moche.
  # }
  
  SEobject@metadata[["correction_batch_method"]] <- "limma (removeBatchEffect)"
  
  # SEobject@metadata[["integration_table"]] <- rbeRes # on ecrase le tableau de resultats
  SummarizedExperiment::assay(SEobject) <- rbeRes
  
  
  return(SEobject)
}
