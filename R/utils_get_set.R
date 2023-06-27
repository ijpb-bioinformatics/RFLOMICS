# ---- Get Factor types : ----
#' @title Get design matrix used for a differential analysis
#'
#' @param object a MAE object (produced by Flomics)
#' @return a dataframe
#' @export
#'
#' @examples 
#' 

getFactorTypes <- function(object){
  if (class(object) == "MultiAssayExperiment") {
    object@metadata$design@Factors.Type
  }else{
    stop("object is not a MultiAssayExperiment.")
  }
}

# TODO add more getters for accessing bio, batch and meta directly.


# ---- Get Design Matrix : ----
#' @title Get design matrix used for a differential analysis
#'
#' @param object a MAE object (produced by Flomics)
#' @return a dataframe
#' @export
#'
#' @examples 
#' 

getDesignMat <- function(object){
  # TODO check if it exists... 
  if (class(object) == "MultiAssayExperiment") {
    object@metadata$design@ExpDesign
  }else{
    stop("object is not a MultiAssayExperiment.")
  }
}

# ---- Get Model Formula : ----
#' @title Get model formula from a Flomics multiassayexperiment. 
#'
#' @param object a MAE object (produced by Flomics)
#' @return a formula
#' @export
#'
#' @examples 
#' 

getModelFormula <- function(object){
  # TODO check if it exists... 
  if (class(object) == "MultiAssayExperiment") {
    object@metadata$design@Model.formula
  }else{
    stop("object is not a MultiAssayExperiment.")
  }
}

# ---- Get possible contrasts : ----
#' @title Get selected contrasts for the differential analysis
#'
#' @param object a MAE object (produced by Flomics) or a SE (expect to find a diffAnalysis slot.)
#' @return a dataTable
#' @export
#'
#' @examples 
#' 

getSelectedContrasts <- function(object){
  # TODO check if it exists... 
  if (class(object) == "MultiAssayExperiment") {
    object@metadata$design@Contrasts.Sel
  }else if (class(object) == "SummarizedExperiment") {
    object@metadata$DiffExpAnal$contrasts
  }else{
    stop("object is not a MultiAssayExperiment or a SummarizedExperiment.")
  }
}

# ---- Set Valid Contrasts : (after differential analysis) ----
#' @title Set Valid Contrasts
#'
#' @param object a SE object or MAE object (produced by Flomics)
#' @param contrasts a vector of contrasts names
#' @return a Flomics SE or MAE
#' @export
#'
#' @examples 
#' 

setValidContrasts <- function(object, 
                              contrasts){
  
  # TODO : check if there are DE entities for each contrasts before really validating them.
  
  if (is.character(contrasts)) {
    if (class(object) == "SummarizedExperiment") {
      object@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName <- contrasts
    }else{
      stop("object is not a SummarizedExperiment.")
    }
  }
  
  return(object)
}

#' @title Get Valid Contrasts
#'
#' @param object a SE object or MAE object (produced by Flomics)
#' @return a list of vectors (if object is a MAE) or a vector of contrasts names (if object is a SE)
#' @export
#'
#' @examples 
#' 
getValidContrasts <- function(object){
  
  if (class(object) == "SummarizedExperiment") {
    return(object@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName) 
  }else if (class(object) == "MultiAssayExperiment") {
    list_res <- lapply(names(object), FUN = function(tableName){
      object[[tableName]]@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName
    })
    names(list_res) <- names(object)
    return(list_res)
  }
}




# ---- Get DE matrix from DiffExpAnalysis ----


#' @title Get DE matrix
#'
#' @param object a SE object (produced by Flomics)
#' @return a matrix of results from the differential analyses. 
#' @export
#'
#' @examples 
#' 
getDEMatrix <- function(object){
  
  if (class(object) == "SummarizedExperiment") {
    if (!is.null(object@metadata$DiffExpAnal$mergeDEF)) return(object@metadata$DiffExpAnal$mergeDEF) 
    else stop("There is no DE matrix in this object.")
  }else stop("object is not a SummarizedExperiment.")
  
}

# ----- Check if character vectors are contrasts Names : -----

#' @title Check if character vectors are contrasts Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expect to find 
#'  a slot of differential analysis. 
#' @param contrastName vector of characters. 
#' @return boolean. TRUE if all of contrastName are indeed contrasts Names. 
#' @export
#'
#' @examples 
#' 
isContrastName <- function(object, contrastName){
  
  df_contrasts <- RFLOMICS::getSelectedContrasts(object)
  
  search_match   <- sapply(contrastName, FUN = function(cn){grep(cn, df_contrasts$contrastName, fixed = TRUE)})
  search_success <- sapply(search_match, identical, integer(0)) # if TRUE, not a success at all. 
  
  if (!any(search_success)) {
    # Congratulations, it's a contrast name!
    return(TRUE)
  }else return(FALSE)
  
}

# ----- Check if character vectors are tags Names : -----

#' @title Check if character vectors are tags Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expect to find 
#'  a slot of differential analysis. 
#' @param tagName vector of characters. 
#' @return boolean. TRUE if all of tagName are indeed tags Names. 
#' @export
#'
#' @examples 
#' 
isTagName <- function(object, tagName){ 
  
  df_contrasts <- RFLOMICS::getSelectedContrasts(object)
  
  search_match   <- sapply(tagName, FUN = function(cn){grep(cn, df_contrasts$tag, fixed = TRUE)})
  search_success <- sapply(search_match, identical, integer(0)) # if TRUE, not a success at all. 
  
  if (!any(search_success)) {
    # Congratulations, it's a tag name!
    return(TRUE)
  }else return(FALSE)
  
}


# ---- convert tag to contrastName ----

#' @title Convert tags names to contrast Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expects to find 
#'  a slot of differential analysis. 
#' @param tagName Vector of characters, expect to be tags (in the form of H1, H2, etc.).
#' @return character vector, contrastNames associated to tags.
#' @export
#'
#' @examples 
#' 
convertTagToContrast <- function(object, tagName){
  
  df_contrasts <- RFLOMICS::getSelectedContrasts(object)
  
  df_contrasts %>% 
    dplyr::filter(tag %in% tagName) %>% 
    dplyr::select(contrastName) %>% 
    unlist(use.names = FALSE)
  
}

# ---- convert contrastName to tag ----

#' @title Convert contrast Names names to tags
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expects to find 
#'  a slot of differential analysis. 
#' @param contrasts Vector of characters, expect to be contrast names. 
#' @return character vector, tags associated to contrast names.
#' @export
#'
#' @examples 
#' 
convertContrastToTag <- function(object, contrasts){
  
  df_contrasts <- RFLOMICS::getSelectedContrasts(object)
  
  df_contrasts %>% 
    dplyr::filter(contrastName %in% contrasts) %>% 
    dplyr::select(tag) %>% 
    unlist(use.names = FALSE)
  
}


# ---- Get union or intersection from list of contrasts ----

# very similar to filter_DE_from_SE but returns a vector instead of a SE. 

#' @title Get union vector of DE entities from list of contrasts
#'
#' @param object a SE object (produced by Flomics). Expects to find a slot with differential analyses results.
#' @param contrasts Vector of characters, expect to be contrast names. Default is null, the operation (union) is performed
#' on every contrasts found.
#' @param operation character. Either union or intersection. 
#' Defines the operation to perform on the DE lists from the contrasts.
#' @return vector of unique DE entities
#' @export
#'
#' @examples 
#' 
opDEList <- function(object, contrasts = NULL, operation = "union"){
  
  # object <- MAE[["RNAseq_norm"]]
  # contrasts <- NULL
  
  if (class(object) != "SummarizedExperiment") stop("Object is not a SummarizedExperiment")
  if (is.null(object@metadata$DiffExpAnal$Validcontrasts)) stop("Please validate your differential analyses first.")
  
  if (is.null(contrasts)) contrasts <- RFLOMICS::getSelectedContrasts(object)[["tag"]]
  if (RFLOMICS::isContrastName(object, contrasts)) contrasts <- RFLOMICS::convertContrastToTag(object, contrasts)
  
  validTags <- RFLOMICS::convertContrastToTag(object, RFLOMICS::getValidContrasts(object))
  
  tagsConcerned <- intersect(contrasts, validTags)
  if (length(tagsConcerned) == 0) stop("It seems there is no contrasts to select DE entities from.")
  
  df_DE <- RFLOMICS::getDEMatrix(object) %>% 
    dplyr::select(c("DEF", tidyselect::any_of(tagsConcerned)))
  
  if (operation == "intersection") {
    
    DETab <- df_DE %>%
      dplyr::mutate(SUMCOL = dplyr::select(., tidyselect::starts_with("H")) %>% 
                      rowSums(na.rm = TRUE))  %>%
      dplyr::filter(SUMCOL == length(validTags))
    
  }else{
    
    DETab <- df_DE %>%
      dplyr::mutate(SUMCOL = dplyr::select(., tidyselect::starts_with("H")) %>% 
                      rowSums(na.rm = TRUE))  %>%
      dplyr::filter(SUMCOL >= 1)
  }
  
  return(DETab$DEF)
  
}


# ---- Get omics experiments and their types (vector) ----

#' @title Get omics experiments and their types
#'
#' @param object a MAE object (produced by Flomics) or a Summarized Experiment object. 
#' @return a named vector with each omics name and its type.
#' @export
#'
#' @examples 
#' 

getOmicsTypes <- function(object){
  
  if (!class(object) %in% c("MultiAssayExperiment", "SummarizedExperiment"))
    stop("Object is not a MultiAssayExperiment nor a SummarizedExperiment")
  
  if (class(object) == "MultiAssayExperiment") {
    sapply(names(object), FUN = function(x){
      object[[x]]@metadata$omicType
    })
  }else{
    object@metadata$omicType
  }
  
}

# ---- Get normalization coefficients ----

#' @title Get normalization coefficients
#'
#' @param object a SE object (produced by Flomics). 
#' @return Normalisation coefficient. If TMM was applied, a list with library size and coefficients.
#' @export
#'
#' @examples 
#' 

getNormCoeff <- function(object){
  
  if (class(object) != "SummarizedExperiment") stop("Object is not a SummarizedExperiment")
  
  return(object@metadata$Normalization$coefNorm)
}



