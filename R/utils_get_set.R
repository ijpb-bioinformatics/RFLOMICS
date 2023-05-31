
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
  if(class(object) == "MultiAssayExperiment"){
    object@metadata$design@ExpDesign
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
  if(class(object) == "MultiAssayExperiment"){
    object@metadata$design@Contrasts.Sel
  }else if(class(object) == "SummarizedExperiment"){
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
  
  if(is.character(contrasts)){
    if(class(object) == "SummarizedExperiment"){
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
  
  if(class(object) == "SummarizedExperiment"){
    return(object@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName) 
  }else if(class(object) == "MultiAssayExperiment"){
    list_res <- lapply(names(object), FUN = function(tableName){
      object[[tableName]]@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName
    })
    names(list_res) <- names(object)
    return(list_res)
  }
}

# ---- Get summary from ORA : ----

#' @title Get summary tables from ORA analyses - once an enrichment has been conducted. 
#'
#' @param object a SE object (produced by Flomics)
#' @param ont either NULL, GO, KEGG or custom. if NULL, all tables are returned in a list. 
#' @param from either DiffExpEnrichAnal or CoExpAnal. 
#' @return a list of tables or a table 
#' @export
#'
#' @examples 
#' 
sumORA <- function(SE, ont = NULL, from = "DiffExpEnrichAnal"){
  
  if(!is.null(ont)){
    return(SE@metadata[[from]][[ont]]$summary)
  }else{
    list_res <- lapply(names(SE@metadata[[from]]), FUN = function(ontres){
      SE@metadata[[from]][[ontres]]$summary
    })
    names(list_res) <- names(SE@metadata[[from]])
    return(list_res)
  }
} 

# ----- Get Summary for diffExpAnalysis : -----

#' @title Get summary table from diffExpAnalysis analysis
#'
#' @param object a SE object (produced by Flomics)
#' @return a table 
#' @export
#'
#' @examples 
#' 

sumDiffExp <- function(object){
  
  pcut <- object@metadata$DiffExpAnal$Adj.pvalue.cutoff
  lcut <- object@metadata$DiffExpAnal$abs.logFC.cutoff
  
  df_sim <- lapply(object@metadata$DiffExpAnal$DEF, FUN = function(tab){

    tab  <- tab %>% dplyr::filter(Adj.pvalue < pcut) %>%
      dplyr::filter(abs(logFC) > lcut)
    
    c("All"  = nrow(tab), 
      "Up"   = nrow(tab %>% dplyr::filter(logFC>0)), 
      "Down" = nrow(tab %>% dplyr::filter(logFC<0))
    )
  })
  return(do.call("rbind", df_sim))
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
  
  if(!any(search_success)){
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
  
  if(!any(search_success)){
    # Congratulations, it's a tag name!
    return(TRUE)
  }else return(FALSE)
  
}


# ---- convert tag to contrastName ----

#' @title Convert tags names to contrast Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expect to find 
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
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expect to find 
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
