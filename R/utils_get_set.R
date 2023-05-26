
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
#' @param object a MAE object (produced by Flomics)
#' @return a dataTable
#' @export
#'
#' @examples 
#' 

getSelectedContrasts <- function(object){
  # TODO check if it exists... 
  if(class(object) == "MultiAssayExperiment"){
    object@metadata$design@Contrasts.Sel
  }else{
    stop("object is not a MultiAssayExperiment.")
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




