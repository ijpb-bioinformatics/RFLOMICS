# ----- Get Summary for diffExpAnalysis : -----

#' @title Get summary table from diffExpAnalysis analysis
#'
#' @param object a SE object (produced by Flomics) or a MAE
#' @return a table 
#' @export
#'
#' @examples 
#' 

sumDiffExp <- function(object, SE.name = NULL){
  
  # TODO valid contrasts only by default
  # TODO if asked by user (all = TRUE), then provide all results even non-validated ones.
  
  if (class(object) == "MultiAssayExperiment") {
    if (!is.null(SE.name)) {
      object <- object[[SE.name]]
    }
  }
  
  pcut <- object@metadata$DiffExpAnal$Adj.pvalue.cutoff
  lcut <- object@metadata$DiffExpAnal$abs.logFC.cutoff
  n_entities <- nrow(SummarizedExperiment::assay(object))
  n_samples  <- ncol(SummarizedExperiment::assay(object))
  
  cat("Parameters:\n|adjusted-pvalue cutoff: ", pcut, 
      "\n|logFC cutoff: ", lcut, 
      "\n|number of features: ", n_entities,
      "\n|number of samples: ", n_samples, "\n")
  
  df_sim <- lapply(object@metadata$DiffExpAnal$DEF, FUN = function(tab){
    
    tab  <- tab %>% dplyr::filter(Adj.pvalue < pcut) %>%
      dplyr::filter(abs(logFC) > lcut)
    
    c("All"  = nrow(tab), 
      "Up"   = nrow(tab %>% dplyr::filter(logFC > 0)), 
      "Down" = nrow(tab %>% dplyr::filter(logFC < 0))
    )
  })
  return(do.call("rbind", df_sim))
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
sumORA <- function(SE,  from = "DiffExpEnrichAnal", ont = NULL, contrast = NULL){
  
  if (!class(SE) %in% "SummarizedExperiment") stop("SE is not a summarizedExperiment.")
  
  if(!is.null(contrast)){
    if(RFLOMICS::isTagName(SE, contrast)) contrast <- RFLOMICS::convertTagToContrast(SE, contrast)
  }
  
  if (!is.null(ont)) {
    toReturn <- SE@metadata[[from]][[ont]]$summary
    if(!is.null(contrast)) toReturn <- toReturn[which(toReturn$Contrast == contrast),]
    return(toReturn)
  }else{
    list_res <- lapply(names(SE@metadata[[from]]), FUN = function(ontres){
      interRes <- SE@metadata[[from]][[ontres]]$summary
      if(!is.null(contrast)) interRes <- interRes[which(interRes$Contrast == contrast),]
      interRes
    })
    names(list_res) <- names(SE@metadata[[from]])
    return(list_res)
  }
} 

# ---- MixOmics summary ----

#' @title Get an overview of MixOmics integration results
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character. Useful if MixOmics was run on several response variable. If NULL, all variables are taken into account. 
#' @return A data frame or a list of dataframe (if selectedResponse is NULL).
#' @export
#'
#' @examples 

sumMixOmics <- function(object, selectedResponse = NULL){
  
  if (class(object) != "MultiAssayExperiment") stop("Object is not a MultiAssayExperiment.")
  if (is.null(object@metadata$mixOmics)) stop("It seems this object has no mixOmics results.")
  
  if (is.null(selectedResponse)) {
    
    res <- lapply(names(object@metadata$mixOmics), FUN = function(selResponse) {
      RFLOMICS:::getOneMORes(object, selectedResponse = selResponse)
    })
    names(res) <- names(object@metadata$mixOmics)
    return(res)
    
  }else{
    RFLOMICS:::getOneMORes(object, selectedResponse = selectedResponse)
  }
}

#' @title Get an overview of MixOmics integration results for a specific response variable. 
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string. 
#' @return A data frame.
#' @keywords internal
#' @NoRd
getOneMORes <- function(object, selectedResponse){
  Data_res <- object@metadata$mixOmics[[selectedResponse]]$MixOmics_results
  
  df <- t(sapply(Data_res$X, dim))
  colnames(df) <- c("Ind", "Features")
  
  if (!is.null(object@metadata$mixOmics[[selectedResponse]]$MixOmics_results$MixOmics_tuning_results)) {
    df <- cbind(df, do.call("rbind", Data_res$keepX))
    colnames(df)[!colnames(df) %in% c("Ind", "Features")] <- paste("Comp", 1:length(Data_res$keepX[[1]]))
  }
  
  return(df)
}