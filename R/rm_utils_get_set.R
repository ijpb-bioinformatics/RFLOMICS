



# ---- Get integration setting ----

#' @title Get MOFA analysis setting parameters
#'
#' @param object a RflomicsMAE object (produced by Flomics).
#' @return List of differential analysis setting parameters
#' @export
#'
getMOFASetting <- function(object) {
  if (!is(object, "RflomicsMAE")) stop("Object is not a RflomicsMAE")
  
  return(object@metadata$IntegrationAnalysis$MOFA$setting)
}

#' @title Get mixOmics analysis setting parameters
#'
#' @param object a RflomicsMAE object (produced by Flomics).
#' @return List of differential analysis setting parameters
#' @export
#'
getMixOmicsSetting <- function(object) {
  if (!is(object, "RflomicsMAE")) stop("Object is not a RflomicsMAE")
  
  return(object@metadata$IntegrationAnalysis$mixOmics$setting)
}

# ----  Get a particular multi-omics result ----
#
#' @title Get a particular multi-omics result
#'
#' @param object a MAE object (produced by Flomics).
#' @param response a character giving the response variable to access specifically. 
#' @param onlyResults default return only the MixOmics or MOFA2 results. If you want to access all information of the integration, 
#' set onlyResuts to FALSE. In MixOmics case, works only when response is specified.
#' @return in getMixOmics, if response is NULL, then all the mixOmics results are returned. 
#' Otherwise, it gives the particular mixOmics result. 
#' @rdname Multi-omics-access
#' @export

getMixOmics <- function(object,
                        response = NULL,
                        onlyResults = TRUE){
  
  toreturn <- metadata(object)[["IntegrationAnalysis"]][["mixOmics"]]
  
  if (is.null(toreturn)) {
    return(toreturn)
  }
  
  if (!is.null(response)) {
    toreturn <- toreturn[[response]]
    if (onlyResults) toreturn <- toreturn$MixOmics_results
    return(toreturn)
  }else{
    return(toreturn)
  }
  
}

#' @rdname Multi-omics-access
#' @export
getMOFA <- function(object, onlyResults = TRUE){
  
  toreturn <- metadata(object)[["IntegrationAnalysis"]][["MOFA"]]
  
  if (onlyResults && !is.null(toreturn)) toreturn <- toreturn[["MOFA_results"]]
  
  return(toreturn)
}

#' @rdname Multi-omics-access
#' @export
setMOFA <- function(object, results = NULL){
  
  metadata(object)[["MOFA"]] <- results
  
  return(object)
}

#' @rdname Multi-omics-access
#' @export
setMixOmics <- function(object, results = NULL){
  
  metadata(object)[["mixOmics"]] <- results
  
  return(object)
}


# ---- INTERNAL - Get a particular enrichment result ----
#
#' @title Get a particular enrichment result
#'
#' @param object a SE object or a MAE object (produced by Flomics).
#' @return enrichment result.
#' @export
#' @noRd
#' @keywords internal

getEnrichRes <- function(object,
                         contrast = NULL,
                         experiment = NULL,
                         from = "DiffExpEnrichAnal",
                         ont = "GO",
                         domain = NULL) {
  
  if (!is(object, "RflomicsMAE") && !is(object, "RflomicsSE")) 
    stop("Object is not a RflomicsSE or a RflomicsMAE")
  
  if (toupper(from) %in% c("DIFFEXP", "DIFFEXPANAL", "DIFFEXPENRICHANAL")) from <- "DiffExpEnrichAnal"
  if (toupper(from) %in% c("COEXP", "COEXPANAL", "COEXPENRICHANAL"))     from <- "CoExpEnrichAnal"
  
  classObj <- NULL
  if (is(object, "RflomicsSE")) classObj <- "SE"
  if (is(object, "RflomicsMAE")) classObj <- "MAE"
  
  res_return <- NULL
  
  switch(classObj,
         "SE" = {
           if (is.null(contrast)) {
             res_return <- object@metadata[[from]][[ont]][["enrichResult"]]
           } else {
             if (isTagName(object, contrast)) contrast <- convertTagToContrast(object, contrast)
             res_return <- object@metadata[[from]][[ont]][["enrichResult"]][[contrast]]
           }
         },
         "MAE" = {
           if (is.null(experiment)) {
             stop("Please indicate from which data you want to extract the enrichment results.")
           }
           
           if (is.null(contrast)) {
             res_return <- object[[experiment]]@metadata[[from]][[ont]][["enrichResult"]]
           } else {
             if (isTagName(object, contrast)) contrast <- convertTagToContrast(object[[experiment]], contrast)
             res_return <- object[[experiment]]@metadata[[from]][[ont]][["enrichResult"]][[contrast]]
           }
         },
         stop("Object is not a RflomicsMAE nor a RflomicsSE")
  )
  
  if (!is.null(domain) && !is.null(contrast)) {
    return(res_return[[domain]])
  } else {
    return(res_return)
  }
}

# ---- INTERNAL - Get a particular enrichment summary ----
# TODO equivalent to sumORA (external)
#
#' @title Get a particular enrichment result
#'
#' @param object a SE object or a MAE object (produced by Flomics).
#' @return enrichment summary
#' @noRd
#' @keywords internal

getEnrichSum <- function(object,
                         experiment = NULL,
                         from = "DiffExpEnrichAnal",
                         dom = "GO") {
  if (!is(object, "RflomicsMAE") && !is(object, "RflomicsSE")) 
    stop("Object is not a RflomicsSE or a RflomicsMAE")
  
  if (is(object, "RflomicsSE")) {
    return(object@metadata[[from]][[dom]][["summary"]])
  } else if (is(object, "RflomicsMAE")) {
    if (is.null(experiment)) {
      stop("Please indicate from which data you want to extract the enrichment results.")
    }
    
    return(object[[experiment]]@metadata[[from]][[dom]][["summary"]])
  }
}


# ---- INTERNAL - Get a pvalue threshold used in enrichment analysis ----
#
#' @title Get a pvalue threshold used in enrichment analysis
#'
#' @param object a SE object
#' @return pvalue
#' @noRd
#' @keywords internal

getEnrichPvalue <- function(object,
                            from = "DiffExpEnrichAnal",
                            dom = "GO") {
  
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  if (!from %in% c("DiffExpEnrichAnal", "CoExpEnrichAnal")) stop(from, " don't existe")
  if (!dom  %in% c("GO", "KEGG", "custom")) stop(from, " not valide value. Choose from c(GO, KEGG, custom)")
  if (is.null(object@metadata[[from]][[dom]]$list_args$pvalueCutoff)) stop("P-value not found")
  
  return(object@metadata[[from]][[dom]]$list_args$pvalueCutoff)
}



