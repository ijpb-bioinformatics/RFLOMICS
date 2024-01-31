






# ---- Get filtering parameters ----

#' @title Get filtering setting
#'
#' @param object a SE object (produced by Flomics).
#' @return List of filtering setting
#' @export
#'
getFilterSetting <- function(object) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  return(object@metadata$DataProcessing$Filtering$setting)
}

# ---- Get filtred features ----

#' @title Get filtering setting
#'
#' @param object a SE object (produced by Flomics).
#' @return List of filtered features
#' @export
#'
getFilteredFeatures <- function(object) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  return(object@metadata$DataProcessing$Filtering$results$filteredFeatures)
}

# ---- Get normalization coefficients ----

#' @title Get normalization coefficients
#'
#' @param object a SE object (produced by Flomics).
#' @return Normalisation coefficient. If TMM was applied, a list with library size and coefficients.
#' @export
#'
getCoeffNorm <- function(object) {
  metadata(object)[["DataProcessing"]][["Normalization"]][["results"]][["coefNorm"]]
}


# ---- Get normalization parameters ----

#' @title Get normalization parameters
#'
#' @param object a SE object (produced by Flomics).
#' @return List of normalization parameters.
#' @export
#'
getNormSetting <- function(object) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  return(object@metadata$DataProcessing$Normalization$setting)
}

# ---- Get transformation parameters ----

#' @title Get transformation parameters
#'
#' @param object a SE object (produced by Flomics).
#' @return List of transformation parameters.
#' @export
#'
getTransSetting <- function(object) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  return(object@metadata$DataProcessing$Transformation$setting)
}


# ---- Get diff setting ----

#' @title Get differential analysis setting parameters
#'
#' @param object a SE object (produced by Flomics).
#' @return List of differential analysis setting parametres.
#' @export
#'
getDiffSetting <- function(object) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  return(object@metadata$DiffExpAnal$setting)
}


# ---- Get coseq setting ----

#' @title Get co-expression analysis setting parameters
#'
#' @param object a SE object (produced by Flomics).
#' @return List of co-expression analysis setting parametres.
#' @export
#'
getCoexpSetting <- function(object) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  return(object@metadata$CoExpAnal$setting)
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



