




# ---- Get DE matrix from DiffExpAnalysis ----


#' @title Get DE matrix
#'
#' @param object a SE object (produced by Flomics)
#' @return a matrix of results from the differential analyses.
#' @export
#' @rdname getDE
#'
getDEMatrix <- function(object) {
  if (is(object, "RflomicsSE")) {
    if (!is.null(object@metadata$DiffExpAnal$mergeDEF)) {
      object@metadata$DiffExpAnal$mergeDEF
    } else {
      stop("There is no DE matrix in this object.")
    }
  } else {
    stop("object is not a RflomicsSE.")
  }
}

#' @param contrast character name (can be a vector of name) for the contrast to select.
#' @param union Boolean value. TRUE : union; FALSE : intersection
#' @export
#' @importFrom tidyselect any_of
#' @importFrom dplyr select
#' @rdname getDE
getDE <- function(object, contrast, union = TRUE) {
  
  if (isContrastName(object, contrast)) 
    contrast <- convertContrastToTag(object, contrast)
  
  DEmat <- getDEMatrix(object)
  DEmat <- DEmat %>% select(any_of(c("DEF", contrast)))
  
  if (union) return(DEmat[rowSums(DEmat[,-1]) >= 1,])
  
  return(DEmat[rowSums(DEmat[,-1]) >= length(contrast),])
}

# ---- Get union or intersection from list of contrasts ----

# very similar to filter_DE_from_SE but returns a vector instead of a SE.

#' @title Operation on differential analyses lists. Get union vector of DE entities from list of contrasts
#'
#' @param object a SE object (produced by Flomics). Expects to find a slot with differential analyses results.
#' @param contrasts Vector of characters, expect to be contrast names. Default is null, the operation (union) is performed
#' on every contrasts found.
#' @param operation character. Either union or intersection.
#' Defines the operation to perform on the DE lists from the contrasts.
#' @return vector of unique DE entities
#' @rdname getDE
#' @export
#' @importFrom tidyselect starts_with any_of
#' @importFrom dplyr select mutate filter
opDEList <- function(object, SE.name = NULL, contrasts = NULL, operation = "union") {
  
  if (!is(object, "RflomicsMAE") && !is(object, "RflomicsSE")) 
    stop("Object is not a RflomicsSE or a RflomicsMAE")
  
  if (is(object, "RflomicsMAE") && is.null(SE.name)) 
    stop("Please provide SE.name argument.")
  
  if (is(object, "RflomicsMAE")) object <- object[[SE.name]] 
  
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
}





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

# ---- INTERNAL - get members of a cluster or coseq clusters ----
#
#' @title get members of a cluster
#'
#' @param object a RflomicsSE, produced by rflomics
#' @param name name of the cluster
#' @return The list of entities inside this cluster.
#' @noRd
#' @importFrom coseq clusters
#' @keywords internal

.getCluster <- function(object, clusterName) {
  
  clusterName <- gsub("cluster[.]", "", clusterName)
  res <- object@metadata$CoExpAnal$coseqResults
  
  if (!is.null(res)) {
    clList <- clusters(res)
    return(names(clList == clusterName))
  } else {
    return(NULL)
  }
  
}

#' @param object a RflomicsSE, produced by rflomics
#' @return all clusters
#' @noRd
#' @importFrom coseq clusters
#' @keywords internal

.getCoseqClusters <- function(object) {
  
  res <- object@metadata$CoExpAnal$coseqResults
  
  if (!is.null(res)) {
    return(clusters(res))
  } else {
    return(NULL)
  }
  
}



