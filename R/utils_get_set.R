# getModelFormula !!!! GetModelFormulae

# ---- Get Factor types : ----
#' @title Get factor types or names of particular factors given a type.
#'
#' @param object a MAE object (produced by Flomics). 
#' @return a named vector (getFactorTypes) or a vector of factor names.
#' @rdname getFactorTypes
#' @export
#' 
getFactorTypes <- function(object) {
  
  if (is(object, "MultiAssayExperiment")) {
    metadata(object)$design@Factors.Type
  } 
  else if(is(object, "SummarizedExperiment")){
    metadata(object)$design$factorType
  }
  else {
    stop("object is not a MultiAssayExperiment or SummarizedExperiment.")
  }
  
}


#' @title Get bio factor.
#'
#' @param object a MAE object or SE object (produced by Flomics). 
#' @return a vector with biological factors
#' @rdname bioFactors
#' @export
#' 
bioFactors <- function(object){
  
  factVect <- toupper(getFactorTypes(object))
  return(names(factVect)[factVect == "BIO"])
  
}


#' @title Get batch factor.
#'
#' @param object a MAE object or SE object (produced by Flomics). 
#' @return a vector with batch factors
#' @rdname batchFactors
#' @export
#' 
batchFactors <- function(object){
  
  factVect <- toupper(getFactorTypes(object))
  return(names(factVect)[factVect == "BATCH"])
  
}


#' @title Get metadata factor.
#'
#' @param object a MAE object or SE object (produced by Flomics). 
#' @return a vector with metadata factors
#' @rdname metaFactors
#' @export
#' 
metaFactors <- function(object){
  
  factVect <- toupper(getFactorTypes(object))
  return(names(factVect)[factVect == "META"])
  
}

# ---- Get Design Matrix : ----
#' @title Get design matrix used for a differential analysis
#'
#' @param object a MAE object (produced by Flomics)
#' @return a dataframe
#' @export
#'
getDesignMat <- function(object) {
  if (is(object, "MultiAssayExperiment") || is(object, "SummarizedExperiment")) {
    return(as.data.frame(object@colData))
    
  } else {
    stop("object is not a MultiAssayExperiment nor a SummarizedExperiment.")
  }
}


# ---- Get Model Formula : ----
#' @title Get model formula from a Flomics multiassayexperiment.
#'
#' @param object a MAE object (produced by Flomics)
#' @return a formula
#' @export
#'
getModelFormula <- function(object) {
  # TODO check if it exists...
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@Model.formula
  } 
  else if(is(object, "SummarizedExperiment")){
    object@metadata$design$Model.formula
  }
  else {
    stop("object is not a MultiAssayExperiment.")
  }
}

# ---- Set Model Formula : ----
#' @title Set model formula from a Flomics multiassayexperiment.
#'
#' @param object a MAE or SE object (produced by Flomics)
#' @param a formula
#' @return a object SE
#' @export
#'
setModelFormula <- function(object, modelFormula=NULL) {
  # TODO check if it exists...
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@Model.formula <- paste(modelFormula, collapse = " ")
    # set modelFormula foreach SE
    for(se in names(object)){
      object[[se]]@metadata$design$Model.formula <- paste(modelFormula, collapse = " ")
    }
  } 
  else if(is(object, "SummarizedExperiment")){
    object@metadata$design$Model.formula <- paste(modelFormula, collapse = " ")
  }
  else {
    stop("object is not a MultiAssayExperiment.")
  }
  return(object)
}

# ---- Get possible contrasts : ----
#' @title Get selected contrasts for the differential analysis
#'
#' @param object a MAE object (produced by Flomics) or a summarized experiment produced by Flomics after a differential analysis.
#' @param typeContrast the type of contrast from which the possible contrasts are extracted. Default is all contrasts types.
#' @param modalities specific levels for the contrast selection
#' @param returnTable return a dataTable with all contrasts information
#' @return a character vector or a dataTable
#' @rdname ContrastsSelection
#' @export
#' @importFrom data.table rbindlist
#'
getPossibleContrasts <- function(object, typeContrast = c("simple", "averaged", "interaction"),
                                 modalities = NULL, returnTable = FALSE) {
  if (is(object, "MultiAssayExperiment")) {
    if (is.null(typeContrast)) typeContrast <- c("simple", "averaged", "interaction")
    
    allContrasts <- MAE@metadata$design@Contrasts.List
    allContrasts <- allContrasts[which(names(allContrasts) %in% typeContrast)]
    allContrastsdt <- data.table::rbindlist(allContrasts, fill = TRUE)
    
    if (!is.null(modalities)) {
      allVarMod <- lapply(getDesignMat(MAE), FUN = function(vect) levels(factor(vect))[levels(factor(vect)) %in% modalities])
      allVarMod <- Filter(length, allVarMod)
      
      allVarMod <- paste0(rep(names(allVarMod), times = sapply(allVarMod, length)), unlist(allVarMod))
      
      allContrastsdt <- allContrastsdt[grep(paste(allVarMod, collapse = "|"), allContrastsdt$contrastName), ]
    }
    
    if (returnTable) {
      return(allContrastsdt)
    } else {
      return(allContrastsdt$contrast)
    }
  } else if (is(object, "SummarizedExperiment")) {
    # expects to find a diff analysis slot
    allContrasts <- metadata(object)$DiffExpAnal$contrasts
    
    if (returnTable) {
      return(allContrasts)
    } else {
      return(allContrasts$contrast)
    }
    
  } else {
    stop("object is not a MultiAssayExperiment or a SummarizedExperiment.")
  }
}

# ---- Get selected contrasts : ----
#' @title Get selected contrasts for the differential analysis
#'
#' @param object a MAE object (produced by Flomics) or a SE (expect to find a diffAnalysis slot.)
#' @return a dataTable
#' @rdname ContrastsSelection
#' @export
#'
getSelectedContrasts <- function(object) {
  # TODO check if it exists...
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@Contrasts.Sel
  } else
    if (is(object, "SummarizedExperiment")) {
    object@metadata$DiffExpAnal$contrasts
  } else {
    stop("object is not a MultiAssayExperiment or a SummarizedExperiment.")
  }
}

# ---- Set Valid Contrasts : (after differential analysis) ----
#' @title Set Valid Contrasts
#'
#' @param object a SE object or MAE object (produced by Flomics)
#' @param contrasts a vector of contrasts names
#' @return a Flomics SE or MAE
#' @rdname ContrastsSelection
#' @export
#'
setValidContrasts <- function(object,
                              contrasts) {
  
  if (isTagName(object, contrasts)) contrasts <- convertTagToContrast(object, contrasts)
  
  # TODO : check if there are DE entities for each contrasts before really validating them.
  
  if (is.character(contrasts)) {
    if (is(object, "SummarizedExperiment")) {
      allTab <- getPossibleContrasts(object, returnTable = TRUE)
      object@metadata$DiffExpAnal[["Validcontrasts"]] <- allTab[allTab$contrastName %in% contrasts,]
    } else {
      stop("object is not a SummarizedExperiment.")
    }
  }
  
  return(object)
}

#' @title Get Valid Contrasts
#'
#' @param object a SE object or MAE object (produced by Flomics)
#' @return a data frame or a list of data frames.
#' @rdname ContrastsSelection
#' @export
#'
getValidContrasts <- function(object) {
  if (is(object, "SummarizedExperiment")) {
    return(object@metadata$DiffExpAnal[["Validcontrasts"]])
  } else if (is(object, "MultiAssayExperiment")) {
    list_res <- lapply(names(object), FUN = function(tableName) {
      object[[tableName]]@metadata$DiffExpAnal[["Validcontrasts"]]
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
#' @rdname getDE
#'
getDEMatrix <- function(object) {
  if (is(object, "SummarizedExperiment")) {
    if (!is.null(object@metadata$DiffExpAnal$mergeDEF)) {
      object@metadata$DiffExpAnal$mergeDEF
    } else {
      stop("There is no DE matrix in this object.")
    }
  } else {
    stop("object is not a SummarizedExperiment.")
  }
}

#' @param contrast character name (can be a vector of name) for the contrast to select.
#' @param union Booleen value. TRUE : union; FALSE : intersection
#' @export
#' @importFrom tidyselect any_of
#' @importFrom dplyr select
#' @rdname getDE
getDE <- function(object, contrast, union = TRUE) {
  
  if (isContrastName(object, contrast)) contrast <- convertContrastToTag(object, contrast)
  
  DEmat <- getDEMatrix(object)
  DEmat <- DEmat %>% dplyr::select(tidyselect::any_of(c("DEF", contrast)))
  
  if(union) return(DEmat[rowSums(DEmat[,-1]) >= 1,])
  
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
#' @importFrom tidyselect starts_with
#'
opDEList <- function(object, SE.name = NULL, contrasts = NULL, operation = "union") {
  
  if (!is(object, "MultiAssayExperiment") && !is(object, "SummarizedExperiment")) 
    stop("Object is not a SummarizedExperiment or a MultiAssayExperiment")
  
  if (is(object, "MultiAssayExperiment") && is.null(SE.name)) 
    stop("Please provide SE.name argument.")
  
  if (is(object, "MultiAssayExperiment")) object <- object[[SE.name]] 
  
  if (is.null(contrasts) || length(contrasts) == 0) contrasts <- getSelectedContrasts(object)[["tag"]]
  if (isContrastName(object, contrasts)) contrasts <- convertContrastToTag(object, contrasts)
  
  if (!is.null(object@metadata$DiffExpAnal$Validcontrasts)) {
    validTags <- convertContrastToTag(object, getValidContrasts(object)$contrastName)
  } else {
    validTags <- contrasts
  }
  
  tagsConcerned <- intersect(contrasts, validTags)
  if (length(tagsConcerned) == 0) stop("It seems there is no contrasts to select DE entities from.")
  
  df_DE <- getDEMatrix(object) %>%
    dplyr::select(c("DEF", tidyselect::any_of(tagsConcerned)))
  
  if (operation == "intersection") {
    DETab <- df_DE %>%
      dplyr::mutate(SUMCOL = dplyr::select(., tidyselect::starts_with("H")) %>%
                      rowSums(na.rm = TRUE)) %>%
      dplyr::filter(SUMCOL == length(validTags))
 
  } else {
    DETab <- df_DE %>%
      dplyr::mutate(SUMCOL = dplyr::select(., tidyselect::starts_with("H")) %>%
                      rowSums(na.rm = TRUE)) %>%
      dplyr::filter(SUMCOL >= 1)
  }
  
  return(unique(DETab$DEF))
}


# ---- Get omics experiments and their types (vector) ----

#' @title Get omics experiments and their types
#'
#' @param object a MAE object (produced by Flomics) or a Summarized Experiment object.
#' @return a named vector with each omics name and its type.
#' @export
#'
getOmicsTypes <- function(object) {
  
  if (!is(object, "MultiAssayExperiment") && !is(object, "SummarizedExperiment"))
    stop("Object is not a SummarizedExperiment or a MultiAssayExperiment")
  
  if (is(object, "MultiAssayExperiment")) {
    return(names(object@metadata$omicList))
    
  } else {
    return(object@metadata$omicType)
  }
}


# ---- Get omics experiments and their types (vector) ----

#' @title Get dataset names
#'
#' @param object a MAE object (produced by RFlomics) or a Summarized Experiment object.
#' @return a named vector with each omics name and its type.
#' @export
#'
getDatasetNames <- function(object) {
  
  if (!is(object, "MultiAssayExperiment")) stop("Object is not a MultiAssayExperiment")
  
  datasetNames <- unlist(object@metadata$omicList)
  names(datasetNames) <- NULL
  
  return(datasetNames)
}


# ---- Get filtering parameters ----

#' @title Get filtering setting
#'
#' @param object a SE object (produced by Flomics).
#' @return List of filtering setting
#' @export
#'
getFilterSetting <- function(object) {
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  
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
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  
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
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  
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
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  
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
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  
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
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  
  return(object@metadata$CoExpAnal$setting)
}

# ---- Get integration setting ----

#' @title Get MOFA analysis setting parameters
#'
#' @param object a MultiAssayExperiment object (produced by Flomics).
#' @return List of differential analysis setting parameters
#' @export
#'
getMOFASetting <- function(object) {
  if (!is(object, "MultiAssayExperiment")) stop("Object is not a MultiAssayExperiment")
  
  return(object@metadata$IntegrationAnalysis$MOFA$setting)
}

#' @title Get mixOmics analysis setting parameters
#'
#' @param object a MultiAssayExperiment object (produced by Flomics).
#' @return List of differential analysis setting parameters
#' @export
#'
getMixOmicsSetting <- function(object) {
  if (!is(object, "MultiAssayExperiment")) stop("Object is not a MultiAssayExperiment")
  
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
  
  if (!is(object, "MultiAssayExperiment") && !is(object, "SummarizedExperiment")) 
    stop("Object is not a SummarizedExperiment or a MultiAssayExperiment")
  
  if (toupper(from) %in% c("DIFFEXP", "DIFFEXPANAL", "DIFFEXPENRICHANAL")) from <- "DiffExpEnrichAnal"
  if (toupper(from) %in% c("COEXP", "COEXPANAL", "COEXPENRICHANAL"))     from <- "CoExpEnrichAnal"
  
  classObj <- NULL
  if (is(object, "SummarizedExperiment")) classObj <- "SE"
  if (is(object, "MultiAssayExperiment")) classObj <- "MAE"
  
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
         stop("Object is not a MultiAssayExperiment nor a SummarizedExperiment")
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
  if (!is(object, "MultiAssayExperiment") && !is(object, "SummarizedExperiment")) 
    stop("Object is not a SummarizedExperiment or a MultiAssayExperiment")
  
  if (is(object, "SummarizedExperiment")) {
    return(object@metadata[[from]][[dom]][["summary"]])
  } else if (is(object, "MultiAssayExperiment")) {
    if (is.null(experiment)) {
      stop("Please indicate from which data you want to extract the enrichment results.")
    }
    
    return(object[[experiment]]@metadata[[from]][[dom]][["summary"]])
  }
}


# ---- INTERNAL - Get a pvalue threshold used in enrichment analysis ----
# TODO equivalent to sumORA (external)
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
  
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")
  if (!from %in% c("DiffExpEnrichAnal", "CoExpEnrichAnal")) stop(from, " don't existe")
  if (!dom  %in% c("GO", "KEGG", "custom")) stop(from, " not valide value. Choose from c(GO, KEGG, custom)")
  if (is.null(object@metadata[[from]][[dom]]$list_args$pvalueCutoff)) stop("P-value not found")
    
  return(object@metadata[[from]][[dom]]$list_args$pvalueCutoff)
}




