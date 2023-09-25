# ---- Get Factor types : ----
#' @title Get design matrix used for a differential analysis
#'
#' @param object a MAE object (produced by Flomics)
#' @return a dataframe
#' @export
#' 
getFactorTypes <- function(object) {
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@Factors.Type
  } else {
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
getDesignMat <- function(object) {
  # TODO check if it exists...
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@ExpDesign
  } else {
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
getModelFormula <- function(object) {
  # TODO check if it exists...
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@Model.formula
  } else {
    stop("object is not a MultiAssayExperiment.")
  }
}

# ---- Get possible contrasts : ----
#' @title Get selected contrasts for the differential analysis
#'
#' @param object a MAE object (produced by Flomics)
#' @param typeContrast the type of contrast from which the possible contrasts are extracted. Default is all contrasts types.
#' @param modalities specific levels for the contrast selection
#' @param returnTable return a dataTable with all contrasts information
#' @return a character vector or a dataTable
#' @export
#'
getPossibleContrasts <- function(object, typeContrast = c("simple", "averaged", "interaction"),
                                 modalities = NULL, returnTable = FALSE) {
  if (is(object, "MultiAssayExperiment")) {
    if (is.null(typeContrast)) typeContrast <- c("simple", "averaged", "interaction")

    allContrasts <- MAE@metadata$design@Contrasts.List
    allContrasts <- allContrasts[which(names(allContrasts) %in% typeContrast)]
    # allContrasts$fill <- TRUE # do.call need this
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
  } else {
    stop("object is not a MultiAssayExperiment or a SummarizedExperiment.")
  }
}

# ---- Get selected contrasts : ----
#' @title Get selected contrasts for the differential analysis
#'
#' @param object a MAE object (produced by Flomics) or a SE (expect to find a diffAnalysis slot.)
#' @return a dataTable
#' @export
#'
getSelectedContrasts <- function(object) {
  # TODO check if it exists...
  if (is(object, "MultiAssayExperiment")) {
    object@metadata$design@Contrasts.Sel
  } else if (is(object, "SummarizedExperiment")) {
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
#' @export
#'
setValidContrasts <- function(object,
                              contrasts) {
  # TODO : check if there are DE entities for each contrasts before really validating them.

  if (is.character(contrasts)) {
    if (is(object, "SummarizedExperiment")) {
      object@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName <- contrasts
    } else {
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
getValidContrasts <- function(object) {
  if (is(object, "SummarizedExperiment")) {
    return(object@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName)
  } else if (is(object, "MultiAssayExperiment")) {
    list_res <- lapply(names(object), FUN = function(tableName) {
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
getDEMatrix <- function(object) {
  if (is(object, "SummarizedExperiment")) {
    if (!is.null(object@metadata$DiffExpAnal$mergeDEF)) {
      return(object@metadata$DiffExpAnal$mergeDEF)
    } else {
      stop("There is no DE matrix in this object.")
    }
  } else {
    stop("object is not a SummarizedExperiment.")
  }
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
#' @export
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
    validTags <- convertContrastToTag(object, getValidContrasts(object))
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

  return(DETab$DEF)
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
    sapply(names(object), FUN = function(x) {
      object[[x]]@metadata$omicType
    })
  } else {
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
getNormCoeff <- function(object) {
  if (!is(object, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment")

  return(object@metadata$Normalization$coefNorm)
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

  if (toupper(from) %in% c("DIFFEXPANAL", "DIFFEXPENRICHANAL")) from <- "DiffExpEnrichAnal"
  if (toupper(from) %in% c("COEXPANAL", "COEXPENRICHANAL"))     from <- "CoExpEnrichAnal"

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


# ----- INTERNAL - Check if character vectors are contrasts Names : -----

#' @title Check if character vectors are contrasts Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expect to find
#'  a slot of differential analysis.
#' @param contrastName vector of characters.
#' @return boolean. TRUE if all of contrastName are indeed contrasts Names.
#' @noRd
#' @keywords internal
isContrastName <- function(object, contrastName) {
  df_contrasts <- getSelectedContrasts(object)
  
  search_match <- sapply(contrastName, FUN = function(cn) {
    grep(cn, df_contrasts$contrastName, fixed = TRUE)
  })
  search_success <- sapply(search_match, identical, integer(0)) # if TRUE, not a success at all.
  
  if (!any(search_success)) {
    # Congratulations, it's a contrast name!
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# ----- INTERNAL - Check if character vectors are tags Names : -----

#' @title Check if character vectors are tags Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expect to find
#'  a slot of differential analysis.
#' @param tagName vector of characters.
#' @return boolean. TRUE if all of tagName are indeed tags Names.
#' @noRd
#' @keywords internal
isTagName <- function(object, tagName) {
  df_contrasts <- getSelectedContrasts(object)
  
  search_match <- sapply(tagName, FUN = function(cn) {
    grep(cn, df_contrasts$tag, fixed = TRUE)
  })
  search_success <- sapply(search_match, identical, integer(0)) # if TRUE, not a success at all.
  
  if (!any(search_success)) {
    # Congratulations, it's a tag name!
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# ---- INTERNAL - convert tag to contrastName ----

#' @title Convert tags names to contrast Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expects to find
#'  a slot of differential analysis.
#' @param tagName Vector of characters, expect to be tags (in the form of H1, H2, etc.).
#' @return character vector, contrastNames associated to tags.
#' @noRd
#' @keywords internal
convertTagToContrast <- function(object, tagName) {
  df_contrasts <- getSelectedContrasts(object)
  
  df_contrasts %>%
    dplyr::filter(tag %in% tagName) %>%
    dplyr::select(contrastName) %>%
    unlist(use.names = FALSE)
}

# ---- INTERNAL - convert contrastName to tag ----

#' @title Convert contrast Names names to tags
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a summarizedExperiment, expects to find
#'  a slot of differential analysis.
#' @param contrasts Vector of characters, expect to be contrast names.
#' @return character vector, tags associated to contrast names.
#' @noRd
#' @keywords internal
convertContrastToTag <- function(object, contrasts) {
  df_contrasts <- getSelectedContrasts(object)
  
  df_contrasts %>%
    dplyr::filter(contrastName %in% contrasts) %>%
    dplyr::select(tag) %>%
    unlist(use.names = FALSE)
}
