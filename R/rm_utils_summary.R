

# ---- Get summary from ORA : ----

#' @title Get summary tables from ORA analyses - once an enrichment has been conducted.
#'
#' @param object a SE object (produced by Flomics)
#' @param ont either NULL, GO, KEGG or custom. if NULL, all tables are returned in a list.
#' @param from either DiffExpEnrichAnal or CoExpAnal.
#' @return a list of tables or a table
#' @export
#'
sumORA <- function(SE, from = "DiffExpEnrichAnal", ont = NULL, contrast = NULL) {
  if (!class(SE) %in% "RflomicsSE") stop("SE is not a RflomicsSE.")

  if (toupper(from) %in% c("DIFFEXPANAL", "DIFFEXPENRICHANAL")) from <- "DiffExpEnrichAnal"
  if (toupper(from) %in% c("COEXPANAL", "COEXPENRICHANAL"))     from <- "CoExpEnrichAnal"
  
  # cat("|From: ", from, "\n")
  
  if (!is.null(contrast)) {
    if (isTagName(SE, contrast)) contrast <- convertTagToContrast(SE, contrast)
  }

  if (!is.null(ont)) {
    toReturn <- SE@metadata[[from]][[ont]]$summary
    if (!is.null(contrast)) toReturn <- toReturn[which(toReturn$Contrast == contrast), ]
    return(toReturn)
  } else {
    list_res <- lapply(names(SE@metadata[[from]]), FUN = function(ontres) {
      interRes <- SE@metadata[[from]][[ontres]]$summary
      if (!is.null(contrast)) interRes <- interRes[which(interRes$Contrast == contrast), ]
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
sumMixOmics <- function(object, selectedResponse = NULL) {
  if (is(object, "RflomicsMAE")) stop("Object is not a RflomicsMAE.")
  if (is.null(object@metadata$mixOmics)) stop("It seems this object has no mixOmics results.")

  if (is.null(selectedResponse)) {
    res <- lapply(names(object@metadata$mixOmics), FUN = function(selResponse) {
      getOneMORes(object, selectedResponse = selResponse)
    })
    names(res) <- names(object@metadata$mixOmics)
    return(res)
  } else {
    getOneMORes(object, selectedResponse = selectedResponse)
  }
}

#' @title Get an overview of MixOmics integration results for a specific response variable.
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string.
#' @return A data frame.
#' @keywords internal
#' @noRd

getOneMORes <- function(object, selectedResponse) {
  Data_res <- object@metadata$mixOmics[[selectedResponse]]$MixOmics_results

  df <- t(sapply(Data_res$X, dim))
  colnames(df) <- c("Ind", "Features")

  if (!is.null(object@metadata$mixOmics[[selectedResponse]]$MixOmics_results$MixOmics_tuning_results)) {
    df <- cbind(df, do.call("rbind", Data_res$keepX))
    colnames(df)[!colnames(df) %in% c("Ind", "Features")] <- paste("Comp", 1:length(Data_res$keepX[[1]]))
  }

  return(df)
}
