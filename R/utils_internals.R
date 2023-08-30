######## INTERNAL - ANNOTATION CLUSTERPROFILER #########

#' @title see_pathview
#' @param ... Possible arguments for pathview function.
#' @return nothing. Plot on the currently opened device.
#' @keywords internal
#' @noRd
#'
# Code from: https://stackoverflow.com/questions/60141841/how-to-get-pathview-plot-displayed-directly-rather-than-saving-as-a-file-in-r
# It deletes every file created by pathview
see_pathview <- function(...) {
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = TRUE)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  img <- png::readPNG(filename)
  grid::grid.raster(img)
  nam <- stringr::str_split(filename, "[.]")
  invisible(file.remove(filename))
  invisible(file.remove(paste0(nam[[1]][1], ".xml")))
  invisible(file.remove(paste0(nam[[1]][1], ".png")))
  return()
}

######## INTERNAL - CHECKS FUNCTIONS ###########

# check_NA: checks if there are NA/nan in the summarizedExperiment assay
#' @title check_NA
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @return boolean. if TRUE, NA/nan are detected in the SE::assay.
#' @keywords internal
#' @noRd
#'
check_NA <- function(object) {
  NA_detect <- ifelse(any(is.na(SummarizedExperiment::assay(object))), TRUE, FALSE)
  return(NA_detect)
}


######## INTERNAL - Transform the Data ###########

# apply_transformation: apply the transformation method stored in object@metadata[["transform_method"]] and modify the assay.
#' @title apply_transformation
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @keywords internal
#' @noRd
#'

apply_transformation <- function(object) {
  if (is.null(object@metadata[["transform"]][["transform_method"]])) {
    stop("Expect transformation method.")
  }

  if (object@metadata[["transform"]][["transformed"]]) {
    message("WARNING: data were already transformed before!")
  }

  transform_method <- object@metadata[["transform"]][["transform_method"]]
  assayTransform <- SummarizedExperiment::assay(object, withDimnames = TRUE)
  validTransform <- TRUE


  switch(transform_method,
    "log1p" = {
      SummarizedExperiment::assay(object) <- log1p(assayTransform)
    },
    "log2" = {
      SummarizedExperiment::assay(object) <- log2(assayTransform + 1)
    },
    "log10" = {
      SummarizedExperiment::assay(object) <- log10(assayTransform + 1)
    },
    "squareroot" = {
      SummarizedExperiment::assay(object) <- sqrt(assayTransform)
    },
    "none" = {
      SummarizedExperiment::assay(object) <- assayTransform
    },
    {
      SummarizedExperiment::assay(object) <- assayTransform
      message("Could not recognize the transformation method. No transformation applied. Please check your parameters.")
      validTransform <- FALSE
    } # default is none
  )

  if (transform_method != "none" && validTransform) object@metadata[["transform"]][["transformed"]] <- TRUE

  return(object)
}

######## INTERNAL - Normalize the Data ###########

# apply_norm: apply the normalization method stored in object@metadata[["Normalization"]] and modify the assay.
#' @title apply_norm
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @description apply the normalization to the assay. Usually, after the transformation,
#' unless in the case of counts RNASeq data (TMM), where log2 is the second step.
#' @keywords internal
#' @noRd
#'

apply_norm <- function(object) {
  if (is.null(object@metadata[["Normalization"]])) {
    stop("Expect normalization method.")
  }

  if (object@metadata[["Normalization"]]$normalized) {
    message("WARNING: data were already normalized before!")
  }

  norm_method <- object@metadata[["Normalization"]]$methode
  coefNorm <- object@metadata[["Normalization"]]$coefNorm
  validNorm <- TRUE

  assayTransform <- SummarizedExperiment::assay(object)

  switch(norm_method,
    "median" = {
      SummarizedExperiment::assay(object) <- sweep(assayTransform, 2, coefNorm, "-")
    },
    "totalSum" = {
      SummarizedExperiment::assay(object) <- sweep(assayTransform, 2, coefNorm, "/")
    },
    "TMM" = {
      scales_factors <- object@metadata[["Normalization"]]$coefNorm$norm.factors * object@metadata[["Normalization"]]$coefNorm$lib.size
      SummarizedExperiment::assay(object) <- scale(assayTransform + 1, center = FALSE, scale = scales_factors)
    },
    "none" = {
      SummarizedExperiment::assay(object) <- assayTransform
    },
    { # default is none
      SummarizedExperiment::assay(object) <- assayTransform
      message("Could not recognize the normalization method. No normalization applied. Please check your parameters.")
      validNorm <- FALSE
    }
  )

  if (norm_method != "none" && validNorm) object@metadata[["Normalization"]][["normalized"]] <- TRUE

  return(object)
}

######## INTERNAL - Check, transform and normalize the data ###########

# checkTransNorm: check the data, transform them and normalize them.
#' @title checkTransNorm
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @description apply the normalization and the transformation stored into the metadata of the SE object.
#'  Applies TMM and log2 transformation for RNAseq data.
#' @keywords internal
#' @noRd
#'

checkTransNorm <- function(object, raw = FALSE) {
  if (class(object) != "SummarizedExperiment") stop("Object is not a SummarizedExperiment")

  # check things
  if (check_NA(object)) stop("NA detected in the assay.")

  # No transformation (except for RNAseq, which expect counts...)
  if (raw) {
    if (object@metadata[["transform"]][["transformed"]]) message("WARNING: your data are not raw (transformed)")
    if (object@metadata[["Normalization"]]$normalized) message("WARNING: your data are not raw (normalized)")

    if (getOmicsTypes(object) == "RNAseq") {
      SummarizedExperiment::assay(object) <- log2(SummarizedExperiment::assay(object) + 1)
    }
  } else {
    # if RNAseq
    if (getOmicsTypes(object) == "RNAseq") {
      # Really depends if TMM is the normalization or not.
      # Make it easier: force TMM and log2.

      if (object@metadata[["transform"]][["transformed"]]) stop("Expect untransformed RNAseq data at this point.")

      if (object@metadata[["Normalization"]][["normalized"]] && object@metadata[["Normalization"]][["methode"]] == "TMM") {
        SummarizedExperiment::assay(object) <- log2(SummarizedExperiment::assay(object))
      } # +1 in the apply_norm function

      if (object@metadata[["Normalization"]][["normalized"]] && object@metadata[["Normalization"]][["methode"]] != "TMM") {
        message("RNAseq counts expects TMM normalization. Data were already normalized with another method.
                Skipping to the end without transforming or normalizing data.")
      }

      if (!object@metadata[["Normalization"]][["normalized"]]) {
        # Force "none" transformation.
        if (object@metadata[["transform"]][["transform_method"]] != "none") {
          message("RNAseq counts expects TMM normalization. Transformation is done after the normalization,
                  using 'none' as transform method. Data will be transformed using log2 after the normalization anyway")

          object@metadata[["transform"]][["transform_method"]] <- "none"
        }

        # Force TMM normalization
        if (object@metadata[["Normalization"]][["methode"]] != "TMM") {
          message("For RNAseq data (counts), only TMM applies for now. Forcing TMM normalization.")
          object <- RunNormalization(object, NormMethod = "TMM")
        }

        # Finally transforming the data.
        object <- apply_transformation(object) # none
        object <- apply_norm(object) # TMM
        SummarizedExperiment::assay(object) <- log2(SummarizedExperiment::assay(object)) # +1 in the apply_norm function
      }
    } else {
      # in case any other omics type (does not expect counts)
      # transform and norm
      if (!object@metadata[["transform"]][["transformed"]]) object <- apply_transformation(object)
      if (!object@metadata[["Normalization"]][["normalized"]]) object <- apply_norm(object)
    }
  }

  # check things
  if (check_NA(object)) stop("NA detected in the assay.")

  return(object)
}

######## INTERNAL - isNorm, isTransform, getNorm, getTransform ###########

#' @title isNorm, isTransform, getNorm, getTransform
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @description get if an assay has been transformed or normalized.
#' @keywords internal
#' @noRd
#'

isTransformed <- function(object) {
  object@metadata[["transform"]][["transformed"]]
}

isNorm <- function(object) {
  object@metadata[["Normalization"]][["normalized"]]
}

getNorm <- function(object) {
  object@metadata[["Normalization"]][["methode"]]
}

getTrans <- function(object) {
  object@metadata[["transform"]][["transform_method"]]
}
