######## INTERNAL - ANNOTATION CLUSTERPROFILER #########

#' @title see_pathview
#' @param ... Possible arguments for pathview function.
#' @return nothing. Plot on the currently opened device.
#' @keywords internal
#' @importFrom grid grid.raster
#' @importFrom stringr str_split
#' @importFrom pathview pathview
#' @importFrom png readPNG
#' @noRd
#'
# Code from: https://stackoverflow.com/questions/60141841/how-to-get-pathview-plot-displayed-directly-rather-than-saving-as-a-file-in-r
# It deletes every file created by pathview
see_pathview <- function(...) {
  if (!exists("bods")) {
    data(bods, package = "pathview")
  }
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
  NA_detect <- ifelse(any(is.na(assay(object))), TRUE, FALSE)
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
  if (is.null(getTrans(object))) {
    stop("Expect transformation method.")
  }
  
  if (isTransformed(object)) {
    warning("Data were already transformed before!")
  }
  
  transform_method <- getTrans(object)
  assayTransform <- assay(object, withDimnames = TRUE)
  validTransform <- TRUE
  
  
  switch(transform_method,
         "log1p" = {
           assay(object) <- log1p(assayTransform)
         },
         "log2" = {
           assay(object) <- log2(assayTransform + 1)
         },
         "log10" = {
           assay(object) <- log10(assayTransform + 1)
         },
         "squareroot" = {
           assay(object) <- sqrt(assayTransform)
         },
         "none" = {
           assay(object) <- assayTransform
         },
         {
           assay(object) <- assayTransform
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
  if (is.null(getNorm(object))) {
    stop("Expects normalization method.")
  }
  
  if (isNorm(object)) {
    warning("Data were already normalized before!")
  }
  
  norm_method <- getNorm(object)
  coefNorm    <- getCoeffNorm(object)
  validNorm   <- TRUE
  
  assayTransform <- assay(object)
  
  switch(norm_method,
         "median" = {
           assay(object) <- sweep(assayTransform, 2, coefNorm, "-")
         },
         "totalSum" = {
           assay(object) <- sweep(assayTransform, 2, coefNorm, "/")
         },
         "TMM" = {
           scales_factors <- coefNorm$norm.factors * coefNorm$lib.size
           assay(object) <- scale(assayTransform + 1, center = FALSE, scale = scales_factors)
         },
         "none" = {
           assay(object) <- assayTransform
         },
         { # default is none
           assay(object) <- assayTransform
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
    if (isTransformed(object)) warning("Your data are not raw (transformed)")
    if (isNorm(object)) warning("Your data are not raw (normalized)")
    
    if (getOmicsTypes(object) == "RNAseq") {
      assay(object) <- log2(assay(object) + 1)
    }
  } else {
    # if RNAseq
    if (getOmicsTypes(object) == "RNAseq") {
      # Really depends if TMM is the normalization or not.
      # Make it easier: force TMM and log2.
      
      if (isTransformed(object)) stop("Expect untransformed RNAseq data at this point.")
      
      if (isNorm(object) && getNorm(object) == "TMM") {
        assay(object) <- log2(assay(object))
      } # +1 in the apply_norm function
      
      if (isNorm(object) && getNorm(object) != "TMM") {
        message("RNAseq counts expects TMM normalization. Data were already normalized with another method.
                Skipping to the end without transforming or normalizing data.")
      }
      
      if (!isNorm(object)) {
        # Force "none" transformation.
        if (getTrans(object) != "none") {
          message("RNAseq counts expects TMM normalization. Transformation is done after the normalization,
                  using 'none' as transform method. Data will be transformed using log2 after the normalization anyway")
          
          object <- setTrans(object, methode = "none")
        }
        
        # Force TMM normalization
        if (object@metadata[["Normalization"]][["methode"]] != "TMM") {
          message("For RNAseq data (counts), only TMM applies for now. Forcing TMM normalization.")
          object <- RunNormalization(object, NormMethod = "TMM")
        }
        
        # Finally transforming the data.
        object <- apply_transformation(object) # none
        object <- apply_norm(object) # TMM
        assay(object) <- log2(assay(object)) # +1 in the apply_norm function
      }
    } else {
      # in case any other omics type (does not expect counts)
      # transform and norm
      if (!isTransformed(object)) object <- apply_transformation(object)
      if (!isNorm(object)) object <- apply_norm(object)
    }
  }
  
  # check things
  if (check_NA(object)) stop("NA detected in the assay.")
  
  return(object)
}

######## INTERNAL - isNorm, isTransform, getNorm, getTransform ###########

#' @title isNorm, isTransform, getNorm, getTransform, setNorm, setTrans
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @description get if an assay has been transformed or normalized.
#' @keywords internal
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @noRd
#'

isTransformed <- function(object) {
  metadata(object)[["transform"]][["transformed"]]
}

isNorm <- function(object) {
  metadata(object)[["Normalization"]][["normalized"]]
}

getNorm <- function(object) {
  metadata(object)[["Normalization"]][["methode"]]
}

getCoeffNorm <- function(object) {
  metadata(object)[["Normalization"]][["coefNorm"]]
}

getTrans <- function(object) {
  metadata(object)[["transform"]][["transform_method"]]
}

setTrans <- function(object, methode = "none") {
  metadata(object)[["transform"]][["transform_method"]] <- methode
  return(object)
}

setNorm <- function(object, methode = "none") {
  metadata(object)[["Normalization"]][["methode"]] <- methode
  return(object)
}

setCoeffNorm <- function(object, coeff = NULL) {
  metadata(object)[["Normalization"]][["coefNorm"]] <- coeff
  return(object)
}


# ---- DO NOT PLOT function ----

#' @title doNotPlot
#' @description
#' Used mainly for the interface to check some conditions before actually plotting said graph.
#' 
#' @param expr An expression, usually producing a plot but not necessarily.
#' @keywords internal
#' @noRd
#' @importFrom utils capture.output
#' 
.doNotPlot <- function(expr){
  pdf(file = NULL)
  out <- tryCatch({
    capture.output(
        suppressMessages(
          eval(expr))
      )
  },
  error = function(e) e,
  warning = function(w) w
  )
  dev.off()
  return(out)
}


#' @title doNotSpeak
#' @description
#' Used mainly for the interface to silence some functions.
#' 
#' @param expr An expression, usually producing a warning.
#' @keywords internal
#' @noRd
#' 
.doNotSpeak <- function(expr){
  out <- tryCatch({
    capture.output(
      suppressWarnings(
        suppressMessages(eval(expr))
      ))
  },
  error = function(e) e,
  warning = function(w) w
  )
  return(out)
}