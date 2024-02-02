






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


