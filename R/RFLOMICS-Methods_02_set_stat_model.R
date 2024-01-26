
# ---- Set Model Formula : ----
#' @title Set model formula from a Flomics RflomicsMAE
#'
#' @param object a MAE or SE object (produced by Flomics)
#' @param a formula
#' @return a object SE
#' @export
#'
setModelFormula <- function(object, modelFormula=NULL) {
  # TODO check if it exists...
  if (is(object, "RflomicsMAE")) {
    object@metadata$design$Model.formula <- paste(modelFormula, collapse = " ")
    # set modelFormula foreach SE
    for(se in names(object)){
      object[[se]]@metadata$design$Model.formula <- paste(modelFormula, collapse = " ")
    }
  } 
  else if(is(object, "RflomicsSE")){
    object@metadata$design$Model.formula <- paste(modelFormula, collapse = " ")
  }
  else {
    stop("object is not a RflomicsMAE")
  }
  return(object)
}


# ---- INTERNAL - get origin of a particular name ----
#
#' @title get origin of a name given a rflomics MAE
#'
#' @param object a RflomicsSE, produced by rflomics
#' @param name name of the parameter to identify. For clusters, please
#' specify cluster.1, cluster.2, etc.
#' @return The origin of the name, one of Contrast, Tag or CoexCluster.
#' @noRd
#' @keywords internal

.getOrigin <- function(object, name) {
  
  if (isContrastName(object, name)) return("Contrast")
  if (isTagName(object, name)) return("Tag")
  if (isClusterName(object, name)) return("CoexCluster")
  
  return("NoOriginFound")
}

