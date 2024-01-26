### ============================================================================
### Load data Rflomics methods
### ----------------------------------------------------------------------------

### All methods and accessors of MAE or SE are inherited by RflomicsMAE and 
### RflomicsSE, respectively.
### MAE : colData, assays, ..., upset, ...
### SE  : colData, assays, ...

# ---- is RflomicsMAE / RflomicsSE ----
# ---- extraire un RflomicsSE d'un RflomicsMAE ----
# ---- reduire un MAE à une list de experiments : elle existe déjà ----
# ---- ... ----

# ---- Get experiment names from RflomicsMAE object ----

#' @title Get Experiment Names
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with experiment names
#' @exportMethod getExperimentNames
methods::setMethod(f          = "getExperimentNames",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     ExperimentNames <- unlist(object@metadata$omicList)
                     names(ExperimentNames) <- NULL
                     
                     return(ExperimentNames)
                   })

# ---- Get omics types of experiments from RflomicsMAE object ----

#' @title Get omics types of experiments
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a data.frame with 
#' @exportMethod getExperimentTypes
methods::setMethod(f          = "getExperimentTypes",
                   signature  = "RflomicsMAE",
                   definition <- function(object, experimentNames = NULL){
                     
                     all_experiment_names <- getExperimentNames(object)
                     
                     if(is.null(all_experiment_names)) stop("RflomicsMAE object is empty.")
                     if(is.null(experimentNames)) experimentNames <- all_experiment_names

                     if(!any(experimentNames %in% all_experiment_names))
                       stop("No experimentNames matching...")
                     
                     # 
                     df <- lapply(names(object@metadata$omicList), function(x){ 
                       
                       data.frame(names = object@metadata$omicList[[x]], 
                                  types = rep(x, length(object@metadata$omicList[[x]]))) 
                     }) |> purrr::reduce(rbind)
                     
                     
                     # Check if all names from experimentNames exist in the object.
                     if(!all(experimentNames %in% all_experiment_names)) 
                       message("!!!Some experimentNames are not present in the object.")
                     
                     return(df[which(df$names %in% experimentNames),])
                     
                   })





# ---- Get Factor types : ----
#' @title Get factor types or names of particular factors given a type.
#'
#' @param object a MAE object (produced by Flomics). 
#' @return a named vector (getFactorTypes) or a vector of factor names.
#' @rdname getFactorTypes
#' @export
#' 
getFactorTypes <- function(object) {
  
  if (is(object, "RflomicsMAE")) {
    metadata(object)$design$Factors.Type
  } 
  else if(is(object, "RflomicsSE")){
    metadata(object)$design$factorType
  }
  else {
    stop("object is not a RflomicsMAE or RflomicsSE.")
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
  if (is(object, "RflomicsMAE") || is(object, "RflomicsSE")) {
    return(as.data.frame(object@colData))
    
  } else {
    stop("object is not a RflomicsMAE nor a RflomicsSE.")
  }
}


# ---- Get Model Formula : ----
#' @title Get model formula from a Flomics RflomicsMAE
#'
#' @param object a MAE object (produced by Flomics)
#' @return a formula
#' @export
#'
getModelFormula <- function(object) {
  # TODO check if it exists...
  if (is(object, "RflomicsMAE")) {
    object@metadata$design$Model.formula
  } 
  else if(is(object, "RflomicsSE")){
    object@metadata$design$Model.formula
  }
  else {
    stop("object is not a RflomicsMAE")
  }
}



# ---- Get omics experiments and their types (vector) ----

#' @title Get omics experiments and their types
#'
#' @param object a MAE object (produced by Flomics) or a Summarized Experiment object.
#' @return a named vector with each omics name and its type.
#' @export
#'
getOmicsTypes <- function(object) {
  
  if (!is(object, "RflomicsMAE") && !is(object, "RflomicsSE"))
    stop("Object is not a RflomicsSE or a RflomicsMAE")
  
  if (is(object, "RflomicsMAE")) {
    return(names(object@metadata$omicList))
    
  } else {
    return(object@metadata$omicType)
  }
}


