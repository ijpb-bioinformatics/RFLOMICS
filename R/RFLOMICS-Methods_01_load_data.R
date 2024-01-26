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

