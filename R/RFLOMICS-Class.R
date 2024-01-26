### ============================================================================
### RflomicsSE object
### ----------------------------------------------------------------------------

#' @title \link{RflomicsSE-class} Class
#' @description
#' A short description...
#' 
#' @seealso \link{SummarizedExperiment}
#' @name RflomicsSE-class
#' @rdname RflomicsSE-class
#' @exportClass RflomicsSE
methods::setClass(
  Class    = "RflomicsSE",
  contains = "SummarizedExperiment",
)

#' @title initialize
#' @description
#' A short description...
#' 
#' @return An object of class \link{RflomicsSE-class}
#' @seealso SummarizedExperiment
#' @noRd
methods::setMethod(
  f          = "initialize",
  signature  = "RflomicsSE",
  definition = function(.Object, 
                        assays   = NULL,
                        colData  = S4Vectors::DataFrame(), 
                        metadata = list()) {
    
    .Object <- methods::callNextMethod(.Object, assays = assays, colData = colData)
    .Object@metadata <- metadata
    
    return(.Object)
  }
)


### ============================================================================
### RflomicsMAE object
### ----------------------------------------------------------------------------

#' @title \link{RflomicsMAE-class} Class
#' @description
#' A short description...
#' 
#' @seealso \link{MultiAssayExperiment}
#' @name RflomicsMAE-class
#' @rdname RflomicsMAE-class
#' @exportClass RflomicsMAE
methods::setClass(
  Class    = "RflomicsMAE",
  contains = "MultiAssayExperiment",
)


#' @title initialize
#' @description xxx
#' @return An object of class \link{RflomicsMAE-class}
#' @seealso MultiAssayExperiment
#' @noRd
methods::setMethod(
  f          = "initialize",
  signature  = "RflomicsMAE",
  definition = function(.Object,
                        ExperimentList = MultiAssayExperiment::ExperimentList(), 
                        colData        = S4Vectors::DataFrame(), 
                        sampleMap      = S4Vectors::DataFrame(assay = factor(), primary = character(), colname = character()), 
                        metadata       = list()){
    
    # 
    .Object <- methods::callNextMethod(.Object, ExperimentList = ExperimentList, colData = colData, sampleMap = sampleMap)
    # initialization of metadata
    .Object@metadata <- metadata
    
    # Retourner l'objet initialisé
    return(.Object)
  })

