### ============================================================================
### RflomicsSE object
### ----------------------------------------------------------------------------

#' @import methods

#' @title \link{RflomicsSE-class} Class
#' @description
#' RflomicsSE is a class that extends the \link{SummarizedExperiment} by imposing a structure
#' on the metadata slot. This class is used by the Rflomics analysis workflow to store the 
#' experimental design, the settings and results of a single omic analysis. 
#' The slot metadata is structured as follows:
#' 
#' \itemize{
#'   \item omicType: the type of omics dataset
#'   \item design: experimental design
#'   \item DataProcessing: a list containing the data processing settings and results 
#'   \item PCAlist: a list containing the PCA settings and results
#'   \item DiffExpAnal: a list containing the Differential Analysis settings and results
#'   \item CoExpAnal: a list containing the Coexpression Analysis settings and results
#'   \item DiffExpEnrichAnal: a list containing the enrichment analysis of the list of DE features settings and results
#'   \item CoExpEnrichAnal: a list containing the enrichment analysis of the list of co-expressed features settings and results
#'  } 
#' 
#' @seealso \link{SummarizedExperiment}
#' @seealso \link{RflomicsSE-accessors}
#' @name RflomicsSE-class
#' @rdname RflomicsSE-class
#' @aliases RflomicsSE
#' @exportClass RflomicsSE
#' @return A \code{RflomicsSE} object.
setClass(
  Class    = "RflomicsSE",
  contains = "SummarizedExperiment",
)

#' @title initialize
#' @description
#' Constructor for the RflomicsSE class.
#' @return An object of class \link{RflomicsSE-class}
#' @seealso SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @noRd
setMethod(
  f          = "initialize",
  signature  = "RflomicsSE",
  definition = function(.Object, 
                        assays   = NULL,
                        colData  = S4Vectors::DataFrame(),
                        omicType = NULL,
                        design   = list() , 
                        DataProcessing = list() , 
                        PCAlist        = list() , 
                        DiffExpAnal    = list() ,      
                        CoExpAnal      = list() , 
                        DiffExpEnrichAnal = list() ,  
                        CoExpEnrichAnal   = list()) {
    
    .Object <-  callNextMethod(.Object, assays = assays, colData = colData)
    .Object@metadata <- list("omicType" = omicType , 
                             "design"   = design , 
                             "DataProcessing" = DataProcessing , 
                             "PCAlist"        = PCAlist , 
                             "DiffExpAnal"    = DiffExpAnal ,      
                             "CoExpAnal"      = CoExpAnal , 
                             "DiffExpEnrichAnal" = DiffExpEnrichAnal ,  
                             "CoExpEnrichAnal"   = CoExpEnrichAnal)
    
    return(.Object)
  }
)


### ============================================================================
### RflomicsMAE object
### ----------------------------------------------------------------------------

#' @title \link{RflomicsMAE-class} Class
#' @description
#' RflomicsMAE is a class that extends the \link{MultiAssayExperiment}  class by imposing a structure
#' to the metadata slot. This class is used by the Rflomics analysis workflow to store the 
#' experimental design, the settings and results of a multi-omics integration analysis.  
#' The slot metadata is structured as follows:
#'  
#'  \itemize{
#'    \item omicList: Contains the list of omics datasets, with the type and name  
#'    \item projectName: The project name
#'    \item design: The experimental design 
#'    \item IntegrationAnalysis: A list containing the multi-omics integration analysis settings and results.   
#'  }
#' 
#' @seealso \link{MultiAssayExperiment}
#' @seealso \link{RflomicsMAE-accessors}
#' @name RflomicsMAE-class
#' @rdname RflomicsMAE-class
#' @aliases RflomicsMAE-class
#' @return A \code{RflomicsMAE} object.
#' @exportClass RflomicsMAE
setClass(
  Class    = "RflomicsMAE",
  contains = "MultiAssayExperiment",
)


#' @title initialize
#' @description Constructor for the RflomicsMAE class.
#' @return An object of class \link{RflomicsMAE-class}
#' @seealso MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
#' @param ExperimentList Similar to MultiAssayExperiment, and experiment list. 
#' @param colData Similar to MultiAssayExperiment, a data.frame or a list 
#' containing all information about the samples
#' @param sampleMap Similar to MultiAssayExperiment
#' @param omicList list of omics data, each given by a dataframe-like object
#' Expect to find sample in rows and variables in columns.
#' @param projectName The name of the project. Is useful for report generation.
#' @param design The experimental design. 
#' @param IntegrationAnalysis usually empty at this point: list of results
#' for the integration. 
#' @noRd
setMethod(
  f          = "initialize",
  signature  = "RflomicsMAE",
  definition = function(.Object,
                        ExperimentList = MultiAssayExperiment::ExperimentList(), 
                        colData        = S4Vectors::DataFrame(), 
                        sampleMap      = S4Vectors::DataFrame(assay = factor(), primary = character(), colname = character()), 
                        omicList       = list(),
                        design         = list(),
                        IntegrationAnalysis = list(),
                        projectName    = NULL
  ) {
    
    .Object <-  callNextMethod(.Object, ExperimentList = ExperimentList, 
                               colData = colData, sampleMap = sampleMap)
    
    Sys.setlocale('LC_TIME', 'C') # change for english
    date <- format(Sys.time(), '%d %B %Y - %H:%M')
    Sys.setlocale('LC_TIME') # return to default system time
    
    # initialization of metadata
    .Object@metadata <- 
      list("omicList"            = omicList,
           "projectName"         = projectName,
           "design"              = design,
           "IntegrationAnalysis" = IntegrationAnalysis,
           "date"                = date,
           "sessionInfo"         = NULL,
           "rflomicsVersion"     = packageVersion('RFLOMICS'))
    
    return(.Object)
  })

