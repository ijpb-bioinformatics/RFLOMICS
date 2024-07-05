### ============================================================================
### [00_common] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# N. Bessoltane,
# D. Charif,
# A. Hulot

#' @import methods

# ---- resetRflomicsMAE ----
#' @title resetRflomicsMAE
#' @description
#' resetRflomicsMAE allows for initializing the object or initializing a 
#' selection of results.
#' @param object An object of class \link{RflomicsMAE-class}
#' @param singleAnalyses vector of single omics analysis results names 
#' (c("DataProcessing", "PCAlist", "DiffExpAnal", 
#' "DiffExpEnrichAnal", "CoExpAnal", "CoExpEnrichAnal"))
#' @param multiAnalyses vector of multi omics analysis results names 
#' (c("IntegrationAnalysis"))
#' @param datasetNames dataset name. 
#' If dataset == NULL, all datasets will be reset
#' @return An object of class \link{RflomicsMAE-class}
#' @noRd
#' @keywords internal
setMethod(f = "resetRflomicsMAE", 
          signature = "RflomicsMAE",
          definition = function(object,
                                singleAnalyses = NULL, 
                                multiAnalyses  = NULL,
                                datasetNames   = NULL) {
            
            all.datasets <- getDatasetNames(object)
            
            if(!is.null(datasetNames)) {
              datasetNames <- intersect(datasetNames, all.datasets)
              if(length(datasetNames) == 0) stop("")
            }
            
            if(is.null(datasetNames) && is.null(singleAnalyses) && is.null(multiAnalyses)) stop("")
            
            # if analyses is null and datasetNames not null -> remove 
            if(is.null(singleAnalyses)){
              
              if(!is.null(datasetNames)){
                
                index <- names(object) %in% datasetNames
                object <- object[,, !index]
                
              }
            }
            else{
              
              # if dataset is null we take all datasets 
              # present in RflomicsMAE object
              if (is.null(datasetNames)) {
                datasetNames <- getDatasetNames(object)
              }
              
              for (data in datasetNames) {
                
                for (analysis in singleAnalyses) {
                  
                  if(!is.null(object[[data]]@metadata[[analysis]])){
                    object[[data]]@metadata[[analysis]] <- list()
                  }
                }
              }
            }
            
            if(!is.null(multiAnalyses)){
              
              for (analysis in multiAnalyses) {
                
                if(!is.null(object@metadata[[analysis]])){
                  object@metadata[[analysis]] <- list()
                }
              }
            }
            
            return(object)
          })

# ---- getAnalyzedDatasetNames ----
#' @rdname generateReport
#' @description
#' \itemize{
#'    \item getAnalyzedDatasetNames: return a list of performed analysis names.}
#' @param analyses vector of list of analysis name
#' @exportMethod getAnalyzedDatasetNames
#' @aliases getAnalyzedDatasetNames,RflomicsMAE-method
#' @name getAnalyzedDatasetNames
#' @examples
#' # See generateReport for an example that includes getAnalyzedDatasetNames
setMethod(
  f          = "getAnalyzedDatasetNames",
  signature  = "RflomicsMAE",
  definition = function(object, analyses = NULL) {
    
    all.analyses <- c("DataProcessing",
                      "DiffExpAnal", "DiffExpEnrichAnal", 
                      "CoExpAnal", "CoExpEnrichAnal")
    
    if(is.null(analyses)) analyses <- all.analyses
    
    df.list <- list()
    for (dataset in getDatasetNames(object)) {
      
      if(is.null(object[[dataset]])) next
      
      for(analysis in analyses){
        
        if(length(object[[dataset]]@metadata[[analysis]]) == 0)
          next
        
        switch (analysis,
                "DataProcessing" = {
                  if(isTRUE(object[[dataset]]@metadata[[analysis]]$done))
                    df.list[[analysis]] <- c(df.list[[analysis]], dataset)
                },
                "DiffExpAnal" = {
                  if(!is.null(object[[dataset]]@metadata[[analysis]]$Validcontrasts))
                    df.list[[analysis]] <- c(df.list[[analysis]], dataset)
                },
                "CoExpAnal" = {
                  if(isTRUE(object[[dataset]]@metadata[[analysis]]$results))
                    df.list[[analysis]] <- c(df.list[[analysis]], dataset)
                },
                {
                  for(db in names(object[[dataset]]@metadata[[analysis]])){
                    df.list[[analysis]][[db]] <- c(df.list[[analysis]][[db]], dataset)
                  }
                }
        )
      }
    }
    
    if(length(df.list) == 0) return(NULL)
    if(length(analyses) == 1) return(df.list[[1]])
    return(df.list)
  })

## ---- set element to metadata slot in rflomicsSE/MAE ----
#' @title setElementToMetadata
#' @description set element to metadata slot
#' @param object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results.. 
#' @param name the name of element to add to metadata slot.
#' @param subName the name of sub element to add to metadata slot.
#' @param content the content of element to add
#' @return An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE-class}.
#' @keywords internal
#' @noRd
setMethod(
  f = "setElementToMetadata",
  signature = "RflomicsMAE",
  definition = function(object, 
                        name = NULL,
                        subName = NULL,
                        content = NULL) {
    
    object <- 
      .setElementToMetadata(object, name, subName, content)
    
    return(object)
  })

#' @keywords internal
#' @noRd
setMethod(
  f = "setElementToMetadata",
  signature = "RflomicsSE",
  definition = function(object, 
                        name = NULL,
                        subName = NULL,
                        content = NULL) {
    
    object <- 
      .setElementToMetadata(object, name, subName, content)
    
    return(object)
  }
)

## ---- get element from metadata slot from rflomicsSE/MAE ----


#' @aliases getAnalysis,RflomicsMAE-method
#' @name getAnalysis
#' @rdname generateReport
#' @description
#' \itemize{
#'    \item getAnalysis: return list of results from a specific analysis.}
#' @param SE.name name of the experiment where the metadata should be added.
#' @param name the name of element to add to metadata slot.
#' @param subName the name of sub element to add to metadata slot.
#' @exportMethod getAnalysis
setMethod(
  f = "getAnalysis",
  signature = "RflomicsMAE",
  definition = function(object,
                        SE.name = NULL,
                        name    = NULL,
                        subName = NULL){
    
    if(!is.null(SE.name)){
      if(!SE.name %in% getDatasetNames(object))
        stop(SE.name, " ?")

      object <- object[[SE.name]]
    }

    results <- 
      .getAnalysis(object, name, subName)
    
    return(results)
  })

#' @rdname generateReport
#' @aliases getAnalysis,RflomicsSE-method
#' @name getAnalysis
#' @exportMethod getAnalysis
setMethod(
  f = "getAnalysis",
  signature = "RflomicsSE",
  definition = function(object, 
                        name = NULL,
                        subName = NULL) {
    
    results <- 
      .getAnalysis(object, name, subName)
    
    return(results)
  }
)


# ---- generateReport ----
#' @title Generate RFLOMICS rmarkdown report
#' @description
#' This function is used to generate a html report from a
#' \link{RflomicsMAE-class} object or archive with results.
#' @param object a object of \link{RflomicsSE} class or 
#' \link{RflomicsMAE-class} class.
#' @param fileName Name of the html report (default: date()_projectName.html).
#' @param archiveName name of archive with all analysis results 
#' (default: date()_projectName.tar.gz).
#' @param export boolean value to create archive (default: FALSE)
#' @param tmpDir temporary directory (default: working directory)
#' @param ... other arguments to pass into the render function.
#' @return An html report or archive (tar.gz)
#' @importFrom rmarkdown render
#' @importFrom ggpubr as_ggplot 
#' @importFrom ggplot2 rel
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr relocate
#' @importFrom DT datatable formatStyle formatSignif
#' @importFrom stats cor
#' @exportMethod generateReport
#' @rdname generateReport
#' @name generateReport
#' @aliases generateReport,RflomicsMAE-method
#' @example inst/examples/generateReport.R
setMethod(
  f          = "generateReport",
  signature  = "RflomicsMAE",
  definition = function(object,
                        fileName = NULL,
                        archiveName = NULL,
                        export = FALSE,
                        tmpDir = getwd(),
                        ...) {
    # Copy the report file to a temporary directory before processing it, in
    # case we don't have write permissions to the current working dir (which
    # can happen when deployed).
    tempReport <-
      file.path(path.package("RFLOMICS"), "/RFLOMICSapp/report.Rmd")
    
    # Check if the object is properly filled.
    ## function
    
    # project name
    projectName <- getProjectName(object)
    RDataName   <- paste0(projectName, ".MAE.RData")
    
    # tmp dir
    if (file.access(tmpDir, 2) != 0)
      stop("No writing access in ", tmpDir)
    
    tmpDir <-
      file.path(tmpDir, 
                paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName))
    # file.path(tmpDir, 
    #           paste0(projectName, "_report"))
    
    
    dir.create(tmpDir, showWarnings = FALSE)
    
    # html name
    if (is.null(fileName))
      fileName <- file.path(tmpDir,
                            paste0(format(Sys.time(), "%Y_%m_%d"), "_", 
                                   projectName, ".html"))
    
    # save FE rflomics.MAE in .Rdata and load it during report execution
    sessionInfo.light <- .writeSessionInfo()
    rflomics.MAE <- setElementToMetadata(object, 
                                         name = "sessionInfo", 
                                         content = sessionInfo.light)
    save(rflomics.MAE, file = file.path(tmpDir, RDataName))
    
    # Set up parameters to pass to Rmd document
    param.list <-
      list(
        FEdata = file.path(tmpDir, RDataName),
        title  = paste0(projectName, " project"),
        rflomicsVersion = object@metadata$rflomicsVersion,
        date = object@metadata$date,
        outDir = tmpDir
      )
    
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the 
    # document from the code in this app).
    render(
      input             = tempReport,
      output_file       = fileName,
      params            = param.list,
      knit_root_dir     = tmpDir,
      intermediates_dir = tmpDir,
      envir = new.env(parent = globalenv()),
      ...
    )
    
    #Export results
    if (isTRUE(export)) {
      if (is.null(archiveName))
        archiveName <- file.path(dirname(tmpDir),
                                 paste0(format(Sys.time(),"%Y_%m_%d"),
                                        "_", projectName, 
                                        ".tar.gz"))
      
      
      
      # cp html in tmpDir
      file.copy(from = fileName, to = tmpDir)
      cmd <- paste0("tar -C ", dirname(tmpDir),
                    " -czf ", archiveName,  " ", basename(tmpDir))
      system(cmd)
      #message(cmd)
      
    } else{
      file.copy(from = fileName, to = dirname(tmpDir))
    }
    unlink(tmpDir, recursive = TRUE)
    
  }
)
