
#' @title Generate RFLOMICS rmarkdown report
#' @description
#' This function is used to generate a html report from a RFLOMICS
#' \link{MultiAssayExperiment}.
#' @param object a \link{MultiAssayExperiment} produced by RFLOMICS.
#' @param projectName name of the project. Title of the html document.
#' @param outdir where to save the results?
#' @param RDataName Name of the RData that will be saved and used by the
#' function to generate results.
#' @param output_file Name of the html document.
#' @importFrom rmarkdown render
#' @importFrom kableExtra kable_styling
#' @importFrom MultiAssayExperiment upsetSamples
#' @export
#'
generateReport <- function(object,
                           projectName = "project",
                           outDir = getwd(),
                           RDataName = "rflomics_MAE",
                           output_file = paste0(outDir, projectName, "_report.html"),
                           ...) {
  tempReport <-  paste0(path.package("RFLOMICS"), "/RFLOMICSapp/","report.Rmd")
  dir.create(path = outDir)
  save(object, file = file.path(outDir, RDataName))
  
  last <- nchar(outDir)
  if (substr(outDir, start = last - 1, stop = last) != "/"){
    outDir <- paste0(outDir, "/")
  }
  
  params <- list(FEdata = file.path(outDir, RDataName),
                 title  = paste0(projectName, " project"),
                 outDir = outDir)
  
  rmarkdown::render(tempReport,
                    output_file = paste0(outDir, projectName, "_report.html"),
                    params = params,
                    knit_root_dir = tempdir(),
                    intermediates_dir = tempdir(),
                    envir = new.env(parent = globalenv()), ...)
  
}