
#' @title Generate RFLOMICS rmarkdown report
#' @description
#' This function is used to generate a html report from a RFLOMICS
#' \link{MultiAssayExperiment}.
#' @param rflomics.MAE a \link{MultiAssayExperiment} produced by RFLOMICS.
#' @param projectName name of the project. Title of the html document.
#' @param outdir where to save the results?
#' @param RDataName Name of the RData that will be saved and used by the
#' function to generate results.
#' @param output_file Name of the html document.
#' @importFrom rmarkdown render
#' @importFrom kableExtra kable_styling
#' @export
#'
generateReport <- function(rflomics.MAE,
                           projectName = "project",
                           outDir = getwd(),
                           RDataName = "rflomics_MAE",
                           output_file = paste0(outDir, projectName, "_report.html"),
                           ...) {
  # Copy the report file to a temporary directory before processing it, in
  # case we don't have write permissions to the current working dir (which
  # can happen when deployed).
  tempReport <-  paste0(path.package("RFLOMICS"), "/RFLOMICSapp/","report.Rmd")
  dir.create(path = outDir)
  
  # save FE rflomics.MAE in .Rdata and load it during report execution
  save(rflomics.MAE, file = file.path(outDir, RDataName))
  
  last <- nchar(outDir)
  if (substr(outDir, start = last - 1, stop = last) != "/") {
    outDir <- paste0(outDir, "/")
  }
  
  # Set up parameters to pass to Rmd document
  params <- list(FEdata = file.path(outDir, RDataName),
                 title  = paste0(projectName, " project"),
                 outDir = outDir)
  
  # Knit the document, passing in the `params` list, and eval it in a
  # child of the global environment (this isolates the code in the document
  # from the code in this app).
  rmarkdown::render(tempReport,
                    output_file = paste0(outDir, projectName, "_report.html"),
                    params = params,
                    knit_root_dir = tempdir(),
                    intermediates_dir = tempdir(),
                    envir = new.env(parent = globalenv()), ...)
  
}