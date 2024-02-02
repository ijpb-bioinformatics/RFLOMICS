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
# Code from: https://stackoverflow.com/questions/60141841/
# how-to-get-pathview-plot-displayed-directly-rather-than-saving-as-a-file-in-r
# It deletes every file created by pathview
.see_pathview <- function(...) {
  if (!exists("bods")) {
    data(bods, package = "pathview")
  }
  msg <- capture.output(pathview(...), type = "message")
  msg <- grep("image file", msg, value = TRUE)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  if (length(filename) > 0 ) {
    img <- readPNG(filename)
    grid.raster(img)
    nam <- str_split(filename, "[.]")
    invisible(file.remove(filename))
    invisible(file.remove(paste0(nam[[1]][1], ".xml")))
    invisible(file.remove(paste0(nam[[1]][1], ".png")))
  }
  return()
}

# ----- Enrichment results ----
#
#' @title Get a particular enrichment result
#' @description
#' Called inside the getEnrichRes method
#' 
#' @param object a SE object or a MAE object (produced by Flomics).
#' @return enrichment result.
#' @noRd
#' @keywords internal
#' 

.getEnrichResIntSE <- function(
    object, 
    contrastName, 
    from ,
    database,  
    domain
){
  if (toupper(from) %in% c("DIFFEXP", "DIFFEXPANAL", "DIFFEXPENRICHANAL")) {
    from <- "DiffExpEnrichAnal"
  }
  if (toupper(from) %in% c("COEXP", "COEXPANAL", "COEXPENRICHANAL")) {
    from <- "CoExpEnrichAnal"
  }
  
  if (is.null(contrastName)) {
    res_return <- object@metadata[[from]][[database]][["enrichResult"]]
  } else {
    if (isTagName(object, contrastName))
      contrastName <- convertTagToContrast(object, contrastName)
    res_return <- object@metadata[[from]][[database]]$enrichResult[[contrastName]]
  }
  
  
}


#
#' @title Get a particular enrichment result
#' @description
#' Called inside the getEnrichSum method
#' 
#' @param object a SE object or a MAE object (produced by Flomics).
#' @return enrichment result summary
#' @noRd
#' @keywords internal
.getEnrichSumIntSE <- function(object,
                               from = "DiffExpEnrichAnal",
                               database = "GO"){
  return(object@metadata[[from]][[database]][["summary"]])
}
