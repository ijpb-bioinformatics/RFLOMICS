######## INTERNAL - ANNOTATION CLUSTERPROFILER #########

#' @title see_pathview
#' @param ... Possible arguments for pathview function.
#' @return nothing. Plot on the currently opened device.
#' @keywords internal
#' @importFrom grid grid.raster
#' @importFrom stringr str_split
#' @importFrom pathview pathview
#' @importFrom png readPNG
#' @importFrom utils data
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
    filename <- sapply(strsplit(msg, " "), function(x)
        x[length(x)])
    if (length(filename) > 0) {
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

.getEnrichResIntSE <- function(object,
                               contrastName,
                               from,
                               database) {
    if (is.null(contrastName)) {
        toret <- metadata(object)[[from]][[database]][["enrichResult"]]
    } else {
        toret <-
            metadata(object)[[from]][[database]]$enrichResult[[contrastName]]
    }
    
    return(toret)
}

#' @title Determine the origin of from argument
#' @description
#' Called inside several enrichment methods
#'
#' @param from a character string, usually containing either diffexp or coexp
#' @return from of the right form
#' @noRd
#' @keywords internal
#'
.determineFrom <- function(from) {
    searchFrom <-
        as.character(c(1, 2)[c(grepl("DIFFEXP", toupper(from)),
                               grepl("COEXP", toupper(from)))])
    if (length(searchFrom) < 1)
        stop(from, " doesn't exist")
    
    from <- switch(searchFrom,
                   "1" = {
                       "DiffExp"
                   },
                   "2" = {
                       "CoExp"
                   },
                   {
                       message(
                           "Argument from is detected to be neither
                          DiffExp nor CoExp, taking DiffExp results."
                       )
                       "DiffExp"
                   })
    
    return(from)
}

#' @title Get a particular enrichment result
#' @description
#' Called inside several enrichment methods
#'
#' @param from a character string, usually containing either diffexp or coexp
#' @return from of the right form
#' @noRd
#' @keywords internal
#'
.determineFromEnrich  <- function(from) {
    searchFrom <-
        as.character(c(1, 2)[c(grepl("DIFFEXP", toupper(from)),
                               grepl("COEXP", toupper(from)))])
    if (length(searchFrom) < 1)
        stop(from, " doesn't exist")
    
    from <- switch(searchFrom,
                   "1" = {
                       "DiffExpEnrichAnal"
                   },
                   "2" = {
                       "CoExpEnrichAnal"
                   },
                   {
                       message(
                           "Argument from is detected to be neither
                          DiffExp nor CoExp, taking DiffExp results."
                       )
                       "DiffExpEnrichAnal"
                   })
    
    return(from)
}

# ---- Valid URL ----
# From https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r

#' @title Test an url 
#' @description
#' Called for pathview before trying to get the map.
#'
#' @param url the url to test
#' @param timeout timeout for open.connection
#' @return boolean. TRUE if accessible, FALSE if not
#' @noRd
#' @keywords internal

.validUrl <- function(url, timeout=2){
    con <- url(url)
    check <- suppressWarnings(try(
        open.connection(con, 
                        open = "rt",
                        timeout = timeout),
        silent = TRUE)[1])
    suppressWarnings(try(close.connection(con), silent = TRUE))
    
    return(ifelse(is.null(check), TRUE, FALSE))
}
