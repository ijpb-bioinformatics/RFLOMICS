



# ---- DO NOT PLOT function ----

#' @title doNotPlot
#' @description
#' Used mainly for the interface to check some conditions before actually 
#' plotting said graph.
#'
#' @param expr An expression, usually producing a plot but not necessarily.
#' @keywords internal
#' @noRd
#' @importFrom utils capture.output
#'
.doNotPlot <- function(expr) {
    pdf(file = NULL)
    out <- tryCatch(
        {
            capture.output(
                suppressMessages(
                    eval(expr)
                )
            )
        },
        error = function(e) e,
        warning = function(w) w
    )
    dev.off()
    return(out)
}


#' @title doNotSpeak
#' @description
#' Used mainly for the interface to silence some functions.
#'
#' @param expr An expression, usually producing a warning.
#' @keywords internal
#' @noRd
#'
.doNotSpeak <- function(expr) {
    capture.output(out <- tryCatch({
        eval(expr)},
        error = function(e) e,
        warning = function(w) w
    ))
    return(out)
}


.addBSpopify <- function(label="", content="", title="", 
                         color="black", placement="right",
                         trigger = "click"){
    
    id <- paste0("id" , paste0(sample(letters, 4, replace = TRUE), 
                               collapse = ""))
    span(label,
         popify(actionLink(id, icon("question-circle")), 
                title = title, content = content,
                trigger = trigger, placement = placement),
         style = paste0("color:", color))
}

# ---- INTERNAL FUNCTIONS ----
# ---- isContrastName : ----
#' @title Check if character vectors are contrasts Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). 
#' If it's a RflomicsSE, expect to find
#'  a slot of differential analysis.
#' @param contrastName vector of characters.
#' @return boolean. TRUE if all of contrastName are indeed contrasts Names.
#' @noRd
#' @keywords internal
.isContrastName <- function(object, contrastName) {
    df_contrasts <- getSelectedContrasts(object)
    
    search_match <- lapply(contrastName, FUN = function(cn) {
        grep(cn, df_contrasts$contrastName, fixed = TRUE)
    })
    search_success <- unlist(lapply(search_match, identical, integer(0))) 
    # if TRUE, not a success at all.
    
    if (!any(search_success)) {
        # Congratulations, it's a contrast name!
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# ---- convertTagToContrast - convert tag to contrastName ----

#' @title Convert tags names to contrast Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). 
#' If it's a RflomicsSE, expects to find
#'  a slot of differential analysis.
#' @param tagName Vector of characters, expect to be tags 
#' (in the form of H1, H2, etc.).
#' @return character vector, contrastNames associated to tags.
#' @importFrom dplyr filter select
#' @noRd
#' @keywords internal
.convertTagToContrast <- function(object, tagName) {
    df_contrasts <- getSelectedContrasts(object)
    
    df_contrasts %>%
        filter(tag %in% tagName) %>%
        select(contrastName) %>%
        unlist(use.names = FALSE)
}

# ---- convertContrastToTag - convert contrastName to tag ----

#' @title Convert contrast Names names to tags
#'
#' @param object a MAE object or a SE object (produced by Flomics). 
#' If it's a RflomicsSE, expects to find
#'  a slot of differential analysis.
#' @param contrasts Vector of characters, expect to be contrast names.
#' @return character vector, tags associated to contrast names.
#' @importFrom dplyr filter select
#' @noRd
#' @keywords internal
.convertContrastToTag <- function(object, contrasts) {
    df_contrasts <- getSelectedContrasts(object)
    
    df_contrasts %>%
        filter(contrastName %in% contrasts) %>%
        select(tag) %>%
        unlist(use.names = FALSE)
}

# ---- isTagName: Check if character vectors are tags Names ----
#' @title Check if character vectors are tags Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). 
#' If it's a RflomicsSE, expect to find
#'  a slot of differential analysis.
#' @param tagName vector of characters.
#' @return boolean. TRUE if all of tagName are indeed tags Names.
#' @noRd
#' @keywords internal
.isTagName <- function(object, tagName) {
    df_contrasts <- getSelectedContrasts(object)
    
    search_match <- lapply(tagName, FUN = function(cn) {
        grep(cn, df_contrasts$tag, fixed = TRUE)
    })
    # if TRUE, not a success at all.
    search_success <- unlist(lapply(search_match, identical, integer(0)))
    
    if (!any(search_success)) {
        # Congratulations, it's a tag name!
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# ---- .getOrigin - get origin of a particular name ----
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
    
    if (.isContrastName(object, name)) return("Contrast")
    if (.isTagName(object, name)) return("Tag")
    if (.isClusterName(object, name)) return("CoexCluster")
    
    return("NoOriginFound")
}

# ---- isClusterName - Check if character vectors are tags Names : -----

#' @title Check if character vectors is a cluster name
#'
#' @param object a SE object (produced by Flomics). Expects to find
#'  a slot of coExpression analysis
#' @param clusterName vector of characters. For clusters, please
#' specify cluster.1, cluster.2, ... although 1,2,3 can work as well.
#' @return boolean. TRUE if all of tagName are indeed tags Names.
#' @noRd
#' @importFrom coseq clusters
#' @keywords internal
.isClusterName <- function(object, clusterName) {
    resClus <- object@metadata$CoExpAnal$coseqResults
    
    if (is.null(resClus)) {
        warning("No coseq results in this object")
        return(FALSE) 
    }
    
    clusterPoss <- unique(clusters(resClus))
    
    if (is.integer(clusterName)){ 
        clusterName <- paste("cluster", clusterName, sep = ".")
    }
    namesClust <- paste("cluster", clusterPoss, sep = ".")
    
    search_match <- unlist(lapply(clusterName, FUN = function(cn) {
        grep(cn, namesClust, fixed = TRUE)
    }))
    search_success <- unlist(lapply(search_match, identical, integer(0))) 
    # if TRUE, not a success at all.
    
    if (!any(search_success)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}



