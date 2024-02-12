



# ---- DO NOT PLOT function ----

#' @title doNotPlot
#' @description
#' Used mainly for the interface to check some conditions before actually plotting said graph.
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
  capture.output(out <- tryCatch(
    {
      suppressWarnings(
        suppressMessages(eval(expr))
      )
      
    },
    error = function(e) e,
    warning = function(w) w
  ))
  return(out)
}


.addBSpopify <- function(label="", content="", title="", 
                         color="black", placement="right", trigger = "click"){
  
  id <- paste0("id" , paste0(sample(letters, 4, replace = TRUE), collapse = ""))
  span(label,
       popify(actionLink(id, icon("question-circle")), title=title, content=content,
              trigger = trigger, placement=placement),
       style=paste0("color:", color))
  #tags$a(icon("question-circle"))
}






