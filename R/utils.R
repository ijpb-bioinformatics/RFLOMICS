#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' GetDesignFromNames
#'
#' @param samples_name a vector of sample names giving the designs factor, each
#' separated by "_"
#'
#' @return a named list of two elements
#' * nb_dfac: the number of design factors
#' * tmpDesign: a data frame with the names of design factors in columns (ie: dFac1) and all
#' the modalities in row
#' @md
#' @export
#'
#' @examples
#'
#'
GetDesignFromNames <- function(samples_name){

  # Get the number of design factor and the factors from the names of the matrix count
  nb_dFac <- stringr::str_count(samples_name,pattern="_")+1
  # Test if the number of factor are the same for all sample names
  try(if(var(nb_dFac) != 0 ) stop("Column names do not have the same level of factor"))
  nb_dFac <- nb_dFac[1]
  names(nb_dFac) <- "n_dFac"

  #  Get all the factors from the names of the count matrix
  tmpDesign <- tibble::tibble(design=samples_name) %>%
    tidyr::separate(.,design,into=paste("dFac",1:nb_dFac["n_dFac"],sep="")) %>%
    dplyr::mutate_all(.,as.factor)

  return(list("nb_dFac"=nb_dFac,"tmpDesign"=tmpDesign))
}

#' GetModelFormulae
#'
#' @param Factors.Name
#' @param Factors.Type
#'
#' @return a formulae
#' @export
#'
#' @examples
#'
GetModelFormulae <- function(Factors.Name,Factors.Type=NULL){

  formulae <- list()
  nFac <- length(Factors.Name)

  getF <- function(x){
    update(as.formula(paste("~ ","(",paste(x,collapse="+"),")^2")),new=~.)
  }
  for(i in 1:nFac){
  formulae[[i]] <- apply(combn(Factors.Name,i),2,getF)
  }
  return(unlist(formulae))
}


