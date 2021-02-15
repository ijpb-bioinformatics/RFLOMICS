#' @title ExpDesign Class
#' This class, is an internal class to the RFLOMICS package, to store all information related to
#' the experimental design set up.
#' @slot List.Factors A vector of character giving the name of each experimental factor
#' @slot Factors.Type A vector of character giving the type of effect for each factor ("Bio" or "batch")
#' @slot Groups design groups
#' @slot Model.formula The model formula that has been choosed.
#' @slot Model.matrix The associated matrix
#' @slot Contrasts.List A list of contrasts.
#' @slot Contrasts.Sel selected contrast
#' @slot Contrasts.Coeff contrast vector
#' @return ExpDesign object
#' @examples
#' @name ExpDesign-class
#' @rdname ExpDesign-class
#' @exportClass ExpDesign
.ExpDesign <- setClass(
  Class="ExpDesign",
  slots=c(ExpDesign="data.frame",
          List.Factors="list",
          Factors.Type="vector",
          Groups="data.frame",
          Model.formula="vector",
          Model.matrix="vector",
          Contrasts.List="list",
          Contrasts.Sel="data.frame",
          Contrasts.Coeff="data.frame")
  )



#' @title DiffAnalysis class
#' @description  Store the results of a differential analysis
#' @slot method The name of the Differential method applied
#' @slot ListOfDEResults A list of data.frame with the following columns:
#' \describe{
#' \item{id}{feature identifier}
#' \item{BaseMean}{}
#' \item{BaseMeanCondA}{}
#' \item{BaseMeanCondB}{}
#' \item{logFC}{}
#' \item{pvalue}{}
#' \item{FDR}{}
#' }
#' @slot ListOfRawMethodResults A list of objects returned by the used method (ex: list of object of class DGELRT)
#'
#' @return
#' @examples
#' @export
#' @name DiffAnalysis-class
#' @rdname DiffAnalysis-class
#' @exportClass DiffAnalysis


.DiffAnalysis <- setClass(
  Class="DiffAnalysis",
  representation=representation(
  method="vector",
  ListOfDEResults="list",
  ListOfRawMethodResults="list"
))



