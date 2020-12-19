#' @title ExpDesign Class
#' @slot List.Factors A list of factor
#' @slot Factors.Type Either 'Biological' or 'batch'
#' @slot Model.formula Modele formula.
#' @slot Contrasts.List A list of contrasts.
#' @slot Model.matrix .
#' @return ExpDesign object
#' @examples
#' @name ExpDesign-class
#' @rdname ExpDesign-class
#' @exportClass ExpDesign
.ExpDesign <- setClass(
  Class="ExpDesign",
  representation=representation(
      ExpDesign="data.frame",
      List.Factors="list",
      Factors.Type="vector",
      Model.formula="vector",
      Model.matrix="vector",
      Contrasts.List="list",
      Contrasts.Sel="data.frame",
      Contrasts.Coeff="data.frame"
      
  ))

setClass("MultiAssayExperiment")

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



