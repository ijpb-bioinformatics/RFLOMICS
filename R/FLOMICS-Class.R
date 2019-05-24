#' @docType class
#'
#' @title FlomicsExperiment-Class
#' @description This class extend the [\code{\link{SummarizedExperiment}}] class by adding a named vector
#' (colDataStruc) giving the number of factors in the RNAseq design and the number of QC parameters.
#' @slot colDataStruc giving the number of factor in the RNAseq design and the number of QC parameters.
#' @slot ExperimentType One of the following type of experiment: "RNAseq", "Proteomic","Metabolomic"
#' @slot design An object of class [\code{\link{ExpDesign}}]
#' @slot LogFilter A list
#' @slot LogInput  A list
#' @slot LogTransform A list
#' @slot Normalization An object of class [\code{\link{Normalization}}]
#' @slot DiffAnalysis An object of class [\code{\link{DiffAnalysis}}]
#' @return FlomicsExperiment object
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @name FlomicsExperiment-class
#' @rdname FlomicsExperiment-class
#' @exportClass FlomicsExperiment

.FlomicsExperiment <- setClass(
  Class="FlomicsExperiment",
  representation=representation(
    ExperimentType="vector",
    # Design
    design="ExpDesign",
    # colData Structure
    colDataStruc="vector",
    # Pre-processing
    LogFilter="list",
    LogInput="list",
    LogTransform="list",
    # Normalization
    Normalization="Normalization",
    # DE
    DiffAnalysis="DiffAnalysis"),
    contains="SummarizedExperiment"
        )


#' @title ExpDesign Class
#'
#' @slot List.Factors A list of factor
#' @slot Factors.Type Either 'Biological' or 'batch'
#' @slot Model.formula Modele formula.
#' @slot Contrasts.List A list of contrasts.
#'
#' @return ExpDesign object
#' @examples
#' @name ExpDesign-class
#' @rdname ExpDesign-class
#' @exportClass ExpDesign

.ExpDesign <- setClass(
  Class="ExpDesign",
  representation=representation(
      List.Factors="list",
      Factors.Type="vector",
      Model.formula="vector",
      Contrasts.List="list"
  ))


#' @title Normalization class
#'
#' @slot method A vector giving the name of the normalization method
#' @slot Norm.factors A named vector giving for each samples, its factor of normalization
#'
#' @return An instance of the Normalization class
#' @examples
#' @name Normalization-class
#' @rdname Normalization-class
#' @exportClass Normalization
#'


.Normalization <- setClass(
  Class="Normalization",
  representation=representation(
    # TMM by default for RNAseq
    Method="vector",
    Norm.factors="vector"
))


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



