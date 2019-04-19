#' @docType class
#'
#' @title FlomicsExperiment-Class
#' @description This class extend the [\code{\link{SummarizedExperiment}}] class by adding a named vector
#' (colDataStruc) giving the number of factors in the RNAseq design and the number of QC parameters.
#' @slot colDataStruc giving the number of factor in the RNAseq design and the number of QC parameters.
#' @return FlomicsExperiment object
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @name FlomicsExperiment-class
#' @rdname FlomicsExperiment-class
#' @exportClass FlomicsExperiment

.FlomicsExperiment <- setClass(
  Class="FlomicsExperiment",
  representation=representation(colDataStruc="vector",design="ExpDesign"),
  contains="SummarizedExperiment"
        )


#'@title ExpDesign Class
#'
#' @slot List.Factors list.
#' @slot Factors.Type vector.
#' @slot Model.formula formula.
#' @slot Contrasts.List list.
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
      Model.formula="formula",
      Contrasts.List="list"
  ))

