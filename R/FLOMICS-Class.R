#' @title [\code{\link{ExpDesign-class}}] Class
#' @description
#' This class, is an internal class to the RFLOMICS package, to store all the information related to
#' the experimental design. An object of this class is instanciated when the shiny UI
#' @slot ExpDesign  An object of class data.frame. Row names give the name of each sample which has been to be construct
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Column names give the name of an experimental
#' factor which is a vector of character storing the factor madality for each sample.
#' @slot List.Factors A named List of factor giving the modality of the experimental factor for each sample.
#' @slot Factors.Type A named vector of character giving the type of effect of each factor. Two types of effect are required
#' ("Bio" or "batch")
#' @slot Groups An object of class data.frame. Row names give the name of each sample which has been to be construct
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Two columns, the first one give the group
#' of the sample by concatenating all factor modality except the replicat. The second one give the name of the sample.
#' @slot Model.formula An object of class formula corresponding to the model that the user choosed  via the shiny UI
#' at the experimental design set up.
#' @slot Contrasts.List A list of data.frame giving all the hypothesis (contrasts) that can be formulated according to a model.
#' This list is generated thanks to the getExpressionContrast() method.
#' @slot Contrasts.Sel A data.frame giving the hypothesis (contrasts) that the user selected through the interface.
#' @slot projectName A vector of character giving the name of the project.
#' @slot Contrasts.Coeff An object of class data.frame giving the coefficient for each contrast.
#' @return ExpDesign object
#' @examples
#' @name ExpDesign-class
#' @rdname ExpDesign-class
#' @exportClass ExpDesign
.ExpDesign <- setClass(
  Class="ExpDesign",
  slots=c(ExpDesign="data.frame",
          projectName="character",
          List.Factors="list",
          Factors.Type="vector",
          Groups="data.frame",
          Model.formula="vector",
          #Model.matrix="vector",
          Contrasts.List="list",
          Contrasts.Sel="data.frame",
          Contrasts.Coeff="data.frame")
  )



