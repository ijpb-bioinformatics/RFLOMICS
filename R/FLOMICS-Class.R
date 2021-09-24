#' @title \link{ExpDesign-class} Class
#' @description
#' This class, is an internal class to the RFLOMICS package. An object of this class stores all the
#' objects related to the experimental design setup and the statistical modeling of the data.
#' This informations are required to run the omics analysis and this object is an input to
#' the \link{FlomicsMultiAssay.constructor.}
#' @slot ExpDesign  An object of class data.frame. Row names must give the name of each sample which has been to be construct
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Column are vector of character storing the factor
#' madality for each sample. Column names give the name of the experimental factor.
#' @slot List.Factors A named List of factor giving the modality of the experimental factor for each sample.
#' @slot Factors.Type A named vector of character giving the type of effect of each factor. Two types of effect are required
#' ("Bio" or "batch")
#' @slot Groups An object of class data.frame. Row names give the name of each sample which has been to be construct
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Two columns, the first one give the group
#' of the sample by concatenating all factor modality except the replicat. The second one give the name of the sample.
#' @slot Model.formula An object of class formula.
#' @slot Contrasts.List A list of data.frame giving all the hypothesis (contrasts) that can be formulated according to a model.
#' This list is generated thanks to the \link{getExpressionContrast} method. Three type of contrasts could be generated: simple,
#' averaged, interaction. For details see the help of the \link{getExpressionContrast} method.
#' @slot Contrasts.Sel A data.frame giving the hypothesis (contrasts) that the user selected through the interface.
#' @slot projectName A vector of character giving the name of the project.
#' @slot Contrasts.Coeff An object of class data.frame giving the coefficient for each contrast.
#' @return ExpDesign
#' @seealso \link{getExpressionContrast}
#' @examples
#' Design.File <- read.table(file= paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/experimental_design.txt",sep=""),
#' header = TRUE,row.names = 1, sep = "\t")
#'
#' # Define the type of each factor
#' Design.typeList <- c("Bio","Bio","batch")
#'
#' # Define the reference modality for each factor
#' Design.refList <- c("WT","control","rep1")
#'
#' # Initialize an object of class ExpDesign
#' Design.obj <- ExpDesign.constructor(ExpDesign = Design.File, projectName = "Design.Name",
#' refList = Design.refList, typeList = Design.typeList)
#' class(Design.obj)
#'
#' @name ExpDesign-class
#' @rdname ExpDesign-class
#' @exportClass ExpDesign
#'
.ExpDesign <- setClass(
  Class="ExpDesign",
  slots=c(ExpDesign="data.frame",
          List.Factors="list",
          Factors.Type="vector",
          Groups="data.frame",
          Model.formula="vector",
          #Model.matrix="vector",
          Contrasts.List="list",
          Contrasts.Sel="data.frame",
          Contrasts.Coeff="data.frame")
  )


# peut-être mettre une formule dans Model.formula

