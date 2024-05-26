### ============================================================================
### [02_set_stat_model] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------

#' @importFrom dplyr full_join mutate arrange select group_by_at count left_join
#' right_join mutate_at if_else 
#' @importFrom ggplot2 ggplot aes element_blank element_text geom_col theme 
#' labs scale_y_continuous geom_tile scale_fill_manual ylab xlab labs
#' facet_grid 
#' @importFrom purrr reduce
#' @importFrom magrittr "%>%" 

# ---- generateModelFormulae ----
#' @title generateModelFormulae
#' @description
#'  From a vector of character giving the name of the factors of an omics
#'  experiment, and their type of effect: biological or batch, it returns all 
#'  models formulae that can be formulated in association with this factors. 
#'  Batch effect factors do not appear in interaction terms with biological 
#'  factor. Model formulae stop in second order interaction.
#' 
#' @param object a RflomicsMAE object
#' @aliases generateModelFormulae
#'
#' @return a named list of object of class formula
#' @exportMethod generateModelFormulae
setMethod(f          = "generateModelFormulae",
          signature  = "RflomicsMAE",
          definition = function(object){
            
            return(.generateModelFormulae(getBioFactors(object),
                                          getBatchFactors(object)))
          })

# ---- setModelFormula : ----
#' @exportMethod setModelFormula
#' @rdname RflomicsMAE-accessors
setMethod(f          = "setModelFormula",
          signature  = "RflomicsMAE",
          definition = function(object, modelFormula=NULL){
            
            object@metadata$design$Model.formula <- paste(modelFormula,
                                                          collapse = " ")
            
            for(name in names(object)){
              object[[name]] <- setModelFormula(object[[name]], 
                                                modelFormula)
            }
            
            return(object)
          })


#' @exportMethod setModelFormula
#' @rdname RflomicsMAE-accessors
setMethod(f          = "setModelFormula",
          signature  = "RflomicsSE",
          definition = function(object, modelFormula=NULL){
            
            object@metadata$design$Model.formula <- paste(modelFormula,
                                                          collapse = " ")
            
            return(object)
          })

# ---- getModelFormula : Get Model Formula : ----
#' @exportMethod getModelFormula
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getModelFormula",
          signature  = "RflomicsMAE",
          definition = function(object){
            
            return(object@metadata$design$Model.formula)
          })

#' @exportMethod getModelFormula
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getModelFormula",
          signature  = "RflomicsSE",
          definition = function(object){
            
            return(object@metadata$design$Model.formula)
          })


# ---- generateExpressionContrast ----
#' @title Get model formula from a Flomics RflomicsSE
#' @description This function allows, from a model formulae, to give the 
#' expression contrast data frames.
#' Three types of contrasts are expressed:
#' \itemize{
#' \item{pairwise comparison}
#' \item{averaged expression}
#' \item{interaction expression}
#' }
#' @param object An object of class [\code{\link{RflomicsSE-class}}] 
#' or class [\code{\link{RflomicsMAE-class}}]
#' @return list of 1 or 3 data.frames of contrast expression
#' @exportMethod generateExpressionContrast
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
#' @rdname generateExpressionContrast
#' @aliases generateExpressionContrast
setMethod(f          = "generateExpressionContrast",
          signature  = "RflomicsSE",
          definition = function(object){
            
            modelFormula <- getModelFormula(object)
            # check
            if (is.null(modelFormula)) stop("model formula is mandatory.")
            if (is(modelFormula, "formula")) {
              modelFormula <- paste(as.character(modelFormula), 
                                    collapse = " ")
            }
            
            # args for getExpressionContrastF()
            factorBio <- getBioFactors(object)
            ExpDesign <- getDesignMat(object)
            
            Contrasts.List  <-  .getExpressionContrastF(ExpDesign, 
                                                        factorBio, 
                                                        modelFormula=modelFormula)
            
            return(Contrasts.List)
          })

#' @title Get model formula from a Flomics RflomicsSE
#' @description This function allows, from a model formulae, 
#' to give the expression contrast data frames.
#' Three types of contrasts are expressed:
#' \itemize{
#' \item{pairwise comparison}
#' \item{averaged expression}
#' \item{interaction expression}
#' }
#' @param object An object of class [\code{\link{RflomicsSE-class}}] 
#' or class [\code{\link{RflomicsMAE-class}}]
#' @return list of 1 or 3 data.frames of contrast expression
#' @exportMethod generateExpressionContrast
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
#' @rdname generateExpressionContrast
#' @aliases generateExpressionContrast
setMethod(f          = "generateExpressionContrast",
          signature  = "RflomicsMAE",
          definition = function(object){
            
            modelFormula <- getModelFormula(object)
            # check
            if (is.null(modelFormula)) stop("model formula is mandatory.")
            if (is(modelFormula, "formula")) {
              modelFormula <- paste(as.character(modelFormula), 
                                    collapse = " ")
            }
            
            # args for getExpressionContrastF()
            factorBio <- getBioFactors(object)
            ExpDesign <- getDesignMat(object)
            
            Contrasts.List  <- .getExpressionContrastF(ExpDesign, 
                                                       factorBio, 
                                                       modelFormula = modelFormula)
            
            return(Contrasts.List)
          })



#' @title Get model formula from a Flomics RflomicsSE
#' @description This function allows, from a model formulae, 
#' to give the expression contrast data frames.
#' Three types of contrasts are expressed:
#' \itemize{
#' \item{pairwise comparison}
#' \item{averaged expression}
#' \item{interaction expression}
#' }
#' @param object An object of class [\code{\link{RflomicsSE-class}}] 
#' or class [\code{\link{RflomicsMAE-class}}]
#' @return list of 1 or 3 data.frames of contrast expression
#' @exportMethod generateExpressionContrast
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
#' @rdname generateExpressionContrast
#' @aliases generateExpressionContrast
setMethod(f          = "generateExpressionContrast",
          signature  = "RflomicsMAE",
          definition = function(object){
            
            modelFormula <- getModelFormula(object)
            # check
            if (is.null(modelFormula)) stop("model formula is mandatory.")
            if (is(modelFormula, "formula")) {
              modelFormula <- paste(as.character(modelFormula), 
                                    collapse = " ")
            }
            
            # args for getExpressionContrastF()
            factorBio <- getBioFactors(object)
            ExpDesign <- getDesignMat(object)
            
            Contrasts.List  <- .getExpressionContrastF(ExpDesign, 
                                                       factorBio, 
                                                       modelFormula = modelFormula)
            
            return(Contrasts.List)
          })


# ---- setSelectedContrasts ----
#' @rdname RflomicsMAE-accessors
#' @exportMethod setSelectedContrasts
setMethod(f          = "setSelectedContrasts",
          signature  = "RflomicsMAE",
          definition = function(object, contrastList=NULL){
            
            # MAE
            object@metadata$design$Contrasts.Sel <- 
              updateSelectedContrasts(object, contrastList)
            
            # for each dataset
            for(SE.name in names(object)){
              
              object[[SE.name]] <- 
                setSelectedContrasts(object[[SE.name]], contrastList)
            }
            
            return(object)
          })


#' @rdname RflomicsSE-accessors
#' @exportMethod setSelectedContrasts
setMethod(f          = "setSelectedContrasts",
          signature  = "RflomicsSE",
          definition = function(object, contrastList=NULL){
            
            object@metadata$design$Contrasts.Sel <- 
              updateSelectedContrasts(object, contrastList)
            
            return(object)
          })


# ---- getSelectedContrasts : ----
#' @rdname RflomicsMAE-accessors
#' @exportMethod getSelectedContrasts
setMethod(f          = "getSelectedContrasts",
          signature  = "RflomicsMAE",
          definition = function(object){
            
            return(object@metadata$design$Contrasts.Sel)
          })


#' @rdname RflomicsSE-accessors
#' @exportMethod getSelectedContrasts
#' @aliases getSelectedContrasts
setMethod(f          = "getSelectedContrasts",
          signature  = "RflomicsSE",
          definition = function(object){
            
            return(object@metadata$design$Contrasts.Sel)
          })

# ---- generateContrastMatrix ----
#' @title generateContrastMatrix
#' @description Defines contrast matrix or contrast list with contrast 
#' name and contrast coefficients
#' @param object An object of class \link{RflomicsMAE-class}
#' @param contrastList a data.frame of contrast
#' @return An object of class \link{RflomicsSE-class}
#' @seealso generateContrastMatrix
#' @exportMethod generateContrastMatrix
#' @importFrom stats formula terms.formula
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
#' @rdname generateContrastMatrix
#' @aliases generateContrastMatrix
setMethod(f          = "generateContrastMatrix",
          signature  = "RflomicsSE",
          definition = function(object, contrastList=NULL){
            
            if(is.null(contrastList)) stop("contrastList arg is mandatory.")
            
            ExpDesign <- getDesignMat(object)
            
            factorBio <- getBioFactors(object)
            
            modelFormula <- getModelFormula(object)
            object@metadata$design$Contrasts.Coeff <- .getContrastMatrixF(ExpDesign = ExpDesign, factorBio = factorBio, contrastList = contrastList$contrast, modelFormula)
            object@metadata$design$Contrasts.Sel   <- contrastList
            
            return(object)
          })

#coefmatrices <- sapply(unique(names(coefvectors)),
#                       function(n) as.matrix(as.data.frame(coefvectors[names(coefvectors)==n])),
#                       simplify=FALSE, USE.NAMES=TRUE)

#Contrasts = list(D1vsD2          = c(1,  1, -1, -1,  0),
#                 C1vsC2          = c(1, -1,  1, -1,  0),
#                 InteractionDC   = c(1, -1, -1,  1,  0),
#                 C1vsC2forD1only = c(1, -1,  0,  0,  0),
#                 C1vsC2forD2only = c(0,  0,  1, -1,  0),
#                 TreatsvsControl = c(1,  1,  1,  1, -4),
#                 T1vsC           = c(1,  0,  0,  0, -1),
#                 T2vsC           = c(0,  1,  0,  0, -1),
#                 T3vsC           = c(0,  0,  1,  0, -1),
#                 T4vsC           = c(0,  0,  0,  1, -1))

# contrast From emmeans v1.3.5 by Russell Lenth 16th Percentile Contrasts and linear functions of EMMs
# coef returns a data.frame containing the object's grid, along with columns named c.1, c.2, ... containing the contrast coefficients.


#' @title generateContrastMatrix
#' @description Defines contrast matrix or contrast list with contrast name 
#' and contrast coefficients
#' @param object An object of class \link{RflomicsMAE-class}
#' @param SE.name the name of the \link{RflomicsSE-class} to access.
#' @param contrastList A data.frame of contrast
#' @return An object of class \link{RflomicsMAE-class}
#' @seealso generateContrastMatrix
#' @exportMethod generateContrastMatrix
#' @importFrom stats formula terms.formula
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
#' @rdname generateContrastMatrix
#' @aliases generateContrastMatrix
setMethod(f          = "generateContrastMatrix",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, contrastList=NULL){
            
            if (is.null(object[[SE.name]])) {
              stop("no Experiment named ", SE.name, " in MAE object")
            }
            if (is.null(modelFormula)) {
              stop("Model.formula arg is mandatory.")
            }
            if (is.null(contrastList)) stop("contrastList is mandatory.")
            if (any(!c("contrast", "contrastName", "groupComparison", "type") %in% names(contrastList))) {
              stop("contrastList data.frame must contain at least these colomn : contrast, contrastName, groupComparison, type")
            }
            
            object <- setModelFormula(object, modelFormula)
            
            object[[SE.name]] <- generateContrastMatrix(object = object[[SE.name]], contrastList = contrastList)
            
            return(object)
          })





# ---- Set Valid Contrasts : (after differential analysis) ----

#' @title setValidContrasts
#' @description Set Valid Contrasts
#' @param object An object of class \link{RflomicsMAE-class}
#' @param contrastList A data.frame of contrast
#' @param omicName a dataset name
#' @return An object of class \link{RflomicsMAE-class}
#' @exportMethod setValidContrasts
#' @rdname setValidContrasts
#' @aliases setValidContrasts
setMethod(f          = "setValidContrasts",
          signature  = "RflomicsMAE",
          definition <- function(object, omicName=NULL, contrastList=NULL){
            
            object[[omicName]] <- 
              setValidContrasts(object[[omicName]], contrastList = contrastList)
            
            return(object)
          })

#' @title setValidContrasts
#' @description Set Valid Contrasts
#' @param object An object of class \link{RflomicsSE-class}
#' @param contrastList A data.frame of contrast
#' @return An object of class \link{RflomicsSE-class}
#' @exportMethod setValidContrasts
#' @rdname setValidContrasts
#' @aliases setValidContrasts
setMethod(f          = "setValidContrasts",
          signature  = "RflomicsSE",
          definition = function(object, contrastList=NULL){
            
            object@metadata$DiffExpAnal[["Validcontrasts"]] <- contrastList
            
            return(object)
          })

#' @title getValidContrasts
#' @description Set Valid Contrasts
#' @param object An object of class \link{RflomicsMAE-class}
#' @param omicName a dataset name
#' @return An object of class \link{RflomicsMAE-class}
#' @exportMethod getValidContrasts
#' @rdname getValidContrasts
#' @aliases getValidContrasts
setMethod(f          = "getValidContrasts",
          signature  = "RflomicsMAE",
          definition = function(object, omicName){
            
            res <- getValidContrasts(object[[omicName]])
            
            return(res)
          })

#' @title getValidContrasts
#' @description Set Valid Contrasts
#' @param object An object of class \link{RflomicsSE-class}
#' @return An object of class \link{RflomicsSE-class}
#' @exportMethod getValidContrasts
#' @rdname getValidContrasts
#' @aliases getValidContrasts
setMethod(f          = "getValidContrasts",
          signature  = "RflomicsSE",
          definition = function(object){
            
            return(object@metadata$DiffExpAnal[["Validcontrasts"]])
          })




