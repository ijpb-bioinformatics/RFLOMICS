
################################## EXPERIMENTAL DESIGN SET UP #################


### ExpDesign CLASS Constructor

#' @title Constructor for the class \link{ExpDesign-class}
#' @description This method initialize an object of class \link{ExpDesign-class}.
#' @param ExpDesign a data.frame. Row names give the name of each sample which has been to be construct
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Column names give the name of
#' an experimental factor which is a vector of character storing the factor modality for each sample.
#' @param refList A list of string giving the reference modality for each factor.
#' @param typeList A vector of string indicating the type of each experimental factor. Two types of effect
#' are required ("Bio" or "batch")
#' @return An object of class \link{ExpDesign-class}
#' @examples
#' Design.File <- read.table(file= paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/experimental_design.txt",sep=""),header = TRUE,row.names = 1, sep = "\t")
#'
#' # Define the type of each factor
#' Design.typeList <- c("Bio","Bio","batch")
#'
#' # Define the reference modality for each factor
#' Design.refList <- c("WT","control","rep1")
#'
#' # Initialize an object of class ExpDesign
#' Design.obj <- ExpDesign.constructor(ExpDesign = Design.File,
#' refList = Design.refList, typeList = Design.typeList)
#' @name ExpDesign-Constructor
#' @rdname ExpDesign-Constructor
#' @export
#' @importFrom stats relevel
#' @importFrom methods new
ExpDesign.constructor <- function(ExpDesign, refList, typeList){

  # Create the List.Factors list with the choosen level of reference for each factor
  dF.List <- lapply(1:dim(ExpDesign)[2], function(i){
    relevel(as.factor(ExpDesign[[i]]), ref=refList[i])
  })
  names(dF.List) <- names(ExpDesign)

  # Factors.Type
  names(typeList) <- names(ExpDesign)

  # Create the groups data.frame
  # groups <- tidyr::unite(as.data.frame(ExpDesign[typeList == "Bio"]), col="groups", sep="_", remove = TRUE) %>%
  #           dplyr::mutate(samples = rownames(.))

  groups <- ExpDesign %>% as.data.frame() %>%
    dplyr::mutate(samples = rownames(.)) %>%
            tidyr::unite(names(typeList[typeList == "Bio"]), col="groups", sep="_", remove = FALSE)


  Design = new(Class = "ExpDesign",
               ExpDesign=as.data.frame(ExpDesign),
               List.Factors=dF.List,
               Factors.Type=typeList,
               Groups=groups,
               Model.formula=vector(),
               Contrasts.List=list(),
               Contrasts.Sel=data.frame(),
               Contrasts.Coeff=data.frame())

  return(Design)
}

#
# Error quand plus de 3 facteurs bio et plus de 1 facteur batch
# TEST(design_nbbio_3)
# TEST(design_nbbatch_1)


###### METHOD to check the completness of the ExpDesign


#' @title CheckExpDesignCompleteness
#' @description This method check some experimental design characteristics.
#' \itemize{
#'  \item{does it have biological factor ?}
#'  \item{does it have replicat/batch factor ?}
#'  \item{does it have enough replicat/batch ? at least 3 are advised}
#'  \item{is the design completed ? presence of all possible combinations of levels for all factors}
#'  \item{is the design balanced ? presence of the same number of replicat for all possible combinations}
#'  }
#'  Completed design and at least on biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{ExpDesign-class}
#' @return a named list of two objects
#' \itemize{
#' \item{"count:"}{ a data.frame with the number of each possible combinations of levels for all factors.}
#' \item{"message:"}{ results of the check}
#' \itemize{
#'  \item{"true" :"}{  false or "The experimental design is complete and balanced}
#'  \item{"lowRep" :"}{  false or "WARNING : 3 biological replicates are needed.}
#'  \item{"noCompl" :"}{ false or "ERROR : The experimental design is not complete.}
#'  \item{"noBalan" :"}{ warning or "WARNING : The experimental design is complete but not balanced.}
#'  \item{"noBio" :"}{ false or "ERROR : no bio factor !}
#'  \item{"noBatch :"}{ false or "ERROR : no replicat}
#'  }
#'  }
#' @exportMethod CheckExpDesignCompleteness
#' @examples
#' Design.File <- read.table(file= paste(path.package("RFLOMICS"), "/ExamplesFiles/TP/experimental_design.txt",sep=""),header = TRUE,row.names = 1, sep = "\t")
#'
#' # Define the type of each factor
#' Design.typeList <- c("Bio","Bio","batch")
#'
#' # Define the reference modality for each factor
#' Design.refList <- c("WT","control","rep1")
#'
#' # Initialize an object of class ExpDesign
#' Design.obj <- ExpDesign.constructor(ExpDesign = Design.File,
#' refList = Design.refList, typeList = Design.typeList)
#' CheckExpDesignCompleteness(Design.obj)

setMethod(f="CheckExpDesignCompleteness",
          signature="ExpDesign",
          definition <- function(object, colnames=NULL){

            # output list
            output <- list()

            # check presence of bio factors
            if (! table(object@Factors.Type)["Bio"] %in% 1:3){

              stop("ERROR : no bio factor ! or nbr of bio factors exeed 3!")
              
              # message <- "noBio"
              # 
              # group_count  <- object@List.Factors[object@Factors.Type == "batch"] %>% as.data.frame() %>% table() %>% as.data.frame()
              # names(group_count)[names(group_count) == "Freq"] <- "Count"
              # output[["count"]]   <- group_count
            }
            if (table(object@Factors.Type)["batch"] == 0){
              
              stop("ERROR : no replicate!")
            }
            
            # count occurence of bio conditions
            if(is.null(colnames)){
              
              ExpDesign <- object@ExpDesign
            }
            else{
              
              ExpDesign <- dplyr::filter(object@ExpDesign, rownames(object@ExpDesign) %in% colnames)
            }
            
           

            group_count <- ExpDesign[object@Factors.Type == "Bio"] %>% dplyr::group_by(temperature, imbibition) %>% dplyr::count(name = "Count")
            #group_count  <- object@List.Factors[object@Factors.Type == "Bio"] %>% as.data.frame() %>% table() %>% as.data.frame()

           

            # check presence of relicat / batch
            # check if design is complete
            # check if design is balanced
            # check nbr of replicats

            message <- dplyr::if_else(min(group_count$Count) == 0 ,           "noCompl",
                       dplyr::if_else(length(unique(group_count$Count)) != 1, "noBalan",
                       dplyr::if_else(max(group_count$Count) < 3,             "lowRep" , "true")))

          

            output[["count"]]   <- group_count
            # switch pour message complet
            output[["message"]] <- switch(message ,
                                          "true"       = { c("true",    "The experimental design is complete and balanced.") },
                                          "lowRep"     = { c("warning", "WARNING : 3 biological replicates are needed.") },
                                          "noCompl"    = { c("false",   "ERROR : The experimental design is not complete.") },
                                          "noBalan"    = { c("warning", "WARNING : The experimental design is complete but not balanced.") })
            return(output)
          })

# print output
# warining -> warning
# false -> stop


###### METHOD which generate the contrasts expression


#' @title getExpressionContrast
#' @description This function allows, from a model formulae, to give the expression contrast data frames.
#' Three types of contrasts are expressed:
#' \itemize{
#' \item{simple}
#' \item{pairwise comparison}
#' \item{averaged expression}
#' }
#' @param model.formula a model formula
#' @return An object of class [\code{\link{ExpDesign-class}}]
#' @exportMethod getExpressionContrast
#'
#' @examples
#' Design.File <- read.table(file= paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/experimental_design.txt",sep=""), header = TRUE,row.names = 1, sep = "\t")
#'
#' # Define the type of each factor
#' Design.Factors.Type <- c("Bio","Bio","batch")
#'
#' # Define the reference modality for each factor
#' Design.Factors.Ref <- c("WT","control","rep1")
#'
#' # Initialize an object of class ExpDesign
#' Design.obj <- ExpDesign.constructor(ExpDesign = Design.File, projectName = "Design.Name", refList = Design.Factors.Ref,
#' typeList = Design.Factors.Type)
#' Design.Factors.Name <- names(Design.File)
#'
#' # Set the model formulae
#' Design.formulae <- GetModelFormulae(Factors.Name = Design.Factors.Name,Factors.Type=Design.Factors.Type)
#' Design.formulae[[1]]
#'
#' # Obtained the Expression of Contrasts
#' Design.obj <- getExpressionContrast(object = Design.obj, model.formula = names(Design.formulae[1]))
#'
#' @author Christine Paysant-Le Roux
#'
setMethod(f="getExpressionContrast",
          signature="ExpDesign",
          definition <- function(object, model.formula){

            # model formula
            modelFormula <- formula(model.formula)

            #Design@Model.formula <- formula(model.formula)
            object@Model.formula <- model.formula

            # bio factor list in formulat
            labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")

            FactorBioInDesign <- intersect(names(object@Factors.Type[object@Factors.Type == "Bio"]), labelsIntoDesign)

            #BioFactors <- object@List.Factors[FactorBioInDesign]

            treatmentFactorsList <- lapply(FactorBioInDesign, function(x){paste(x, unique(object@List.Factors[[x]]), sep="")})
            names(treatmentFactorsList) <- FactorBioInDesign

            interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
            # define all simple contrasts pairwise comparisons
            allSimpleContrast_df <- defineAllSimpleContrasts(treatmentFactorsList)
            listOfContrastsDF <- list(simple = allSimpleContrast_df)

            # define all simples contrast means
            # exists("allSimpleContrast_df", inherits = FALSE)
            if(length(treatmentFactorsList) != 1){
              allAveragedContrasts_df <- define_averaged_contrasts (allSimpleContrast_df)
              listOfContrastsDF[["averaged"]] <- allAveragedContrasts_df
            }

            # define all interaction contrasts
            if(length(treatmentFactorsList) != 1){
              if(interactionPresent){
                labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
                labelOrder <- attr(terms.formula(modelFormula), "order")
                twoWayInteractionInDesign <- labelsIntoDesign[which(labelOrder==2)]
                groupInteractionToKeep <- gsub(":", " vs ", twoWayInteractionInDesign)
                allInteractionsContrasts_df <- defineAllInteractionContrasts(treatmentFactorsList, groupInteractionToKeep)

                listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
              }
              #allInteractionsContrasts_df <- defineAllInteractionContrasts(treatmentFactorsList)
              #listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
            }
            # choose the contrasts and rbind data frames of contrasts
            #selectedContrasts <- returnSelectedContrasts(listOfContrastsDF)

            # replace interactive selection of contrasts by return all contrasts -> shiny
            object@Contrasts.List  <- listOfContrastsDF
            object@Contrasts.Coeff <- data.frame()
            object@Contrasts.Sel   <- data.frame()

            return(object)
          })


###### METHOD to obtain the Matrix of contrast with their names and coefficients

#
#

#' @title getContrastMatrix
#' @description Define contrast matrix or contrast list with contrast name and contrast coefficients
#' @param An object of class \link{ExpDesign-class}
#' @param contrastList A vector of character of contrast
#' @return An object of class \link{ExpDesign-class}
#' @seealso getExpressionContrast
#' @exportMethod getContrastMatrix
#' @importFrom stats formula terms.formula
#'
#' @author Christine Paysant-Le Roux
setMethod(f="getContrastMatrix",
          signature="ExpDesign",
          definition <- function(object, contrastList){

            contrast <- contrastName <- type <- groupComparison <- NULL

            contrast.sel.list <- list()
            contrast.sel.list <- lapply(names(object@Contrasts.List), function(contrastType) {

              tmp <- object@Contrasts.List[[contrastType]] %>% dplyr::filter(contrast %in% contrastList) %>%
                dplyr::select(contrast, contrastName, type, groupComparison)
              return(tmp)
            })
            object@Contrasts.Sel <- contrast.sel.list %>% purrr::reduce(rbind) %>% dplyr::mutate(tag = paste("H", 1:dim(.)[1], sep=""))


            sampleData <-  object@ExpDesign
            selectedContrasts <- object@Contrasts.Sel$contrast

            modelFormula <- formula(object@Model.formula)
            # bio factor list in formulat
            labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
            FactorBioInDesign <- intersect(names(object@Factors.Type[object@Factors.Type == "Bio"]), labelsIntoDesign)

            #BioFactors <- object@List.Factors[FactorBioInDesign]

            treatmentFactorsList <- lapply(FactorBioInDesign, function(x){paste(x, unique(object@List.Factors[[x]]), sep="")})
            names(treatmentFactorsList) <- FactorBioInDesign

            treatmentCondenv <- new.env()

            interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
            isThreeOrderInteraction <- any(attr(terms.formula(modelFormula),"order") == 3)

            # get model matrix
            modelMatrix <- stats::model.matrix(modelFormula, data = object@List.Factors %>% as.data.frame())
            colnames(modelMatrix)[colnames(modelMatrix) == "(Intercept)"] <- "Intercept"
            # assign treatment conditions(group) to boolean vectors according to the design model matrix
            #treatmentCondenv <- new.env()
            assignVectorToGroups(treatmentFactorsList = treatmentFactorsList,
                                 modelMatrix = modelMatrix,
                                 interactionPresent = interactionPresent,
                                 isThreeOrderInteraction = isThreeOrderInteraction,
                                 treatmentCondenv = treatmentCondenv)
            # get the coefficient vector associated with each selected contrast
            # contrast <- allSimpleContrast_df$contrast[1]
            colnamesGLMdesign <- colnames(modelMatrix)


            #coefficientsMatrix <- sapply(selectedContrasts$contrast, function(x) returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
            coefficientsMatrix <- sapply(selectedContrasts, function(x) returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))

            #coefficientsMatrix <- MASS::as.fractions(coefficientsMatrix)
            colnames(coefficientsMatrix) <- selectedContrasts

            rownames(coefficientsMatrix) <- colnamesGLMdesign
            contrastMatrix <- as.data.frame(t(coefficientsMatrix))
            #contrastMatrix <- as_tibble(t(coefficientsMatrix)) %>%
            #dplyr::mutate(contrast = selectedContrasts, .before = "Intercept") %>%
            #dplyr::mutate(type = selectedContrasts$type, .after = "contrast")
            #contrastMatrix <- MASS::as.fractions(contrastMatrix)
            #contrastMatrix
            # contrastList <- as.list(as.data.frame(coefficientsMatrix))

            object@Contrasts.Coeff <- contrastMatrix
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



################################################### OMICS DATA MANAGMENT AND ANALYSIS #################


###### FlomicsMultiAssay CLASS Constructor for managing omics DATA and RESULTS


#' @title FlomicsMultiAssay.constructor Constructor for the class \link{MultiAssayExperiment-class}
#' @description This function initializes an object of class \link{MultiAssayExperiment-class}
#' from a list of omics data and an object of class \link{ExpDesign-class}.
#' @param inputs A named list of omic dataset. Names must refer to the name of the omic dataset.
#' An omics dataset must be itself a list of three objects:
#' \itemize{
#' \item{dataFile:}{the path to the omic data}
#' \item{qcFile:}{the path to an optional quality check file}
#' \item{omicType:}{Type of omic data type "None", "RNAseq", "proteomics" or "Metabolomics".}
#' }
#' @param Design An object of class \link{ExpDesign-class}
#' @return An object of class \link{MultiAssayExperiment-class}
#' @examples
#' Design.File <- read.table(file= paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/experimental_design.txt",sep=""),header = TRUE,row.names = 1, sep = "\t")
#'
#' # Define the type of each factor
#' Design.Factors.Type <- c("Bio","Bio","batch")
#'
#' # Define the reference modality for each factor
#' Design.Factors.Ref <- c("WT","control","rep1")
#'
#' # Initialize an object of class ExpDesign
#' Design.obj <- ExpDesign.constructor(ExpDesign = Design.File, projectName = "Design.Name", refList = Design.Factors.Ref,
#'  typeList = Design.Factors.Type)
#' Design.Factors.Name <- names(Design.File)
#' Design.formulae <- GetModelFormulae(Factors.Name = Design.Factors.Name,Factors.Type=Design.Factors.Type)
#' Design.formulae[[1]]
#' Design.obj <- getExpressionContrast(object = Design.obj, model.formula = names(Design.formulae[1]))
#' Design.contrastList <- lapply(Design.obj@Contrasts.List, function(x) {
#' return(x[1:2]$contrast)
#' })
#' Design.obj <- getContrastMatrix(object = Design.obj, contrastList = unlist(Design.contrastList))
#'
#'  # Create a list of datasets
#' ListofData <- list("RNAseq1"=list("dataFile"=paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/rnaseq_gene_counts.txt",sep=""),
#' "qcFile"=paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/rnaseq_bioinfo_QC.txt",sep=""), "omicType"="RNAseq"))
#' FlomicsMultiAssay.constructor(inputs = ListofData, Design=Design.obj)
#'
#' @name FlomicsMultiAssay.constructor
#' @rdname FlomicsMultiAssay.constructor
#' @export
#'
#'

# FlomicsMultiAssay.constructor <- function(inputs, Design){
# 
#   # if input == NULL
# 
#   SummarizedExperimentList <- list()
#   listmap  <- list()
#   omicList <- list()
#   k <- 0
#   for (dataName in names(inputs)){
# 
#     k <- k+1
# 
#     ## construct SummarizedExperiment for each data
#     abundance <- read.table(inputs[[dataName]][["dataFile"]], header = TRUE, row.names = 1)
# 
#     if(!is.null(inputs[[dataName]][["qcFile"]])){
#       print("# ... metadata QC...")
#       QCmat <- read.table(inputs[[dataName]][["qcFile"]], header = TRUE)
#     }
#     else{
#       QCmat <- data.frame(primary = colnames(abundance),
#                           colname = colnames(abundance),
#                           stringsAsFactors = FALSE)
#     }
# 
#     # groups
# 
#     SummarizedExperimentList[[dataName]] <- SummarizedExperiment::SummarizedExperiment(assays   = S4Vectors::SimpleList(abundance=as.matrix(abundance)),
#                                                                                        colData  = QCmat,
#                                                                                        metadata = list(omicType = inputs[[dataName]][["omicType"]],
#                                                                                                        Groups   = Design@Groups))
# 
#     # metadata for sampleMap for MultiAssayExperiment
#     listmap[[dataName]] <- data.frame(primary = as.vector(SummarizedExperimentList[[dataName]]@colData$primary),
#                                       colname = as.vector(SummarizedExperimentList[[dataName]]@colData$colname),
#                                       stringsAsFactors = FALSE)
# 
#     #
#     omicType <- inputs[[dataName]][["omicType"]]
# 
#     colnames <- c(names(omicList[[omicType]]), k)
#     omicList[[omicType]] <- c(omicList[[omicType]] ,dataName)
#     names(omicList[[omicType]]) <- colnames
# 
#   }
# 
#   FlomicsMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(experiments = SummarizedExperimentList,
#                                                                   colData     = Design@ExpDesign,
#                                                                   sampleMap   = MultiAssayExperiment::listToMap(listmap),
#                                                                   metadata    = list(design = Design,
#                                                                                      colDataStruc = c(n_dFac = dim(Design@ExpDesign)[2], n_qcFac = 0),
#                                                                                      omicList = omicList))
# 
#   return(FlomicsMultiAssay)
# }
# 
# #
# #


#' @title FlomicsMultiAssay.constructor Constructor for the class \link{MultiAssayExperiment-class}
#' @description This function initializes an object of class \link{MultiAssayExperiment-class}
#' from a list of omics data and an object of class \link{ExpDesign-class}.
#' @param inputs A named list of omic dataset. Names must refer to the name of the omic dataset.
#' An omics dataset must be itself a list of three objects:
#' \itemize{
#' \item{data:}{matrix of omic data}
#' \item{meta:}{an optional quality check data}
#' \item{omicType:}{Type of omic data type "None", "RNAseq", "proteomics" or "Metabolomics".}
#' }
#' @param Design An object of class \link{ExpDesign-class}
#' @param projectName Project name
#' @return An object of class \link{MultiAssayExperiment-class}
#' @examples
#' Design.File <- read.table(file= paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/experimental_design.txt",sep=""),header = TRUE,row.names = 1, sep = "\t")
#'
#' # Define the type of each factor
#' Design.Factors.Type <- c("Bio","Bio","batch")
#'
#' # Define the reference modality for each factor
#' Design.Factors.Ref <- c("WT","control","rep1")
#'
#' # Initialize an object of class ExpDesign
#' Design.obj <- ExpDesign.constructor(ExpDesign = Design.File, projectName = "Design.Name", refList = Design.Factors.Ref,
#'  typeList = Design.Factors.Type)
#' Design.Factors.Name <- names(Design.File)
#' Design.formulae <- GetModelFormulae(Factors.Name = Design.Factors.Name,Factors.Type=Design.Factors.Type)
#' Design.formulae[[1]]
#' Design.obj <- getExpressionContrast(object = Design.obj, model.formula = names(Design.formulae[1]))
#' Design.contrastList <- lapply(Design.obj@Contrasts.List, function(x) {
#' return(x[1:2]$contrast)
#' })
#' Design.obj <- getContrastMatrix(object = Design.obj, contrastList = unlist(Design.contrastList))
#'
#'  # Create a list of datasets
#' ListofData <- list("RNAseq1"=list("dataFile"=paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/rnaseq_gene_counts.txt",sep=""),
#' "qcFile"=paste(path.package("RFLOMICS"),"/ExamplesFiles/TP/rnaseq_bioinfo_QC.txt",sep=""), "omicType"="RNAseq"))
#' FlomicsMultiAssay.constructor(inputs = ListofData, Design=Design.obj)
#'
#' @name FlomicsMultiAssay.constructor
#' @rdname FlomicsMultiAssay.constructor
#' @export
#'
#'

FlomicsMultiAssay.constructor <- function(inputs, Design, projectName){
  
  # if input == NULL
  
  SummarizedExperimentList <- list()
  listmap  <- list()
  omicList <- list()
  k <- 0
  for (dataName in names(inputs)){
    
    k <- k+1
    
    ## construct SummarizedExperiment for each data
    abundance <- inputs[[dataName]][["data"]]
    
    if(!is.null(inputs[[dataName]][["meta"]])){
      QCmat <- inputs[[dataName]][["meta"]]
    }
    else{
      QCmat <- data.frame(primary = colnames(abundance),
                          colname = colnames(abundance),
                          stringsAsFactors = FALSE)
    }
    
    
    ###### remove row with sum == 0
    matrix <- as.matrix(abundance)
    ## nbr of genes with 0 count
    genes_flt0  <- rownames(matrix[rowSums(matrix) <= 0, ])
    ## remove 0 count
    matrix.filt  <- matrix[rowSums(matrix)  > 0, ]
    
    ##### create groups for SE
    
    # groups <- Design@Groups %>%
    #   dplyr::mutate(samples = rownames(.)) %>%
    #   tidyr::unite(names(typeList[typeList == "Bio"]), col="groups", sep="_", remove = FALSE)
    
    ###### create SE object
    SummarizedExperimentList[[dataName]] <- SummarizedExperiment::SummarizedExperiment(assays   = S4Vectors::SimpleList(abundance=as.matrix(matrix.filt)),
                                                                                       colData  = QCmat,
                                                                                       metadata = list(omicType = inputs[[dataName]][["omicType"]],
                                                                                                         Groups = Design@Groups, 
                                                                                                       rowSums.zero = genes_flt0))
    #names(assays(SummarizedExperimentList[[dataName]])) <- c(dataName)
    
    # metadata for sampleMap for MultiAssayExperiment
    listmap[[dataName]] <- data.frame(primary = as.vector(SummarizedExperimentList[[dataName]]@colData$primary),
                                      colname = as.vector(SummarizedExperimentList[[dataName]]@colData$colname),
                                      stringsAsFactors = FALSE)
    
    #
    omicType <- inputs[[dataName]][["omicType"]]
    
    colnames <- c(names(omicList[[omicType]]), k)
    omicList[[omicType]] <- c(omicList[[omicType]] ,dataName)
    names(omicList[[omicType]]) <- colnames
    
  }
  
  

  
  prepFlomicsMultiAssay <- MultiAssayExperiment::prepMultiAssay( ExperimentList = SummarizedExperimentList, 
                                           sampleMap      = MultiAssayExperiment::listToMap(listmap), 
                                           colData        = Design@ExpDesign, outFile = stdout())
  
  
  FlomicsMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(experiments = prepFlomicsMultiAssay$experiments,
                                                                  colData     = prepFlomicsMultiAssay$colData,
                                                                  sampleMap   = prepFlomicsMultiAssay$sampleMap,
                                                                  metadata    = list(colDataStruc = c(n_dFac = dim(prepFlomicsMultiAssay$colData)[2], n_qcFac = 0),
                                                                                     omicList = omicList, projectName = projectName, design = Design))
  
  
  
  return(FlomicsMultiAssay)
}



################################# EXPLORATION OF BIOLOGICAL AND TECHICAL VARIABILITY ##################################


##### Statistical METHODS for exploring biological and technical variability


#' @title RunPCA
#' @description This function performed a scaled principal component analysis on omic data stored in an object of class [\code{\link{MultiAssayExperiment-class}]
#' Results are stored in the metadata slot.
#' @param object An object of class \link{SummarizedExperiment-class}.
#' @return An object of class \link{SummarizedExperiment}
#' @exportMethod RunPCA
#' @examples
#'
setMethod(f="RunPCA",
          signature="SummarizedExperiment",
          definition <- function(object){

            # if the data has undergone a normalization (RNAseq data) 
            if(! is.null(object@metadata[["Normalization"]]$coefNorm)){
              pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE,
                                   scale=object@metadata[["Normalization"]]$coefNorm$norm.factors)+1)
              object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo),ncp = 5,graph=F)
              
            }
            # if the data has undergone a transformation (meta or prot data) 
            else if(! is.null(object@metadata$transform_method)){
              pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE) + 1)
              object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = 5,graph=F)
            }
            
            # if no transformation
            else{
              pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE) + 1)
              object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = 5,graph=F)
            }
            return(object)
          }
)

##### Graphical METHODS for exploring biological and technical variability


#' Library_size_barplot.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @return plot
#' @export
#' @importFrom ggplot2 ggplot geom_bar xlab ylab element_text ggtitle
#'
#' @examples
setMethod(f="Library_size_barplot.plot",
          signature="SummarizedExperiment",
          definition <- function(object){
            
            value    <- NULL
            warnning <- NULL
            
            if (object@metadata$omicType != "RNAseq"){
              warnning <- "WARNING : data are not RNAseq!"
            }
            
            abundances <- SummarizedExperiment::assay(object)
            samples    <- colnames(abundances)
            
            # normalized data
            if(! is.null(object@metadata$Normalization)){
              pseudo  <- scale(SummarizedExperiment::assay(object), center=FALSE,
                               scale=object@metadata$Normalization$coefNorm$norm.factors) %>% colSums(., na.rm = TRUE)
              ylab <- "Total normalized read count per sample"
              title <- "Filtered and normalized (TMM) data"
            }
            # raw data
            else{
              pseudo  <- SummarizedExperiment::assay(object) %>% colSums(., na.rm = TRUE)
              ylab <- "Total read count per sample (method : TMM)"
              title <- "Raw data"
            }
            
            libSizeNorm <- data.frame ( "value" = pseudo , "samples"=names(pseudo)) %>% dplyr::full_join(object@metadata$Groups, by="samples")
            
            libSizeNorm$samples <- factor(libSizeNorm$samples, levels = libSizeNorm$samples)
            
            p <- ggplot(libSizeNorm, aes(x=samples,y=value, fill=groups)) + geom_bar( stat="identity" ) + ylab(ylab) +
              theme(axis.text.x      = element_text(angle = 45, hjust = 1), legend.position  = "none") + labs(x = "") + ggtitle(title)
            #axis.text.x     = element_blank(),
            #axis.ticks      = element_blank())
            #legend.key.size = unit(0.3, "cm"))
            #legend.text     = element_text(size=5))
            print(p)
            
          })





#' @title Data_Distribution_Density.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @export
#' @importFrom ggplot2 geom_density xlab

setMethod(f="Data_Distribution_Density.plot",
          signature="SummarizedExperiment",
          definition <- function(object){
            
            switch (object@metadata$omicType,
                    "RNAseq" = {
                      
                      # before normalization
                      if(is.null(object@metadata[["Normalization"]]$coefNorm)){
                        pseudo <- log2(SummarizedExperiment::assay(object) + 1) 
                        y_lab  <- "log2(gene counts)"
                        title  <- "Raw data"
                        
                      }
                      # after normalization
                      else{
                        pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE,
                                             scale=object@metadata[["Normalization"]]$coefNorm$norm.factors)+1)
                        y_lab  <- "log2(normalized gene counts)"
                        title  <- "Filtered and normalized (TMM) data"
                      }
                      # 
                      pseudo.gg <- pseudo %>% reshape2::melt()
                      colnames(pseudo.gg) <- c("features", "samples", "value")
                      
                      pseudo.gg <- pseudo.gg %>% dplyr::full_join(object@metadata$Groups, by="samples") 
                      
                      p <- ggplot2::ggplot(pseudo.gg) + geom_density(aes(x=value, group = samples, color=groups), trim=FALSE) + xlab(y_lab) +
                        theme(legend.position='none') + ggtitle(title)
                    },
                    "proteomics" = {
                      # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
                      if(is.null(object@metadata$transform_method)){
                        pseudo <- SummarizedExperiment::assay(object)
                        x_lab  <- "Protein abundance (?)"
                        title  <- "Raw data"
                        
                      }
                      # after transformation
                      else{
                        # if log2 transformation was chosen
                        switch (object@metadata$transform_method,
                                "log2" = {
                                  pseudo <- SummarizedExperiment::assay(object)
                                  x_lab  <- "Transformed protein abundance" 
                                  title  <- "Transformed data (method : log2)" },
                                
                                "none" = {
                                  pseudo <- log2(SummarizedExperiment::assay(object) + 1)
                                  x_lab  <- "Transformed protein abundance" 
                                  title  <- "Transformed data (method : ?)"
                                  
                                  }
                        )
                        
                        
                        
                      }
                      # 
                      pseudo.gg <- pseudo %>% reshape2::melt()
                      colnames(pseudo.gg) <- c("features", "samples", "value")
                      
                      pseudo.gg <- pseudo.gg %>% dplyr::full_join(object@metadata$Groups, by="samples")
                      
                      p <- ggplot2::ggplot(pseudo.gg) + geom_density(aes(x=value, group = samples, color=groups), trim=FALSE) + xlab(x_lab) +
                        theme(legend.position='none') + ggtitle(title)
                    },
                    "metabolomics" = {
                      # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
                      if(is.null(object@metadata$transform_method)){
                        pseudo <- SummarizedExperiment::assay(object)
                        x_lab  <- "Metabolite abundance (?)"
                        title  <- "Raw data"
                      }
                      # after transformation
                      else{
                        # if log2 transformation was chosen
                        switch( object@metadata$transform_method,
                                "log2" = {
                                  pseudo <- SummarizedExperiment::assay(object)
                                  x_lab  <- "Transformed protein abundance"
                                  title  <- "Transformed data (method : log2)"
                                },
                                "none" = {
                                  pseudo <- log2(SummarizedExperiment::assay(object) + 1)
                                  x_lab  <- "Transformed protein abundance"
                                  title  <- "Transformed data (method : ?)"
                                }
                        )
                      }
                      # 
                      pseudo.gg <- pseudo %>% reshape2::melt()
                      colnames(pseudo.gg) <- c("features", "samples", "value")
                      
                      pseudo.gg <- pseudo.gg %>% dplyr::full_join(object@metadata$Groups, by="samples")
                      
                      p <- ggplot2::ggplot(pseudo.gg) + geom_density(aes(x=value, group = samples, color=groups), trim=FALSE) + xlab(x_lab) +
                        theme(legend.position='none') + ggtitle(title)
                    }
            )
            print(p)
          }
)





#' #' @title Data_Distribution.plot
#' #'
#' #' @param object An object of class \link{SummarizedExperiment}
#' #' @param plot Type of plot : boxplot ("boxplot") or density plot ("density")
#' #' @exportMethod Data_Distribution.plot
#' #' @export
#' #' @return plot
#' #' @rdname Data_Distribution.plot
#' #' @examples
#' #' @importFrom ggplot2 geom_density boxplot xlab
#' 
#' setMethod(f="Data_Distribution.plot",
#'           signature="SummarizedExperiment",
#'           definition <- function(object, plot = "boxplot"){
#'             
#'             switch (object@metadata$omicType,
#'                     "RNAseq" = {
#'                       
#'                       # before normalization
#'                       if(is.null(object@metadata[["Normalization"]]$coefNorm)){
#'                         pseudo <- log2(SummarizedExperiment::assay(object) + 1) %>% reshape2::melt()
#'                         y_lab  <- "log2(gene counts)"
#'                         title  <- "Raw data"
#'                         
#'                       }
#'                       # after normalization
#'                       else{
#'                         pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE,
#'                                              scale=object@metadata[["Normalization"]]$coefNorm$norm.factors)+1) %>% reshape2::melt()
#'                         y_lab  <- "log2(normalized gene counts)"
#'                         title  <- "Filtered and normalized (TMM) data"
#'                       }
#' 
#'                     },
#'                     "proteomics" = {
#'                       # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
#'                       if(is.null(object@metadata$transform_method)){
#'                         pseudo <- SummarizedExperiment::assay(object) %>% reshape2::melt()
#'                         x_lab  <- "Protein abundance (?)"
#'                         title  <- "Raw data"
#'                         
#'                       }
#'                       # after transformation
#'                       else{
#'                         # if log2 transformation was chosen
#'                         switch (object@metadata$transform_method,
#'                                 "log2" = {
#'                                   pseudo <- SummarizedExperiment::assay(object) %>% reshape2::melt()
#'                                   x_lab  <- "Transformed protein abundance" 
#'                                   title  <- "Transformed data (method : log2)" },
#'                                 
#'                                 "none" = {
#'                                   pseudo <- log2(SummarizedExperiment::assay(object) + 1) %>% reshape2::melt()
#'                                   x_lab  <- "Transformed protein abundance" 
#'                                   title  <- "Transformed data (method : ?)" } )
#'                       }
#'                     },
#'                     "metabolomics" = {
#'                       # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
#'                       if(is.null(object@metadata$transform_method)){
#'                         pseudo <- SummarizedExperiment::assay(object) %>% reshape2::melt()
#'                         x_lab  <- "Metabolite abundance (?)"
#'                         title  <- "Raw data"
#'                       }
#'                       # after transformation
#'                       else{
#'                         # if log2 transformation was chosen
#'                         switch( object@metadata$transform_method,
#'                                 "log2" = {
#'                                   pseudo <- SummarizedExperiment::assay(object) %>% reshape2::melt()
#'                                   x_lab  <- "Transformed protein abundance"
#'                                   title  <- "Transformed data (method : log2)" },
#'                                 "none" = {
#'                                   pseudo <- log2(SummarizedExperiment::assay(object) + 1) %>% reshape2::melt()
#'                                   x_lab  <- "Transformed protein abundance"
#'                                   title  <- "Transformed data (method : ?)" } )
#'                       }
#'                     }
#'             )
#'             
#'             colnames(pseudo) <- c("features", "samples", "value")
#'             pseudo.gg <- dplyr::full_join(pseudo, object@metadata$Groups, by="samples") 
#'             
#'             switch (plot,
#'               "boxplot" = {
#'                 #
#'                 pseudo_bis$samples <- factor(pseudo_bis$samples, levels = unique(pseudo_bis$samples))
#' 
#'                 p <- ggplot(pseudo_bis, aes(x=samples, y=value)) + ggplot2::geom_boxplot(aes(fill=groups)) +
#'                   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + xlab("") + ylab(y_lab) + ggtitle(title)
#'                 
#'                 
#'               },
#'               "density" = {
#'                   # 
#'                   p <- ggplot2::ggplot(pseudo.gg) + geom_density(aes(x=value, group = samples, color=groups), trim=FALSE) + xlab(y_lab) +
#'                     theme(legend.position='none') + ggtitle(title)
#'               }
#'             )
#'             
#'             print(p)
#'           }
#' )




#' @title abundanceBoxplot
#' @description This function produces boxplots from raw and normalized data matrix. One color by level
#' of combination factor. It allows to detect outlier samples and to see the normalization
#' effect.
#' @param object An object of class \link{SummarizedExperiment-class}
#' @exportMethod abundanceBoxplot
#' @rdname abundanceBoxplot
#'
setMethod(f= "abundanceBoxplot",
          signature = "SummarizedExperiment",
          definition <- function(object){
            
            samples <- value <- NULL
            
            # normalized data
            if(! is.null(object@metadata$Normalization)){
              pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE,
                                   scale=object@metadata$Normalization$coefNorm$norm.factors)+1) %>% reshape2::melt()
              y_lab  <- "log2(normalized gene counts)"
              title  <- "Raw data"
            }
            else{
              pseudo <- log2(SummarizedExperiment::assay(object)+1) %>% reshape2::melt()
              y_lab  <- "log2( gene counts)"
              title  <- "Filtred and normalized (TMM) data"
              
            }
            
            colnames(pseudo) <- c("feature", "samples", "value")
            
            pseudo_bis <- dplyr::full_join(pseudo, object@metadata$Groups, by="samples")
            
            pseudo_bis$samples <- factor(pseudo_bis$samples, levels = unique(pseudo_bis$samples))
            
            # boxplot
            p <- ggplot(pseudo_bis, aes(x=samples, y=value)) + ggplot2::geom_boxplot(aes(fill=groups)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + xlab("") + ylab(y_lab) + ggtitle(title)
            
            print(p)
            
          }
)

#' @title plotPCA
#' @description This function plot the factorial map from a PCA object stored
#' in a \link{SummarizedExperiment-class} object. By default, samples are
#' colored by groups (all combinations of level's factor)
#' @param object An object of class \link{SummarizedExperiment-class}
#' @param PCA This argument indicates whether the scaled PCA has to be performed on raw [\sQuote{raw}] or normalized [\sQuote{norm}] data.
#' @param PCs A vector giving the two axis that have to be drawn for the factorial map
#' @param condition All combination of level's factor
#' @return
#' @exportMethod plotPCA
#' @examples
#' @rdname plotPCA
setMethod(f= "plotPCA",
          signature = "SummarizedExperiment",
          definition <- function(object, PCA, PCs=c(1,2), condition="groups"){
            
            #
            PC1 <- paste("Dim.",PCs[1], sep="")
            PC2 <- paste("Dim.",PCs[2], sep="")
            
            score     <- object@metadata$PCAlist[[PCA]]$ind$coord[, PCs] %>% as.data.frame() %>%
              dplyr::mutate(samples=row.names(.)) %>% dplyr::full_join(., object@metadata$Groups, by="samples")
            
            var1 <- round(object@metadata$PCAlist[[PCA]]$eig[PCs,2][1], digits=3)
            var2 <- round(object@metadata$PCAlist[[PCA]]$eig[PCs,2][2], digits=3)
            
            switch (PCA,
                    "raw"  = {title <- paste0("Raw ", object@metadata$omicType, " data")},
                    "norm" = {title <- switch (object@metadata$omicType,
                                               "RNAseq" = {"Filtred and normalized (TMM) RNAseq data"},
                                               "proteomics" = {"Transformed proteomics data (log2)"},
                                               "metabolomics" = {"Transformed proteomics data (log2)"}
                    )}
            )
            
            
            p <- ggplot(score, aes_string(x=PC1, y=PC2, color=condition))  +
              ggplot2::geom_point(size=3) +
              ggplot2::geom_text(aes(label=samples), size=3, vjust = 0) +
              xlab(paste(PC1, " (",var1,"%)", sep="")) +
              ylab(paste(PC2, " (",var2,"%)", sep="")) +
              ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "red") +
              ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "red") +
              theme(strip.text.x = element_text(size=8, face="bold.italic"),
                    strip.text.y = element_text(size=8, face="bold.italic")) +
              ggtitle(title)
            
            print(p)
            
          })




#' @title mvQCdesign
#' @description mvQCdesign is for multivariate quality check of design. For each design factor (one color for each),
#' and each PCA axis this function plot the coordinates of the sample in a PCA axis (y-axis) in an
#' increasing order along the x-axis. It allows to have a quick view of the variability associated to each factor.
#' @param object An object of class \link{MultiAssayExperiment}
#' @param data The name of the omic data for which the PCA plot has to be drawn
#' @param axis The number of PCA axis to keep
#' @param PCA This argument indicates which PCA results to plot: raw ("raw") or normalised ("normalised")
#' @param pngFile The name of the png file for saving the plot.
#' @examples
#' @exportMethod mvQCdesign
#'
#' @rdname mvQCdesign

setMethod(f="mvQCdesign",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA=c("raw","norm"), axis=5, pngFile=NULL){

            resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <-  object@metadata$colDataStruc["n_dFac"]

            bigdf <- list()
            for(i in 1:n_dFac){
              Factor <-  object@colData[,i]
              nameFac <- names(object@colData)[i]
              df <- list()
              for(j in 1:axis){
                qc = as.vector(resPCA$ind$coord[,j])
                o = order(qc)[order(Factor[order(qc)])]
                col = cc[Factor][o]
                y=qc[o]
                Axis = rep(paste("PCA",j, "\n(Var=",round(resPCA$eig[j,2],1),")",sep=""),length(qc))
                Fac = Factor[o]
                df[[j]] <- data.frame(y,col,Axis,dfac=names(object@colData)[1:n_dFac][i],
                                      Levels=Fac,x=1:length(y))
              }
              bigdf[[i]] <- dplyr::bind_rows(df)
            }
            big <- dplyr::bind_rows(bigdf)
            out <- by(data = big, INDICES = big$dfac, FUN = function(m) {
              m <- droplevels(m)
              names(m) <- c("y", "col", "Axis", "dfac", unique(m$dfac), "x")
              m <- ggplot(m,aes(y=y,x=x),aes_string(color=unique(m$dfac)))+
                ggplot2::geom_bar(stat = "identity",position = ggplot2::position_dodge(),aes_string(fill=unique(m$dfac)))+
                facet_grid(as.factor(dfac)~Axis) +
                labs(x = "Samples", y="Coordinates on \n the PCA axis")+
                theme(axis.title.y = element_text(size = 5),
                      axis.title.x=element_text(size = 5),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
            })
            p <- do.call(gridExtra::grid.arrange, out)
            print(p)

            if(! is.null(pngFile)){
              ggplot2::ggsave(filename = pngFile,  plot = p)
            }
})

# copier coldataStruct dans metadata
# tester que PCA existe
#

#' @title mvQCdata
#' @description mvQCdata is for multivariate quality check of metadata.
#' This function helps to control if some experimental parameters given as metadata (numeric one) in input
#' explain much variability than expected in the data or if their effect could be confused with biological one.
#' This function correlates quantitative variable describing technical aspect for each sample with
#' their coordinate on the PCA axis.
#' \itemize{
#' \item{Technical parameters from sample preparation as the day of the RNAseq library preparation}
#' \item{Statistics results after the bioinformatics workflow as the percent of sequences with primers or % of rrna in the library}
#'  }
#' @param object An object of class \link{MultiAssayExperiment-class}
#' @param data The name of the omic data for which the PCA plot has to be drawn
#' @param axis The number of PCA axis to keep
#' @param PCA This argument indicates which type of PCA results to take: on raw ("raw") or normalized ("norm") data.
#' @param pngFile The name of the png file for saving the plot.
#' @exportMethod mvQCdata
#' @rdname mvQCdata
setMethod(f="mvQCdata",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA=c("raw","norm"),axis=3, pngFile=NULL){

            resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <- object@metadata$colDataStruc["n_dFac"]
            n_qcFac <- dim(object[[data]]@colData[,-c(1:2)])[2]
            var <- names(object[[data]]@colData[,-c(1:2)])

            corA=list()
            VarAxis = list()
            for(i in 1:axis){
            corA[[i]] <- cor(resPCA$ind$coord[,i],
                             as.data.frame(object[[data]]@colData[,-c(1:2)]),
                             method="spearman")
            VarAxis[[i]] <- paste("\n(Var=",round(resPCA$eig[i,2],1),")",sep="")
            }

            df = data.frame("QCparam"=rep(var,axis),
                            "Spearman"= unlist(corA),
                            "Axis"=rep(paste(rep("Axis",axis),1:axis,VarAxis,sep=""), each=n_qcFac))

            p <- ggplot(df,aes(x=Axis, y=abs(Spearman),fill=QCparam))+
              geom_bar(stat="identity",position=ggplot2::position_dodge(),width=0.7)+ggplot2::ylim(0,1)+
              ggplot2::labs(x = "Axis number", y="Cor(Coord_dFactor_PCA,QCparam)")

            print(p)

            if(! is.null(pngFile)){
              ggplot2::ggsave(filename = pngFile, plot = p)
            }

          })



########################################## TRANSFORM DATA #################

#### METHOD to transform data

#' @title TransformData
#'
#' @param An object of class \link{SummarizedExperiment}
#'
#' @return An object of class \link{SummarizedExperiment}
#'
#' @exportMethod TransformData
#'
#' @examples
setMethod(f= "TransformData",
          signature = "SummarizedExperiment",
          definition <- function(object, transform_method = "log2"){

            objectTransform <- object

            assayTransform  <- SummarizedExperiment::assay(objectTransform)

            switch(transform_method,
                   "log2" = {
                     SummarizedExperiment::assay(objectTransform) <- log2(assayTransform)
                     objectTransform@metadata[["transform_method"]] <- transform_method
                     },
                   "none"= {
                     SummarizedExperiment::assay(objectTransform) <- assayTransform
                     objectTransform@metadata[["transform_method"]] <- transform_method
                   }
            )
            return(objectTransform)
          })



########################################## FILTER DATA #################

#### METHOD to filter data

# Cette method est propre au RNASEQ => Est-ce que c'est vraiment ce que l'on souhaite ?
# Plutot qu'une fonction interface pour tous les omics ?
# Pourquoi ne pas avoir utilisée directement la fonction de edgeR ?

#' @title FilterLowAbundance
#' @description This function aims at removing genes/transcript from the count data matrix of an omic of type "RNAseq".
#' by applying filtering criterion described in reference.
#' By default, gene/transcript with 0 count are removed from the data. The function then
#' computes the count per million or read (CPM) for each gene in each sample and gives by
#' genes the number of sample(s) which are over the CPM_cutoff (NbOfsample_over_cpm).
#' Then Two filtering strategies are proposed:
#' \itemize{
#' \item{NbConditions: }{keep gene if the NbOfsample_over_cpm >= NbConditions}
#' \item{NbReplicates: }{keep gene if the NbOfsample_over_cpm >= min(NbReplicat)}
#' \item{filterByExpr:} {the default filtering method implemented in the edgeR filterByExpr() function.}
#' }
#' @param object An object of class \link{SummarizedExperiment}
#' @param Filter_Strategy The filtering strategy ("NbConditions" or "NbReplicates")
#' @param CPM_Cutoff The CPM cutoff.
#' @return An object of class \link{SummarizedExperiment}
#' @details
#' Filtered dataset is stored in the ExperimentList slot of the \link{SummarizedExperiment} object
#' as a List named (DataName.filtred).
#' List of filtered features are stored as a named list ("FilteredFeatures") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{SummarizedExperiment} object.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod FilterLowAbundance
#' @seealso edgeR::filterByExpr
#' @examples
setMethod(f= "FilterLowAbundance",
          signature = "SummarizedExperiment",
          definition <- function(object, Filter_Strategy = "NbConditions", CPM_Cutoff = 5){

            objectFilt <- object

            assayFilt  <- MultiAssayExperiment::assay(objectFilt)

            ## nbr of genes with 0 count
            genes_flt0  <- objectFilt[rowSums(assayFilt) <= 0, ]@NAMES

            ## remove 0 count
            objectFilt  <- objectFilt[rowSums(assayFilt)  > 0, ]
            assayFilt   <- MultiAssayExperiment::assay(objectFilt)

            ## filter cpm
            NbReplicate  <- table(object@metadata$Groups$groups)
            NbConditions <- length(unique(object@metadata$Groups$groups))

            switch(Filter_Strategy,
                   "NbConditions" = { keep <- rowSums(edgeR::cpm(assayFilt) >= CPM_Cutoff) >=  NbConditions },
                   "NbReplicates" = { keep <- rowSums(edgeR::cpm(assayFilt) >= CPM_Cutoff) >=  min(NbReplicate) },
                   "filterByExpr" = { dge  <- edgeR::DGEList(counts = assayFilt, genes = rownames(assayFilt))
                                      #keep <- filterByExpr(dge, GLM_Model)
                                      keep <- edgeR::filterByExpr(dge)
                                      }
                   )

            ## nbr of genes filtered
            genes_flt1  <- objectFilt[!keep]@NAMES

            objectFilt@metadata[["FilteredFeatures"]] <-  c(genes_flt0, genes_flt1)

            object <- objectFilt[keep]

            object@metadata$FilteredOptions <- list()

            object@metadata$FilteringOptions[["Filter_Strategy"]] <- Filter_Strategy
            object@metadata$FilteringOptions[["CPM_Cutoff"]] <- CPM_Cutoff

            return(object)

          })


######### NORMALIZATION #################

#### METHOD to normalize data

# Function non generique pour les autres data

#' @title RunNormalization
#' @description This function applied a normalization method on an omic data sets stored in an object of
#' class \link{SummarizedExperiment}.
#' \itemize{
#' \item{For RNAseq data:}{the TMM function of edgeR is proposed by default, see the ref}
#' \item{For Proteomic data:}{}
#' \item{For Metabolomic data:}{}
#' }
#' @param object An object of class \link{SummarizedExperiment}
#' @param data The name of the data set for which the normalization has to be performed.
#' @param NormMethod Normalization method
#' @return An object of class \link{SummarizedExperiment}
#' The applied normalization method and computed scaling factors (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "methode" and "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{SummarizedExperiment} object.
#' @exportMethod RunNormalization
#' @seealso TMM.Normalization
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @examples
#'
setMethod(f="RunNormalization",
          signature="SummarizedExperiment",
          definition <- function(object, NormMethod){

            coefNorm  = switch(NormMethod,
                               "TMM"=TMM.Normalization(SummarizedExperiment::assay(object), object@metadata$Groups$groups)
            )
            object@metadata[["Normalization"]] <- list(methode = NormMethod, coefNorm = coefNorm)
            return(object)
          })



################################### DIFF-ANALYSIS #############################


###### Statistical METHOD

## METHOD to perform differential analysis

#' @title RunDiffAnalysis
#' @description This is an interface method which run a differential analysis method on
#' omic datasets stored in an object of class \link{SummarizedExperiment}.
#' According to the type of omic and to a list of contrasts,
#' a differential analysis method is applied to each contrasts (or hypothesis).
#' Three methods are available according to the type of object:
#' \itemize{
#' \item{For RNAseq data: }{the \code{glmFit} function of the \code{edgeR} package}
#' \item{For proteomic and metabolomic data: }{the \code{lmFit} function of the \code{limma} package}
#' }
#' Parameters used for RNAseq are those recommended in DiCoExpress workflow (see the paper in reference)
#' @return
#' All the results are stored as a named list \code{DiffExpAnal} in the metadata slot of a
#' given \code{SummarizedExperiment} object.
#' Objects are:
#' \itemize{
#' \item{contrasts: }{The selected contrasts for which the differential analysis has been conducted}
#' \item{method: }{The method used for the differential analysis}
#' \item{Adj.pvalue.method: The method applied for the pvalue adjustment}
#' \item{Adj.pvalue.cutoff: The threshold applied for the pvalue adjustment}
#' \item{FDR: }{The false discovery rate given in input}
#' \item{RawDEFres: }{a list giving for each contrast the raw results of the differential analysis method}
#' \item{DEF: }{a list giving for each contrast a data.frame of non filtered differential expressed features}}
#' \item{TopDEF: }{a list giving for each contrast a data.frame of differential expressed features by Adj.pvalue.cutoff}
#' \item{mergeDEF: }{A data frame indicating for each features in row, if it is DE in a given contrasts in column}
#' }
#' @param object an object of class [\code\link{SummarizedExperiment}]
#' @param design an object of class [\code{\link{ExpDesign-class}]
#' @param DiffAnalysisMethod A character vector giving the name of the differential analysis method
#' to run. Either "edgeRglmfit", "limmalmFit", ...
#' @param contrastList The list of contrast to test
#' @param Adj.pvalue.method The method choosen to adjust pvalue.
#' @param Adj.pvalue.cutoff The adjusted pvalue cut-off
#' @param clustermq A boolean indicating whether the constrasts have to be computed in local or in a distant machine
#' @param filter_only A boolean indicating whether only filter on DE results have to be applied (\code{filter_only=TRUE}). FALSE by default.
#' @return An object of class [\code\link{SummarizedExperiment}]
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod RunDiffAnalysis
#' @examples
#'
#'
setMethod(f="RunDiffAnalysis",
          signature="SummarizedExperiment",
          definition <- function(object, design, Adj.pvalue.method="BH",
                                 contrastList, DiffAnalysisMethod, clustermq=FALSE){



            contrastName <- NULL
            Contrasts.Sel <- dplyr::filter(design@Contrasts.Sel, contrastName %in% contrastList)

            object@metadata$DiffExpAnal <- list()
            object@metadata$DiffExpAnal[["contrasts"]] <- Contrasts.Sel
            object@metadata$DiffExpAnal[["method"]]    <- DiffAnalysisMethod
            object@metadata$DiffExpAnal[["Adj.pvalue.method"]]  <- Adj.pvalue.method
            #object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]  <- Adj.pvalue.cutoff

            # move in ExpDesign Constructor
            model_matrix <- model.matrix(as.formula(design@Model.formula), data=as.data.frame(design@List.Factors))


            ListRes <- switch(DiffAnalysisMethod,
                                        "edgeRglmfit"=edgeR.AnaDiff(count_matrix    = SummarizedExperiment::assay(object),
                                                                    model_matrix    = model_matrix,
                                                                    group           = object@metadata$Normalization$coefNorm$group,
                                                                    lib.size        = object@metadata$Normalization$coefNorm$lib.size,
                                                                    norm.factors    = object@metadata$Normalization$coefNorm$norm.factors,
                                                                    Contrasts.Sel   = object@metadata$DiffExpAnal[["contrasts"]],
                                                                    Contrasts.Coeff = design@Contrasts.Coeff,
                                                                    FDR             = 1,
                                                                    clustermq=clustermq),
                                          "limmalmFit"=limma.AnaDiff(count_matrix      = SummarizedExperiment::assay(object),
                                                                     model_matrix      = model_matrix,
                                                                     Contrasts.Sel     = object@metadata$DiffExpAnal[["contrasts"]],
                                                                     Contrasts.Coeff   = design@Contrasts.Coeff,
                                                                     Adj.pvalue.cutoff = 1,
                                                                     Adj.pvalue.method = Adj.pvalue.method,
                                                                     clustermq=clustermq))


            ### RawDEFres: Raw results from the given diff method
            object@metadata$DiffExpAnal[["RawDEFres"]] <- ListRes[["RawDEFres"]]
            #names(ListRes[["TopDEF"]]) <- names( ListRes[["RawDEFres"]])
            object@metadata$DiffExpAnal[["DEF"]] <- ListRes[["TopDEF"]]

            return(object)
          })

# limma
# Warning quand pas de F DE
# Recuperer les messages d'erreurs de limma ou

## METHOD to filter differential analysis

#' Title
#'
#' @param SummarizedExperiment
#'
#' @return
#' @exportMethod FilterDiffAnalysis
#'
#' @examples
#'
setMethod(f="FilterDiffAnalysis",
          signature="SummarizedExperiment",
          definition <- function(object, Adj.pvalue.cutoff = 0.05){

            if(is.null(object@metadata$DiffExpAnal)){
              stop("can't filter the DiffExpAnal object because it doesn't exist")
            }

            object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]  <- Adj.pvalue.cutoff

            ## TopDEF: Top differential expressed features
            DEF_filtred <- lapply(1:length(object@metadata$DiffExpAnal[["DEF"]]),function(x){
              res <- object@metadata$DiffExpAnal[["DEF"]][[x]]
              res <- res[res$Adj.pvalue <= Adj.pvalue.cutoff,]
              return(res)
            })
            names(DEF_filtred) <- names(object@metadata$DiffExpAnal[["RawDEFres"]])
            object@metadata$DiffExpAnal[["TopDEF"]] <- DEF_filtred

            ## stats
            stats_list <- lapply(1:length(object@metadata$DiffExpAnal[["TopDEF"]]), function(x){
              gN = dim(object@metadata$DiffExpAnal[["DEF"]][[x]])[1]
              gDE =  dim(object@metadata$DiffExpAnal[["TopDEF"]][[x]])[1]
              pgDE =   round((gDE/gN)*100,0)
              gDEup =  dim(dplyr::filter(object@metadata$DiffExpAnal[["TopDEF"]][[x]],logFC > 0))[1]
              pgDEup =  round((gDEup/gDE)*100,0)
              gDEdown =  dim(dplyr::filter(object@metadata$DiffExpAnal[["TopDEF"]][[x]],logFC < 0))[1]
              pgDEdown =  round((gDEdown/gDE)*100,0)
              list(
              "gN" = gN,
              "gDE" =  gDE,
              "pgDE" =  pgDE,
              "gDEup" =  gDEup,
              "pgDEup" =  pgDEup,
              "gDEdown" =  gDEdown,
              "pgDEdown" =  pgDEdown)
            })
            names(stats_list) <- names(object@metadata$DiffExpAnal[["RawDEFres"]])
            object@metadata$DiffExpAnal[["stats"]] <- stats_list


            ## merge results in bin matrix
            DEF_list <- lapply(names(object@metadata$DiffExpAnal[["TopDEF"]]), function(x){
              print(x)
              res <- object@metadata$DiffExpAnal[["TopDEF"]][[x]]
              tmp <- data.frame(DEF = rownames(res), bin = rep(1,length(rownames(res))))
              colnames(tmp) <- c("DEF", dplyr::filter(object@metadata$DiffExpAnal$contrasts, contrastName == x)$tag)
              print(dplyr::filter(object@metadata$DiffExpAnal$contrasts, contrastName == x)$tag)
              return(tmp)
            })
            names(DEF_list) <- names(object@metadata$DiffExpAnal[["TopDEF"]])

            object@metadata$DiffExpAnal[["mergeDEF"]] <- DEF_list %>% purrr::reduce(dplyr::full_join, by="DEF") %>%
              dplyr::mutate_at(.vars = 2:(length(DEF_list)+1), 
                               .funs = function(x){
                                  dplyr::if_else(is.na(x), 0, 1)}) %>% 
              data.table::data.table()

            return(object)
          })


###### Graphical METHOD

## Method to plot results of a differential analysis

#' @title DiffAnal.plot
#' @description
#' This is an interface method which draw a MAplot from the results of a differential analysis
#' performed on omic datasets stored in an object of class \link{SummarizedExperiment}
#' @param object An object of class \link{SummarizedExperiment}
#' @param data The name of the omic data for which the MAplot has to be drawn
#' @param hypothesis The hypothesis for which the MAplot has to be drawn
#' @return plot
#' @exportMethod DiffAnal.plot
#' @export
#'
#' @examples
setMethod(f="DiffAnal.plot",
          signature="SummarizedExperiment",
          definition <- function(object, hypothesis,Adj.pvalue.cutoff = 0.05){

            plots <- list()

            res      <- object@metadata$DiffExpAnal[["RawDEFres"]][[hypothesis]]
            resTable <- object@metadata$DiffExpAnal[["DEF"]][[hypothesis]]

            plots[["MA.plot"]]     <- MA.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff)
            plots[["Pvalue.hist"]] <- pvalue.plot(data =resTable)

            return(plots)
          })




################################### CO-EXPRESSION #############################


#' @title runCoExpression
#' @description This is an interface method which performed co-expression/co-abundance analysis
#' of omic-data.
#' @details For instance, only the coseq function of the package coseq is proposed.
#' For RNAseq data, parameters used are those recommended in DiCoExpress workflow (see the reference).
#' This parameters are: \code{model="normal"}, \code{transformation="arcsin"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="TMM"}, \code{meanFilterCutoff = 50}
#' For proteomic or metabolomic, data are scaled by protein or metabolite to groups them by expression
#' profiles rather than by expression intensity.
#' After data scaling, recommended parameters (from \code{coseq} developers) for co-expression analysis are:
#' \code{model="normal"}, \code{transformation="none"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="none",  \code{meanFilterCutoff = NULL}
#'
#' @return
#' All the results are stored as a named list \code{CoExpAnal} in the metadata slot of a
#' given \code{SummarizedExperiment} object. Objects are:
#' The runCoExpression method return several results, for \link{coseq} method, objects are:
#' \itemize{
#' \item{\code{model:} }{see model params description}
#' \item{\code{transformation:} }{see transformation params description}
#' \item{\code{normFactors:} }{see normFactors params description}
#' \item{\code{meanFilterCutoff:} }{set to 50 for RNA and to NULL for others}
#' \item{\code{gene.list.names:} }{see nameList in Arguments description}
#' \item{\code{merge.type:} }{see merge params description}
#' \item{\code{coseqResults:} }{the raw results of \code{coseq}}
#' \item{\code{clusters:} }{a List of clusters}
#' \item{\code{cluster.nb:} }{The number of cluster}
#' \item{\code{plots:} }{The plots of \code{coseq} results}
#' }
#' @param object An object of class \link{SummarizedExperiment}
#' @param geneList A list of genes
#' @param K Number of clusters (a single value or a vector of values)
#' @param replicates The number of iteration for each K.
#' @param model Type of mixture model to use \code{"Poisson"} or \code{"normal"}. By default, it is the normal.
#' @param GaussianModel Type of \code{GaussianModel} to be used for the Normal mixture model only. This parameters
#' is set to \code{"Gaussian_pk_Lk_Ck"} by default and doesn't have to be changed except if an error message proposed
#' to try another model like \code{"Gaussian_pk_Lk_Bk"}.
#' @param transformation The transformation type to be used. By default, it is the "arcsin" one.
#' @param normFactors The type of estimator to be used to normalize for differences in library size.
#' By default, it is the "TMM" one.
#' @param merge \code{"union"} or \code{"intersection"}
#' @return An S4 object of class \link{SummarizedExperiment}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runCoExpression
#' @seealso \code{\link{coseq::coseq}}
setMethod(f="runCoExpression",
          signature="SummarizedExperiment",
          definition <- function(object, geneList, K=2:20, replicates=5, nameList, merge="union",
                                 model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation, normFactors, clustermq=FALSE){

            
            CoExpAnal <- list()
            
            CoExpAnal[["tools"]]            <- "CoSeq"
            CoExpAnal[["gene.list.names"]]  <- nameList
            CoExpAnal[["merge.type"]]       <- merge
            
            counts = SummarizedExperiment::assay(object)[geneList,]
  
            
            # set default parameters based on data type
            param.list <- list()
            switch (object@metadata$omicType,
                    
              "RNAseq" = {
                param.list[["model"]]            <- model
                param.list[["transformation"]]   <- "arcsin"
                param.list[["normFactors"]]      <- "TMM"
                param.list[["meanFilterCutoff"]] <- 50
                param.list[["GaussianModel"]]    <- GaussianModel

              },
              "proteomics" = {
                # Print the selected GaussianModel
                print(paste("Use ",GaussianModel,sep=""))
                print("Scale each protein (center=TRUE,scale = TRUE)")
                CoExpAnal[["transformation.prot"]] <- "scaleProt"
                counts[] <- t(apply(counts,1,function(x){ scale(x, center=TRUE,scale = TRUE) }))
                
                # param
                param.list[["model"]]            <- model
                param.list[["transformation"]]   <- "none"
                param.list[["normFactors"]]      <- "none"
                param.list[["meanFilterCutoff"]] <- NULL
                param.list[["GaussianModel"]]    <- GaussianModel
              },
              "metabolomics" = {
                # Print the selected GaussianModel
                print(paste("Use ",GaussianModel,sep=""))
                print("Scale each metabolite (center=TRUE,scale = TRUE)")
                CoExpAnal[["transformation.metabo"]] <- "scaleMetabo"
                counts[] <- t(apply(counts,1,function(x){ scale(x, center=TRUE,scale = TRUE) }))
                
                # param
                param.list[["model"]]            <- model
                param.list[["transformation"]]   <- "none"
                param.list[["normFactors"]]      <- "none"
                param.list[["meanFilterCutoff"]] <- NULL
                param.list[["GaussianModel"]]    <- GaussianModel
              }
            )

            CoExpAnal[["param"]] <- param.list
          
            # run coseq : on local machine or remote cluster
            coseq.res.list <- list()
            coseq.res.list <- try_rflomics(
                    runCoseq(counts, K=K, replicates=replicates, param.list=param.list , clustermq=clustermq))
            

            if( ! is.null(coseq.res.list$value) ){
              
              
                    # ICL plot
                    ICL.vec <- sapply(1:(length(K)*replicates), function(x){ coseq::ICL(coseq.res.list[["value"]][[x]]) }) 
                    ICL.tab <- data.frame(K=stringr::str_replace(names(ICL.vec), "K=", ""), ICL=ICL.vec)
                    ICL.p   <- ggplot(data = ICL.tab) + geom_boxplot(aes(x=K, y=ICL))
                    
                    #CoExpAnal[["plots"]][["ICL"]] <- ICL.p
                    
                    # logLike plot
                    logLike.vec <- sapply(1:(length(K)*replicates), function(x){ coseq::likelihood(coseq.res.list[["value"]][[x]]) }) 
                    logLike.tab <- data.frame(K=stringr::str_replace(names(logLike.vec), "K=", ""), logLike=logLike.vec)
                    logLike.p   <- ggplot(data = logLike.tab) + geom_boxplot(aes(x=K, y=logLike))
                    
                    #CoExpAnal[["plots"]][["logLike"]] <- logLike.p
                    
                    # results with 
                    coseq.res <- coseq.res.list[["value"]][[which.min(ICL.vec)]]
              
                    CoExpAnal[["results"]]      <- TRUE
                    CoExpAnal[["warning"]]      <- coseq.res.list$warning
                    CoExpAnal[["coseqResults"]] <- coseq.res
                    #CoExpAnal[["coseqResults"]] <- coseq.res.list$value
                    #coseq.res <- coseq.res.list$value
      
                    # list of genes per cluster
                    clusters <- lapply(1:length(table(coseq::clusters(coseq.res))), function(i){
                      names(coseq::clusters(coseq.res)[coseq::clusters(coseq.res) == i])
                      })
                    CoExpAnal[["clusters"]] <- clusters
                    names(CoExpAnal[["clusters"]]) <- paste("cluster", 1:length(table(coseq::clusters(coseq.res))), sep = ".")
      
                    # nbr of cluster
                    nb_cluster <- coseq.res@metadata$nbCluster[min(coseq.res@metadata$ICL) == coseq.res@metadata$ICL]
                    CoExpAnal[["cluster.nb"]] <- nb_cluster
      
                    # plot
                    plot.coseq.res <- coseq::plot(coseq.res, conds = object@metadata$Groups$groups,
                                                  graphs = c("profiles", "boxplots", "probapost_boxplots",
                                                             "probapost_barplots", "probapost_histogram")) # , collapse_reps = "average"
                    CoExpAnal[["plots"]] <- plot.coseq.res
                    CoExpAnal[["plots"]][["ICL"]]     <- ICL.p
                    CoExpAnal[["plots"]][["logLike"]] <- logLike.p
            }
            # Réinitialisation de l'objet CoExpAnal
            else{
                    CoExpAnal[["results"]] <- FALSE
                    CoExpAnal[["warning"]] <- coseq.res.list$warning
                    CoExpAnal[["error"]]   <- coseq.res.list$error
            }


            object@metadata$CoExpAnal <- CoExpAnal
            return(object)
})


# Pour utiliser la fonction repeatable(), "seed"  pourrait être ajouté en paramètre.

################################### ANNOTATION #############################

#' @title runAnnotationEnrichment
#' @description This function performs enrichment test from functional gene annotation data. This data could be
#' GO, KEGG or other... For instance, the hypergeometric test is applied. Parameters used are those
#' recommended in DiCoExpress workflow (see the paper in reference)
#' @param object An object of class \link{SummarizedExperiment}
#' @param CoExpListNames A list of clusters names.
#' @param DiffListNames A list of differential expressed genes given as constrast names.
#' @param alpha The pvalue cut-off
#' @param probaMethod The probabilistic method to use.
#' @param annotation The gene annotation file.
#' @return An object of class \link{SummarizedExperiment}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runAnnotationEnrichment
#'
setMethod(f="runAnnotationEnrichment",
          signature="SummarizedExperiment",

          definition <- function(object, DiffListNames = NULL, CoExpListNames = NULL, annotation,
                                 alpha = 0.01, probaMethod = "hypergeometric"){

            EnrichAnal <- list()
            EnrichAnal[["list.names"]] <- c(CoExpListNames, DiffListNames)
            EnrichAnal[["alpha"]]      <- alpha
            EnrichAnal[["proba.test"]] <- probaMethod


            ## list of gene list to annotate
            geneLists <- list()
            geneLists.diff <- list()
            geneLists.diff <- lapply(DiffListNames, function(listname){

              row.names(object@metadata$DiffExpAnal[["TopDEF"]][[listname]])
            })
            names(geneLists.diff) <- DiffListNames

            geneLists.coseq <- list()
            #geneLists.coseq <- lapply(CoExpListNames, function(listname){

            geneLists.coseq <-  object@metadata[["CoExpAnal"]][["clusters"]][CoExpListNames]
            #})
            #names(geneLists.coseq) <- CoExpListNames

            geneLists <- c(geneLists.diff, geneLists.coseq)


            Results <- list()

            for(geneList in names(geneLists)){

              Results[[geneList]] <- switch(probaMethod,
                     "hypergeometric"=EnrichmentHyperG(annotation, geneLists[[geneList]], alpha = 0.01)
                     )
            }

            EnrichAnal[["results"]] <- Results

            object@metadata$EnrichAnal <- EnrichAnal

            return(object)
          })



