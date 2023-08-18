
################################## EXPERIMENTAL DESIGN SET UP #################


### ExpDesign CLASS Constructor

#' @title Constructor for the class \link{ExpDesign-class}
#' @description This method initializes an object of class \link{ExpDesign-class}.
#' @param ExpDesign a data.frame. Row names give the name of each sample which has to be constructed
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Column names give the name of
#' an experimental factor which is a vector of character storing the factor modality for each sample.
#' @param refList A list of string giving the reference modality for each factor.
#' @param typeList A vector of string indicating the type of each experimental factor. Two types of effects
#' are required ("Bio" or "batch"). A third one ("meta") is allowed but is not necessary.
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
#' @noRd
#' @export
#' @importFrom stats relevel
#' @importFrom methods new
ExpDesign.constructor <- function(ExpDesign, refList, typeList){
  
  # check ExpDesign dimension
  if(dim(ExpDesign)[1] == 0 || dim(ExpDesign)[2] == 0){
    stop("Error: ExpDesign matrix is empty!")
  }
  
  # check refList length
  if(length(refList) != length(names(ExpDesign))){
    stop("Error: refList length is different from the dimension of ExpDesign matrix!")
  }
  
  # check typeList length
  if(length(typeList) != length(names(ExpDesign))){
    stop("Error: typeList length is different from the dimension of ExpDesign matrix!")
  }
  
  # Create the List.Factors list with the choosen level of reference for each factor
  names(refList)  <- names(ExpDesign)
  names(typeList) <- names(ExpDesign)
  
  # for(i in c(names(typeList[typeList == "batch"]), names(typeList[typeList == "Bio"]))){
  #   ExpDesign      <- dplyr::arrange(ExpDesign, get(i))
  # }
  
  
  dF.List <- list()
  for(i in names(ExpDesign)){
    ExpDesign[[i]] <- relevel(as.factor(ExpDesign[[i]]), ref=refList[i])
    dF.List[[i]]   <- ExpDesign[[i]]
  }
  # dF.List <- lapply(1:dim(ExpDesign)[2], function(i){
  #
  #   relevel(as.factor(ExpDesign[[i]]), ref=refList[i])
  #
  # })
  names(dF.List) <- names(ExpDesign)
  
  # Create the groups data.frame
  # groups <- tidyr::unite(as.data.frame(ExpDesign[typeList == "Bio"]), col="groups", sep="_", remove = TRUE) %>%
  #           dplyr::mutate(samples = rownames(.))
  
  groups <- ExpDesign %>% dplyr::mutate(samples = rownames(.)) %>%
    tidyr::unite(names(typeList[typeList == "Bio"]), col="groups", sep="_", remove = FALSE)
  
  groups$samples <- factor(groups$samples, levels = unique(groups$samples))
  groups$groups  <- factor(groups$groups,  levels = unique(groups$groups))
  
  Design = new(Class           = "ExpDesign",
               ExpDesign       = as.data.frame(ExpDesign),
               List.Factors    = dF.List,
               Factors.Type    = typeList,
               Groups          = groups,
               Model.formula   = vector(),
               Contrasts.List  = list(),
               Contrasts.Sel   = data.frame(),
               Contrasts.Coeff = data.frame())
  
  return(Design)
}

#
# Error quand plus de 3 facteurs bio et plus de 1 facteur batch
# TEST(design_nbbio_3)
# TEST(design_nbbatch_1)


###### METHOD to check the completness of the ExpDesign


#' @title CheckExpDesignCompleteness
#' @description This method checks some experimental design characteristics.
#'  A complete design and at least one biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{MultiAssayExperiment-class}
#' @param sampleList list of samples to check.
#' @return a named list of two objects
#' \itemize{
#'  \item{"count:"}{ a data.frame with the number of each possible combinations of levels for all factors.}
#'  \item{"plot:"}{ plot of count data.frame.}
#'  \item{"warning:"}{ warning message if design is not balanced}
#'  \item{"error:"}{error message if:}
#'  \itemize{
#'   \item{number of biological factor is not between 1 and 3}
#'   \item{number of batch factor is not between 1 and 2}
#'   \item{design is uncomplete}
#'   }
#'  }
#'  
#' @exportMethod CheckExpDesignCompleteness
#' @examples
#' @noRd

methods::setMethod(f         = "CheckExpDesignCompleteness",
                   signature = "MultiAssayExperiment",
                   definition <- function(object, sampleList=NULL){
                     
                     Design <- object@metadata$design
                     
                     # output list
                     output <- list()
                     output[["error"]] <- NULL
                     output[["warning"]] <- NULL
                     
                     
                     # check presence of bio factors
                     if (! table(Design@Factors.Type)["Bio"] %in% 1:3){
                       
                       output[["error"]] <- "ERROR: no bio factor! or nbr of bio factors exeed 3!"
                       
                     }
                     if (table(Design@Factors.Type)["batch"] == 0){
                       
                       output[["error"]] <- "ERROR: no replicate!"
                     }
                     
                     
                     # count occurence of bio conditions
                     if(is.null(sampleList)){
                       # tmp <- sampleMap(object) %>% data.frame()
                       # sampleList <- lapply(unique(tmp$assay), function(dataset){filter(tmp, assay == dataset)$primary }) %>%
                       #   purrr::reduce(dplyr::union)
                       
                       #sampleList <- sampleMap(object)$primary
                       sampleList.tmp <- dplyr::group_by(data.frame(MultiAssayExperiment::sampleMap(object)), primary) %>% dplyr::count() %>% 
                         dplyr::ungroup() %>% dplyr::mutate(max=max(n)) %>% dplyr::filter(n==max)
                       sampleList <- sampleList.tmp$primary
                     }
                     
                     
                     # Only work with bio and batch factors for the rest of the function
                     namFact <- names(Design@Factors.Type)[Design@Factors.Type %in% c("Bio", "batch")]
                     expDesign_mod <- Design@ExpDesign %>% dplyr::select(tidyselect::any_of(namFact))
                     
                     dF.List <- lapply(1:ncol(expDesign_mod), function(i){
                       factor(expDesign_mod[[i]], levels = unique(expDesign_mod[[i]]))
                     })
                     names(dF.List) <- names(expDesign_mod)
                     
                     ExpDesign <- dplyr::filter(expDesign_mod, rownames(expDesign_mod) %in% sampleList)
                     
                     bio.fact <- names(dF.List[Design@Factors.Type == "Bio"])
                     tmp <- ExpDesign %>% dplyr::mutate(samples=row.names(.))
                     
                     group_count <- as.data.frame(dF.List) %>% 
                       table() %>% 
                       as.data.frame() %>% 
                       dplyr::full_join(tmp, by=names(dF.List)) %>% 
                       dplyr::mutate_at(.vars = "samples", .funs = function(x) dplyr::if_else(is.na(x), 0, 1)) %>%
                       dplyr::group_by_at((bio.fact)) %>% 
                       dplyr::summarise(Count=sum(samples), .groups = "keep")
                     
                     # check presence of relicat / batch
                     # check if design is complete
                     # check if design is balanced
                     # check nbr of replicats
                     if(min(group_count$Count) == 0){
                       message <- "ERROR: The experimental design is not complete."
                       output[["error"]] <- message
                     }
                     else if(min(group_count$Count) == 1){
                       message <- "ERROR: You need at least 2 biological replicates."
                       output[["error"]] <- message
                     }
                     else if(length(unique(group_count$Count)) != 1){
                       message <- "WARNING: The experimental design is complete but not balanced."
                       output[["warning"]] <- message
                     }
                     else{
                       message <- "The experimental design is complete and balanced."
                     }
                     
                     ### plot
                     output[["plot"]] <- RFLOMICS::plotExperimentalDesign(group_count, message=message)
                     output[["counts"]] <- group_count
                     return(output)
                   })

# print output
# warining -> warning
# false -> stop



#' @title Datasets overview plot
#' @description This function plot overview of loaded datasets aligned per sample 
#' (n=number of entities (genes/metabolites/proteins); k=number of samples)
#' @param An object of class \link{MultiAssayExperiment-class}
#' @exportMethod Datasets_overview_plot
#' @return plot

methods::setMethod(f         = "Datasets_overview_plot",
                   signature = "MultiAssayExperiment",
                   definition <- function(object){
                     
                     if (class(object) != "MultiAssayExperiment") stop("ERROR: object is not MultiAssayExperiment class.")
                     if (length(object@ExperimentList) == 0) stop("ERROR: object@ExperimentList is NULL")
                     
                     nb_entities <- lapply(object@ExperimentList, function(SE){ dim(SE)[1] }) %>% unlist()
                     
                     data <- data.frame(nb_entities = nb_entities, assay = names(nb_entities)) %>%
                       dplyr::full_join(data.frame(MultiAssayExperiment::sampleMap(object)), by="assay") %>%
                       dplyr::mutate(y.axis = paste0(assay, "\n", "n=", nb_entities)) %>% dplyr::arrange(primary)
                     
                     p <- ggplot2::ggplot(data, ggplot2::aes(x=primary, y=y.axis)) +
                       ggplot2::geom_tile(aes(fill = y.axis), colour = "grey50") +
                       ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                             panel.background = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position="none",
                             axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
                       ggplot2::xlab(paste0("Samples (k=", length(unique(MultiAssayExperiment::sampleMap(object)$primary)), ")")) +
                       ggplot2::ylab("")
                     
                     print(p)
                     
                   })


###### METHOD which generate the contrasts expression


#' @title getExpressionContrast
#' @description This function allows, from a model formulae, to give the expression contrast data frames.
#' Three types of contrasts are expressed:
#' \itemize{
#' \item{simple}
#' \item{pairwise comparison}
#' \item{averaged expression}
#' }
#' @param model.formula a model formula (characters or formula)
#' @return An object of class [\code{\link{MultiAssayExperiment-class}}]
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
#' # Set the model formula
#' Design.formulae <- GetModelFormulae(Factors.Name = Design.Factors.Name,Factors.Type=Design.Factors.Type)
#' Design.formulae[[1]]
#'
#' # Obtained the Expression of Contrasts
#' Design.obj <- getExpressionContrast(object = Design.obj, model.formula = names(Design.formulae[1]))
#'
#' @author Christine Paysant-Le Roux
#' @noRd
methods::setMethod(f          = "getExpressionContrast",
                   signature  = "MultiAssayExperiment",
                   definition <- function(object, model.formula){
                     
                     Design <- object@metadata$design
                     
                     # model formula
                     if (class(model.formula) == "formula") model.formula <- as.character(model.formula)
                     # modelFormula <- formula(model.formula) # deprecated? 
                     modelFormula <- formula(paste(model.formula, collapse = " "))
                     
                     #Design@Model.formula <- formula(model.formula)
                     Design@Model.formula <- model.formula
                     
                     # bio factor list in formulat
                     labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
                     
                     FactorBioInDesign <- intersect(names(Design@Factors.Type[Design@Factors.Type == "Bio"]), labelsIntoDesign)
                     
                     #BioFactors <- Design@List.Factors[FactorBioInDesign]
                     
                     treatmentFactorsList <- lapply(FactorBioInDesign, function(x){(paste(x, levels(Design@List.Factors[[x]]), sep=""))})
                     names(treatmentFactorsList) <- FactorBioInDesign
                     
                     interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
                     
                     listOfContrastsDF <- list()
                     # define all simple contrasts pairwise comparisons
                     
                     allSimpleContrast_df <- defineAllSimpleContrasts(treatmentFactorsList)
                     # if 1 factor or more than 1 + interaction
                     if(length(treatmentFactorsList) == 1 || !isFALSE(interactionPresent)){
                       
                       listOfContrastsDF[["simple"]] <- allSimpleContrast_df
                     }
                     
                     # define all simples contrast means
                     # exists("allSimpleContrast_df", inherits = FALSE)
                     if(length(treatmentFactorsList) != 1){
                       allAveragedContrasts_df <- define_averaged_contrasts (allSimpleContrast_df)
                       listOfContrastsDF[["averaged"]] <- allAveragedContrasts_df
                     }
                     
                     # define all interaction contrasts
                     if(length(treatmentFactorsList) != 1){
                       if(interactionPresent){
                         labelsIntoDesign            <- attr(terms.formula(modelFormula),"term.labels")
                         labelOrder                  <- attr(terms.formula(modelFormula), "order")
                         twoWayInteractionInDesign   <- labelsIntoDesign[which(labelOrder == 2)]
                         groupInteractionToKeep      <- gsub(":", " vs ", twoWayInteractionInDesign)
                         allInteractionsContrasts_df <- RFLOMICS::defineAllInteractionContrasts(treatmentFactorsList, groupInteractionToKeep)
                         
                         listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
                       }
                       #allInteractionsContrasts_df <- defineAllInteractionContrasts(treatmentFactorsList)
                       #listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
                     }
                     # choose the contrasts and rbind data frames of contrasts
                     #selectedContrasts <- returnSelectedContrasts(listOfContrastsDF)
                     
                     # replace interactive selection of contrasts by return all contrasts -> shiny
                     Design@Contrasts.List  <- listOfContrastsDF
                     Design@Contrasts.Coeff <- data.frame()
                     Design@Contrasts.Sel   <- data.frame()
                     
                     object@metadata$design <- Design
                     
                     return(object)
                   })


###### METHOD to obtain the Matrix of contrast with their names and coefficients

#
#

#' @title getContrastMatrix
#' @description Defines contrast matrix or contrast list with contrast name and contrast coefficients
#' @param An object of class \link{MultiAssayExperiment-class}
#' @param contrastList A vector of character of contrast
#' @return An object of class \link{MultiAssayExperiment-class}
#' @seealso getExpressionContrast
#' @exportMethod getContrastMatrix
#' @importFrom stats formula terms.formula
#' @noRd
#' @author Christine Paysant-Le Roux
methods::setMethod(f          = "getContrastMatrix",
                   signature  = "MultiAssayExperiment",
                   definition <- function(object, contrastList){
                     
                     Design <- object@metadata$design
                     
                     contrast <- contrastName <- type <- groupComparison <- NULL
                     
                     contrast.sel.list <- list()
                     contrast.sel.list <- lapply(names(Design@Contrasts.List), function(contrastType) {
                       
                       tmp <- Design@Contrasts.List[[contrastType]] %>% dplyr::filter(contrast %in% contrastList) %>%
                         dplyr::select(contrast, contrastName, type, groupComparison)
                       return(tmp)
                     })
                     Design@Contrasts.Sel <- contrast.sel.list %>% purrr::reduce(rbind) %>% dplyr::mutate(tag = paste("H", 1:dim(.)[1], sep=""))
                     
                     
                     sampleData <-  Design@ExpDesign
                     selectedContrasts <- Design@Contrasts.Sel$contrast
                     
                     # modelFormula <- formula(Design@Model.formula) # deprecated?
                     modelFormula <- formula(paste(Design@Model.formula, collapse = " ")) 
                     
                     # bio factor list in formula
                     labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
                     FactorBioInDesign <- intersect(names(Design@Factors.Type[Design@Factors.Type == "Bio"]), labelsIntoDesign)
                     
                     #BioFactors <- Design@List.Factors[FactorBioInDesign]
                     
                     treatmentFactorsList <- lapply(FactorBioInDesign, function(x){paste(x, unique(Design@List.Factors[[x]]), sep="")})
                     names(treatmentFactorsList) <- FactorBioInDesign
                     
                     treatmentCondenv <- new.env()
                     
                     interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
                     isThreeOrderInteraction <- any(attr(terms.formula(modelFormula),"order") == 3)
                     
                     # get model matrix
                     modelMatrix <- stats::model.matrix(modelFormula, data = Design@List.Factors %>% as.data.frame())
                     colnames(modelMatrix)[colnames(modelMatrix) == "(Intercept)"] <- "Intercept"
                     # assign treatment conditions(group) to boolean vectors according to the design model matrix
                     #treatmentCondenv <- new.env()
                     RFLOMICS::assignVectorToGroups(treatmentFactorsList    = treatmentFactorsList,
                                                    modelMatrix             = modelMatrix,
                                                    interactionPresent      = interactionPresent,
                                                    isThreeOrderInteraction = isThreeOrderInteraction,
                                                    treatmentCondenv        = treatmentCondenv)
                     # get the coefficient vector associated with each selected contrast
                     # contrast <- allSimpleContrast_df$contrast[1]
                     colnamesGLMdesign <- colnames(modelMatrix)
                     
                     #coefficientsMatrix <- sapply(selectedContrasts$contrast, function(x) returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
                     coefficientsMatrix <- sapply(selectedContrasts, function(x) RFLOMICS::returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
                     
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
                     
                     Design@Contrasts.Coeff <- contrastMatrix
                     
                     object@metadata$design <- Design
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
#' \item{data:}{matrix of omic data}
#' \item{meta:}{an optional quality check data}
#' \item{omicType:}{Type of omic data type "None", "RNAseq", "proteomics" or "metabolomics".}
#' }
#' @param Design An object of class \link{ExpDesign-class}
#' @param projectName Project name
#' @param ExpDesign a data.frame. Row names give the name of each sample which has been to be construct
#' by combining factor's modality separated by a "_" (EX: WT_treated_rep1). Column names give the name of
#' an experimental factor which is a vector of character storing the factor modality for each sample.
#' @param refList A list of string giving the reference modality for each factor.
#' @param typeList A vector of string indicating the type of each experimental factor. Two types of effect
#' are required ("Bio" or "batch")
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
#' @noRd
#'

FlomicsMultiAssay.constructor <- function(inputs, projectName, ExpDesign , refList , typeList){
  
  # consctuct ExpDesign object
  for (i in names(ExpDesign)){
    ExpDesign      <- dplyr::arrange(ExpDesign, get(i))
  }
  
  Design <- RFLOMICS::ExpDesign.constructor(ExpDesign = ExpDesign, refList = refList, typeList = typeList)
  
  
  ## creat SE object of each dataset
  SummarizedExperimentList <- list()
  listmap  <- list()
  omicList <- list()
  k <- 0
  for (dataName in names(inputs)){
    
    k <- k+1
    
    ## construct SummarizedExperiment for each data
    abundance <- inputs[[dataName]][["data"]]
    
    # check overlap between design and data
    sample.intersect <- intersect(colnames(abundance), row.names(ExpDesign))
    if(length(sample.intersect) == 0){
        stop("samples in omics data should match the names in experimental design")
    }
    
    ##### Comment : 01/03/2023 : 
    # TODO : uncomment !!!! 
    # if(length(sample.intersect) < length(row.names(ExpDesign))){
    #   
    #   message("more than half of samples don't match to experimental design")
    # }
    
    # select abundance from design table
    abundance <- dplyr::select(abundance, tidyselect::all_of(sample.intersect))
    
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
    
    # => creat colData
    if(!is.null(inputs[[dataName]][["meta"]])){
      QCmat <- inputs[[dataName]][["meta"]]
    }
    else{
      
      #sample.intersect <- intersect(colnames(matrix.filt), row.names(ExpDesign))
      
      QCmat <- data.frame(primary = sample.intersect,
                          colname = sample.intersect,
                          stringsAsFactors = FALSE)
    }
    
    omicType <- inputs[[dataName]][["omicType"]]
    Groups  <- dplyr::filter(Design@Groups, samples %in% colnames(as.matrix(matrix.filt)))
    
    SE <- SummarizedExperiment::SummarizedExperiment(assays   = S4Vectors::SimpleList(abundance = as.matrix(matrix.filt)),
                                                     colData  = QCmat,
                                                     metadata = list(omicType      = omicType, 
                                                                     Groups        = Groups, 
                                                                     rowSums.zero  = genes_flt0,
                                                                     Normalization = list(methode = "none", coefNorm = NULL, normalized = FALSE),
                                                                     transform     = list(transform_method = "none", transformed = FALSE)
                                                                     ))
    
    #SummarizedExperimentList[[dataName]] <- SE[, SE$primary %in% row.names(ExpDesign)]
    
    #### run PCA for raw count
    SummarizedExperimentList[[dataName]] <- RFLOMICS::RunPCA(SE, raw = TRUE)
    
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



################################# EXPLORATION OF BIOLOGICAL AND TECHNICAL VARIABILITY ##################################


##### Statistical METHODS for exploring biological and technical variability


#' @title RunPCA
#' @description This function performs a principal component analysis on omic data stored in an object of class [\code{\link{SummarizedExperiment-class}]
#' Results are stored in the metadata slot of the same object. If a "Normalization" slot is present in the metadata slot, then data are normalized before running the PCA according to the indicated transform method.
#' @param object An object of class \link{SummarizedExperiment-class}.
#' @param nbcp Number of components to compute. Default is 5.
#' @param raw boolean. Does the pca have to be ran on raw data or transformed and normalized data? Default is FALSE, pca is ran on transformed and normalized data.
#' @return An object of class \link{SummarizedExperiment}
#' @exportMethod RunPCA
#' @examples
#'
methods::setMethod(f          = "RunPCA",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, nbcp = 5, raw = FALSE){

                     # Check for NA/nan
                     if(RFLOMICS:::check_NA(object)){
                       message("STOP: NA or nan detected in your data")
                       return(object)
                     }
                     
                     # Compute PCA on raw data
                     # Assume assay(object) is untransformed and un-normalized
                     if (raw) {
                       
                       if (object@metadata[["transform"]][["transformed"]])
                         message("WARNING: your data are not raw (transformed)")
                       
                       if (object@metadata[["Normalization"]]$normalized) 
                         message("WARNING: your data are not raw (normalized)")
                       
                       pseudo <- SummarizedExperiment::assay(object)
                       
                       if (object@metadata$omicType == "RNAseq") {
                         pseudo <- log2(pseudo + 1) 
                         object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                       } else {
                         object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                       }
                       
                       return(object)
                     }
                     
                     # Compute PCA, transform and normalize the data if needed.
                     else{
                       
                      objectPCA <- object

                      if (!objectPCA@metadata[["transform"]][["transformed"]]) objectPCA <- RFLOMICS:::apply_transformation(objectPCA)
                      if (!objectPCA@metadata[["Normalization"]]$normalized)   objectPCA <- RFLOMICS:::apply_norm(objectPCA)
                
                      pseudo <- SummarizedExperiment::assay(objectPCA)
                       
                      if (objectPCA@metadata$omicType == "RNAseq")  pseudo <- log2(pseudo) # + 1 in apply norm
                       
                       object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                       
                       return(object)
                     }
                     
                     
                     
                     ### OLD CODE
                     # # Transformation of the data (log2, log10, etc.)
                     # if(transformData){
                     #   if(!is.null(object@metadata[["transform_method"]])){
                     #     print("PCA: transforming data")
                     #     objectPCA <- RFLOMICS::TransformData(object, transform_method = transformMethod)
                     # 
                     #     # Check for NA/nan
                     #     if(RFLOMICS::check_NA(objectPCA)){
                     #       message("STOP: NA or nan detected in your data")
                     #       return(object)
                     #     }
                     # 
                     #     pseudo <- SummarizedExperiment::assay(objectPCA)
                     #   }else{
                     #     message("PCA: asking for transforming the data but no transform method in the object. Keeping untransformed data")
                     #     pseudo <- SummarizedExperiment::assay(object)
                     #   }
                     # }else{
                     #   pseudo <- SummarizedExperiment::assay(object)
                     # } # end if transform
                     # 
                     # # if the data has undergone a transformation
                     # if(!is.null(object@metadata[["Normalization"]]$methode)){
                     #   
                     #   # RNASeq data
                     #   if(object@metadata[["Normalization"]]$methode == "TMM"  && ! is.null(object@metadata[["Normalization"]]$coefNorm)){
                     #     pseudo <- log2(scale(pseudo+1, center=FALSE,
                     #                          scale=object@metadata[["Normalization"]]$coefNorm$norm.factors*object@metadata[["Normalization"]]$coefNorm$lib.size))
                     #     object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     #     
                     #   }
                     #   # Proteo and metabo median
                     #   else if(object@metadata[["Normalization"]]$methode == "median"){
                     #     pseudo <- apply(pseudo, 2, FUN = function(sample_vect) sample_vect - median(sample_vect)) 
                     #     
                     #     object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     #   }
                     #   # Proteo and metabo totalSum
                     #   else if(object@metadata[["Normalization"]]$methode == "totalSum"){
                     #     pseudo <- apply(pseudo, 2, FUN = function(sample_vect) sample_vect/sum(sample_vect^2)) 
                     #     
                     #     object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     #   }else{ # method = "none"
                     #     pseudo <- pseudo
                     #     
                     #     object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     #   }
                     # }
                     # 
                     # # if no transformation: differentiate RNASeq from the rest
                     # else{
                     #   if(object@metadata$omicType == "RNAseq"){
                     #     pseudo <- log2(pseudo + 1)
                     #     object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     #   }else{
                     #     pseudo <- pseudo # do nothing and compute the PCA
                     #     object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     #   }
                     #  }
                     
                     
                     ### OLD CODE (older)
                     # # if the data has undergone a transformation (meta or prot data)
                     # else if(! is.null(object@metadata$transform_method)){
                     #   pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE) + 1) # pourquoi du log2 ?! 
                     #   object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = 5,graph=F)
                     # }
                     # 
                     # # if no transformation
                     # else{
                     #   pseudo <- log2(scale(SummarizedExperiment::assay(object), center=FALSE) + 1)
                     #   object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = 5,graph=F)
                     # }
                     # return(object)
                     # 
                   }
)


methods::setMethod(f          = "RunPCA",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, nbcp = 5, raw = FALSE){
                     
                     object[[SE.name]] <- RFLOMICS::RunPCA(object[[SE.name]], 
                                                           nbcp = nbcp, 
                                                           raw  = raw)
                     
                     return(object)
                     
                   })

##### Graphical METHODS for exploring biological and technical variability


#' Library_size_barplot.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @return plot
#' @export
#' @importFrom ggplot2 ggplot geom_bar xlab ylab element_text ggtitle
#' @noRd
#' @examples
methods::setMethod(f          = "Library_size_barplot.plot",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, raw = FALSE){
                     
                     value    <- NULL
                     warnning <- NULL
                     
                     if (object@metadata$omicType != "RNAseq") stop("WARNING: data are not RNAseq!")
                     
                     abundances <- SummarizedExperiment::assay(object)
                     samples    <- colnames(abundances)
                     
                     if (raw) {
                       
                       pseudo <- SummarizedExperiment::assay(object) %>% colSums(., na.rm = TRUE)
                       title  <- "Raw data"
                       
                     }else {
                      
                       # RNAseq, not expected to find any transformation method.
                       if (!object@metadata[["Normalization"]]$normalized) pseudo <- SummarizedExperiment::assay(RFLOMICS:::apply_norm(object)) 
  
                       pseudo <- pseudo %>% colSums(., na.rm = TRUE)
                       title <- paste0("Filtered and normalized (", object@metadata$Normalization$methode, ") data")
                     }
                     
                     # # normalized data
                     # if(! is.null(object@metadata$Normalization)) {
                     #   pseudo  <- scale(SummarizedExperiment::assay(object), center=FALSE,
                     #                    scale=object@metadata$Normalization$coefNorm$norm.factors*object@metadata$Normalization$coefNorm$lib.size) %>% colSums(., na.rm = TRUE)
                     #   title <- paste0("Filtered and normalized (", object@metadata$Normalization$methode, ") data")
                     # }
                     # 
                     # # raw data
                     # else{
                     #   pseudo  <- SummarizedExperiment::assay(object) %>% colSums(., na.rm = TRUE)
                     #   title <- "Raw data"
                     # }
                     
                     libSizeNorm <-  dplyr::full_join(object@metadata$Groups, data.frame("value" = pseudo , "samples" = names(pseudo)), by = "samples") %>%
                       dplyr::arrange(groups)
                     
                     libSizeNorm$samples <- factor(libSizeNorm$samples, levels = unique(libSizeNorm$samples))
                     
                     p <- ggplot2::ggplot(libSizeNorm, ggplot2::aes(x = samples, y = value, fill = groups)) + 
                       ggplot2::geom_bar(stat = "identity" ) + 
                       ggplot2::ylab(ylab) + # useful?
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position  = "none") + 
                       ggplot2::labs(x = "", y = "Total read count per sample") + 
                       ggplot2::ggtitle(title)
                     #axis.text.x     = element_blank(),
                     #axis.ticks      = element_blank())
                     #legend.key.size = unit(0.3, "cm"))
                     #legend.text     = element_text(size=5))
                     print(p)
                     
                   })


methods::setMethod(f          = "Library_size_barplot.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, raw = FALSE){
                     
                     if (RFLOMICS::getOmicsTypes(object[[SE.name]]) == "RNAseq") {
                       return(RFLOMICS::Library_size_barplot.plot(object[[SE.name]], raw = raw))
                     }else{
                       stop("This function only applies to RNAseq data")
                     }

                   })


#' @title Data_Distribution_plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param plot plot type ("boxplot" or "density")
#' @export
#' @exportMethod Data_Distribution_plot
#' @importFrom ggplot2 geom_density xlab
#' @noRd

methods::setMethod(f          = "Data_Distribution_plot",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, plot = "boxplot", raw = FALSE){

                     # Raw data 
                     if (raw) {
                       
                       pseudo <- SummarizedExperiment::assay(object)
                       
                       if (object@metadata$omicType == "RNAseq") {
                         pseudo <- log2(pseudo + 1) 
                         x_lab  <- paste0("log2(", object@metadata$omicType, " data)")
                         title  <- paste0(object@metadata$omicType, " raw data")
                       } else {
                         x_lab  <- paste0(object@metadata$omicType, " data")
                         title  <- paste0(object@metadata$omicType, " raw data")
                       }
                       
                     } else {
                      # Already normalized or transformed Data
                       
                       x_lab  <- paste0(object@metadata$omicType, " data")
                       title  <- paste0(object@metadata$omicType, " data")
                       
                       # object <- MAE[["Metabolites"]]
                       
                       if (!object@metadata[["transform"]][["transformed"]] && object@metadata[["transform"]][["transform_method"]] != "none") {
                         object <- RFLOMICS:::apply_transformation(object)
                         title  <- paste0("Transformed (", object@metadata[["transform"]][["transform_method"]], ") ", title)
                       }
                       if (!object@metadata[["Normalization"]]$normalized && object@metadata[["Normalization"]][["methode"]] != "none") {
                         object <- RFLOMICS:::apply_norm(object)
                         title <- paste0(title, " - normalization: ", object@metadata[["Normalization"]]$methode)
                       }                     
                       
                       pseudo <- SummarizedExperiment::assay(object)
                       
                       if (object@metadata$omicType == "RNAseq") {
                         pseudo <- log2(pseudo) # +1 in NormMethod
                         x_lab  <- paste0("log2(", object@metadata$omicType, " data)")
                        }
                       
                       
                        
                      
                     }
                     
                     # switch(object@metadata$omicType,
                     #        "RNAseq" = {
                     #          
                     #          # before normalization
                     #          if (is.null(object@metadata[["Normalization"]]$coefNorm)) {
                     #            pseudo <- log2(SummarizedExperiment::assay(object) + 1)
                     #            x_lab  <- paste0("log2(", object@metadata$omicType, " data)")
                     #            title  <- paste0(object@metadata$omicType, " raw data")
                     #            
                     #          }
                     #          # after normalization
                     #          else{
                     #            pseudo <- log2(scale(SummarizedExperiment::assay(object) + 1, center = FALSE,
                     #                                 scale = object@metadata[["Normalization"]]$coefNorm$norm.factors*object@metadata[["Normalization"]]$coefNorm$lib.size))
                     #            x_lab  <- paste0("log2(", object@metadata$omicType, " data)")
                     #            title  <- paste0("Filtered and normalized ", object@metadata$omicType, " (",
                     #                             object@metadata$Normalization$methode, ") data")
                     #          }
                     #        },
                     #        "proteomics" = {
                     #          # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
                     #          if (is.null(object@metadata$transform_method) || 
                     #              object@metadata$transform_method == "none" || 
                     #              is.null(object@metadata[["Normalization"]]$methode)) {
                     #            
                     #            pseudo <- SummarizedExperiment::assay(object)
                     #            x_lab  <- paste0(object@metadata$omicType, " data")
                     #            title  <- paste0(object@metadata$omicType, " raw data")
                     #            
                     #          }
                     #          # after transformation
                     #          else{
                     #            
                     #            x_lab  <- paste0(object@metadata$omicType, " data")
                     #            title  <- paste0("Transformed ", object@metadata$omicType, " (" , object@metadata$transform_method, ") data")
                     #            if (is.null(object@metadata[["Normalization"]]$methode)) title <- paste0(title, " - No Normalization")
                     #            else if (object@metadata[["Normalization"]]$methode != "none") title <- paste0(title, " - normalization: ", object@metadata[["Normalization"]]$methode)
                     #            
                     #            switch(object@metadata$transform_method,
                     #                   "log1p" = {
                     #                     pseudo <- log1p(SummarizedExperiment::assay(object)) },
                     #                   "squareroot" = {
                     #                     pseudo <- sqrt(SummarizedExperiment::assay(object)) },
                     #                   "log2" = {
                     #                     pseudo <- log2(SummarizedExperiment::assay(object) + 1 )},
                     #                   "log10" = {
                     #                     pseudo <- log10(SummarizedExperiment::assay(object) + 1)}
                     #            )
                     #            
                     #          }
                     #          
                     #          # Normalization is after the transformation 
                     #          if (!is.null(object@metadata[["Normalization"]]$methode)) {
                     #            
                     #            switch(object@metadata[["Normalization"]]$methode,
                     #                   "none" = {
                     #                     pseudo <- pseudo },
                     #                   "median" = {
                     #                     pseudo <- apply(pseudo, 2, FUN = function(sample_vect) sample_vect - median(sample_vect)) },
                     #                   "totalSum" = {
                     #                     pseudo <- apply(pseudo, 2, FUN = function(sample_vect) sample_vect/sum(sample_vect^2))}
                     #            )
                     #          }
                     #        },
                     #        "metabolomics" = {
                     #          # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
                     #          if (is.null(object@metadata$transform_method) || 
                     #              object@metadata$transform_method == "none" || 
                     #              is.null(object@metadata[["Normalization"]]$methode)) {
                     #            pseudo <- SummarizedExperiment::assay(object)
                     #            x_lab  <- paste0(object@metadata$omicType, " data")
                     #            title  <- paste0(object@metadata$omicType, " raw data")
                     #            
                     #          }
                     #          # after transformation
                     #          else{
                     #            
                     #            x_lab  <- paste0(object@metadata$omicType, " data")
                     #            title  <- paste0("Transformed ", object@metadata$omicType, " (" , object@metadata$transform_method, ") data")
                     #            if (is.null(object@metadata[["Normalization"]]$methode)) title <- paste0(title, " - No Normalization")
                     #            else if (object@metadata[["Normalization"]]$methode != "none") title <- paste0(title, " - normalization: ", object@metadata[["Normalization"]]$methode)
                     #            
                     #            switch(object@metadata$transform_method,
                     #                   
                     #                   "log1p" = {
                     #                     pseudo <- log1p(SummarizedExperiment::assay(object)) },
                     #                   "squareroot" = {
                     #                     pseudo <- sqrt(SummarizedExperiment::assay(object)) },
                     #                   "log2" = {
                     #                     pseudo <- log2(SummarizedExperiment::assay(object) + 1 )},
                     #                   "log10" = {
                     #                     pseudo <- log10(SummarizedExperiment::assay(object) + 1 )}
                     #            )
                     #            
                     #            
                     #          }
                     #          
                     #          # Normalization is after the transformation 
                     #          if (!is.null(object@metadata[["Normalization"]]$methode)) {
                     #            
                     #            switch(object@metadata[["Normalization"]]$methode,
                     #                   "none" = {
                     #                     pseudo <- pseudo },
                     #                   "median" = {
                     #                     pseudo <- apply(pseudo, 2, FUN = function(sample_vect) sample_vect - median(sample_vect)) },
                     #                   "totalSum" = {
                     #                     pseudo <- apply(pseudo, 2, FUN = function(sample_vect) sample_vect/sum(sample_vect^2))}
                     #            )
                     #          }
                     #        }
                     # )
                     # 
                     pseudo.gg <- pseudo %>% reshape2::melt()
                     colnames(pseudo.gg) <- c("features", "samples", "value")
                     
                     pseudo.gg <- pseudo.gg %>% dplyr::full_join(object@metadata$Groups, by = "samples") %>%
                       dplyr::arrange(groups)
                     
                     pseudo.gg$samples <- factor(pseudo.gg$samples, levels = unique(pseudo.gg$samples))
                     
                     switch(plot,
                            "density" = {
                              p <- ggplot2::ggplot(pseudo.gg) + 
                                ggplot2::geom_density(ggplot2::aes(x = value, group = samples, color = groups), trim = FALSE) +
                                ggplot2::xlab(x_lab) + 
                                ggplot2::theme(legend.position = 'none') + 
                                ggplot2::ggtitle(title)
                            },
                            "boxplot" = {
                              p <- ggplot2::ggplot(pseudo.gg, ggplot2::aes(x = samples, y = value,label = features)) +
                                ggplot2::geom_boxplot(ggplot2::aes(fill = groups), outlier.size = 0.3) +
                                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
                                ggplot2::xlab("") + ggplot2::ylab(x_lab) + ggplot2::ggtitle(title) #+
                              #geom_point(alpha = 1/100,size=0)
                            })
                     
                     print(p)
                     
                   }
)


methods::setMethod(f          = "Data_Distribution_plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, plot = "boxplot", raw = FALSE){
                     
                     RFLOMICS::Data_Distribution_plot(object = object[[SE.name]],
                                                      plot = plot,
                                                      raw = raw)
                     
                   })

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
methods::setMethod(f= "plotPCA",
                   signature = "SummarizedExperiment",
                   definition <- function(object, PCA, PCs=c(1,2), condition="groups"){
                     
                     #
                     PC1 <- paste("Dim.",PCs[1], sep="")
                     PC2 <- paste("Dim.",PCs[2], sep="")
                     
                     if(PC1 == PC2){
                       stop("PC1 and PC2 must be different")
                     }
                     
                     score     <- object@metadata$PCAlist[[PCA]]$ind$coord[, PCs] %>% as.data.frame() %>%
                       dplyr::mutate(samples=row.names(.)) %>% dplyr::full_join(., object@metadata$Groups, by="samples")
                     
                     var1 <- round(object@metadata$PCAlist[[PCA]]$eig[PCs,2][1], digits=3)
                     var2 <- round(object@metadata$PCAlist[[PCA]]$eig[PCs,2][2], digits=3)
                     
                     switch (PCA,
                             "raw"  = {title <- paste0("Raw ", object@metadata$omicType, " data")},
                             "norm" = {title <- switch (object@metadata$omicType,
                                                        "RNAseq" = { paste0("Filtred and normalized ", object@metadata$omicType," data (", object@metadata$Normalization$methode, ")")  },
                                                        "proteomics" = {paste0("Transformed and normalized ", object@metadata$omicType," data (", object@metadata$transform_method, 
                                                                               " - norm: ", object@metadata$Normalization$methode, ")")},
                                                        "metabolomics" = {paste0("Transformed and normalized ", object@metadata$omicType," data (", object@metadata$transform_method, 
                                                                                 " - norm: ", object@metadata$Normalization$methode, ")")}
                             )}
                     )
                     
                     
                     
                     p <- ggplot2::ggplot(score, ggplot2::aes_string(x=PC1, y=PC2, color=condition))  +
                       ggplot2::geom_point(size=2) +
                       ggplot2::geom_text(ggplot2::aes(label=samples), size=2, vjust = 0) +
                       ggplot2::xlab(paste(PC1, " (",var1,"%)", sep="")) +
                       ggplot2::ylab(paste(PC2, " (",var2,"%)", sep="")) +
                       ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "red") +
                       ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "red") +
                       ggplot2::theme(strip.text.x = ggplot2::element_text(size=8, face="bold.italic"),
                             strip.text.y = ggplot2::element_text(size=8, face="bold.italic")) +
                       ggplot2::ggtitle(title)
                     
                     # ellipse corr
                     aa <- dplyr::select(score, tidyselect::all_of(condition), PC1, PC2)
                     bb <- FactoMineR::coord.ellipse(aa, bary = TRUE)
                     p <- p + ggplot2::geom_polygon(data = bb$res, aes_string(x=PC1, y=PC2, fill = condition),
                                           show.legend = FALSE,
                                           alpha = 0.1)
                     
                     
                     # if(condition != "groups"){
                     #   p <- p + stat_ellipse(geom="polygon", aes_string(fill = condition),
                     #                         alpha = 0.01,
                     #                         show.legend = FALSE,
                     #                         level = 0.95)
                     # }
                     
                     print(p)
                     
                   })

methods::setMethod(f          = "plotPCA",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, PCA, PCs=c(1,2), condition="groups"){
                     
                     RFLOMICS::plotPCA(object[[SE.name]], PCA, PCs, condition)
                     
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

methods::setMethod(f="mvQCdesign",
                   signature="MultiAssayExperiment",
                   definition <- function(object, data, PCA=c("raw","norm"), axis=5, pngFile=NULL){
                     
                     # Stop if the PCA object does not exist
                     resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
                     
                     if(is.null(resPCA)){
                       stop(paste0(PCA,"PCA does not exist for the ",data," data"))
                     }
                     
                     n_dFac <-  object@metadata$colDataStruc["n_dFac"]
                     
                     # corespondance between coordinate and factor modalities thanks to the sample's name
                     tab_tmp <- merge(as.data.frame(resPCA$ind$coord),as.data.frame(object@colData),by='row.names',all=TRUE)
                     
                     df <- list()
                     bigdf <- list()
                     
                     for(i in 1:n_dFac){
                       
                       Factor <-  object@colData[,i]
                       FactorName <- names(object@colData)[i]
                       nF <- length(Factor)
                       
                       for(j in 1:axis){
                         Dim=paste0("Dim.",j)
                         # select the coordinates and the factor columns
                         df[[j]] <- dplyr::select(tab_tmp, all_of(FactorName),all_of(Dim)) %>%
                           # rename
                           dplyr::rename(.,"Levels"=FactorName,"y"=starts_with("Dim.")) %>%
                           # sort by factor modalities then by coordinate
                           dplyr::arrange(Levels,y) %>%
                           # add column
                           dplyr::mutate(.,"x"=1:nF,
                                         "Axis"=rep(paste("PCA",j, "\n(",round(resPCA$eig[j,2],1),"%)",sep=""),nF),
                                         "FactorN"=rep(FactorName,nF))
                         
                       }
                       bigdf[[i]] <- dplyr::bind_rows(df)
                     }
                     big <- dplyr::bind_rows(bigdf)
                     
                     out <- by(data = big, INDICES = big$FactorN, FUN = function(m) {
                       
                       m <- droplevels(m)
                       
                       m <- ggplot(m,aes(y=y,x=x,colour=Levels))+
                         ggplot2::geom_bar(stat = "identity",position = ggplot2::position_dodge(),aes(fill=Levels))+
                         facet_grid(as.factor(FactorN)~Axis) +
                         labs(x = "Samples", y="Coordinates \n on the PCA axis")+
                         theme(axis.title.y = element_text(size = 10),
                               axis.title.x=element_text(size = 10),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
                     })
                     p <- do.call(gridExtra::grid.arrange, c(out,ncol=1))
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
methods::setMethod(f="mvQCdata",
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
                     
                     if (!is.null(pngFile)){
                       ggplot2::ggsave(filename = pngFile, plot = p)
                     }
                     
                   })



########################################## TRANSFORM DATA #################

#### METHOD to transform data

#' @title TransformData
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param transform_method The transformation to store in the metadata or to store and apply if modify_assay is TRUE.
#' @param modify_assay Boolean. Do the transformation need to be applied on the data? The raw data will be replaced by the transformed ones.
#'
#' @return An object of class \link{SummarizedExperiment}
#'
#' @exportMethod TransformData
#'
#' @examples
methods::setMethod(f= "TransformData",
                   signature = "SummarizedExperiment",
                   definition <- function(object, transform_method = "log2", modify_assay = FALSE){
                     
                     objectTransform <- object
                     
                     # assayTransform  <- SummarizedExperiment::assay(objectTransform)
                     objectTransform@metadata[["transform"]][["transform_method"]] <- transform_method  
                     
                     if (modify_assay){
                       objectTransform <- RFLOMICS:::apply_transformation(object)
                       objectTransform@metadata[["transform"]][["transformed"]] <- TRUE
                     }
                     
                     return(objectTransform)
                     
                   })

methods::setMethod(f          = "TransformData",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, transform_method = "log2", modify_assay = FALSE){
                     
                     object[[SE.name]] <- RFLOMICS::TransformData(object[[SE.name]], 
                                                                  transform_method = transform_method, 
                                                                  modify_assay = modify_assay)
                     
                     return(object)
                     
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
methods::setMethod(f         = "FilterLowAbundance",
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


methods::setMethod(f          = "FilterLowAbundance",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, Filter_Strategy = "NbConditions", CPM_Cutoff = 5){
                     
                     if(RFLOMICS::getOmicsTypes(object[[SE.name]]) == "RNAseq"){
                       object[[SE.name]] <- RFLOMICS::FilterLowAbundance(object          = object[[SE.name]], 
                                                                         Filter_Strategy = Filter_Strategy, 
                                                                         CPM_Cutoff      = CPM_Cutoff)
                       return(object)
                     }else{
                       message("Can't apply this method to omics types other than RNAseq.")
                     }
                     
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
methods::setMethod(f="RunNormalization",
                   signature="SummarizedExperiment",
                   definition <- function(object, NormMethod){
                     
                     coefNorm  = switch(NormMethod,
                                        "TMM"        = TMM.Normalization(SummarizedExperiment::assay(object), object@metadata$Groups$groups),
                                        "median"     = apply(SummarizedExperiment::assay(object), 2, FUN = function(sample_vect) {median(sample_vect)}),
                                        "totalSum"   = apply(SummarizedExperiment::assay(object), 2, FUN = function(sample_vect) {sum(sample_vect^2)})
                     )
                     object@metadata[["Normalization"]]$methode <- NormMethod
                     object@metadata[["Normalization"]]$coefNorm <- coefNorm
                     
                     return(object)
                   })

methods::setMethod(f          = "RunNormalization",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, NormMethod){
                     
                     object[[SE.name]] <- RFLOMICS::RunNormalization(object     = object[[SE.name]],
                                                                     NormMethod = NormMethod)
                     
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
#' \item{method: }{The method used for the differential analysis. }
#' \item{Adj.pvalue.method: The method applied for the pvalue adjustment.}
#' \item{Adj.pvalue.cutoff: The threshold applied for the pvalue adjustment}
#' \item{FDR: }{The false discovery rate given in input}
#' \item{RawDEFres: }{a list giving for each contrast the raw results of the differential analysis method}
#' \item{DEF: }{a list giving for each contrast a data.frame of non filtered differential expressed features}
#' \item{TopDEF: }{a list giving for each contrast a data.frame of differential expressed features by Adj.pvalue.cutoff}
#' \item{mergeDEF: }{A data frame indicating for each features in row, if it is DE in a given contrasts in column}
#' }
#' @param object an object of class [\code\link{SummarizedExperiment}]
#' @param design an object of class [\code{\link{ExpDesign-class}]
#' @param DiffAnalysisMethod A character vector giving the name of the differential analysis method
#' to run. Either "edgeRglmfit" or "limmalmFit".
#' @param contrastList The list of contrast to test
#' @param Adj.pvalue.method The method choosen to adjust pvalue. Takes the same values as the ones of adj.p.adjust method.
#' @param Adj.pvalue.cutoff The adjusted pvalue cut-off
#' @param clustermq A boolean indicating whether the constrasts have to be computed in local or in a distant machine
#' @return An object of class [\code\link{SummarizedExperiment}]
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod RunDiffAnalysis
#' @examples
#'
#'
methods::setMethod(f         = "RunDiffAnalysis",
                   signature = "SummarizedExperiment",
                   definition <- function(object, design, Adj.pvalue.method="BH",
                                          contrastList, DiffAnalysisMethod, Adj.pvalue.cutoff=0.05, logFC.cutoff=0, clustermq=FALSE, parallel = FALSE, nworkers = 1){
                     
                     contrastName <- NULL
                     Contrasts.Sel <- dplyr::filter(design@Contrasts.Sel, contrastName %in% contrastList)
                     
                     object@metadata$DiffExpAnal <- list()
                     object@metadata$DiffExpAnal[["contrasts"]] <- Contrasts.Sel
                     object@metadata$DiffExpAnal[["method"]]    <- DiffAnalysisMethod
                     object@metadata$DiffExpAnal[["Adj.pvalue.method"]]  <- Adj.pvalue.method
                     #object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]  <- Adj.pvalue.cutoff
                     
                     # transform and norm if needed
                     if (DiffAnalysisMethod == "limmalmFit") {
                       
                       object2 <- object
                       
                       if (!object2@metadata[["transform"]][["transformed"]]) object2 <- RFLOMICS:::apply_transformation(object2)
                       if (!object2@metadata[["Normalization"]]$normalized)   object2 <- RFLOMICS:::apply_norm(object2)

                     }
                     
                     # move in ExpDesign Constructor
                     model_matrix <- model.matrix(as.formula(paste(design@Model.formula, collapse = " ")), data = as.data.frame(design@List.Factors))
                     rownames(model_matrix) <- rownames(design@ExpDesign)
                     
                     ListRes <- switch(DiffAnalysisMethod,
                                       "edgeRglmfit" = try_rflomics(edgeR.AnaDiff(count_matrix  = SummarizedExperiment::assay(object),
                                                                                model_matrix    = model_matrix[colnames(object),],
                                                                                group           = object@metadata$Normalization$coefNorm$group,
                                                                                lib.size        = object@metadata$Normalization$coefNorm$lib.size,
                                                                                norm.factors    = object@metadata$Normalization$coefNorm$norm.factors,
                                                                                Contrasts.Sel   = object@metadata$DiffExpAnal[["contrasts"]],
                                                                                Contrasts.Coeff = design@Contrasts.Coeff,
                                                                                FDR             = 1,
                                                                                clustermq       = clustermq,
                                                                                parallel        = parallel,
                                                                                nworkers        = nworkers)),
                                       "limmalmFit" = try_rflomics(limma.AnaDiff(count_matrix    = SummarizedExperiment::assay(object2),
                                                                               model_matrix      = model_matrix[colnames(object2),],
                                                                               Contrasts.Sel     = object2@metadata$DiffExpAnal[["contrasts"]],
                                                                               Contrasts.Coeff   = design@Contrasts.Coeff,
                                                                               Adj.pvalue.cutoff = 1,
                                                                               Adj.pvalue.method = Adj.pvalue.method,
                                                                               clustermq         = clustermq)))
                     
                     if(! is.null(ListRes$value)){
                       if(! is.null(ListRes$value[["RawDEFres"]])){
                         object@metadata$DiffExpAnal[["results"]] <- TRUE
                         object@metadata$DiffExpAnal[["RawDEFres"]] <- ListRes$value[["RawDEFres"]]
                         object@metadata$DiffExpAnal[["DEF"]] <- ListRes$value[["TopDEF"]]
                       }else{
                         object@metadata$DiffExpAnal[["results"]]    <- FALSE
                         object@metadata$DiffExpAnal[["ErrorStats"]] <- ListRes$value[["ErrorTab"]]
                       }
                     }else{
                       object@metadata$DiffExpAnal[["results"]]    <- FALSE
                       object@metadata$DiffExpAnal[["Error"]]      <- ListRes$error
                       object@metadata$DiffExpAnal[["ErrorStats"]] <- NULL
                     }
                     
                     ## filtering
                     object <- RFLOMICS::FilterDiffAnalysis(object = object, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff)
                     
                     return(object)
                   })



methods::setMethod(f          = "RunDiffAnalysis",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, design = NULL, Adj.pvalue.method="BH",
                                          contrastList = NULL, DiffAnalysisMethod = NULL, Adj.pvalue.cutoff=0.05, logFC.cutoff=0, clustermq=FALSE, parallel = FALSE, nworkers = 1){
                     
                     # Check for design existence inside MAE object
                     if (is.null(design)) {
                       if (!is.null(object@metadata$design))  design <- object@metadata$design
                       else stop("Argument design is missing and does not exist in the object.")
                     }
                     
                     # Check for contrastList existence inside MAE object
                     if (is.null(contrastList)) {
                       if (!is.null(object@metadata$design@Contrasts.Sel)) contrastList <- object@metadata$design@Contrasts.Sel$contrastName
                       else stop("Argument contrastList is missing and no selected contrasts have been found.")
                     }
                     
                     # Type of DiffAnalysis automatically decided if argument is null
                     if (is.null(DiffAnalysisMethod)) {
                       
                       switch(RFLOMICS::getOmicsTypes(object[[SE.name]]),
                              "RNAseq" = {DiffAnalysisMethod = "edgeRglmfit"
                              message("DiffAnalyseMethod was missing. Detected omic type is RNASeq, using edgeRglmFit for differential analysis.")},
                              "proteomics" = {DiffAnalysisMethod = "limmalmFit"
                              message("DiffAnalyseMethod was missing. Detected omic type is proteomics, using limmalmFit for differential analysis.")},
                              "metabolomics" = {DiffAnalysisMethod = "limmalmFit"
                              message("DiffAnalyseMethod was missing. Detected omic type is metabolomics, using limmalmFit for differential analysis.")}
                       )
                       
                     }
                     
                     object[[SE.name]] <- RFLOMICS::RunDiffAnalysis(object = object[[SE.name]],
                                                                    design = design, 
                                                                    Adj.pvalue.method = Adj.pvalue.method,
                                                                    contrastList = contrastList,
                                                                    DiffAnalysisMethod = DiffAnalysisMethod,
                                                                    Adj.pvalue.cutoff = Adj.pvalue.cutoff,
                                                                    logFC.cutoff = logFC.cutoff,
                                                                    clustermq = clustermq,
                                                                    parallel = parallel,
                                                                    nworkers = nworkers
                                                                    )
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
methods::setMethod(f          = "FilterDiffAnalysis",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, Adj.pvalue.cutoff = 0.05, logFC.cutoff = 0){
                     
                     if(is.null(object@metadata$DiffExpAnal[["RawDEFres"]])){
                       stop("can't filter the DiffExpAnal object because it doesn't exist")
                     }
                     
                     object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]  <- Adj.pvalue.cutoff
                     object@metadata$DiffExpAnal[["abs.logFC.cutoff"]]  <- logFC.cutoff
                     
                     ## TopDEF: Top differential expressed features
                     DEF_filtred <- lapply(1:length(object@metadata$DiffExpAnal[["DEF"]]),function(x){
                       res <- object@metadata$DiffExpAnal[["DEF"]][[x]]
                       keep <- (res$Adj.pvalue <= Adj.pvalue.cutoff) & (abs(res$logFC) > logFC.cutoff)
                       res <- res[keep,]
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
                     DEF_list <- list()
                     for(x in names(object@metadata$DiffExpAnal[["TopDEF"]])){
                       res <- object@metadata$DiffExpAnal[["TopDEF"]][[x]]
                       tmp <- data.frame(DEF = rownames(res), bin = rep(1,length(rownames(res))))
                       colnames(tmp) <- c("DEF", dplyr::filter(object@metadata$DiffExpAnal$contrasts, contrastName == x)$tag)
                       
                       if(dim(tmp)[1] != 0){ DEF_list[[x]] <- tmp }
                     }
                     
                     object@metadata$DiffExpAnal[["mergeDEF"]] <- NULL
                     
                     if(length(DEF_list) != 0){
                       
                       object@metadata$DiffExpAnal[["mergeDEF"]] <- DEF_list %>% purrr::reduce(dplyr::full_join, by="DEF") %>%
                         dplyr::mutate_at(.vars = 2:(length(DEF_list)+1),
                                          .funs = function(x){
                                            dplyr::if_else(is.na(x), 0, 1)}) %>%
                         data.table::data.table()
                     }
                     
                     return(object)
                   })


methods::setMethod(f          = "FilterDiffAnalysis",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, Adj.pvalue.cutoff = 0.05, logFC.cutoff = 0){
                     
                     object[[SE.name]] <- RFLOMICS::FilterDiffAnalysis(object = object[[SE.name]],
                                                                       Adj.pvalue.cutoff = Adj.pvalue.cutoff,
                                                                       logFC.cutoff = logFC.cutoff)
                     
                     return(object)
                     
                   })

###### Graphical METHOD

## Method to plot results of a differential analysis

#' @title DiffAnal.plot
#' @description
#' This is an interface method which draw a MAplot, a volcano plot and the pvalues distribution from the results of a differential analysis
#' performed on omic datasets stored in an object of class \link{SummarizedExperiment}
#' @param object An object of class \link{SummarizedExperiment}
#' @param hypothesis The hypothesis for which the plots has to be drawn
#' @param typeofplots The plots you want to return. Default is all possible plots: MA plot, Volcano plot and non adjusted pvalues histogram.
#' @return plot
#' @exportMethod DiffAnal.plot
#' @export
#'
#' @examples
methods::setMethod(f="DiffAnal.plot",
                   signature="SummarizedExperiment",
                   
                   definition <- function(object, hypothesis, typeofplots = c("MA.plot", "volcano", "histogram")){
                     
                     if (RFLOMICS::isTagName(object, hypothesis)) hypothesis <- RFLOMICS::convertTagToContrast(object, hypothesis)
                     
                     plots <- list()
                     
                     res      <- object@metadata$DiffExpAnal[["RawDEFres"]][[hypothesis]]
                     resTable <- object@metadata$DiffExpAnal[["DEF"]][[hypothesis]]
                     
                     logFC.cutoff      <- object@metadata$DiffExpAnal[["abs.logFC.cutoff"]]
                     Adj.pvalue.cutoff <- object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]
                     
                     if ("MA.plot" %in% typeofplots) plots[["MA.plot"]]        <- RFLOMICS::MA.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff, hypothesis=hypothesis)
                     if ("volcano" %in% typeofplots) plots[["Volcano.plot"]]   <- RFLOMICS::Volcano.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff, hypothesis=hypothesis)
                     if ("histogram" %in% typeofplots) plots[["Pvalue.hist"]]  <- RFLOMICS::pvalue.plot(data =resTable, hypothesis=hypothesis)
                     return(plots)
                   })


methods::setMethod(f          = "DiffAnal.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, hypothesis, typeofplots = c("MA.plot", "volcano", "histogram")){
                     
                     if (RFLOMICS::isTagName(object, hypothesis)) hypothesis <- RFLOMICS::convertTagToContrast(object, hypothesis)
                     
                     return(RFLOMICS::DiffAnal.plot(object      = object[[SE.name]],
                                                    hypothesis  = hypothesis,
                                                    typeofplots = typeofplots))
                     
                   })


#' @title heatmap.plot
#' @description
#' This is an interface method which draw a heatmap from the results of a differential analysis
#' performed on omic datasets stored in an object of class \link{SummarizedExperiment}
#' @param object An object of class \link{SummarizedExperiment}
#' @param hypothesis The hypothesis for which the MAplot has to be drawn
#' @param condition characters. Default to none. Name of a feature in the design matrix, splits the samples on the heatmap according to its modalities.  
#' @param title characters. Title of the heatmap. 
#' @param annot_to_show vector. Names of the annotations to keep in the Heatmap. Default takes all available information.
#' @param subset_list named list of vectors of modalities to subset and print on the heatmap. 
#' @param draw_args,heatmap_args  named lists. Any additional parameter passed to ComplexHeatmap::Heatmap or ComplexHeatmap::draw
#' @return plot
#' @exportMethod heatmap.plot
#' @export
#'
#' @examples
methods::setMethod(f          = "heatmap.plot",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, hypothesis, condition="none", title = "", annot_to_show = NULL, subset_list = NULL, draw_args = list(), heatmap_args = list()){

                     if(is.null(object@metadata$DiffExpAnal[["TopDEF"]][[hypothesis]])){
                       stop("no DE variables")
                     }
                     
                     resTable <- dplyr::arrange(object@metadata$DiffExpAnal[["TopDEF"]][[hypothesis]], Adj.pvalue)
                     
                     if(dim(resTable)[1] == 0){
                       stop("no differentially expressed variables...")
                     }
                     
                     
                     if(dim(resTable)[1] > 2000){
                       message("differentially expressed variables exceeding 2000 variables")
                       resTable <- resTable[1:2000,]
                       title = ifelse(title == "", paste0(title, "plot only 2000 TOP DE variables"),
                                                  paste0(title, "\nplot only 2000 TOP DE variables"))
                     }
                     
                     # m.def <- assays(object)[[1]][,object@metadata$Groups$samples]
                     
                     object2 <- object
                     
                     if (!object2@metadata[["transform"]][["transformed"]]) object2 <- RFLOMICS:::apply_transformation(object2)
                     if (!object2@metadata[["Normalization"]]$normalized)   object2 <- RFLOMICS:::apply_norm(object2)
                     
                     m.def <- SummarizedExperiment::assays(object2)[[1]]
                     if (object2@metadata$omicType == "RNAseq")  m.def <- log2(m.def) 
                     m.def <- as.data.frame(m.def) %>% dplyr::select(tidyselect::any_of(object2@metadata$Groups$samples))
                     
                     # switch (object@metadata$omicType,
                     #         "RNAseq" = {
                     #           
                     #           m.def <- log2(scale(m.def+1, center = FALSE,
                     #                               scale = object@metadata$Normalization$coefNorm$lib.size*object@metadata$Normalization$coefNorm$lib.size))
                     #         },
                     #         "proteomics" = {
                     #           m.def <- m.def # data already transformed at this point (applied on filtred)
                     #          
                     #           if(!is.null(object@metadata[["Normalization"]]$methode)){
                     #             
                     #             switch (object@metadata[["Normalization"]]$methode,
                     #                     "none" = {
                     #                       m.def <- m.def },
                     #                     "median" = {
                     #                       m.def <- t(t(m.def) - object@metadata[["Normalization"]]$coefNorm)},
                     #                     "totalSum" = {
                     #                       m.def <- t(t(m.def)/object@metadata[["Normalization"]]$coefNorm)}
                     #             )
                     #           }
                     #           
                     #            # switch (object@metadata$transform_method,
                     #           #         "log1p"      = { m.def <- log1p(m.def) },
                     #           #         "squareroot" = { m.def <- sqrt(m.def) },
                     #           #         "log2"       = { m.def <- log2(m.def + 1) },
                     #           #         "log10"      = { m.def <- log10(m.def + 1) },
                     #           #         "none"       = { m.def <- m.def }
                     #           # )
                     #           
                     #          
                     #         },
                     #         "metabolomics" = {
                     #           
                     #           m.def <- m.def # data already transformed at this point (applied on filtred)
                     #           
                     #           if(!is.null(object@metadata[["Normalization"]]$methode)){
                     #             
                     #             switch (object@metadata[["Normalization"]]$methode,
                     #                     "none" = {
                     #                       m.def <- m.def },
                     #                     "median" = {
                     #                       m.def <- t(t(m.def) - object@metadata[["Normalization"]]$coefNorm)},
                     #                     "totalSum" = {
                     #                       m.def <- t(t(m.def)/object@metadata[["Normalization"]]$coefNorm)}
                     #             )
                     #           }
                               
                               # switch (object@metadata$transform_method,
                               #         "log1p"      = { m.def <- log1p(m.def) },
                               #         "squareroot" = { m.def <- sqrt(m.def) },
                               #         "log2"       = { m.def <- log2(m.def + 1) },
                               #         "log10"      = { m.def <- log10(m.def + 1) },
                               #         "none"       = { m.def <- m.def }
                               # )
                             # }
                     # )
                     
     
                     
                     # filter by DE
                     m.def.filter <- subset(m.def, rownames(m.def) %in% row.names(resTable))
                     
                     # normalize count
                     
                     # Center
                     # m.def.filter.center <- scale(m.def.filter,center=TRUE,scale=FALSE)
                     m.def.filter.center <- t(scale(t(m.def.filter), center = TRUE, scale = FALSE)) # Modified 221123 : centered by genes and not by samples
         
  
                     # Annotations datatable
                     df_annotation <- object@metadata$Groups %>% dplyr::select(!samples & !groups)  
                     df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),] 
                     
                     # Subset the dataset to print only interesting modalities
                     if (!is.null(subset_list)) {
                       if (is.null(names(subset_list))) {
                         message("In plot.heatmap, subset_list argument needs a named list. Not subsetting")
                       }else{ 
                         samplesToKeep <- Reduce("intersect", lapply(
                           1:length(subset_list),
                           FUN = function(i){
                             col_nam <- names(subset_list)[i]
                             rownames(df_annotation[which(df_annotation[[col_nam]] %in% subset_list[[i]]),])
                           }
                         ))
                         
                         # print(samplesToKeep)
                         
                         df_annotation <- df_annotation[which(rownames(df_annotation) %in% samplesToKeep),]
                         m.def.filter.center <- m.def.filter.center[, which(colnames(m.def.filter.center) %in% samplesToKeep)]
                         
                         # print(dim(df_annotation))
                         # print(m.def.filter.center)
                         
                         df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),]
                       }
                     }
                     
                     # Select the right columns
                     if (!is.null(annot_to_show)) {
                       df_annotation <- df_annotation %>% dplyr::select(tidyselect::any_of(annot_to_show))
                     }
                     
                     # Split management
                     column_split.value <- if (condition != "none") { df_annotation[, condition] }else{NULL}
                     
                     # Color annotations
                     set.seed(10000) ; selectPal <- sample(rownames(RColorBrewer::brewer.pal.info),  size = ncol(df_annotation), replace = FALSE)
                     
                     color_list <- lapply(1:ncol(df_annotation), FUN = function(i){
                       annot_vect <- unique(df_annotation[,i])
                       
                        col_vect <-  grDevices::colorRampPalette(
                          suppressWarnings({ RColorBrewer::brewer.pal(n = min(length(annot_vect), 8), name = selectPal[i])}))(length(annot_vect)) 
                       names(col_vect) <- annot_vect 
                       col_vect[!is.na(names(col_vect))] # RcolorBrewer::brewer.pal n is minimum 3, remove NA names if only 2 levels
                     })
                     names(color_list) <- colnames(df_annotation)
                     
                     column_ha <- ComplexHeatmap::HeatmapAnnotation(df = df_annotation, col = color_list)
                     
                     
                     # names(formals(ComplexHeatmap::Heatmap))
                     
                     namArg <- ifelse(RFLOMICS::getOmicsTypes(object) == "RNAseq", "normalized counts", "XIC")
                     
                     # Arguments for Heatmap
                     heatmap_args <- c(
                       list(matrix = m.def.filter.center,
                            name = namArg,
                            show_row_names = ifelse( dim(m.def.filter.center)[1] > 50, FALSE, TRUE),
                            row_names_gp = grid::gpar(fontsize = 8),
                            column_names_gp = grid::gpar(fontsize = 12),
                            row_title_rot = 0 ,
                            clustering_method_columns = "ward.D2",
                            cluster_column_slice = FALSE,
                            column_split = column_split.value,
                            top_annotation = column_ha,
                            column_title = title),
                       heatmap_args)
                     
                     # Arguments for drawing the heatmap
                     draw_args <- c(list(merge_legend = TRUE),
                                    draw_args)            

                     # Drawing heatmap in a null file to not plot it
                     pdf(file = NULL)
                     ha <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
                     
                     draw_args$object <- ha
                     ha <- do.call(ComplexHeatmap::draw, draw_args)
                     
                     # ha <- ComplexHeatmap::draw(ComplexHeatmap::Heatmap(m.def.filter.center, name = "normalized counts\nor XIC",
                     #                               show_row_names= ifelse( dim(m.def.filter.center)[1] > 50, FALSE, TRUE),
                     #                               row_names_gp = grid::gpar(fontsize = 8),
                     #                               column_names_gp = grid::gpar(fontsize = 12),
                     #                               row_title_rot = 0 ,
                     #                               clustering_method_columns = "ward.D2",
                     #                               cluster_column_slice = FALSE,
                     #                               column_split = column_split.value,
                     #                               top_annotation = column_ha,
                     #                               column_title = title, match.arg(..., names(formals(ComplexHeatmap::Heatmap)))),  merge_legend = TRUE, match.arg(..., names(formals(ComplexHeatmap::draw))))
                     dev.off()
                     
                     return(ha)
                   })


methods::setMethod(f          = "heatmap.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, hypothesis, condition="none", title = "", annot_to_show = NULL, subset_list = NULL, draw_args = list(), heatmap_args = list()){
                     
                     
                     if (RFLOMICS::isTagName(object, hypothesis)) hypothesis <- RFLOMICS::convertTagToContrast(object, hypothesis)
                     
                     return(RFLOMICS::heatmap.plot(object        = object[[SE.name]],
                                                   hypothesis    = hypothesis,
                                                   condition     = condition,
                                                   title         = title,
                                                   annot_to_show = annot_to_show,
                                                   subset_list   = subset_list,
                                                   draw_args     = draw_args,
                                                   heatmap_args  = heatmap_args))
                     
                   })
                   
                   
                   

#' @title boxplot.DE.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param DE variable name (gene/protein/metabolite name)
#' @export
#' @exportMethod boxplot.DE.plot
#' @importFrom ggplot2 geom_density xlab
#' @noRd

methods::setMethod(f          = "boxplot.DE.plot",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, DE = NULL, condition="groups", raw = FALSE){
                     
                     # check variable name
                     if (is.null(DE) | DE == "" | length(DE) != 1) {
                       message("set variable name")
                       
                       p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("set variable name") 
                       
                       return(p)
                     }
                     
                     if (!raw) {
                       if (!object@metadata[["transform"]][["transformed"]]) object <- RFLOMICS:::apply_transformation(object)
                       if (!object@metadata[["Normalization"]]$normalized)   object <- RFLOMICS:::apply_norm(object)
                     }
                     
                     # check presence of variable in SE
                     object.DE <- tryCatch(object[DE], error = function(e) e)
                     if (!is.null(object.DE$message)) {
                       message(object.DE$message)
                       
                       p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(object.DE$message) 
                       
                       return(p)
                     }
                     
                     if (raw) {
                       if (object.DE@metadata$omicType != "RNAseq") {
                         
                         pseudo <- SummarizedExperiment::assay(object.DE)
                         x_lab  <- DE
                         title  <- DE
                         
                       } else {
                         
                         pseudo <- log2(SummarizedExperiment::assay(object.DE) + 1)
                         
                         x_lab  <- paste0("log2(", DE, " data)")
                         title  <- DE
                         
                       }
                       
                     } else{ 
                       if (object.DE@metadata$omicType != "RNAseq") {
                         
                         title <- DE
                         
                         if (!object.DE@metadata[["transform"]][["transformed"]] && object.DE@metadata[["transform"]][["transform_method"]] != "none") {
                           title  <- paste0("Transformed (", object.DE@metadata[["transform"]][["transform_method"]], ") ", title)
                         }
                         if (!object.DE@metadata[["Normalization"]]$normalized && object.DE@metadata[["Normalization"]][["methode"]] != "none") {
                           title <- paste0(title, " - normalization: ", object.DE@metadata[["Normalization"]]$methode)
                         }                     
                         
                         pseudo <- SummarizedExperiment::assay(object.DE)
                         x_lab  <- paste0(DE, " data")
                         
                       } else {
                         
                         pseudo <- log2(SummarizedExperiment::assay(object.DE)) # +1 inside the norm method
                         title <- DE
                         x_lab  <- paste0("log2(", DE, " data)") 
                         
                       }
                     }
                     
                     
                     
                     # switch(object.DE@metadata$omicType,
                     #         "RNAseq" = {
                     #           
                     #           # before normalization
                     #           if(is.null(object.DE@metadata[["Normalization"]]$coefNorm)){
                     #             pseudo <- log2(SummarizedExperiment::assay(object.DE))
                     #             x_lab  <- paste0("log2(", DE, " data)")
                     #             title  <- paste0(DE)
                     #             
                     #           }
                     #           # after normalization
                     #           else{
                     #             pseudo <- log2(scale(SummarizedExperiment::assay(object.DE)+1, center=FALSE,
                     #                                  scale=object.DE@metadata[["Normalization"]]$coefNorm$norm.factors*object.DE@metadata[["Normalization"]]$coefNorm$lib.size))
                     #             x_lab  <- paste0("log2(",DE, " data)")
                     #             title  <- paste0("Filtered and normalized (",object.DE@metadata$Normalization$methode, ") : " , DE)
                     #           }
                     #         },
                     #         "proteomics" = {
                     #           # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
                     #           if(object.DE@metadata$transform_method == "none" || is.null(object.DE@metadata$transform_method)){
                     #             pseudo <- SummarizedExperiment::assay(object.DE)
                     #             x_lab  <- paste0(DE, " data")
                     #             title  <- paste0(DE)
                     #             
                     #           }
                     #           # after transformation
                     #           else{
                     #             pseudo <- SummarizedExperiment::assay(object.DE)
                     #             x_lab  <- paste0(DE, " data")
                     #             title  <- paste0("Transformed (" , object.DE@metadata$transform_method, ") : ", DE)
                     #             
                     #             if(!is.null(object.DE@metadata[["Normalization"]]$methode)){
                     #               
                     #               switch (object.DE@metadata[["Normalization"]]$methode,
                     #                       "none" = {
                     #                         pseudo <- pseudo },
                     #                       "median" = {
                     #                         pseudo <- t(t(pseudo) - object.DE@metadata[["Normalization"]]$coefNorm)},
                     #                       "totalSum" = {
                     #                         pseudo <- t(t(pseudo)/object.DE@metadata[["Normalization"]]$coefNorm)}
                     #               )
                     #             }
                     #             
                     #           }
                     #         },
                     #         "metabolomics" = {
                     #           # before rflomics transformation (plot without log2; because we don't know if input prot/meta are transformed or not)
                     #           if(object.DE@metadata$transform_method == "none" || is.null(object.DE@metadata$transform_method)){
                     #             pseudo <- SummarizedExperiment::assay(object.DE)
                     #             x_lab  <- paste0(DE, " data")
                     #             title  <- paste0(DE)
                     #             
                     #           }
                     #           # after transformation
                     #           else{
                     #             
                     #             pseudo <- SummarizedExperiment::assay(object.DE)
                     #             x_lab  <- paste0(DE, " data")
                     #             title  <- paste0("Transformed (" , object.DE@metadata$transform_method, ") : ", DE)
                     #             
                     #             if(!is.null(object.DE@metadata[["Normalization"]]$methode)){
                     #               
                     #               switch (object.DE@metadata[["Normalization"]]$methode,
                     #                       "none" = {
                     #                         pseudo <- pseudo },
                     #                       "median" = {
                     #                         pseudo <- t(t(pseudo) - object.DE@metadata[["Normalization"]]$coefNorm)},
                     #                       "totalSum" = {
                     #                         pseudo <- t(t(pseudo)/object.DE@metadata[["Normalization"]]$coefNorm)}
                     #               )
                     #             }
                     #           }
                     #         }
                     # )
                     
                     pseudo.gg <- pseudo %>% reshape2::melt()
                     colnames(pseudo.gg) <- c("features", "samples", "value")
                     
                     pseudo.gg <- pseudo.gg %>% dplyr::full_join(object@metadata$Groups, by="samples") %>%
                       dplyr::arrange(groups)
                     
                     pseudo.gg <- dplyr::arrange(pseudo.gg, get(condition))
                     
                     pseudo.gg$groups <- factor(pseudo.gg$groups, levels = unique(pseudo.gg$groups))
                     
                     p <- ggplot(pseudo.gg, aes(x=groups, y=value, label = features)) +
                       ggplot2::geom_boxplot(aes(fill=get(condition))) +
                       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                       guides(fill=guide_legend(title="condition")) + 
                       xlab("") + ylab(x_lab) + ggtitle(title) #+
                     #geom_point(alpha = 1/100,size=0)
                     
                     
                     return(p)
                     
                   }
)


methods::setMethod(f          = "boxplot.DE.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, DE = NULL, condition="groups"){
                   
                     RFLOMICS::boxplot.DE.plot(object = object[[SE.name]], 
                                               DE = DE,
                                               condition = condition)
                     
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
methods::setMethod(f="runCoExpression",
                   signature="SummarizedExperiment",
                   definition <- function(object, K=2:20, replicates=5, nameList, merge="union",
                                          model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation, normFactors, clustermq=FALSE){
                     
                     
                     CoExpAnal <- list()
                     
                     CoExpAnal[["tools"]]            <- "CoSeq"
                     CoExpAnal[["gene.list.names"]]  <- nameList
                     CoExpAnal[["merge.type"]]       <- merge
                     CoExpAnal[["replicates.nb"]]    <- replicates
                     CoExpAnal[["K.range"]]          <- K
                     
                     geneList <- dplyr::select(object@metadata$DiffExpAnal[["mergeDEF"]], DEF, tidyselect::all_of(nameList)) %>% 
                       dplyr::mutate(intersection = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) == length(nameList), "YES", "NO"), 
                                     union = dplyr::if_else(rowSums(dplyr::select(., tidyselect::contains(nameList))) != 0 , "YES", "NO")) %>% 
                       dplyr::filter(union != "NO", get(merge) == "YES") 
                     geneList <- geneList$DEF
                     
                     # if (!object@metadata[["transform"]][["transformed"]]) object <- RFLOMICS:::apply_transformation(object)
                     # if (!object@metadata[["Normalization"]]$normalized)   object <- RFLOMICS:::apply_norm(object)
                     # counts = SummarizedExperiment::assay(object)[geneList,]
                     
                     
                     # set default parameters based on data type
                     param.list <- list("meanFilterCutoff"=NULL)
                     switch (object@metadata$omicType,
                             
                             "RNAseq" = {
                               counts = SummarizedExperiment::assay(object)[geneList,]
                               
                               param.list[["model"]]            <- model
                               param.list[["transformation"]]   <- "arcsin"
                               param.list[["normFactors"]]      <- "TMM"
                               param.list[["meanFilterCutoff"]] <- 50
                               param.list[["GaussianModel"]]    <- GaussianModel
                               
                             },
                             "proteomics" = {
                               if (!object@metadata[["transform"]][["transformed"]]) object <- RFLOMICS:::apply_transformation(object)
                               if (!object@metadata[["Normalization"]]$normalized)   object <- RFLOMICS:::apply_norm(object)
                               
                               counts = SummarizedExperiment::assay(object)[geneList,]
                               
                               # Print the selected GaussianModel
                               print(paste("Use ", GaussianModel, sep = ""))
                               print("Scale each protein (center=TRUE,scale = TRUE)")
                               CoExpAnal[["transformation.prot"]] <- "scaleProt"
                               counts[] <- t(apply(counts,1,function(x){ scale(x, center = TRUE, scale = TRUE) }))
                               
                               # param
                               param.list[["model"]]            <- model
                               param.list[["transformation"]]   <- "none"
                               param.list[["normFactors"]]      <- "none"
                               #param.list[["meanFilterCutoff"]] <- NULL
                               param.list[["GaussianModel"]]    <- GaussianModel
                             },
                             "metabolomics" = {
                               if (!object@metadata[["transform"]][["transformed"]]) object <- RFLOMICS:::apply_transformation(object)
                               if (!object@metadata[["Normalization"]]$normalized)   object <- RFLOMICS:::apply_norm(object)
                               
                               counts = SummarizedExperiment::assay(object)[geneList,]
                               
                               # Print the selected GaussianModel
                               print(paste("Use ",GaussianModel,sep=""))
                               print("Scale each metabolite (center=TRUE,scale = TRUE)")
                               CoExpAnal[["transformation.metabo"]] <- "scaleMetabo"
                               counts[] <- t(apply(counts,1,function(x){ scale(x, center = TRUE, scale = TRUE) }))
                               
                               # param
                               param.list[["model"]]            <- model
                               param.list[["transformation"]]   <- "none"
                               param.list[["normFactors"]]      <- "none"
                               #param.list[["meanFilterCutoff"]] <- NULL
                               param.list[["GaussianModel"]]    <- GaussianModel
                             }
                     )
                     
                     CoExpAnal[["param"]] <- param.list
                     
                     # run coseq : on local machine or remote cluster
                     
                     print("#     => coseq... ")
                     
                     coseq.res.list <- list()
                     
                     coseq.res.list <- switch (as.character(clustermq),
                                               `FALSE` = {
                                                 
                                                 try_rflomics(
                                                   runCoseq_local(counts, conds = object@metadata$Groups$groups, K=K, replicates=replicates, param.list=param.list))
                                                 
                                               },
                                               `TRUE` = {
                                                 
                                                 try_rflomics(
                                                   runCoseq_clustermq(counts, conds = object@metadata$Groups$groups, K=K, replicates=replicates, param.list=param.list))
                                                 
                                               })
                     
                     # If coseq could run (no problem with SSH connexion in case of clustermq=TRUE)
                     
                     if(! is.null(coseq.res.list$value)){
                       
                       CoExpAnal <- c(CoExpAnal, coseq.res.list$value)
                     }
                     else{
                       CoExpAnal[["results"]] <- FALSE
                       CoExpAnal[["stats"]]   <- NULL
                       CoExpAnal[["error"]]   <- coseq.res.list$error
                     }
                     
                     object@metadata$CoExpAnal <- CoExpAnal
                     return(object)
                   })


methods::setMethod(f          = "runCoExpression",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, K=2:20, replicates=5, nameList, merge="union",
                                         model = "Normal", GaussianModel = "Gaussian_pk_Lk_Ck", transformation, normFactors, clustermq=FALSE){
                     
                     
                     object[[SE.name]] <- RFLOMICS::runCoExpression(object = object[[SE.name]],
                                                                    K = K,
                                                                    replicates = replicates,
                                                                    nameList = nameList, 
                                                                    merge = merge, 
                                                                    model = model,
                                                                    GaussianModel = GaussianModel,
                                                                    transformation = transformation,
                                                                    normFactors = normFactors,
                                                                    clustermq = clustermq)
                     
                     return(object)
                     
                   })

# Pour utiliser la fonction repeatable(), "seed"  pourrait être ajouté en paramètre.

#' @title coseq.profile.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param numCluster cluster number
#' @param condition 
#' @param observation 
#' @export
#' @exportMethod coseq.profile.plot
#' @noRd

methods::setMethod(f="coseq.profile.plot",
                   signature="SummarizedExperiment",
                   definition <- function(object, numCluster = 1, condition="groups", observation=NULL){
                     
                     coseq.res  <- object@metadata$CoExpAnal[["coseqResults"]]
                     assays.data <- dplyr::filter(as.data.frame(coseq.res@assays@data[[1]]), get(paste0("Cluster_",numCluster)) > 0.8)
                     
                     y_profiles.gg <- coseq.res@y_profiles[rownames(assays.data),] %>% data.frame() %>% dplyr::mutate(observations=rownames(.)) %>% 
                       reshape2::melt(id="observations", value.name = "y_profiles") %>%  dplyr::rename(samples = variable) %>%
                       dplyr::full_join(object@metadata$Groups , by = "samples")
                     
                     #y_profiles.gg$samples <- factor(y_profiles.gg$samples, levels = unique(conds$samples))
                     
                     y_profiles.gg <- dplyr::arrange(y_profiles.gg, get(condition))
                     y_profiles.gg$groups <- factor(y_profiles.gg$groups, levels = unique(y_profiles.gg$groups))
                     
                     p <- ggplot2::ggplot(data = y_profiles.gg, ggplot2::aes(x = groups, y = y_profiles)) +
                       ggplot2::geom_boxplot(ggplot2::aes_string(fill = condition), outlier.size = 0.3) +
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +  
                       ggplot2::xlab("Conditions") + ggplot2::ylab("Expression profiles") +
                       ggplot2::ggtitle(paste0("Cluster: ",numCluster, "; nb_observations : ", dim(assays.data)[1]))
                     
                     if(!is.null(observation)){
                       
                       df <- dplyr::filter(y_profiles.gg, observations == observation) %>% 
                         dplyr::group_by(groups) %>% 
                         dplyr::summarise(mean.y_profiles=mean(y_profiles))
                       p <- p + 
                         geom_point(data = df, aes(x=groups, y=mean.y_profiles), color="red", size = 2) +
                         geom_line( data = df, aes(x=groups, y=mean.y_profiles), color="red", group = 1) +
                         ggtitle(paste0("Cluster: ",numCluster, "; nb_observations : ", dim(assays.data)[1], "; red : ", observation))
                     }
                     -                     
                       return(p)
                   })


methods::setMethod(f          = "coseq.profile.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, numCluster = 1, condition="groups", observation=NULL){
                     
                     RFLOMICS::coseq.profile.plot(object = object[[SE.name]],
                                                  numCluster = numCluster,
                                                  condition = condition,
                                                  observation = observation)
                     
                   })


##### RESETTING ####

#' resetFlomicsMultiAssay
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param results vector of results names
#' @param dataset dataset name. If dataset == NULL, all datasets will be reset
#' @return An object of class \link{MultiAssayExperiment}
#' @export
#' @exportMethod resetFlomicsMultiAssay
#' @examples
#' @noRd
#'
methods::setMethod(f="resetFlomicsMultiAssay", signature="MultiAssayExperiment",
                   
                   definition <- function(object, results, datasets = NULL){
                     
                     # if dataset is null we take all datasets presente in MultiAssayExperiment object
                     if(is.null(datasets)){
                       datasets <- paste0(unlist(object@metadata$omicList), ".filtred")
                     }
                     else{
                       # check if given dataset name include in datasets presente in MultiAssayExperiment object
                       if(!datasets %in% paste0(unlist(object@metadata$omicList), ".filtred")){
                         print("Warning : The given dataset name is not present in MultiAssayExperiment object")
                         return(object)
                       }
                     }
                     
                     for(data in datasets){
                       
                       if(!is.null(object[[data]])){
                         
                         dataset <- object[[data]]
                         
                         for(res in results){
                           if(!is.null(dataset@metadata[[res]])){ dataset@metadata[[res]] <- NULL }
                         }
                         
                         object[[data]] <- dataset
                       }
                       
                     }
                     
                     return(object)
                   })


