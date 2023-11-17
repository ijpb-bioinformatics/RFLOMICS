
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
    stop("ExpDesign matrix is empty!")
  }
  
  # check refList length
  if(length(refList) != length(names(ExpDesign))){
    stop("refList length is different from the dimension of ExpDesign matrix!")
  }
  
  # check typeList length
  if(length(typeList) != length(names(ExpDesign))){
    stop("typeList length is different from the dimension of ExpDesign matrix!")
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



#' @title CheckExpDesign
#' @description This method checks some experimental design characteristics.
#'  A complete design and at least one biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{MultiAssayExperiment-class}
#' @return a named list of two objects
#' \itemize{
#'  \item{"plot:"}{ plot of count data.frame.}
#'  }
#'  
#' @exportMethod CheckExpDesign
#' @noRd

methods::setMethod(f         = "CheckExpDesign",
                   signature = "MultiAssayExperiment",
                   definition <- function(object){
                     
                     Design <- object@metadata$design
                     
                     # check presence of bio factors
                     if (!table(Design@Factors.Type)["Bio"] %in% 1:3){ stop("No bio factor! or nbr of bio factors exceed 3!") }
                     if ( table(Design@Factors.Type)["batch"] == 0){ stop("No replicates found!") }
                     
                     ####################
                     
                     BioFact <- bioFactors(object)
                     coldata <- MultiAssayExperiment::colData(object)
                     coldata[["samples"]] <- rownames(coldata)
                     coldata <- tibble::as_tibble(coldata)
                     coldata <- MultiAssayExperiment::sampleMap(object) %>% tibble::as_tibble() %>% 
                       dplyr::left_join(coldata, by = c("primary" = "samples"))
                     
                     all_combin_cond <- lapply(BioFact, function(x){ 
                       df <- unique(Design@List.Factors[[x]]) %>% as.data.frame()
                       names(df) <- x
                       return(df) 
                     }) %>% purrr::reduce(merge)
                     
                     counts <- coldata %>% dplyr::select(assay, all_of(BioFact)) %>% unique() %>% 
                       dplyr::group_by_at(BioFact) %>% dplyr::count(name = "Count") %>% 
                       dplyr::right_join(all_combin_cond, by = BioFact) %>% 
                       dplyr::mutate_at(.vars = "Count", .funs = function(x) dplyr::if_else(is.na(x), 0, x))
                     
                     ####################                 
                     
                     counts <- counts %>% dplyr::mutate(status = dplyr::if_else(Count == length(object) , "all_data", dplyr::if_else(Count == 0 , "no_data", "some_data")))
                     
                     #list of factor names
                     factors <- names(counts)[1:(dim(counts)[2]-2)]
                     
                     col.panel <- c("all_data", "some_data", "no_data")
                     names(col.panel) <- c("#00BA38", "orange", "red")
                     
                     col.panel.u <- col.panel[col.panel %in% unique(counts$status)]
                     
                     switch (length(factors),
                             "1" = { p <- ggplot2::ggplot(counts ,ggplot2::aes_string(x = factors[1], y = 1)) + 
                               ggplot2::theme(axis.text.y = ggplot2::element_blank()) + ggplot2::ylab("") },
                             "2" = { p <- ggplot2::ggplot(counts ,ggplot2::aes_string(x = factors[1], y = factors[2])) },
                             "3" = {
                               #get factor with min conditions -> to select for "facet_grid"
                               factors.l <- lapply(factors, function(x){ length(unique(counts[[x]])) }) %>% unlist()
                               names(factors.l) <- factors
                               factor.min <- names(factors.l[factors.l == min(factors.l)][1])
                               
                               factors <- factors[factors != factor.min]
                               
                               #add column to rename facet_grid
                               counts <- counts %>% dplyr::mutate(grid = paste0(factor.min, "=",get(factor.min)))
                               
                               p <- ggplot2::ggplot(counts ,ggplot2::aes_string(x = factors[1], y = factors[2])) +
                                 ggplot2::facet_grid(grid~.) })
                     
                     p <- p + ggplot2::geom_tile(ggplot2::aes(fill = status), color = "white", size = 1, width = 1, height = 1)  + ggplot2::geom_text(ggplot2::aes(label = Count)) +
                       ggplot2::scale_fill_manual(values = names(col.panel.u), breaks = col.panel.u) +
                       ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                      axis.ticks = ggplot2::element_blank(), axis.text.x=ggplot2::element_text(angle=90, hjust=1))
                     return(p)
                     
                   })


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
#' @noRd

methods::setMethod(f         = "CheckExpDesignCompleteness",
                   signature = "MultiAssayExperiment",
                   definition <- function(object, datasetList=NULL, sampleList=NULL){
                     
                     # test object type
                     # test if object exist
                     
                     Design <- object@metadata$design
                     
                     if(is.null(datasetList)){ datasetList <- names(object) }
                     
                     
                     # check presence of bio factors
                     if (!table(Design@Factors.Type)["Bio"] %in% 1:3){ stop("no bio factor! or nbr of bio factors exceed 3!") }
                     if ( table(Design@Factors.Type)["batch"] == 0){ stop("No replicate!") }
                     
                     
                     
                     # # count occurence of bio conditions
                     # if(is.null(sampleList)){
                     #   # tmp <- sampleMap(object) %>% data.frame()
                     #   # sampleList <- lapply(unique(tmp$assay), function(dataset){filter(tmp, assay == dataset)$primary }) %>%
                     #   #   purrr::reduce(dplyr::union)
                     #   
                     #   #sampleList <- sampleMap(object)$primary
                     #   sampleList.tmp <- dplyr::group_by(data.frame(MultiAssayExperiment::sampleMap(object)), primary) %>% 
                     #     dplyr::count() %>% 
                     #     dplyr::ungroup() %>% 
                     #     dplyr::mutate(max = max(n)) %>% 
                     #     dplyr::filter(n == max)
                     #   sampleList <- sampleList.tmp$primary
                     # }
                     
                     # Only works with bio and batch factors for the rest of the function
                     namFact <- c(bioFactors(object), batchFactors(object))
                     expDesign_mod <- Design@ExpDesign %>% dplyr::select(tidyselect::any_of(namFact))
                     
                     dF.List <- lapply(1:ncol(expDesign_mod), function(i){
                       factor(expDesign_mod[[i]], levels = unique(expDesign_mod[[i]]))
                     })
                     names(dF.List) <- names(expDesign_mod)
                     
                     # output list
                     output <- list()
                     output[["summary"]] <- data.frame()
                     
                     for(dataset in datasetList){
                       
                       if(is.null(object[[dataset]])) stop(paste(dataset, "not exist."))
                       
                       if(is.null(sampleList)){
                         
                         sampleList_bis <- colnames(object[[dataset]])
                       }
                       else if(length(intersect(sampleList, colnames(object[[dataset]]))) == 0){
                         
                         stop(paste("sampleList values not exist in ", dataset))
                       }
                       else{ sampleList_bis <- sampleList }
                       
                       ExpDesign <- dplyr::filter(expDesign_mod, rownames(expDesign_mod) %in% sampleList_bis)
                       
                       bio.fact <- names(Design@Factors.Type[Design@Factors.Type == "Bio"])
                       #bio.fact <- names(dF.List[Design@Factors.Type == "Bio"])
                       tmp <- ExpDesign %>% dplyr::mutate(samples=row.names(.))
                       
                       group_count <- as.data.frame(dF.List) %>%
                         table() %>%
                         as.data.frame() %>%
                         dplyr::full_join(tmp, by=names(dF.List)) %>%
                         dplyr::mutate_at(.vars = "samples", .funs = function(x) dplyr::if_else(is.na(x), 0, 1)) %>%
                         dplyr::group_by_at((bio.fact)) %>%
                         dplyr::summarise(Count=sum(samples), .groups = "keep")
                       
                       # remplacer le code ci-dessus par celui en bas                       
                       # group_count <- ExpDesign %>% 
                       #   dplyr::group_by_at((bio.fact)) %>% dplyr::count(name = "Count")
                       
                       # check presence of relicat / batch
                       # check if design is complete
                       # check if design is balanced
                       # check nbr of replicats
                       if(min(group_count$Count) == 0){
                         
                         output[["summary"]] <- rbind(output[["summary"]], c(dataset, "error", "The experimental design is not complete."))
                         output[["error"]]   <- TRUE
                       }
                       else if(min(group_count$Count) == 1){
                         
                         output[["summary"]] <- rbind(output[["summary"]], c(dataset, "error", "You need at least 2 biological replicates."))
                         output[["error"]]   <- TRUE
                       }
                       else if(length(unique(group_count$Count)) != 1){
                         
                         output[["summary"]] <- rbind(output[["summary"]], c(dataset, "warning", "The experimental design is complete but not balanced."))
                       }
                       else{
                         output[["summary"]] <- rbind(output[["summary"]], c(dataset, "pass", "The experimental design is complete and balanced."))
                       }
                       
                       ### plot
                       output[[dataset]][["counts"]] <- group_count
                     }
                     names(output[["summary"]]) <- c('dataset', "status", "message")
                     
                     if(length(datasetList) == 1){
                       output[["plot"]]   <- plotExperimentalDesign(counts = output[[datasetList[1]]][["counts"]], 
                                                                    message= output[["summary"]][1,3])
                     }
                     
                     return(output)
                   })



#' @title Datasets overview plot
#' @description This function plot overview of loaded datasets aligned per sample 
#' (n=number of entities (genes/metabolites/proteins); k=number of samples)
#' @param object An object of class \link{MultiAssayExperiment-class}
#' @exportMethod Datasets_overview_plot
#' @return plot

methods::setMethod(f         = "Datasets_overview_plot",
                   signature = "MultiAssayExperiment",
                   definition <- function(object, dataset.list=NULL, real.size=FALSE){
                     
                     if(length(object@ExperimentList) == 0) stop("object@ExperimentList is NULL")
                     if(!is.null(dataset.list)){
                       if(length(intersect(dataset.list, names(object))) == 0){
                         stop(paste0(paste0(dataset.list, collapse = ","), " is not part of dataset names"))
                       }
                       else{
                         object <- object[,, dataset.list]
                       }
                     }
                     
                     nb_entities <- lapply(object@ExperimentList, function(SE){ dim(SE)[1] }) %>% unlist()
                     
                     data <- data.frame(nb_entities = nb_entities, assay = names(nb_entities)) %>%
                       dplyr::full_join(data.frame(MultiAssayExperiment::sampleMap(object)), by="assay") %>%
                       dplyr::mutate(y.axis = paste0(assay, "\n", "n=", nb_entities)) %>% dplyr::arrange(primary)
                     
                     nb_entities_ord <- dplyr::select(data, y.axis, nb_entities) %>% unique() %>% dplyr::arrange(desc(nb_entities))
                     nb_entities_ord$nb_entities <- log(nb_entities_ord$nb_entities)
                     tmp.vec <- c(0)
                     breaks  <- vector()
                     for(i in 1:length(nb_entities_ord$nb_entities)){ 
                       tmp.vec[i+1] <- tmp.vec[i] + nb_entities_ord$nb_entities[i]
                       breaks[i] <- tmp.vec[i] + nb_entities_ord$nb_entities[i]/2 
                     } 
                     
                     switch (as.character(real.size),
                       "TRUE"  = {
                         p <- ggplot2::ggplot(data, ggplot2::aes(x=primary, y=log(nb_entities))) +
                           ggplot2::geom_col(ggplot2::aes(fill = y.axis)) + 
                           ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
                                          panel.background = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
                                          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), legend.position="none",
                                          axis.text.y = ggplot2::element_text(hjust = 0)) +  
                           ggplot2::labs(x=paste0("Samples (k=", length(unique(MultiAssayExperiment::sampleMap(object)$primary)), ")"), y="") +
                           ggplot2::scale_y_continuous(breaks = (breaks), labels = nb_entities_ord$y.axis)
                         
                       },
                       "FALSE" = {
                         p <- ggplot2::ggplot(data, ggplot2::aes(x=primary, y=y.axis)) +
                           ggplot2::geom_tile(ggplot2::aes(fill = y.axis), colour = "grey50") +
                           ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                          panel.background = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position="none",
                                          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
                           ggplot2::labs(x=paste0("Samples (k=", length(unique(MultiAssayExperiment::sampleMap(object)$primary)), ")"), y="")
                       }
                     )
                     return(p)
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
                     if (is(model.formula, "formula")) model.formula <- paste(as.character(model.formula), collapse = " ")
                     modelFormula <- formula(model.formula) 
                     # modelFormula <- formula(paste(model.formula, collapse = " "))
                     # modelFormula <- model.formula
                     
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
                         allInteractionsContrasts_df <- defineAllInteractionContrasts(treatmentFactorsList, groupInteractionToKeep)
                         
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
#' @param object An object of class \link{MultiAssayExperiment-class}
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
                     assignVectorToGroups(treatmentFactorsList    = treatmentFactorsList,
                                          modelMatrix             = modelMatrix,
                                          interactionPresent      = interactionPresent,
                                          isThreeOrderInteraction = isThreeOrderInteraction,
                                          treatmentCondenv        = treatmentCondenv)
                     # get the coefficient vector associated with each selected contrast
                     # contrast <- allSimpleContrast_df$contrast[1]
                     colnamesGLMdesign <- colnames(modelMatrix)
                     
                     #coefficientsMatrix <- sapply(selectedContrasts$contrast, function(x) returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
                     coefficientsMatrix <- sapply(selectedContrasts, function(x)  returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
                     
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
#'
FlomicsMultiAssay.constructor <- function(inputs, projectName, ExpDesign , refList , typeList){
  
  if (sum(duplicated(names(inputs))) > 0) {
    warning("Some names of input are duplicated. Adding suffix.")
    
    names(inputs)[duplicated(names(inputs))] <- paste0(names(inputs)[duplicated(names(inputs))], 1:sum(duplicated(names(inputs))))
    
  }
  
  # consctuct ExpDesign object
  for (i in names(ExpDesign)){
    ExpDesign      <- dplyr::arrange(ExpDesign, get(i))
  }
  
  Design <-  ExpDesign.constructor(ExpDesign = ExpDesign, refList = refList, typeList = typeList)
  
  
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
                                                                     DataProcessing = list(rowSumsZero = genes_flt0,
                                                                                           Filtering = NULL, 
                                                                                           Normalization =  list(setting = list(method = "none"), results = NULL,  normalized = FALSE), 
                                                                                           Transformation = list(setting = list(method = "none"), results = NULL,  transformed = FALSE))
                                                     ))
    
    #SummarizedExperimentList[[dataName]] <- SE[, SE$primary %in% row.names(ExpDesign)]
    
    #### run PCA for raw count
    SummarizedExperimentList[[dataName]] <-  RunPCA(SE, raw = TRUE)
    
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
#' @description This function performs a principal component analysis on omic data stored in an object of class \link{SummarizedExperiment-class}
#' Results are stored in the metadata slot of the same object. If a "Normalization" slot is present in the metadata slot, then data are normalized before running the PCA according to the indicated transform method.
#' @param object An object of class \link{SummarizedExperiment-class}.
#' @param nbcp Number of components to compute. Default is 5.
#' @param raw boolean. Does the pca have to be ran on raw data or transformed and normalized data? Default is FALSE, pca is ran on transformed and normalized data.
#' @return An object of class \link{SummarizedExperiment}
#' @exportMethod RunPCA
#' @importFrom FactoMineR PCA
#' @rdname RunPCA
#' 
methods::setMethod(f          = "RunPCA",
                   signature  = "SummarizedExperiment",
                   definition = function(object, nbcp = 5, raw = FALSE){
                     
                     object2 <- checkTransNorm(object, raw = raw)
                     pseudo  <- SummarizedExperiment::assay(object2)
                     
                     if (raw) object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     else    object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = nbcp, graph = FALSE)
                     
                     return(object)
                     
                   }
)

#' @rdname RunPCA
#' @title RunPCA
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod RunPCA
methods::setMethod(f          = "RunPCA",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, nbcp = 5, raw = FALSE){
                     
                     object[[SE.name]] <-  RunPCA(object[[SE.name]], 
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
#' @rdname Library_size_barplot.plot
#' @noRd
methods::setMethod(f          = "Library_size_barplot.plot",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, raw = FALSE){
                     
                     value    <- NULL
                     warnning <- NULL
                     
                     if (getOmicsTypes(object) != "RNAseq") stop("WARNING: data are not RNAseq!")
                     
                     abundances <- SummarizedExperiment::assay(object)
                     samples    <- colnames(abundances)
                     
                     if (raw) {
                       
                       pseudo <- SummarizedExperiment::assay(object) %>% colSums(., na.rm = TRUE)
                       title  <- "Raw data"
                       
                     }else {
                       
                       # RNAseq, not expected to find any transformation method.
                       if (!isNorm(object)) pseudo <- SummarizedExperiment::assay(apply_norm(object)) 
                       
                       pseudo <- pseudo %>% colSums(., na.rm = TRUE)
                       title <- paste0("Filtered and normalized (", getNormSetting(object)$method, ") data")
                     }
                     
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

#' @rdname Library_size_barplot.plot
#' @noRd
#' @title Library_size_barplot.plot
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod Library_size_barplot.plot
methods::setMethod(f          = "Library_size_barplot.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, raw = FALSE){
                     
                     if (getOmicsTypes(object[[SE.name]]) == "RNAseq") {
                       return( Library_size_barplot.plot(object[[SE.name]], raw = raw))
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
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join arrange
#' @rdname Data_Distribution_plot
#' @noRd

methods::setMethod(
  f = "Data_Distribution_plot",
  signature = "SummarizedExperiment",
  definition = function(object, plot = "boxplot", raw = FALSE) {
    
    object2 <- checkTransNorm(object, raw = raw)
    pseudo <- SummarizedExperiment::assay(object2)
    
    omicsType <- getOmicsTypes(object2)
    
    x_lab <- paste0(omicsType, " data")
    if (omicsType == "RNAseq") {
      x_lab <- paste0("log2(", omicsType, " data)")
    }
    
    # Raw data
    if (raw) {
      title <- paste0(omicsType, " raw data")
    } else {
      # Already normalized or transformed Data
      title <- paste0(omicsType, " data")
      
      if (isTransformed(object2)) {
        title <- paste0("Transformed (", getTransSetting(object2)$method, ") ", title)
      }
      
      if (isNorm(object2)) {
        title <- paste0(title, " - normalization: ", getNormSetting(object2)$method)
      }
      
      if (omicsType == "RNAseq") {
        x_lab <- paste0("log2(", omicsType, " data)")
      }
    }
    
    pseudo.gg <- pseudo %>% reshape2::melt()
    colnames(pseudo.gg) <- c("features", "samples", "value")
    
    pseudo.gg <- pseudo.gg %>%
      dplyr::full_join(object@metadata$Groups, by = "samples") %>%
      dplyr::arrange(groups)
    
    pseudo.gg$samples <- factor(pseudo.gg$samples, levels = unique(pseudo.gg$samples))
    
    switch(plot,
           "density" = {
             p <- ggplot2::ggplot(pseudo.gg) +
               ggplot2::geom_density(ggplot2::aes(x = value, group = samples, color = groups), trim = FALSE) +
               ggplot2::xlab(x_lab) +
               ggplot2::theme(legend.position = "none") +
               ggplot2::ggtitle(title)
           },
           "boxplot" = {
             p <- ggplot2::ggplot(pseudo.gg, ggplot2::aes(x = samples, y = value, label = features)) +
               ggplot2::geom_boxplot(ggplot2::aes(fill = groups), outlier.size = 0.3) +
               ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
               ggplot2::xlab("") +
               ggplot2::ylab(x_lab) +
               ggplot2::ggtitle(title) #+
             # geom_point(alpha = 1/100,size=0)
           }
    )
    
    print(p)
  }
)

#' @rdname Data_Distribution_plot
#' @noRd
#' @title Data_Distribution_plot
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod Data_Distribution_plot
methods::setMethod(
  f = "Data_Distribution_plot",
  signature = "MultiAssayExperiment",
  definition = function(object, SE.name, plot = "boxplot", raw = FALSE) {
    Data_Distribution_plot(
      object = object[[SE.name]],
      plot = plot,
      raw = raw
    )
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
#' @return PCA plot (ggplot2 object)
#' @importFrom dplyr mutate full_join select
#' @importFrom FactoMineR coord.ellipse
#' @importFrom ggplot2 geom_polygon
#' @exportMethod plotPCA
#' @rdname plotPCA
#' 
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
                     
                     omicsType <- getOmicsTypes(object)
                     
                     switch (PCA,
                             "raw"  = {title <- paste0("Raw ", omicsType, " data")},
                             "norm" = {title <- switch (omicsType,
                                                        "RNAseq" = { paste0("Filtred and normalized ", omicsType, 
                                                                            " data (", getNormSetting(object)$method, ")")  },
                                                        "proteomics" = {paste0("Transformed and normalized ", omicsType, 
                                                                               " data (", getTransSetting(object)$method, 
                                                                               " - norm: ", getNormSetting(object)$method, ")")},
                                                        "metabolomics" = {paste0("Transformed and normalized ", omicsType, 
                                                                                 " data (", getTransSetting(object)$method, 
                                                                                 " - norm: ", getNormSetting(object)$method, ")")}
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
                     
                     print(p)
                     
                   })

#' @rdname plotPCA
#' @title plotPCA
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod plotPCA
methods::setMethod(f          = "plotPCA",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, PCA, PCs=c(1,2), condition="groups"){
                     
                     plotPCA(object[[SE.name]], PCA, PCs, condition)
                     
                   })

########################################## TRANSFORM DATA #################

#### METHOD to transform data

#' @title TransformData
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param transformMethod The transformation to store in the metadata or to store and apply if modify_assay is TRUE.
#' @param modify_assay Boolean. Do the transformation need to be applied on the data? The raw data will be replaced by the transformed ones.
#'
#' @return An object of class \link{SummarizedExperiment}
#'
#' @exportMethod TransformData
#' @rdname TransformData
#' 
methods::setMethod(f          = "TransformData",
                   signature  = "SummarizedExperiment",
                   definition = function(object, transformMethod = NULL, modify_assay = FALSE){
                     
                     if (is.null(transformMethod)) {
                       if (!modify_assay && getOmicsTypes(object) == "RNAseq") {
                         message("No transform method indicated, no assay modification asked, RNAseq detected -> transformation_method set to \"none\"")
                         transformMethod <- "none"
                       } else {
                         message("No transform method indicated, using log2 transformation")
                         transformMethod <- "log2"
                       }
                     }
                     
                     if ( getOmicsTypes(object) == "RNAseq" && transformMethod != "none") {
                       message("Transformation other than 'none' are not allowed for RNAseq for now. Forcing none transformation. 
                               Data will be transformed using log2 after the normalization is ran.")
                       transformMethod <- "none"
                     }
                     
                     object@metadata[["DataProcessing"]][["Transformation"]][["setting"]][["method"]] <- transformMethod
                     
                     if (modify_assay) {
                       object <-  apply_transformation(object)
                     }
                     
                     return(object)
                     
                   })

#' @rdname TransformData
#' @title TransformData
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod TransformData
methods::setMethod(f          = "TransformData",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, transformMethod = NULL, modify_assay = FALSE){
                     
                     object[[SE.name]] <-  TransformData(object[[SE.name]], 
                                                         transformMethod = transformMethod, 
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
#' genes the number of sample(s) which are over the cpmCutoff (NbOfsample_over_cpm).
#' Then Two filtering strategies are proposed:
#' \itemize{
#' \item{NbConditions: }{keep gene if the NbOfsample_over_cpm >= NbConditions}
#' \item{NbReplicates: }{keep gene if the NbOfsample_over_cpm >= min(NbReplicat)}
#' \item{filterByExpr:} {the default filtering method implemented in the edgeR filterByExpr() function.}
#' }
#' @param object An object of class \link{SummarizedExperiment}
#' @param filterMethod The filtering model ("CPM")
#' @param filterStrategy The filtering strategy ("NbConditions" or "NbReplicates")
#' @param cpmCutoff The CPM cutoff.
#' @return An object of class \link{SummarizedExperiment}
#' @details
#' Filtered dataset is stored in the ExperimentList slot of the \link{SummarizedExperiment} object
#' as a List named (DataName.filtred).
#' List of filtered features are stored as a named list ("FilteredFeatures") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{SummarizedExperiment} object.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod FilterLowAbundance
#' @importFrom edgeR DGEList filterByExpr cpm
#' @seealso edgeR::filterByExpr#' 
#' @rdname FilterLowAbundance


methods::setMethod(f         = "FilterLowAbundance",
                   signature = "SummarizedExperiment",
                   definition <- function(object, filterMethod= "CPM", filterStrategy = "NbConditions", cpmCutoff = 5){

                     objectFilt <- object
                     assayFilt  <- assay(objectFilt)
                     
                     ## nbr of genes with 0 count
                     genes_flt0  <- objectFilt[rowSums(assayFilt) <= 0, ]@NAMES
                     
                     ## remove 0 count
                     objectFilt  <- objectFilt[rowSums(assayFilt)  > 0, ]
                     assayFilt   <- assay(objectFilt)
                     
                     ## filter cpm
                     NbReplicate  <- table(object@metadata$Groups$groups)
                     NbConditions <- length(unique(object@metadata$Groups$groups))
                     
                     switch(filterStrategy,
                            "NbConditions" = { keep <- rowSums(edgeR::cpm(assayFilt) >= cpmCutoff) >=  NbConditions },
                            "NbReplicates" = { keep <- rowSums(edgeR::cpm(assayFilt) >= cpmCutoff) >=  min(NbReplicate) },
                            "filterByExpr" = { dge  <- edgeR::DGEList(counts = assayFilt, genes = rownames(assayFilt))
                            #keep <- filterByExpr(dge, GLM_Model)
                            keep <- edgeR::filterByExpr(dge)
                            }
                     )
                     
                     ## nbr of genes filtered
                     genes_flt1  <- objectFilt[!keep]@NAMES
                     
                     #objectFilt@metadata[["FilteredFeatures"]] <-  c(genes_flt0, genes_flt1)
                     
                     object <- objectFilt[keep]
                     # object@metadata$FilteredOptions <- list()
                     # object@metadata$FilteringOptions[["filterStrategy"]] <- filterStrategy
                     # object@metadata$FilteringOptions[["cpmCutoff"]] <- cpmCutoff
                     
                     Filtering <- list()
                     Filtering[["setting"]][["method"]]           <- filterMethod
                     Filtering[["setting"]][["filterStrategy"]]   <- filterStrategy
                     Filtering[["setting"]][["cpmCutoff"]]        <- cpmCutoff
                     Filtering[["results"]][["filteredFeatures"]] <- c(genes_flt0, genes_flt1)
                     #Filtering[["results"]][["filteredSamples"]] 
                     
                     if(is.null(object@metadata$DataProcessing)) object@metadata$DataProcessing <- list()
                     object@metadata$DataProcessing$Filtering <- Filtering
                     
                     return(object)

                   })


#' @rdname FilterLowAbundance
#' @title FilterLowAbundance
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod FilterLowAbundance
methods::setMethod(f          = "FilterLowAbundance",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, filterStrategy = "NbConditions", cpmCutoff = 5){
                     
                     if ( getOmicsTypes(object[[SE.name]]) == "RNAseq") {
                       object[[SE.name]] <-  FilterLowAbundance(object          = object[[SE.name]], 
                                                                filterStrategy = filterStrategy, 
                                                                cpmCutoff      = cpmCutoff)
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
#' @param NormMethod Normalization method
#' @param modify_assay Does the normalization have to be applied or just stored for later? Recommended it stays FALSE.
#' @return An object of class \link{SummarizedExperiment}
#' The applied normalization method and computed scaling factors (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "methode" and "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{SummarizedExperiment} object.
#' @exportMethod RunNormalization
#' @seealso TMM.Normalization
#' @rdname RunNormalization
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).

methods::setMethod(f          = "RunNormalization",
                   signature  = "SummarizedExperiment",
                   definition = function(object, NormMethod = NULL, modify_assay = FALSE){
                     
                     if (isNorm(object)) 
                       warning("Data were already normalized before!")
                     
                     if (isTransformed(object)) 
                       warning("Data were transformed before!")
                     
                     if (is.null(NormMethod)) {
                       if (getOmicsTypes(object) == "RNAseq" && getTransSetting(object)$method == "none") {
                         message("Using TMM normalization for RNAseq (counts) data")
                         NormMethod <- "TMM"
                       } else {
                         NormMethod <- "none"
                       }
                     }
                     
                     if (getOmicsTypes(object) == "RNAseq" && NormMethod != "TMM") {
                       message("Forcing TMM normalization for RNAseq (counts) data")
                       NormMethod <- "TMM"
                     }
                     
                     object2 <- object
                     
                     # Run normalization on transformed data (except for RNAseq data, transformation is expected to be "none" anyway)
                     if (!isTransformed(object) && getTransSetting(object)$method != "none")
                       object2 <-  apply_transformation(object2)
                     
                     switch(NormMethod,
                            "TMM"        = {coefNorm  <- TMM.Normalization(assay(object2), object2@metadata$Groups$groups) },
                            "median"     = {coefNorm  <- apply(assay(object2), 2, FUN = function(sample_vect) {median(sample_vect)}) },
                            "totalSum"   = {coefNorm  <- apply(assay(object2), 2, FUN = function(sample_vect) {sum(sample_vect^2)}) },
                            "none"       = {coefNorm  <- rep(1, ncol(assay(object2))) },
                            { message("Could not recognize the normalization method, applying 'none'. Please check your parameters.")
                              NormMethod <- "none"
                              coefNorm  <-  rep(1, ncol(assay(object2)))
                            }
                     )
                     
                     # object <- setNorm(object, NormMethod)
                     # object <- setCoeffNorm(object, coefNorm)
                     
                     Normalization <- list()
                     Normalization[["setting"]][["method"]]   <- NormMethod
                     Normalization[["results"]][["coefNorm"]] <- coefNorm
                     Normalization[["normalized"]]            <- FALSE

                     if(is.null(object@metadata$DataProcessing)) object@metadata$DataProcessing <- list()
                     object@metadata$DataProcessing[["Normalization"]] <- Normalization
                     
                     if (modify_assay) object <-  apply_norm(object)

                     return(object)
                   })

#' @rdname RunNormalization
#' @title RunNormalization
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod RunNormalization
methods::setMethod(f          = "RunNormalization",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, NormMethod, modify_assay = FALSE){
                     
                     object[[SE.name]] <-  RunNormalization(object       = object[[SE.name]],
                                                            NormMethod   = NormMethod,
                                                            modify_assay = modify_assay)
                     
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
#' @param object an object of class \link{SummarizedExperiment} or \link{MultiAssayExperiment} 
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment 
#' @param design an object of class \link{ExpDesign-class}
#' @param DiffAnalysisMethod A character vector giving the name of the differential analysis method
#' to run. Either "edgeRglmfit" or "limmalmFit".
#' @param contrastList The list of contrast to test
#' @param Adj.pvalue.method The method choosen to adjust pvalue. Takes the same values as the ones of adj.p.adjust method.
#' @param Adj.pvalue.cutoff The adjusted pvalue cut-off
#' @param clustermq A boolean indicating whether the constrasts have to be computed in local or in a distant machine
#' @return An object of class \link{SummarizedExperiment}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod RunDiffAnalysis
#' @importFrom dplyr filter
#' @rdname RunDiffAnalysis
#' 
methods::setMethod(f         = "RunDiffAnalysis",
                   signature = "SummarizedExperiment",
                   definition <- function(object, design, Adj.pvalue.method="BH",
                                          contrastList, DiffAnalysisMethod, 
                                          Adj.pvalue.cutoff=0.05, logFC.cutoff=0, clustermq=FALSE, parallel = FALSE, nworkers = 1,
                                          cmd = FALSE){
                     
                     contrastName <- NULL
                     
                     # TODO if contrastList missing, takes contrasts.Sel
                     Contrasts.Sel <- dplyr::filter(design@Contrasts.Sel, contrastName %in% contrastList)
                     
                     object@metadata$DiffExpAnal <- list()
                     object@metadata$DiffExpAnal[["contrasts"]] <- Contrasts.Sel
                     # object@metadata$DiffExpAnal[["method"]]    <- DiffAnalysisMethod
                     # object@metadata$DiffExpAnal[["Adj.pvalue.method"]]  <- Adj.pvalue.method
                     # #object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]  <- Adj.pvalue.cutoff
                     
                     # remplacera à terme les lignes ci-dessus
                     object@metadata$DiffExpAnal[["setting"]][["method"]] <- DiffAnalysisMethod
                     object@metadata$DiffExpAnal[["setting"]][["Adj.pvalue.method"]] <- Adj.pvalue.method
                     object@metadata$DiffExpAnal[["setting"]][["Adj.pvalue.cutoff"]] <- Adj.pvalue.cutoff
                     object@metadata$DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
                     
                     # transform and norm if needed
                     if (DiffAnalysisMethod == "limmalmFit") {
                       
                       object2 <- object
                       
                       if (!isTransformed(object2)) object2 <-  apply_transformation(object2)
                       if (!isNorm(object2))        object2 <-  apply_norm(object2)
                       
                     }
                     
                     # move in ExpDesign Constructor
                     model_matrix <- model.matrix(as.formula(paste(design@Model.formula, collapse = " ")), data = as.data.frame(design@List.Factors))
                     # model_matrix <- model.matrix(as.formula(paste(design@Model.formula, collapse = " ")), data = as.data.frame(design@ExpDesign))
                     rownames(model_matrix) <- rownames(design@ExpDesign)
                     
                     ListRes <- switch(DiffAnalysisMethod,
                                       "edgeRglmfit" = try_rflomics(edgeR.AnaDiff(count_matrix  = SummarizedExperiment::assay(object),
                                                                                  model_matrix    = model_matrix[colnames(object),],
                                                                                  group           = getCoeffNorm(object)$group,
                                                                                  lib.size        = getCoeffNorm(object)$lib.size,
                                                                                  norm.factors    = getCoeffNorm(object)$norm.factors,
                                                                                  Contrasts.Sel   = object@metadata$DiffExpAnal[["contrasts"]],
                                                                                  Contrasts.Coeff = design@Contrasts.Coeff,
                                                                                  FDR             = 1,
                                                                                  clustermq       = clustermq,
                                                                                  parallel        = parallel,
                                                                                  nworkers        = nworkers,
                                                                                  cmd             = cmd)),
                                       "limmalmFit" = try_rflomics(limma.AnaDiff(count_matrix    = SummarizedExperiment::assay(object2),
                                                                                 model_matrix      = model_matrix[colnames(object2),],
                                                                                 Contrasts.Sel     = object2@metadata$DiffExpAnal[["contrasts"]],
                                                                                 Contrasts.Coeff   = design@Contrasts.Coeff,
                                                                                 Adj.pvalue.cutoff = 1,
                                                                                 Adj.pvalue.method = Adj.pvalue.method,
                                                                                 clustermq         = clustermq,
                                                                                 cmd               = cmd)))
                     
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
                     object <-  FilterDiffAnalysis(object = object, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff)
                     
                     return(object)
                   })


#' @rdname RunDiffAnalysis
#' @title RunDiffAnalysis
#' @exportMethod RunDiffAnalysis
methods::setMethod(f          = "RunDiffAnalysis",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, design = NULL, Adj.pvalue.method="BH",
                                         contrastList = NULL, DiffAnalysisMethod = NULL,
                                         Adj.pvalue.cutoff=0.05, logFC.cutoff=0, clustermq=FALSE, 
                                         parallel = FALSE, nworkers = 1, cmd = FALSE){
                     
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
                       
                       switch( getOmicsTypes(object[[SE.name]]),
                               "RNAseq" = {DiffAnalysisMethod = "edgeRglmfit"
                               message("DiffAnalyseMethod was missing. Detected omic type is RNASeq, using edgeRglmFit for differential analysis.")},
                               "proteomics" = {DiffAnalysisMethod = "limmalmFit"
                               message("DiffAnalyseMethod was missing. Detected omic type is proteomics, using limmalmFit for differential analysis.")},
                               "metabolomics" = {DiffAnalysisMethod = "limmalmFit"
                               message("DiffAnalyseMethod was missing. Detected omic type is metabolomics, using limmalmFit for differential analysis.")}
                       )
                       
                     }
                     
                     object[[SE.name]] <-  RunDiffAnalysis(object = object[[SE.name]],
                                                           design = design, 
                                                           Adj.pvalue.method = Adj.pvalue.method,
                                                           contrastList = contrastList,
                                                           DiffAnalysisMethod = DiffAnalysisMethod,
                                                           Adj.pvalue.cutoff = Adj.pvalue.cutoff,
                                                           logFC.cutoff = logFC.cutoff,
                                                           clustermq = clustermq,
                                                           parallel = parallel,
                                                           nworkers = nworkers,
                                                           cmd = cmd
                     )
                     return(object)
                   })

## METHOD to filter differential analysis

#' Filter differential analysis
#'
#' @param object A SummarizedExperiment object
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment
#' @param Adj.pvalue.cutoff adjusted pvalue cutoff. Default is the parameter from the differential analysis.
#' @param logFC.cutoff cutoff for absolute value of log2FC. Default is the parameter from the differential analysis. 
#'
#' @return A SummarizedExperiment object or a MultiAssayExperiment, depending on the object type, 
#' where the differential analysis results have been actualized with the new parameters.
#' @exportMethod FilterDiffAnalysis
#' @rdname FilterDiffAnalysis
#' @importFrom dplyr filter if_else mutate_at
#' @importFrom data.table data.table
#' @importFrom purrr reduce
#'
methods::setMethod(f          = "FilterDiffAnalysis",
                   signature  = "SummarizedExperiment",
                   definition <- function(object, Adj.pvalue.cutoff = NULL, logFC.cutoff = NULL){

                     if(is.null(object@metadata$DiffExpAnal[["RawDEFres"]])){
                       stop("can't filter the DiffExpAnal object because it doesn't exist")
                     }
                     
                     if (is.null(Adj.pvalue.cutoff)) 
                       Adj.pvalue.cutoff <- getDiffSetting(object)$Adj.pvalue.cutoff
                     
                     if (is.null(logFC.cutoff))
                       logFC.cutoff <- getDiffSetting(object)$abs.logFC.cutoff
                     
                     # object@metadata$DiffExpAnal[["Adj.pvalue.cutoff"]]  <- Adj.pvalue.cutoff
                     # object@metadata$DiffExpAnal[["abs.logFC.cutoff"]]  <- logFC.cutoff
                     
                     # remplacera à terme les lignes ci-dessus
                     object@metadata$DiffExpAnal[["setting"]][["Adj.pvalue.cutoff"]] <- Adj.pvalue.cutoff
                     object@metadata$DiffExpAnal[["setting"]][["abs.logFC.cutoff"]]  <- logFC.cutoff
                     
                     ## TopDEF: Top differential expressed features
                     DEF_filtred <- lapply(1:length(object@metadata$DiffExpAnal[["DEF"]]), function(x){
                       res <- object@metadata$DiffExpAnal[["DEF"]][[x]]
                       keep <- (res$Adj.pvalue < Adj.pvalue.cutoff) & (abs(res$logFC) > logFC.cutoff)
                       res <- res[keep,]
                       return(res)
                     })
                     names(DEF_filtred) <- names(object@metadata$DiffExpAnal[["RawDEFres"]])
                     object@metadata$DiffExpAnal[["TopDEF"]] <- DEF_filtred
                     
                     ## stats
                     object@metadata$DiffExpAnal[["stats"]] <- sumDiffExp(object)
                     
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

#' @rdname FilterDiffAnalysis
#' @title FilterDiffAnalysis
#' @exportMethod FilterDiffAnalysis
methods::setMethod(f          = "FilterDiffAnalysis",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, 
                                         Adj.pvalue.cutoff = NULL, logFC.cutoff = NULL){
                     
                     if(!SE.name %in% names(object))
                       stop(SE.name, " isn't the name of an experiment in ", object)
                     
                     if (is.null(Adj.pvalue.cutoff) && is.null(logFC.cutoff)) {
                       
                       message("Parameter Adj.pvalue.cutoff and logFC.cutoff are both NULL. Not changing anything")
                       return(object)
                       
                     }else{
                       
                       if (is.null(Adj.pvalue.cutoff)) 
                         Adj.pvalue.cutoff <- getDiffSetting(object[[SE.name]])$Adj.pvalue.cutoff
                       
                       if (is.null(logFC.cutoff))
                         logFC.cutoff <- getDiffSetting(object[[SE.name]])$abs.logFC.cutoff
                       
                       object[[SE.name]] <-  FilterDiffAnalysis(object = object[[SE.name]],
                                                                Adj.pvalue.cutoff = Adj.pvalue.cutoff,
                                                                logFC.cutoff = logFC.cutoff)
                       
                       return(object)
                     }
                     
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
#' @rdname DiffAnal.plot
#' @export
#' 
methods::setMethod(f="DiffAnal.plot",
                   signature="SummarizedExperiment",
                   
                   definition <- function(object, hypothesis, typeofplots = c("MA.plot", "volcano", "histogram")){
                     
                     if ( isTagName(object, hypothesis)) hypothesis <-  convertTagToContrast(object, hypothesis)
                     
                     plots <- list()
                     
                     res      <- object@metadata$DiffExpAnal[["RawDEFres"]][[hypothesis]]
                     resTable <- object@metadata$DiffExpAnal[["DEF"]][[hypothesis]]
                     
                     logFC.cutoff      <- getDiffSetting(object)[["abs.logFC.cutoff"]]
                     Adj.pvalue.cutoff <- getDiffSetting(object)[["Adj.pvalue.cutoff"]]
                     
                     if ("MA.plot" %in% typeofplots) plots[["MA.plot"]]        <-  MA.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff, hypothesis=hypothesis)
                     if ("volcano" %in% typeofplots) plots[["Volcano.plot"]]   <-  Volcano.plot(data = resTable, Adj.pvalue.cutoff = Adj.pvalue.cutoff, logFC.cutoff = logFC.cutoff, hypothesis=hypothesis)
                     if ("histogram" %in% typeofplots) plots[["Pvalue.hist"]]  <-  pvalue.plot(data =resTable, hypothesis=hypothesis)
                     return(plots)
                   })

#' @rdname DiffAnal.plot
#' @title DiffAnal.plot
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment
#' @exportMethod DiffAnal.plot
methods::setMethod(f          = "DiffAnal.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, hypothesis, typeofplots = c("MA.plot", "volcano", "histogram")){
                     
                     if (isTagName(object, hypothesis)) hypothesis <-  convertTagToContrast(object, hypothesis)
                     
                     return(DiffAnal.plot(object      = object[[SE.name]],
                                          hypothesis  = hypothesis,
                                          typeofplots = typeofplots))
                     
                   })


#' @title heatmapPlot
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
#' @exportMethod heatmapPlot
#' @export
#' @importFrom dplyr arrange select
#' @importFrom tidyselect any_of
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importClassesFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importMethodsFrom ComplexHeatmap draw
#' @importFrom grid gpar
#' @rdname heatmapPlot
#' 
methods::setMethod(f          = "heatmapPlot",
                   signature  = "SummarizedExperiment",
                   definition = function(object, 
                                         hypothesis, 
                                         condition="none", 
                                         title = "", 
                                         annot_to_show = NULL, 
                                         subset_list = NULL, 
                                         draw_args = list(), 
                                         heatmap_args = list()){
                     
                     if (is.null(object@metadata$DiffExpAnal[["TopDEF"]][[hypothesis]])) {
                       stop("no DE variables")
                     }
                     
                     resTable <- dplyr::arrange(object@metadata$DiffExpAnal[["TopDEF"]][[hypothesis]], Adj.pvalue)
                     
                     if (dim(resTable)[1] == 0) {
                       stop("no differentially expressed variables...")
                     }
                     
                     
                     if (dim(resTable)[1] > 2000) {
                       message("differentially expressed variables exceeding 2000 variables, only the first 2000 will be displayed")
                       resTable <- resTable[1:2000,]
                       title = ifelse(title == "", paste0(title, "plot only 2000 TOP DE variables"),
                                      paste0(title, "\nplot only 2000 TOP DE variables"))
                     }
                     
                     object2 <- checkTransNorm(object, raw = FALSE)
                     m.def  <- assay(object2)
                     
                     m.def <- as.data.frame(m.def) %>%
                       dplyr::select(tidyselect::any_of(object2@metadata$Groups$samples))
                     
                     # filter by DE
                     m.def.filter <- subset(m.def, rownames(m.def) %in% row.names(resTable))
                     
                     # normalize count
                     
                     # Center
                     m.def.filter.center <- t(scale(t(m.def.filter), center = TRUE, scale = FALSE))
                     
                     # Annotations datatable
                     df_annotation <- object@metadata$Groups %>% dplyr::select(!samples & !groups)  
                     df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),] 
                     
                     # Subset the dataset to print only interesting modalities
                     if (!is.null(subset_list)) {
                       if (is.null(names(subset_list))) {
                         message("In heatmapPlot, subset_list argument needs a named list. Not subsetting")
                       }else{ 
                         samplesToKeep <- Reduce("intersect", lapply(
                           1:length(subset_list),
                           FUN = function(i){
                             col_nam <- names(subset_list)[i]
                             rownames(df_annotation[which(df_annotation[[col_nam]] %in% subset_list[[i]]),])
                           }
                         ))
                         
                         df_annotation <- df_annotation[which(rownames(df_annotation) %in% samplesToKeep),]
                         m.def.filter.center <- m.def.filter.center[, which(colnames(m.def.filter.center) %in% samplesToKeep)]
                         
                         df_annotation <- df_annotation[match(colnames(m.def.filter.center), rownames(df_annotation)),]
                       }
                     }
                     
                     # Split management
                     column_split.value <- if (condition != "none") { df_annotation[, condition] } else { NULL }
                     
                     # Select the right columns
                     if (!is.null(annot_to_show)) {
                       df_annotation <- df_annotation %>% dplyr::select(tidyselect::any_of(annot_to_show))
                     }
                     
                     # Color annotations
                     set.seed(10000) # seed for chosing palette. Ensure running twice the function get the same annotations colors.
                     selectPal <- sample(rownames(RColorBrewer::brewer.pal.info),  size = ncol(df_annotation), replace = FALSE)
                     
                     color_list <- lapply(1:ncol(df_annotation), FUN = function(i){
                       annot_vect <- unique(df_annotation[,i])
                       
                       col_vect <-  grDevices::colorRampPalette(
                         suppressWarnings({ RColorBrewer::brewer.pal(n = min(length(annot_vect), 8), name = selectPal[i])}))(length(annot_vect)) 
                       names(col_vect) <- annot_vect 
                       col_vect[!is.na(names(col_vect))] # RcolorBrewer::brewer.pal n is minimum 3, remove NA names if only 2 levels
                     })
                     names(color_list) <- colnames(df_annotation)
                     
                     column_ha <- ComplexHeatmap::HeatmapAnnotation(df = df_annotation, col = color_list)
                     
                     namArg <- ifelse(getOmicsTypes(object) == "RNAseq", "normalized counts", "XIC")
                     
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
                     
                     dev.off()
                     
                     return(ha)
                   })


#' @rdname heatmapPlot
#' @title heatmapPlot
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment
#' @exportMethod heatmapPlot
methods::setMethod(f          = "heatmapPlot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, hypothesis, condition="none", title = "", annot_to_show = NULL, subset_list = NULL, draw_args = list(), heatmap_args = list()){
                     
                     
                     if (isTagName(object, hypothesis)) hypothesis <- convertTagToContrast(object, hypothesis)
                     
                     return(heatmapPlot(object        = object[[SE.name]],
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
#' @importFrom ggplot2 geom_density xlab theme_void ggtitle ggplot geom_boxplot guide_legend guides
#' @importFrom dplyr full_join  arrange
#' @noRd

methods::setMethod(f          = "boxplot.DE.plot",
                   signature  = "SummarizedExperiment",
                   definition = function(object, DE = NULL, condition="groups", raw = FALSE){
                     
                     # check variable name
                     if (is.null(DE) | DE == "" | length(DE) != 1) {
                       message("set variable name")
                       
                       p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("set variable name") 
                       
                       return(p)
                     }
                     
                     object <- checkTransNorm(object, raw = raw)
                     
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
                         pseudo <- SummarizedExperiment::assay(object.DE)
                         x_lab  <- paste0(DE, " data")
                         
                         if (isTransformed(object.DE) && getTransSetting(object.DE)$method != "none") {
                           title  <- paste0("Transformed (", getTransSetting(object.DE)$method, ") ", title)
                         }
                         if (isNorm(object.DE) && getNormSetting(object.DE)$method != "none") {
                           title <- paste0(title, " - normalization: ", getNormSetting(object.DE)$method)
                         }  
                       } else {
                         
                         pseudo <- SummarizedExperiment::assay(object.DE) 
                         title  <- DE
                         x_lab  <- paste0("log2(", DE, " data)") 
                         
                       }
                     }
                     
                     pseudo.gg <- pseudo %>% reshape2::melt()
                     colnames(pseudo.gg) <- c("features", "samples", "value")
                     
                     pseudo.gg <- pseudo.gg %>% dplyr::full_join(object@metadata$Groups, by="samples") %>%
                       dplyr::arrange(groups)
                     
                     pseudo.gg <- dplyr::arrange(pseudo.gg, get(condition))
                     
                     pseudo.gg$groups <- factor(pseudo.gg$groups, levels = unique(pseudo.gg$groups))
                     
                     p <- ggplot2::ggplot(pseudo.gg, ggplot2::aes(x=groups, y=value, label = features)) +
                       ggplot2::geom_boxplot(ggplot2::aes(fill=get(condition))) +
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
                       ggplot2::guides(fill=guide_legend(title="condition")) + 
                       ggplot2::xlab("") + 
                       ggplot2::ylab(x_lab) + 
                       ggtitle(title) #+
                     #geom_point(alpha = 1/100,size=0)
                     
                     return(p)
                     
                   }
)

#' @rdname boxplot.DE.plot
#' @title boxplot.DE.plot
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment
#' @exportMethod boxplot.DE.plot
methods::setMethod(f          = "boxplot.DE.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, DE = NULL, condition="groups"){
                     
                     boxplot.DE.plot(object = object[[SE.name]], 
                                     DE = DE,
                                     condition = condition)
                     
                   })


################################### CO-EXPRESSION #############################



#' @title runCoExpression
#' @description This is an interface method which performs co-expression/co-abundance analysis
#' of omic-data.
#' @details For now, only the coseq function of the coseq package is used.
#' For RNAseq data, parameters used are those recommended in DiCoExpress workflow (see the reference).
#' This parameters are: \code{model="normal"}, \code{transformation="arcsin"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="TMM"}, \code{meanFilterCutoff = 50}
#' For proteomic or metabolomic, data are scaled by protein or metabolite to group them by expression
#' profiles rather than by expression intensity.
#' After data scaling, recommended parameters (from \code{coseq} developers) for co-expression analysis are:
#' \code{model="normal"}, \code{transformation="none"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="none"},  \code{meanFilterCutoff = NULL}.
#'
#' @return
#' An S4 object of class \link{SummarizedExperiment}
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
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment
#' @param nameList names of the contrasts from which the DE entities are taken. Can be NULL, in that case every contrasts from the differential analysis is taken into consideration.
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
#' @param clustermq_arg boolean. Does the computation need to be executed on a distant server?
#' @param silent if TRUE, coseq run silently (without any console print or message)
#' @param cmd if TRUE, print steps of the analysis. Used inside the coseq module in the shiny interface.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runCoExpression
#' @seealso \code{\link{coseq::coseq}}
#' @rdname runCoExpression
#' 
methods::setMethod(f = "runCoExpression",
                   signature = "SummarizedExperiment",
                   definition = function(object,
                                          K=2:20, 
                                          replicates=5, 
                                          nameList = NULL, 
                                          merge="union",
                                          model = "Normal",
                                          GaussianModel = NULL, 
                                          transformation = NULL, 
                                          normFactors = NULL, 
                                          clustermq_arg = FALSE,
                                          meanFilterCutoff = NULL, 
                                          scale = NULL,
                                          silent = TRUE, 
                                          cmd = FALSE){
                    
                     if (is.null(object@metadata$DiffExpAnal[["mergeDEF"]]))
                       stop("Please run a differential analysis. runCoExpression uses these results.")
                     
                     if (is.null(nameList) && !is.null(getValidContrasts(object)[["tag"]])) 
                       nameList <- getValidContrasts(object)[["tag"]]
                     else if (is.null(nameList) && is.null(getValidContrasts(object)[["tag"]])) 
                       nameList <- colnames(object@metadata$DiffExpAnal[["mergeDEF"]])[-1]
                     
                     CoExpAnal <- list()
                     
                     CoExpAnal[["setting"]][["method"]]           <- "coseq"
                     CoExpAnal[["setting"]][["gene.list.names"]]  <- nameList
                     names(CoExpAnal[["setting"]][["gene.list.names"]])  <- dplyr::filter(object@metadata$DiffExpAnal$contrasts, tag %in% nameList)$contrastName
                     CoExpAnal[["setting"]][["merge.type"]]       <- merge
                     CoExpAnal[["setting"]][["replicates.nb"]]    <- replicates
                     CoExpAnal[["setting"]][["K.range"]]          <- K
                     CoExpAnal[["setting"]][["scale"]]            <- scale
                     
                     geneList <- opDEList(object = object, contrasts = nameList, operation = merge)
                     
                     # set default parameters based on data type
                     param.list <- list("model" = model)
                     
                     switch(object@metadata$omicType,
                            
                            "RNAseq" = {
                              counts <- SummarizedExperiment::assay(object)[geneList,]
                              
                              param.list[["transformation"]]   <- ifelse(is.null(transformation), "arcsin", transformation)
                              param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "TMM", normFactors)
                              param.list[["meanFilterCutoff"]] <- ifelse(is.null(meanFilterCutoff), 50, meanFilterCutoff)
                              param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Ck", GaussianModel)
                              
                            },
                            "proteomics" = {
                              object <- checkTransNorm(object)
                              counts <- SummarizedExperiment::assay(object)[geneList,]
                              
                              # Print the selected GaussianModel
                              if (cmd) print(paste("Use ", GaussianModel, sep = ""))
                              if (cmd) print("Scale each protein (center = TRUE, scale = TRUE)")
                              CoExpAnal[["transformation.prot"]] <- "scaleProt"
                              counts[] <- t(apply(counts,1,function(x){scale(x, center = TRUE, scale = TRUE) }))
                              
                              # param
                              param.list[["transformation"]]   <- ifelse(is.null(transformation), "none", transformation)
                              param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "none", normFactors)
                              param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Bk", GaussianModel)
                            },
                            "metabolomics" = {
                              object <- checkTransNorm(object)
                              counts <- SummarizedExperiment::assay(object)[geneList,]
                              
                              # Print the selected GaussianModel
                              if (cmd) print(paste("Use ", GaussianModel, sep= ""))
                              if (cmd) print("Scale each metabolite (center = TRUE,scale = TRUE)")
                              CoExpAnal[["transformation.metabo"]] <- "scaleMetabo"
                              counts[] <- t(apply(counts,1,function(x){ scale(x, center = TRUE, scale = TRUE) }))
                              
                              # param
                              param.list[["transformation"]]   <- ifelse(is.null(transformation), "none", transformation)
                              param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "none", normFactors)
                              param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Bk", GaussianModel)
                            }
                     )
                     
                     CoExpAnal[["setting"]] <- c(CoExpAnal[["setting"]], param.list)
                     
                     # run coseq : on local machine or remote cluster
                     
                     if (cmd) print("#     => coseq... ")
                     
                     conds <- object@metadata$Groups
                     counts <- counts[, match(rownames(conds), colnames(counts))]
                     if (!identical(colnames(counts), rownames(conds), attrib.as.set = FALSE)) {
                       stop("colnames counts and rownames conds don't match!")
                     }
                     conds <- conds$groups
                     
                     coseq.res.list <- list()
                     
                     coseq.res.list <- switch(as.character(clustermq_arg),
                                              `FALSE` = {
                                                try_rflomics(
                                                  runCoseq_local(counts, 
                                                                 conds = object@metadata$Groups$groups,
                                                                 K = K, 
                                                                 replicates = replicates, 
                                                                 param.list = param.list,
                                                                 silent = silent,
                                                                 cmd = cmd))
                                              },
                                              `TRUE` = {
                                                try_rflomics(
                                                  runCoseq_clustermq(counts, 
                                                                     conds = object@metadata$Groups$groups,
                                                                     K = K, 
                                                                     replicates = replicates, 
                                                                     param.list = param.list,
                                                                     silent = silent, 
                                                                     cmd = cmd))
                                                
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

#' @rdname runCoExpression
#' @title runCoExpression
#' @exportMethod runCoExpression
methods::setMethod(f          = "runCoExpression",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, K=2:20, replicates=5, nameList, merge="union",
                                         model = "Normal", GaussianModel = NULL, 
                                         transformation = NULL, normFactors = NULL, clustermq_arg = FALSE,
                                         meanFilterCutoff = NULL, scale = NULL, silent = TRUE, cmd = FALSE){
                     
                     
                     if (!SE.name %in% names(object)) 
                       stop(paste0(SE.name, " is not part of ", object))
                     
                     object[[SE.name]] <-  runCoExpression(object = object[[SE.name]],
                                                           K = K,
                                                           replicates = replicates,
                                                           nameList = nameList, 
                                                           merge = merge, 
                                                           model = model,
                                                           GaussianModel = GaussianModel,
                                                           transformation = transformation,
                                                           normFactors = normFactors,
                                                           clustermq_arg = clustermq_arg,
                                                           meanFilterCutoff = meanFilterCutoff,
                                                           scale = scale,
                                                           silent = silent,
                                                           cmd = cmd)
                     
                     return(object)
                     
                   })

# Pour utiliser la fonction repeatable(), "seed"  pourrait être ajouté en paramètre.




#' @title CoExpressionPlots
#' 
#' @param object An object of class \link{SummarizedExperiment}
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @importFrom coseq plot
#' @importFrom ggplot2 ggplot geom_boxplot geom_text
#' @export
#' @exportMethod CoExpressionPlots
#' @noRd
#' 
coExpressionPlots <- methods::setMethod(f="CoExpressionPlots",
                                        signature="SummarizedExperiment",
                                        definition <- function(object){
                                          
                                          if(is.null(object@metadata$CoExpAnal) || length(object@metadata$CoExpAnal) == 0) stop("No co-expression results!")
                                          CoExpAnal <- object@metadata$CoExpAnal
                                          
                                          coseq.res     <- CoExpAnal[["coseqResults"]]
                                          ICL.list      <- CoExpAnal[["plots"]][["ICL"]] 
                                          logLike.list  <- CoExpAnal[["plots"]][["logLike"]]
                                          K             <- CoExpAnal[["K.range"]]
                                          conds         <- object@metadata$Groups$groups
                                          
                                          #### Plots
                                          ### plot ICL
                                          ICL.p   <- ggplot2::ggplot(data = ICL.list[["ICL.tab"]]) +
                                            ggplot2::geom_boxplot(ggplot2::aes(x = as.factor(K), y = ICL, group = K)) +
                                            ggplot2::geom_text(data = ICL.list[["ICL.n"]], ggplot2::aes(x = 1:length(K), y = max(ICL.list[["ICL.vec"]], na.rm = TRUE),
                                                                                                        label = paste0("n=", n)), col = 'red', size = 4) +
                                            ggplot2::ylim(min(ICL.list[["ICL.vec"]], na.rm = TRUE), max(ICL.list[["ICL.vec"]], na.rm = TRUE)) +
                                            ggplot2::xlab("K")
                                          
                                          ### plot logLike
                                          logLike.p   <- ggplot2::ggplot(data = logLike.list[["logLike.tab"]]) +
                                            ggplot2::geom_boxplot(ggplot2::aes(x = as.factor(K), y = logLike, group = K)) +
                                            ggplot2::xlab("K") +
                                            ggplot2::geom_text(data = logLike.list[["logLike.n"]], ggplot2::aes(x = 1:length(K), y = max(logLike.list[["logLike.vec"]], na.rm = TRUE),
                                                                                                                label = paste0("n=", n)), col = 'red', size = 4)
                                          
                                          ### coseq plots
                                          plot.coseq.res <- coseq::plot(coseq.res, conds = conds, collapse_reps = "average",
                                                                        graphs = c("profiles", "boxplots", "probapost_boxplots",
                                                                                   "probapost_barplots", "probapost_histogram"))
                                          
                                          # CoExpAnal[["plots"]] <- plot.coseq.res
                                          # CoExpAnal[["plots"]][["ICL"]]     <- ICL.p
                                          # CoExpAnal[["plots"]][["logLike"]] <- logLike.p
                                          
                                          return(c(plot.coseq.res, list("ICL" = ICL.p, "logLike" = logLike.p)))
                                        })


#' @title coseq.profile.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param SE.name the name of the data to fetch in the object if the object is a MultiAssayExperiment
#' @param numCluster cluster number
#' @param condition 
#' @param observation 
#' @export
#' @exportMethod coseq.profile.plot
#' @importFrom dplyr filter mutate rename full_join arrange group_by summarise
#' @importFrom reshape2 melt 
#' @noRd

methods::setMethod(f="coseq.profile.plot",
                   signature="SummarizedExperiment",
                   definition <- function(object, numCluster = 1, condition="groups", observation=NULL){
                     
                     coseq.res  <- object@metadata$CoExpAnal[["coseqResults"]]
                     assays.data <- dplyr::filter(as.data.frame(coseq.res@assays@data[[1]]), get(paste0("Cluster_",numCluster)) > 0.8)
                     
                     y_profiles.gg <- coseq.res@y_profiles[rownames(assays.data),] %>% 
                       data.frame() %>% 
                       dplyr::mutate(observations=rownames(.)) %>% 
                       reshape2::melt(id="observations", value.name = "y_profiles") %>%  
                       dplyr::rename(samples = variable) %>%
                       dplyr::full_join(object@metadata$Groups , by = "samples")
                     
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
                         ggplot2::geom_point(data = df, ggplot2::aes(x = groups, y = mean.y_profiles), color = "red", size = 2) +
                         ggplot2::geom_line( data = df, ggplot2::aes(x = groups, y = mean.y_profiles), color = "red", group = 1) +
                         ggplot2::ggtitle(paste0("Cluster: ",numCluster, "; nb_observations : ", dim(assays.data)[1], "; red : ", observation))
                     }
                     -                     
                       return(p)
                   })

#' @rdname coseq.profile.plot
#' @title coseq.profile.plot
#' @exportMethod coseq.profile.plot
methods::setMethod(f          = "coseq.profile.plot",
                   signature  = "MultiAssayExperiment",
                   definition = function(object, SE.name, numCluster = 1, condition="groups", observation=NULL){
                     
                     coseq.profile.plot(object = object[[SE.name]],
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
#' @noRd
#'
methods::setMethod(f="resetFlomicsMultiAssay", signature="MultiAssayExperiment",
                   
                   definition <- function(object, results, datasets = NULL){
                     
                     # if dataset is null we take all datasets present in MultiAssayExperiment object
                     if(is.null(datasets)){
                       datasets <- unlist(object@metadata$omicList)
                     }
                     else{
                       # check if given dataset name include in datasets presente in MultiAssayExperiment object
                       if(!datasets %in% unlist(object@metadata$omicList)){
                         warning("The given dataset name is not present in MultiAssayExperiment object")
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


