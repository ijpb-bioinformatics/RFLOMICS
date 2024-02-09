


###### METHOD to check the completness of the ExpDesign



#' @title CheckExpDesign
#' @description This method checks some experimental design characteristics.
#'  A complete design and at least one biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{RflomicsMAE-class}
#' @return a gg plot object
#' \itemize{
#'  \item{"plot:"}{ plot of count data.frame.}
#'  }
#'  
#' @exportMethod CheckExpDesign
#' @noRd

methods::setMethod(f         = "CheckExpDesign",
                   signature = "RflomicsMAE",
                   definition <- function(object){
                     
                     # check presence of bio factors
                     if (!length(bioFactors(object)) %in% 1:3){ stop("No bio factor! or nbr of bio factors exceed 3!") }
                     if (!length(batchFactors(object)) %in% 1:2){ stop("No replicates found!") }
                     
                     ####################
                     
                     BioFact <- bioFactors(object)
                     coldata <- getDesignMat(object) %>%
                       dplyr::mutate(samples=rownames(.))
                     #coldata <- tibble::as_tibble(coldata)
                     coldata <- MultiAssayExperiment::sampleMap(object) %>% as.data.frame() %>% 
                       dplyr::left_join(coldata, by = c("primary" = "samples"))
                     
                     all_combin_cond <- lapply(BioFact, function(x){ 
                       df <- unique(coldata[x])
                       rownames(df) <- 1:nrow(df)
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
#' @param An object of class \link{RflomicsSE-class}
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
                   signature = "RflomicsSE",
                   definition <- function(object, sampleList=NULL){
                     
                     object <- runSampleFiltering(object, samples = sampleList)
                     output <- list()
                     
                     # check presence of bio factors
                     # check presence of bio factors
                     if (!length(bioFactors(object)) %in% 1:3){ 
                       output[["messages"]] <-  "Error : You need at least 1 biological factor with at least 2 modalities."
                       output[["error"]]    <- TRUE
                       return(output)
                     }
                     if (!length(batchFactors(object)) %in% 1:2){ 
                       output[["messages"]] <-  "Error : You need at least 1 batch factor with at least 2 replicats."
                       output[["error"]]    <- TRUE 
                       return(output)
                     }
                     
                     # Only works with bio and batch factors for the rest of the function
                     ExpDesign <- getDesignMat(object)
                     
                     bio.fact <- bioFactors(object)
                     
                     # group_count <- as.data.frame(dF.List) %>%
                     #   table() %>% as.data.frame() %>%
                     #   dplyr::full_join(ExpDesign, by=names(dF.List)) %>%
                     #   dplyr::mutate_at(.vars = "samples", .funs = function(x) dplyr::if_else(is.na(x), 0, 1)) %>%
                     #   dplyr::group_by_at((bio.fact)) %>%
                     #   dplyr::summarise(Count=sum(samples), .groups = "keep")
                     
                     #remplacer le code ci-dessus par celui en bas
                     group_count <- ExpDesign %>% dplyr::group_by_at((bio.fact)) %>% dplyr::count(name = "Count")
                     
                     
                     output[["counts"]] <- group_count
                     output[["plot"]]   <- plotExperimentalDesign(counts = output[["counts"]], message= output[["messages"]])
                     
                     # check presence of relicat / batch
                     # check if design is complete
                     # check if design is balanced
                     # check nbr of replicats
                     if(min(group_count$Count) == 0){
                       
                       output[["messages"]] <- "Error : The experimental design is not complete."
                       output[["error"]]   <- TRUE
                     }
                     else if(min(group_count$Count) == 1){
                       
                       output[["messages"]] <-  "Error : You need at least 2 biological replicates."
                       output[["error"]]   <- TRUE
                     }
                     else if(length(unique(group_count$Count)) != 1){
                       
                       output[["messages"]] <- "The experimental design is complete but not balanced."
                       output[["error"]]   <- FALSE
                     }
                     else{
                       output[["messages"]] <- "The experimental design is complete and balanced."
                       output[["error"]]   <- FALSE
                     }
                     
                     ### plot
                     
                     
                     return(output)
                   })


#' @title CheckExpDesignCompleteness
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod CheckExpDesignCompleteness
#' @noRd
methods::setMethod(f         = "CheckExpDesignCompleteness",
                   signature = "RflomicsMAE",
                   definition <- function(object, SE.name, sampleList=NULL){
                     
                     SEObject <- object[[SE.name]]
                     
                     CheckExpDesignCompleteness(SEObject, 
                                                sampleList = sampleList)
                     
                   })

#' @title Datasets overview plot
#' @description This function plot overview of loaded datasets aligned per sample 
#' (n=number of entities (genes/metabolites/proteins); k=number of samples)
#' @param object An object of class \link{RflomicsMAE-class}
#' @importFrom dplyr full_join mutate arrange select
#' @importFrom ggplot2 ggplot aes element_blank element_text geom_col theme 
#' labs scale_y_continuous geom_tile
#' 
#' @exportMethod Datasets_overview_plot
#' @return plot

methods::setMethod(f         = "Datasets_overview_plot",
                   signature = "RflomicsMAE",
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
                     
                     Groups <- getDesignMat(object)
                     
                     nb_entities <- lapply(object@ExperimentList, function(SE){ dim(SE)[1] }) %>% unlist()
                     
                     data <- data.frame(nb_entities = nb_entities, assay = names(nb_entities)) %>%
                       full_join(data.frame(sampleMap(object)), by="assay") %>%
                       mutate(y.axis = paste0(assay, "\n", "n=", nb_entities)) %>% arrange(primary)
                     
                     data$primary <- factor(data$primary, levels = levels(Groups$samples)) 
                     
                     nb_entities_ord <- select(data, y.axis, nb_entities) %>% unique() %>% arrange(desc(nb_entities))
                     nb_entities_ord$nb_entities <- log(nb_entities_ord$nb_entities)
                     tmp.vec <- c(0)
                     breaks  <- vector()
                     for(i in 1:length(nb_entities_ord$nb_entities)){ 
                       tmp.vec[i+1] <- tmp.vec[i] + nb_entities_ord$nb_entities[i]
                       breaks[i] <- tmp.vec[i] + nb_entities_ord$nb_entities[i]/2 
                     } 
                     
                     switch (as.character(real.size),
                             "TRUE"  = {
                               p <- ggplot(data, aes(x=primary, y=log(nb_entities))) +
                                 geom_col(aes(fill = y.axis)) + 
                                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.ticks = element_blank(), 
                                                axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none",
                                                axis.text.y =  element_text(hjust = 0)) +  
                                  labs(x=paste0("Samples (k=", length(unique(sampleMap(object)$primary)), ")"), y="") +
                                 scale_y_continuous(breaks = (breaks), labels = nb_entities_ord$y.axis)
                               
                             },
                             "FALSE" = {
                               p <-  ggplot(data,  aes(x=primary, y=y.axis)) +
                                  geom_tile( aes(fill = y.axis), colour = "grey50") +
                                  theme(panel.grid.major =  element_blank(), panel.grid.minor = element_blank(),
                                                panel.background =  element_blank(), axis.ticks = element_blank(), legend.position="none",
                                                axis.text.x =  element_text(angle = 90, hjust = 1)) +
                                 labs(x=paste0("Samples (k=", length(unique(sampleMap(object)$primary)), ")"), y="")
                             }
                     )
                     return(p)
                   })

###### METHOD to obtain the Matrix of contrast with their names and coefficients

#
#

#' @title getContrastMatrix
#' @description Defines contrast matrix or contrast list with contrast name and contrast coefficients
#' @param object An object of class \link{RflomicsMAE-class}
#' @param contrastList A data.frame of contrast
#' @return An object of class \link{RflomicsMAE-class}
#' @seealso getExpressionContrast
#' @exportMethod getContrastMatrix
#' @importFrom stats formula terms.formula
#' @noRd
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
methods::setMethod(f          = "getContrastMatrix",
                   signature  = "RflomicsMAE",
                   definition <- function(object, SE.name, modelFormula = NULL, contrastList=NULL){
                     
                     if (is.null(object[[SE.name]])) stop("no Experiment named ", SE.name, " in MAE object")
                     if (is.null(modelFormula)) stop("Model.formula arg is mandatory.")
                     if (is.null(contrastList)) stop("contrastList is mandatory.")
                     if (any(!c("contrast", "contrastName", "groupComparison", "type") %in% names(contrastList))) 
                       stop("contrastList data.frame must contain at least these colomn : contrast, contrastName, groupComparison, type")
                     
                     object <- setModelFormula(object, modelFormula)
                    
                     object[[SE.name]] <- getContrastMatrix(object = object[[SE.name]], contrastList = contrastList, modelFormula = modelFormula)
                     
                     return(object)
                   })

#' @title getContrastMatrix
#' @description Defines contrast matrix or contrast list with contrast name and contrast coefficients
#' @param object An object of class \link{RflomicsMAE-class}
#' @param contrastList a data.frame of contrast
#' @param modelFormula a model formula
#' @return An object of class \link{RflomicsSE-class}
#' @seealso getExpressionContrast
#' @exportMethod getContrastMatrix
#' @importFrom stats formula terms.formula
#' @noRd
#' @author Christine Paysant-Le Roux, adapted by Nadia Bessoltane
methods::setMethod(f          = "getContrastMatrix",
                   signature  = "RflomicsSE",
                   definition <- function(object, modelFormula = NULL, contrastList=NULL){

                     
                     if(is.null(modelFormula)) stop("Model.formula arg is mandatory.")
                     if(is.null(contrastList)) stop("contrastList arg is mandatory.")
                     
                     ExpDesign <- getDesignMat(object)
                     
                     factorBio <- bioFactors(object)
                      
                     object <- setModelFormula(object, modelFormula)
                     object@metadata$design$Contrasts.Coeff <- getContrastMatrixF(ExpDesign = ExpDesign, factorBio = factorBio, contrastList = contrastList$contrast, modelFormula)
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


##### RESETTING ####

#' resetFlomicsMultiAssay
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param results vector of results names
#' @param dataset dataset name. If dataset == NULL, all datasets will be reset
#' @return An object of class \link{RflomicsMAE}
#' @export
#' @exportMethod resetFlomicsMultiAssay
#' @noRd
#'
methods::setMethod(f="resetFlomicsMultiAssay", signature="RflomicsMAE",
                   
                   definition <- function(object, results, datasets = NULL){
                     
                     # if dataset is null we take all datasets present in RflomicsMAE object
                     if(is.null(datasets)){
                       datasets <- unlist(object@metadata$omicList)
                     }
                     else{
                       # check if given dataset name include in datasets presente in RflomicsMAE object
                       if(!datasets %in% unlist(object@metadata$omicList)){
                         warning("The given dataset name is not present in RflomicsMAE object")
                         return(object)
                       }
                     }
                     
                     for(res in results){
                       
                       for(data in datasets){
                         if(!is.null(object[[data]])){
                           
                           if(!is.null(object[[data]]@metadata[[res]])){ object[[data]]@metadata[[res]] <- NULL }
                         }
                       }
                       
                       object@metadata[[res]] <- NULL
                     }
                     
                     
                     # for(data in datasets){
                     #   
                     #   if(!is.null(object[[data]])){
                     #     
                     #     dataset <- object[[data]]
                     #     
                     #     for(res in results){
                     #       if(!is.null(dataset@metadata[[res]])){ dataset@metadata[[res]] <- NULL }
                     #     }
                     #     
                     #     object[[data]] <- dataset
                     #   }
                     #   
                     # }
                     return(object)
                   })


