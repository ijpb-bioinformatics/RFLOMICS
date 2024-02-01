


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
                       dplyr::full_join(data.frame(MultiAssayExperiment::sampleMap(object)), by="assay") %>%
                       dplyr::mutate(y.axis = paste0(assay, "\n", "n=", nb_entities)) %>% dplyr::arrange(primary)
                     
                     data$primary <- factor(data$primary, levels = levels(Groups$samples)) 
                     
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


################################# EXPLORATION OF BIOLOGICAL AND TECHNICAL VARIABILITY ##################################


##### Statistical METHODS for exploring biological and technical variability


#' @title RunPCA
#' @description This function performs a principal component analysis on omic data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a "Normalization" slot is present in the metadata slot, then data are normalized before running the PCA according to the indicated transform method.
#' @param object An object of class \link{RflomicsSE-class}.
#' @param nbcp Number of components to compute. Default is 5.
#' @param raw boolean. Does the pca have to be ran on raw data or transformed and normalized data? Default is FALSE, pca is ran on transformed and normalized data.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod RunPCA
#' @importFrom FactoMineR PCA
#' @rdname RunPCA
#' 
methods::setMethod(f          = "RunPCA",
                   signature  = "RflomicsSE",
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
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, nbcp = 5, raw = FALSE){
                     
                     object[[SE.name]] <-  RunPCA(object[[SE.name]], 
                                                  nbcp = nbcp, 
                                                  raw  = raw)
                     
                     return(object)
                     
                   })

##### Graphical METHODS for exploring biological and technical variability


#' Library_size_barplot.plot
#'
#' @param object An object of class \link{RflomicsSE}
#' @return plot
#' @export
#' @importFrom ggplot2 ggplot geom_bar xlab ylab element_text ggtitle
#' @rdname Library_size_barplot.plot
#' @noRd
methods::setMethod(f          = "Library_size_barplot.plot",
                   signature  = "RflomicsSE",
                   definition <- function(object, raw = FALSE){
                     
                     if (getOmicsTypes(object) != "RNAseq") stop("WARNING: data are not RNAseq!")
                     
                     abundances <- SummarizedExperiment::assay(object)
                     Groups     <- getDesignMat(object)
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
                     
                     libSizeNorm <-  dplyr::full_join(Groups, data.frame("value" = pseudo , "samples" = names(pseudo)), by = "samples") %>%
                       dplyr::arrange(groups)
                     
                     libSizeNorm$samples <- factor(libSizeNorm$samples, levels = levels(Groups$samples))
                     
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
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, raw = FALSE){
                     
                     if (getOmicsTypes(object[[SE.name]]) == "RNAseq") {
                       return( Library_size_barplot.plot(object[[SE.name]], raw = raw))
                     }else{
                       stop("This function only applies to RNAseq data")
                     }
                     
                   })


#' @title Data_Distribution_plot
#'
#' @param object An object of class \link{RflomicsSE}
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
  signature = "RflomicsSE",
  definition = function(object, plot = "boxplot", raw = FALSE) {
    
    object2 <- checkTransNorm(object, raw = raw)
    pseudo <- SummarizedExperiment::assay(object2)
    Groups <- getDesignMat(object)
    
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
      dplyr::full_join(Groups, by = "samples") %>%
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
               ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none",
                              plot.margin=ggplot2::margin(0.5,0.5,0.5,1,"cm")) +
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
  signature = "RflomicsMAE",
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
#' in a \link{RflomicsSE-class} object. By default, samples are
#' colored by groups (all combinations of level's factor)
#' @param object An object of class \link{RflomicsSE-class}
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
                   signature = "RflomicsSE",
                   definition <- function(object, PCA, PCs=c(1,2), condition="groups"){
                     
                     ExpDesign <- getDesignMat(object)
                     
                     #
                     PC1 <- paste("Dim.",PCs[1], sep="")
                     PC2 <- paste("Dim.",PCs[2], sep="")
                     
                     if(PC1 == PC2) PC2 <- PC1+1
                     
                     score     <- object@metadata$PCAlist[[PCA]]$ind$coord[, PCs] %>% as.data.frame() %>%
                       dplyr::mutate(samples=row.names(.)) %>% dplyr::right_join(., ExpDesign, by="samples")
                     
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
                       ggplot2::geom_text(ggplot2::aes(label=samples), size=2, vjust="inward",hjust="inward") +
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
                     p <- p + ggplot2::geom_polygon(data = bb$res, ggplot2::aes_string(x=PC1, y=PC2, fill = condition),
                                                    show.legend = FALSE,
                                                    alpha = 0.1)
                     
                     print(p)
                     
                   })

#' @rdname plotPCA
#' @title plotPCA
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod plotPCA
methods::setMethod(f          = "plotPCA",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, PCA, PCs=c(1,2), condition="groups"){
                     
                     plotPCA(object[[SE.name]], PCA, PCs, condition)
                     
                   })

########################################## TRANSFORM DATA #################

#### METHOD to transform data

#' @title TransformData
#'
#' @param object An object of class \link{RflomicsSE}
#' @param transformMethod The transformation to store in the metadata or to store and apply if modify_assay is TRUE.
#' @param modify_assay Boolean. Do the transformation need to be applied on the data? The raw data will be replaced by the transformed ones.
#'
#' @return An object of class \link{RflomicsSE}
#'
#' @exportMethod TransformData
#' @rdname TransformData
#' 
methods::setMethod(f          = "TransformData",
                   signature  = "RflomicsSE",
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
                   signature  = "RflomicsMAE",
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
#' @param object An object of class \link{RflomicsSE}
#' @param filterMethod The filtering model ("CPM")
#' @param filterStrategy The filtering strategy ("NbConditions" or "NbReplicates")
#' @param cpmCutoff The CPM cutoff.
#' @return An object of class \link{RflomicsSE}
#' @details
#' Filtered dataset is stored in the ExperimentList slot of the \link{RflomicsSE} object
#' as a List named (DataName.filtred).
#' List of filtered features are stored as a named list ("FilteredFeatures") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{RflomicsSE} object.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod FilterLowAbundance
#' @importFrom edgeR DGEList filterByExpr cpm
#' @seealso edgeR::filterByExpr#' 
#' @rdname FilterLowAbundance


methods::setMethod(f         = "FilterLowAbundance",
                   signature = "RflomicsSE",
                   definition <- function(object, filterMethod= "CPM", filterStrategy = "NbConditions", cpmCutoff = 5){
                     
                     if(isFALSE(filterStrategy %in% c("NbReplicates","NbConditions"))) 
                       stop("filterStrategy argument must be one of this tow options : NbReplicates or NbConditions")
                     
                     
                     objectFilt <- object
                     assayFilt  <- SummarizedExperiment::assay(objectFilt)
                     Groups     <- getDesignMat(object)
                     
                     ## nbr of genes with 0 count
                     genes_flt0  <- objectFilt[rowSums(assayFilt) <= 0, ]@NAMES
                     
                     ## remove 0 count
                     objectFilt  <- objectFilt[rowSums(assayFilt)  > 0, ]
                     assayFilt   <- assay(objectFilt)
                     
                     ## filter cpm
                     NbReplicate  <- table(Groups$groups)
                     NbConditions <- length(unique(Groups$groups))
                     
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
                   signature  = "RflomicsMAE",
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
#' class \link{RflomicsSE}.
#' \itemize{
#' \item{For RNAseq data:}{the TMM function of edgeR is proposed by default, see the ref}
#' \item{For Proteomic data:}{}
#' \item{For Metabolomic data:}{}
#' }
#' @param object An object of class \link{RflomicsSE}
#' @param NormMethod Normalization method
#' @param modify_assay Does the normalization have to be applied or just stored for later? Recommended it stays FALSE.
#' @return An object of class \link{RflomicsSE}
#' The applied normalization method and computed scaling factors (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "methode" and "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{RflomicsSE} object.
#' @exportMethod RunNormalization
#' @seealso TMM.Normalization
#' @rdname RunNormalization
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).

methods::setMethod(f          = "RunNormalization",
                   signature  = "RflomicsSE",
                   definition = function(object, NormMethod = NULL, modify_assay = FALSE){
                     
                     Groups     <- getDesignMat(object)
                     
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
                            "TMM"        = {coefNorm  <- TMM.Normalization(assay(object2), Groups$groups) },
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
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, NormMethod, modify_assay = FALSE){
                     
                     object[[SE.name]] <-  RunNormalization(object       = object[[SE.name]],
                                                            NormMethod   = NormMethod,
                                                            modify_assay = modify_assay)
                     
                     return(object)
                     
                   })



# ------ process data -----

# Function non generique pour les autres data

#' @title runDataProcessing
#' @description This function applied a processing (filtering, normalization and/or transformation, PCA) on an omic data sets stored in an object of
#' class \link{RflomicsSE}.
#' \itemize{
#' \item{For RNAseq data:}{}
#' \item{For Proteomic data:}{}
#' \item{For Metabolomic data:}{}
#' }
#' @param object An object of class \link{RflomicsSE}
#' @param samples samples to keep.
#' @param lowCountFiltering_strategy strategy of RNAseq low count filtering. Mandatory for RNAseq data. Default value : "NbReplicates".
#' @param lowCountFiltering_CPM_Cutoff CPM cutoff for RNAseq low count filtering. Mandatory for RNAseq data. Default value : 1.
#' @param normalisation_method method of normalisation. Mandatory for RNAseq data. Default value : RNAseq = TMM.
#' @param transformation_method method of transformation.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runDataProcessing
#' @seealso runSampleFiltering
#' @seealso FilterLowAbundance
#' @seealso RunNormalization
#' @seealso TransformData
#' @rdname runDataProcessing

methods::setMethod(f          = "runDataProcessing",
                   signature  = "RflomicsSE",
                   definition = function(object, samples=NULL, lowCountFiltering_strategy = "NbReplicates", lowCountFiltering_CPM_Cutoff = 1, 
                                         normalisation_method = "none", transformation_method = "none")
                   {
                     
                     # keep selected samples
                     print("#    => select samples...")
                     object <- runSampleFiltering(object, samples)
                     
                     if(nrow(getDesignMat(object)) == 0) stop("no samples in object!")
                     
                     # spported values:
                     lowCountFiltering_strategy.sup     <- c("NbReplicates","NbConditions")
                     transformation_method.sup          <- c("log1p", "squareroot", "log2", "log10", "none")
                     normalisation_method.abundance.sup <- c("median", "totalSum", "none")
                     normalisation_method.count.sup     <- c("TMM")
                     
                     # apply data processing
                     switch(object@metadata$omicType,
                            "RNAseq" = {
                              
                              # Filter low abundance
                              print("#    => Low counts Filtering...")
                              if(is.null(lowCountFiltering_strategy)   || !lowCountFiltering_strategy %in% lowCountFiltering_strategy.sup) 
                                stop("the low count filtering strategy : ", lowCountFiltering_strategy, " isn't supported by rflomics package. Supported values : ",  paste(lowCountFiltering_strategy.sup, collapse = ", "))
                              
                              if(is.null(lowCountFiltering_CPM_Cutoff) || !is.numeric(lowCountFiltering_CPM_Cutoff)) 
                                stop(lowCountFiltering_CPM_Cutoff, " must be a integer value > 1")
                              
                              SE.processed <- FilterLowAbundance(object = object, filterStrategy = lowCountFiltering_strategy, cpmCutoff = lowCountFiltering_CPM_Cutoff)
                              
                              # Run Normalisation 
                              print("#    => Counts normalization...")
                              if(is.null(normalisation_method) || normalisation_method != "TMM"){
                                normalisation_method <- "TMM"
                                warning("only ", normalisation_method.count.sup, " method is supported for ", object@metadata$omicType, " normalisation.")
                              }

                              SE.processed <- RunNormalization(SE.processed, NormMethod = normalisation_method)
                            },
                            {
                              print("#    => transformation data...")
                              if(is.null(transformation_method)) transformation_method <- "none"
                              if(!transformation_method %in% transformation_method.sup) 
                                stop("the transformation method ", transformation_method," is not support in rflomics package.Supported values : ", paste(transformation_method.sup, collapse = ", "))
                              SE.processed <- TransformData(object, transformMethod = transformation_method)
                              
                              print("#    => Run normalization...")
                              if(is.null(normalisation_method)) normalisation_method <- "none"
                              if(!normalisation_method %in% normalisation_method.abundance.sup) 
                                stop("the normalisation method ", normalisation_method," is not support in rflomics package. Supported values : ", paste(normalisation_method.abundance.sup, collapse = ", "))
                              SE.processed <- RunNormalization(SE.processed, NormMethod = normalisation_method)
                            }
                     )
                     
                     #### Run PCA for filtred & normalized data ####
                     print("#    => Compute PCA ")
                     SE.processed <- RunPCA(SE.processed)  
                     
                     SE.processed@metadata$DataProcessing[["done"]] <- TRUE
                     
                     return(SE.processed)
                   })


#' @rdname runDataProcessing
#' @title runDataProcessing
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod runDataProcessing
methods::setMethod(f          = "runDataProcessing",
                   signature  = "RflomicsMAE",
                   definition = function(object, samples=NULL, lowCountFiltering_strategy = "NbReplicates", lowCountFiltering_CPM_Cutoff = 1, 
                                         normalisation_method = "none", transformation_method = "none", SE.name){
                     
                     # if paste0(SE.name, ".raw") exist
                     if (!paste0(SE.name, ".raw") %in% names(object)){
                       
                       if (!SE.name %in% names(object)) stop("no ", SE.name, " SE object.")
                       stop("raw SE must be tagged by .raw")
                     }
                       
                     # if paste0(SE.name, ".raw") exist
                     if(is.null(object[[paste0(SE.name, ".raw")]])) stop("raw SE must be tagged by .raw")
                     
                     SE.raw       <- object[[paste0(SE.name, ".raw")]]
                     SE.processed <- runDataProcessing(object = SE.raw,
                                                       samples = samples, 
                                                       lowCountFiltering_strategy = lowCountFiltering_strategy, 
                                                       lowCountFiltering_CPM_Cutoff = lowCountFiltering_CPM_Cutoff,
                                                       normalisation_method = normalisation_method, 
                                                       transformation_method = transformation_method)
                     
                     
                     # remove SE processed if exist
                     if (SE.name %in% names(object)) object <- object[,, -which(names(object) == SE.name)]

                     # add new SE with processed data
                     object <- eval(parse(text = paste0('c( object ,', SE.name, ' = SE.processed )')))

                     return(object)
                   })

# ------ filtering per sample -----

#' @title runSampleFiltering
#' @description This function applied sample filtering on an omic data sets stored in an object of
#' class \link{RflomicsSE}.
#' \itemize{
#' \item{For RNAseq data:}{}
#' \item{For Proteomic data:}{}
#' \item{For Metabolomic data:}{}
#' }
#' @param object An object of class \link{RflomicsSE}
#' @param samples samples to keep.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runSampleFiltering
#' @rdname runSampleFiltering

methods::setMethod(f          = "runSampleFiltering",
                   signature  = "RflomicsSE",
                   definition = function(object, samples=NULL) {
                     
                     # if no samples to filter
                     if(is.null(samples)) return(object) 
                     
                     # check if samples overlap with
                     samples <- intersect(colnames(object), samples)
                     
                     # if no iverlap with data colnames
                     if(length(intersect(samples, colnames(object))) == 0) {
                       warning("No overlap between samples and colnames of data")
                       return(object)
                       }
                     
                     # if 100% ovelap!
                     if(length(samples) == length(colnames(object))) return(object)
                     
                     # keep selected samples
                     # keep only samples in data matrix, and colData
                     SE.new <- object[, object$samples %in% samples]
                     
                     # if we remove all samples :
                     if(nrow(SE.new@colData) == 0) SE.new@colData <- SE.new@colData[c("samples", "groups")]
                     
                     for (factor in c(bioFactors(object), batchFactors(object))){
                       
                       # if only one category remains after the filter, it's will be removed
                       if(length(unique(SE.new@colData[[factor]])) <= 1 ) {
                         SE.new@colData[[factor]] <- NULL
                         factor.types <- getFactorTypes(SE.new)
                         SE.new@metadata$design$factorType <- factor.types[which(names(factor.types) != factor)]
                         # replace with setFactorTypes
                       }
                       else{
                         F.levels <- levels(SE.new@colData[[factor]])
                         SE.new@colData[[factor]] <- factor(SE.new@colData[[factor]], levels = intersect(F.levels, unique(SE.new@colData[[factor]])))
                       }
                     }
                     colData <- as.data.frame(SE.new@colData)
                     order_levels <- with(colData, do.call(order, colData[c(bioFactors(SE.new), batchFactors(SE.new))]))
                     SE.new$samples <- factor(SE.new$samples, levels = unique(SE.new$samples[order_levels]))
                     SE.new$groups  <- factor(SE.new$groups,  levels = unique(SE.new$groups[order_levels]))
                     
                     # à retirer dès que je remplace Groups par colData
                     SE.new@metadata$Groups <- dplyr::filter(SE.new@metadata$Groups, samples %in% SE.new$samples)
                     
                     SE.new@metadata$DataProcessing$filteredSamples <- setdiff(colnames(object), samples)
                     
                     return(SE.new)
                   })




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


