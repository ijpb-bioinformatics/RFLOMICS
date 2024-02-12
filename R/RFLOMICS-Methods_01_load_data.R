### ============================================================================
### [01_Load_Data] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------

#' @importFrom dplyr full_join mutate arrange select group_by_at count left_join
#' right_join mutate_at if_else 
#' @importFrom ggplot2 ggplot aes element_blank element_text geom_col theme 
#' labs scale_y_continuous geom_tile scale_fill_manual ylab xlab labs
#' facet_grid 
#' @importFrom purrr reduce

# ---- ACCESSORS ----
# ---- getDesignMat:    get colData ----
#' @title Get design matrix / colData
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a data frame with design matrix
#' @exportMethod getDesignMat
methods::setMethod(f          = "getDesignMat",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     return(as.data.frame(object@colData))
                   })


#' @title Get design matrix / colData
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a data frame with design matrix
#' @exportMethod getDesignMat
methods::setMethod(f          = "getDesignMat",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     return(as.data.frame(object@colData))
                   })

# ---- getDatasetNames: get experiment names from RflomicsMAE object ----

#' @title Get Experiment Names
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with omic dataset names
#' @exportMethod getDatasetNames
methods::setMethod(f          = "getDatasetNames",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     ExperimentNames <- unlist(object@metadata$omicList)
                     names(ExperimentNames) <- NULL
                     
                     return(ExperimentNames)
                   })

#' @title Get omic dataset Names
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with omic dataset names
#' @exportMethod getDatasetNames
methods::setMethod(f          = "getDatasetNames",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     return(names(object@metadata$omicType))
                   })

# ---- getOmicsTypes:   get omics types of experiments from RflomicsMAE object ----

#' @title Get omics types of experiments
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @importFrom purrr reduce
#' @return names vector
#' @exportMethod getOmicsTypes
methods::setMethod(f          = "getOmicsTypes",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     OmicsTypes <- lapply(names(object@metadata$omicList), function(x){ 
                       
                       datasetNames <- object@metadata$omicList[[x]]
                       datasetType  <- rep(x, length(datasetNames))
                       names(datasetType) <- datasetNames
                       
                       return(datasetType)
                       
                     }) |> reduce(c)
                     
                     return(OmicsTypes)
                   })

#' @title Get omics types of experiments
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return names vector
#' @exportMethod getOmicsTypes
methods::setMethod(f          = "getOmicsTypes",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     return(object@metadata$omicType)
                   })


# ---- getFactorNames:  get Factor names ----
#' @title Get factor names
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with names of design factors
#' @exportMethod getFactorNames
methods::setMethod(f          = "getFactorNames",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     return(names(object@metadata$design$Factors.Type))
                   })


#' @title Get factor names
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with names of design factors
#' @exportMethod getFactorNames
methods::setMethod(f          = "getFactorNames",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     return(names(object@metadata$design$factorType))
                   })


# ---- getFactorTypes:  get Factor types ----
#' @title Get factor types
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a named vector
#' @exportMethod getFactorTypes
methods::setMethod(f          = "getFactorTypes",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     return(object@metadata$design$Factors.Type)
                   })


#' @title Get factor types
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a named vector
#' @exportMethod getFactorTypes
methods::setMethod(f          = "getFactorTypes",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     return(object@metadata$design$factorType)
                   })


# ---- getBioFactors:   get bio factor ----
#' @title Get bio factor.
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with biological factors
#' @exportMethod getBioFactors
methods::setMethod(f          = "getBioFactors",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BIO"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })


#' @title Get bio factor.
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with biological factors
#' @exportMethod getBioFactors
methods::setMethod(f          = "getBioFactors",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BIO"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

#' @title Get bio factor.
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with biological factors
#' @exportMethod bioFactors
methods::setMethod(f          = "bioFactors",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BIO"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })


#' @title Get bio factor.
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with biological factors
#' @exportMethod bioFactors
methods::setMethod(f          = "bioFactors",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BIO"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

# ---- getBatchFactors: get batch factor names ----
#' @title Get batch factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with batch factors
#' @exportMethod getBatchFactors
methods::setMethod(f          = "getBatchFactors",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BATCH"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

#' @title Get batch factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with batch factors
#' @exportMethod getBatchFactors
methods::setMethod(f          = "getBatchFactors",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BATCH"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })


#' @title Get batch factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with batch factors
#' @exportMethod batchFactors
methods::setMethod(f          = "batchFactors",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BATCH"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

#' @title Get batch factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with batch factors
#' @exportMethod batchFactors
methods::setMethod(f          = "batchFactors",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "BATCH"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

# ---- getMetaFactors:  get meta Factor names ----
#' @title Get meta factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with batch factors
#' @exportMethod getMetaFactors
methods::setMethod(f          = "getMetaFactors",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "META"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

#' @title Get meta factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with batch factors
#' @exportMethod getMetaFactors
methods::setMethod(f          = "getMetaFactors",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "META"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })


#' @title Get meta factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @return a vector with batch factors
#' @exportMethod metaFactors
methods::setMethod(f          = "metaFactors",
                   signature  = "RflomicsMAE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "META"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

#' @title Get meta factor
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @return a vector with batch factors
#' @exportMethod metaFactors
methods::setMethod(f          = "metaFactors",
                   signature  = "RflomicsSE",
                   definition <- function(object){
                     
                     factVect <- toupper(getFactorTypes(object))
                     res <- names(factVect)[factVect == "META"]
                     
                     if(length(res) == 0) return(NULL)
                     return(res)
                   })

# ---- getRflomicsSE:   get RflomicsSE object of one omic dataset ----
#' @title Get RflomicsSE object of one omic dataset
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @param datasetName name of omic dataset
#' @return An object of class \link{RflomicsSE-class}
#' @exportMethod getRflomicsSE
methods::setMethod(f          = "getRflomicsSE",
                   signature  = "RflomicsMAE",
                   definition <- function(object, datasetName = NULL){
                     
                     if(is.null(datasetName)) return(NULL)
                     
                     return(object[[datasetName]])
                   })

# ---- getFactorModalities: ----
#' @title getFactorModalities
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @param factorName a name of factor
#' @return a vector
#' @exportMethod getFactorModalities
methods::setMethod(f          = "getFactorModalities",
                   signature  = "RflomicsMAE",
                   definition <- function(object, factorName){
                     
                     if(is.null(factorName)) return(NULL)
                     if(!factorName %in% getFactorNames(object)) return(NULL) 
                     
                     return(levels(getDesignMat(object)[[factorName]]))
                   })


#' @title getFactorModalities
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsSE-class}
#' @param factorName a name of factor
#' @return a vector
#' @exportMethod getFactorModalities
methods::setMethod(f          = "getFactorModalities",
                   signature  = "RflomicsSE",
                   definition <- function(object, factorName){
                     
                     if(is.null(factorName)) return(NULL)
                     if(!factorName %in% getFactorNames(object)) return(NULL) 
                     
                     return(levels(getDesignMat(object)[[factorName]]))
                   })

# ---- get names of executed analysis ----
# ---- METHODS ----
# ---- subRflomicsMAE:  subset a RflomicsMAE from ----
#' @title extract sub RflomicsMAE
#' @description
#' A short description...
#' 
#' @param object An object of class \link{RflomicsMAE-class}
#' @param omicNames a vector with the dataset names
#' @return object An object of class \link{RflomicsMAE-class}
#' @exportMethod subRflomicsMAE
methods::setMethod(f          = "subRflomicsMAE",
                   signature  = "RflomicsMAE",
                   definition <- function(object, omicNames = NULL){
                     
                     if(is.null(omicNames)) return(object)
                     dataset.names <- names(object)
                     
                     if(!all(omicNames %in% dataset.names)) return(NULL)
                     
                     return(object[,, omicNames])
                   })

# ---- checkExpDesignCompleteness ----
#' @title checkExpDesignCompleteness
#' @description This method checks some experimental design characteristics.
#'  A complete design and at least one biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{RflomicsSE-class}
#' @param sampleList list of samples to check.
#' @return a string with message
#' @exportMethod checkExpDesignCompleteness
methods::setMethod(f         = "checkExpDesignCompleteness",
                   signature = "RflomicsSE",
                   definition <- function(object, sampleList=NULL){
                     
                     object <- RFLOMICS::runSampleFiltering(object, samples = sampleList)
                     
                     output <- list()
                     output[["error"]] <- FALSE
                     
                     # Only works with bio and batch factors for the rest of the function
                     ExpDesign <- getDesignMat(object)
                     bio.fact <- getBioFactors(object)
                     
                     # check presence of bio factors
                     if (!length(getBioFactors(object)) %in% 1:3){ 
                       output[["messages"]] <-  "Error : You need at least 1 biological factor with at least 2 modalities."
                       output[["error"]]    <- TRUE
                       return(output)
                     }
                     # check presence of bash factors
                     if (!length(getBatchFactors(object)) %in% 1:2){ 
                       output[["messages"]] <-  "Error : You need at least 1 batch factor with at least 2 replicats."
                       output[["error"]]    <- TRUE 
                       return(output)
                     }
                     
                     #remplacer le code ci-dessus par celui en bas
                     group_count <- .countSamplesPerCondition(ExpDesign, bio.fact)
                     
                     # check presence of relicat / batch
                     # check if design is complete
                     # check if design is balanced
                     # check nbr of replicats
                     if(min(group_count$Count) == 0){
                       
                       output[["messages"]] <- "Error : The experimental design is not complete."
                       output[["error"]]    <- TRUE
                     }
                     else if(min(group_count$Count) == 1){
                       
                       output[["messages"]] <-  "Error : You need at least 2 biological replicates."
                       output[["error"]]    <- TRUE
                     }
                     else if(length(unique(group_count$Count)) != 1){
                       
                       output[["messages"]] <- "The experimental design is complete but not balanced."
                     }
                     else{
                       output[["messages"]] <- "The experimental design is complete and balanced."
                     }
                     
                     return(output)
                   })

#' @title checkExpDesignCompleteness
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod checkExpDesignCompleteness
#' @noRd
methods::setMethod(f         = "checkExpDesignCompleteness",
                   signature = "RflomicsMAE",
                   definition <- function(object, omicName, sampleList=NULL){
                     
                     if(is.null(omicName)) stop("Arg omicName missed")
                     
                     SEObject <- getRflomicsSE(object, omicName)
                     
                     checkExpDesignCompleteness(SEObject, sampleList = sampleList)
                   })



# ---- PLOTS ----
# ---- plotDataOverview ----
#' @title Dataset overview plot
#' @description This function plot overview of loaded datasets aligned per sample 
#' (n=number of entities (genes/metabolites/proteins); k=number of samples)
#' @param object An object of class \link{RflomicsMAE-class}
#' @param omicNames a vector with dataset names
#' @importFrom MultiAssayExperiment sampleMap
#' @param realSize booleen value
#' @exportMethod plotDataOverview
#' @return plot
methods::setMethod(f         = "plotDataOverview",
                   signature = "RflomicsMAE",
                   definition <- function(object, omicNames=NULL, realSize=FALSE){
                     
                     if(length(object) == 0) stop("object is NULL")
                     
                     object <- subRflomicsMAE(object, omicNames)
                     
                     if(is.null(object)) return(NULL)
                     
                     Groups <- getDesignMat(object)
                     
                     nb_entities <- lapply(names(object), function(SE){ dim(object[[SE]])[1] }) %>% unlist()
                     names(nb_entities) <- names(object)
                     
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
                     
                     switch (as.character(realSize),
                             "TRUE"  = {
                               p <- ggplot(data, aes(x=primary, y=log(nb_entities))) +
                                 geom_col(aes(fill = y.axis)) + 
                                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.ticks = element_blank(), 
                                                axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none",
                                                axis.text.y = element_text(hjust = 0)) +  
                                 labs(x=paste0("Samples (k=", length(unique(sampleMap(object)$primary)), ")"), y="") +
                                 scale_y_continuous(breaks = (breaks), labels = nb_entities_ord$y.axis)
                               
                             },
                             "FALSE" = {
                               p <- ggplot(data, aes(x=primary, y=y.axis)) +
                                 geom_tile(aes(fill = y.axis), colour = "grey50") +
                                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                panel.background = element_blank(), axis.ticks = element_blank(), legend.position="none",
                                                axis.text.x = element_text(angle = 90, hjust = 1)) +
                                 labs(x=paste0("Samples (k=", length(unique(sampleMap(object)$primary)), ")"), y="")
                             }
                     )
                     return(p)
                   })



# ---- plotConditionsOverview ----
#' @title plotConditionsOverview
#' @description this method 
#'  A complete design and at least one biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{RflomicsMAE-class}
#' @param omicNames a vector with list of dataset names
#' @importFrom MultiAssayExperiment sampleMap
#' @return a gg plot object
#' \itemize{
#'  \item{"plot:"}{ plot of count data.frame.}
#'  }
#'  
#' @exportMethod plotConditionsOverview
#' @noRd

methods::setMethod(f         = "plotConditionsOverview",
                   signature = "RflomicsMAE",
                   definition <- function(object, omicNames = NULL){
                     
                     # check presence of bio factors
                     if (!length(getBioFactors(object)) %in% 1:3){ stop("No bio factor! or nbr of bio factors exceed 3!") }
                     if (!length(getBatchFactors(object)) %in% 1:2){ stop("No replicates found!") }
                     
                     ####################
                     
                     BioFact <- getBioFactors(object)
                     coldata <- getDesignMat(object) %>%
                       mutate(samples=rownames(.))
                     #coldata <- tibble::as_tibble(coldata)
                     coldata <- sampleMap(object) %>% as.data.frame() %>% 
                       left_join(coldata, by = c("primary" = "samples"))
                     
                     all_combin_cond <- lapply(BioFact, function(x){ 
                       df <- unique(coldata[x])
                       rownames(df) <- 1:nrow(df)
                       return(df) 
                     }) %>% reduce(merge)
                     
                     counts <- coldata %>% select(assay, all_of(BioFact)) %>% unique() %>% 
                       group_by_at(BioFact) %>% count(name = "Count") %>% 
                       right_join(all_combin_cond, by = BioFact) %>% 
                       mutate_at(.vars = "Count", .funs = function(x) if_else(is.na(x), 0, x))
                     
                     ####################                 
                     
                     counts <- counts %>% mutate(status = if_else(Count == length(object) , "all_data", if_else(Count == 0 , "no_data", "some_data")))
                     
                     #list of factor names
                     factors <- names(counts)[1:(dim(counts)[2]-2)]
                     
                     col.panel <- c("all_data", "some_data", "no_data")
                     names(col.panel) <- c("#00BA38", "orange", "red")
                     
                     col.panel.u <- col.panel[col.panel %in% unique(counts$status)]
                     
                     switch (length(factors),
                             "1" = { p <- ggplot(counts ,aes_string(x = factors[1], y = 1)) + 
                               theme(axis.text.y = element_blank()) + ylab("") },
                             "2" = { p <- ggplot(counts ,aes_string(x = factors[1], y = factors[2])) },
                             "3" = {
                               #get factor with min conditions -> to select for "facet_grid"
                               factors.l <- lapply(factors, function(x){ length(unique(counts[[x]])) }) %>% unlist()
                               names(factors.l) <- factors
                               factor.min <- names(factors.l[factors.l == min(factors.l)][1])
                               
                               factors <- factors[factors != factor.min]
                               
                               #add column to rename facet_grid
                               counts <- counts %>% mutate(grid = paste0(factor.min, "=",get(factor.min)))
                               
                               p <- ggplot(counts ,aes_string(x = factors[1], y = factors[2])) +
                                 facet_grid(grid~.) })
                     
                     p <- p + geom_tile(aes(fill = status), color = "white", size = 1, width = 1, height = 1)  + 
                       geom_text(aes(label = Count)) +
                       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             axis.ticks = element_blank(), axis.text.x=element_text(angle=90, hjust=1),
                             legend.position = "none")
                     return(p)
                     
                   })

# ---- plotExpDesignCompleteness ----
#' @title plotExpDesignCompleteness
#' @description This method checks some experimental design characteristics.
#'  A complete design and at least one biological and one batch factors are required for using RFLOMICS workflow.
#' @param An object of class \link{RflomicsSE-class}
#' @param sampleList list of samples to check.
#' @return plot
#' @exportMethod plotExpDesignCompleteness

methods::setMethod(f         = "plotExpDesignCompleteness",
                   signature = "RflomicsSE",
                   definition <- function(object, sampleList=NULL){
                     
                     # reduce object to sample list
                     object <- RFLOMICS::runSampleFiltering(object, samples = sampleList)
                     
                     check <- checkExpDesignCompleteness(object)
                     
                     if(isTRUE(check$error)) return(NULL)
                     
                     # Only works with bio and batch factors for the rest of the function
                     ExpDesign <- getDesignMat(object)
                     bio.fact <- getBioFactors(object)
                     
                     group_count <- .countSamplesPerCondition(ExpDesign, bio.fact)

                     plot <- .plotExperimentalDesign(counts = group_count, message= check$message)
                     
                     return(plot)
                   })

#' @title plotExpDesignCompleteness
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod plotExpDesignCompleteness
#' @noRd
methods::setMethod(f         = "plotExpDesignCompleteness",
                   signature = "RflomicsMAE",
                   definition <- function(object, omicName, sampleList=NULL){
                     
                     SEObject <- getRflomicsSE(object, omicName)
                     
                     plotExpDesignCompleteness(SEObject, sampleList = sampleList)
                   })

