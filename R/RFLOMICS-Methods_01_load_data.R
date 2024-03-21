### ============================================================================
### [01_Load_Data] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------




#' @importFrom dplyr full_join mutate arrange select group_by_at count left_join
#' right_join mutate_at if_else 
#' @importFrom ggplot2 ggplot aes element_blank element_text geom_col theme 
#' labs scale_y_continuous geom_tile scale_fill_manual ylab xlab labs
#' facet_grid 
#' @importFrom purrr reduce


#' @name RflomicsMAE-methods
#' @title Access or modify the metadata slot of an object of class RflomicsMAE
#'
#' @description A set of accessor and setter generic functions to access and modify 
#' objects of the slot metadata of a \link{RflomicsMAE} object. 
#'
#' @section Accessors:
#' \itemize{
#'    \item getProjectName: Access to the name of the omic analysis.
#'    \item getDesignMat: Design matrix of the omic analysis.
#'    \item getDatasetNames: Names of the analyzed omics datasets.
#'    \item getOmicsTypes: Access to the type of omics datasets.
#'    \item getFactorNames: Access to the experimental factor names containing in the design list of the metadata.
#'    \item getFactorTypes: List the experimental factor types ("bio", "batch" or "meta") containing in the design list of the metadata.
#'    \item getBioFactors: List the biological factors.
#'    \item getBatchFactors: List the batch factors.
#'    \item getMetaFactors: List the metadata factors
#'    \item getFactorModalities: Access to the factor modalities 
#'    \item getRflomicsSE: Access to a RflomicsSE object 
#'    \item subRflomics: Extract a subset of a RflomicsMAE object.
#'    \item getModelFormula: Access to the model formula of the experimental design
#'    \item getSelectedContrasts: List the selected contrasts 
#'    \item getValidContrasts: List the valid contrasts
#'    \item getContrastMatrix: Get the contrast matrix 
#'    \item getTransSettings: Access to the transformation settings of a given omic dataset.
#'    \item getFilterSettings: Access to the filtering settings of a given omic dataset.
#'    \item getFilteredFeatures: Access to the filtered features of a given omic dataset.
#'    \item getNormSettings: Access to the normalization settings of a given omic dataset.
#'    \item getCoeffNorm: Access to the normalization coefficients of a given omic dataset
#'    \item getDiffSettings: Access to the differential expression analysis settings of a given omic dataset
#'    \item getCoexpSettings : Access to the co-expression analysis settings of a given omic dataset
#'    \item getClusterEntities: Access to the cluster entities of a given omic dataset
#'    \item getCoseqClusters: Access to the co-expression clusters of a given omic dataset
#' }
#'
#' @section Setters:
#' Setter method values (i.e., '\code{function(x) <- value}'):
#' \itemize{
#'     \item setModelFormula: Set the model formula stored in \code{metadata} slot.
#'     \item setSelectedContrasts: Set the selected contrasts stored in \code{metadata} slot.
#'     \item setValidContrasts: Set the valid contrasts stored in \code{metadata} slot.
#'     \item setTrans: Set the transformation settings stored in \code{metadata} slot.
#'     \item setNorm: Set the normalization settings stored in \code{metadata} slot.
#' }
#' @param object A \code{RflomicsMAE} object
#'
#'
#' @return Accessors: Either a \code{sampleMap}, \code{ExperimentList}, or
#' \code{DataFrame} object
#' @return Setters: A \code{RflomicsMAE}  object
NULL


#' @name RflomicsSE-methods
#' @title Accessing and modifying information in RflomicsSE
#'
#' @description A set of accessor and setter generic functions to modify 
#' objects from the slot metadata of a \code{RflomicsSE} object. 
#' The slot metadata contains: 
#'
#' \itemize{
#'   \item omicType: the type of omics dataset
#'   \item design: a list of design
#'   \item DataProcessing: a list of data processing settings and results 
#'   \item PCAlist: a list of PCA settings and results
#'   \item DiffExpAnal: a list of DiffExpAnal settings and results
#'   \item CoExpAnal: a list of CoExpAnal settings and results
#'   \item DiffExpEnrichAnal: a list of DiffExpEnrichAnal settings and results
#'   \item CoExpEnrichAnal: a list of CoExpEnrichAnal settings and results
#'  } 
#'
#' @section Accessors:
#' \itemize{
#'    \item getDesignMat: Design matrix of the omic analysis.
#'    \item getDatasetNames: Names of the analyzed omics datasets.
#'    \item getOmicsTypes: Access to the type of omics datasets.
#'    \item getFactorNames: Access to the experimental factor names containing in the design list of the metadata.
#'    \item getFactorTypes: List the experimental factor types ("bio", "batch" or "meta") containing in the design list of the metadata.
#'    \item getBioFactors: List the biological factors.
#'    \item getBatchFactors: List the batch factors.
#'    \item getMetaFactors: List the metadata factors
#'    \item getFactorModalities: Access to the factor modalities 
#'    \item getRflomicsSE: Access to a RflomicsSE object 
#'    \item subRflomics: Extract a subset of a RflomicsMAE object.
#'    \item getModelFormula: Access to the model formula of the experimental design
#'    \item getSelectedContrasts: List the selected contrasts 
#'    \item getValidContrasts: List the valid contrasts
#'    \item getContrastMatrix: Get the contrast matrix 
#'    \item getTransSettings: Access to the transformation settings of a given omic dataset.
#'    \item getFilterSettings: Access to the filtering settings of a given omic dataset.
#'    \item getFilteredFeatures: Access to the filtered features of a given omic dataset.
#'    \item getNormSettings: Access to the normalization settings of a given omic dataset.
#'    \item getCoeffNorm: Access to the normalization coefficients of a given omic dataset
#'    \item getDiffSettings: Access to the differential expression analysis settings of a given omic dataset
#'    \item getCoexpSettings : Access to the co-expression analysis settings of a given omic dataset
#'    \item getClusterEntities: Access to the cluster entities of a given omic dataset
#'    \item getCoseqClusters: Access to the co-expression clusters of a given omic dataset
#' }
#'
#' @section Setters:
#' Setter method values (i.e., '\code{function(x) <- value}'):
#' \itemize{
#'     \item setValid
#' }
#' @param object,x A \code{RflomicsMAE} object
#'
#' @return Accessors: an object stored in the metadata slot of the \code{RflomicsSE} object
#' @return Setters: A \code{RflomicsSE} object
NULL


# ---- ACCESSORS ----
# ---- getProjectName ----
# ---- getDesignMat:    get colData ----
#' @exportMethod getProjectName
#' @return a character string with the project name
#' @rdname RflomicsMAE-methods
setMethod(f          = "getProjectName",
          signature  = "RflomicsMAE",
          definition <- function(object){
              return(object@metadata$projectName)
          })


# ---- getDesignMat:    get colData ----
#' @return a data frame with design matrix
#' @exportMethod getDesignMat
#' @rdname RflomicsMAE-methods
setMethod(f          = "getDesignMat",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              return(as.data.frame(object@colData))
          })


#' @param object An object of class \link{RflomicsSE-class}
#' @return a data frame with design matrix
#' @exportMethod getDesignMat
#' @rdname RflomicsSE-methods
setMethod(f          = "getDesignMat",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(as.data.frame(object@colData))
          })

# ---- getDatasetNames: get experiment names from RflomicsMAE object ----

#' @return a vector with omic dataset names
#' @exportMethod getDatasetNames
#' @rdname RflomicsMAE-methods
setMethod(f          = "getDatasetNames",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              ExperimentNames <- unlist(object@metadata$omicList)
              names(ExperimentNames) <- NULL
              
              return(ExperimentNames)
          })


#' @return a vector with the omic dataset name
#' @exportMethod getDatasetNames
#' @rdname RflomicsSE-methods
setMethod(f          = "getDatasetNames",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(names(object@metadata$omicType))
          })

# ---- getOmicsTypes:   get omics types of experiments from RflomicsMAE object ----


#' @importFrom purrr reduce
#' @return names vector
#' @exportMethod getOmicsTypes
#' @rdname RflomicsMAE-methods
setMethod(f          = "getOmicsTypes",
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


#' @return named vector
#' @exportMethod getOmicsTypes
#' @rdname RflomicsSE-methods
setMethod(f          = "getOmicsTypes",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(object@metadata$omicType)
          })


# ---- getFactorNames:  get Factor names ----
#' @return a vector with the names of the design factors
#' @exportMethod getFactorNames
#' @rdname RflomicsMAE-methods
setMethod(f          = "getFactorNames",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              return(names(object@metadata$design$Factors.Type))
          })


#' @return a vector with names of design factors
#' @exportMethod getFactorNames
#' @rdname RflomicsSE-methods
setMethod(f          = "getFactorNames",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(names(object@metadata$design$factorType))
          })


# ---- getFactorTypes:  get Factor types ----
#' @return a named vector
#' @exportMethod getFactorTypes
#' @rdname RflomicsMAE-methods
setMethod(f          = "getFactorTypes",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              return(object@metadata$design$Factors.Type)
          })



#' @return a named vector
#' @exportMethod getFactorTypes
#' @rdname RflomicsSE-methods
setMethod(f          = "getFactorTypes",
          signature  = "RflomicsSE",
          definition <- function(object){
              return(object@metadata$design$factorType)
          })


# ---- getBioFactors:   get bio factor ----
#' @return a vector with biological factors
#' @exportMethod getBioFactors
#' @rdname RflomicsMAE-methods
setMethod(f          = "getBioFactors",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BIO"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

#' @return a vector with biological factors
#' @exportMethod getBioFactors
#' @rdname RflomicsSE-methods
setMethod(f          = "getBioFactors",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BIO"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

# ---- getBatchFactors: get batch factor names ----
#' @return a vector with batch factors
#' @exportMethod getBatchFactors
#' @rdname RflomicsMAE-methods
setMethod(f          = "getBatchFactors",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BATCH"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })


#' @return a vector with batch factors
#' @exportMethod getBatchFactors
#' @rdname RflomicsSE-methods
setMethod(f          = "getBatchFactors",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BATCH"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

# ---- getMetaFactors:  get meta Factor names ----
#' @return a vector with batch factors
#' @exportMethod getMetaFactors
#' @rdname RflomicsMAE-methods
setMethod(f          = "getMetaFactors",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "META"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })


#' @return a vector with batch factors
#' @exportMethod getMetaFactors
#' @rdname RflomicsSE-methods
setMethod(f          = "getMetaFactors",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "META"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

# ---- getRflomicsSE:   get RflomicsSE object of one omic dataset ----
#' @param datasetName name of omic dataset
#' @return An object of class \link{RflomicsSE-class}
#' @exportMethod getRflomicsSE
#' @rdname RflomicsMAE-methods
setMethod(f          = "getRflomicsSE",
          signature  = "RflomicsMAE",
          definition <- function(object, datasetName = NULL){
              
              if(is.null(datasetName)) return(NULL)
              
              return(object[[datasetName]])
          })

# ---- getFactorModalities: ----
#' @param factorName a name of factor
#' @return a vector
#' @exportMethod getFactorModalities
#' @rdname RflomicsMAE-methods
setMethod(f          = "getFactorModalities",
          signature  = "RflomicsMAE",
          definition <- function(object, factorName){
              
              if(is.null(factorName)) return(NULL)
              if(!factorName %in% getFactorNames(object)) return(NULL) 
              
              return(levels(getDesignMat(object)[[factorName]]))
          })


#' @param factorName a name of factor
#' @return a vector
#' @exportMethod getFactorModalities
#' @rdname RflomicsSE-methods
setMethod(f          = "getFactorModalities",
          signature  = "RflomicsSE",
          definition <- function(object, factorName){
              
              if(is.null(factorName)) return(NULL)
              if(!factorName %in% getFactorNames(object)) return(NULL) 
              
              return(levels(getDesignMat(object)[[factorName]]))
          })

# ---- get names of executed analysis ----
# ---- METHODS ----
# ---- subRflomicsMAE:  subset a RflomicsMAE from ----
#' @param omicNames a vector with the dataset names
#' @return object An object of class \link{RflomicsMAE-class}
#' @exportMethod subRflomicsMAE
#' @rdname RflomicsMAE-methods
setMethod(f          = "subRflomicsMAE",
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
#'  A complete design (all combinations of factor modalities with at least 2 replicates for each have to be present) with
#'  at least one biological and one batch factors are required to use the RFLOMICS workflow.
#' @param An object of class \link{RflomicsSE-class}
#' @param sampleList list of samples to check.
#' @return a string with message
#' @exportMethod checkExpDesignCompleteness
#' @rdname checkExpDesignCompleteness
setMethod(f         = "checkExpDesignCompleteness",
          signature = "RflomicsSE",
          definition <- function(object, sampleList=NULL){
              
              object <- runSampleFiltering(object, samples = sampleList)
              
              output <- list()
              output[["error"]] <- FALSE
              
              # Only works with bio and batch factors for the rest of the function
              ExpDesign <- getDesignMat(object)
              bio.fact  <- getBioFactors(object)
              
              # check presence of bio factors
              if (!length(getBioFactors(object)) %in% seq_len(3)){ 
                  output[["messages"]] <-  "Error: You need at least 1 biological factor with at least 2 modalities."
                  output[["error"]]    <- TRUE
                  return(output)
              }
              # check presence of bash factors
              if (!length(getBatchFactors(object)) %in% c(1,2)){ 
                  output[["messages"]] <-  "Error: You need at least 1 batch factor with at least 2 replicats."
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
                  
                  output[["messages"]] <- "Error: The experimental design is not complete."
                  output[["error"]]    <- TRUE
              }
              else if(min(group_count$Count) == 1){
                  
                  output[["messages"]] <-  "Error: You need at least 2 biological replicates."
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
setMethod(f         = "checkExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){
              
              if(is.null(omicName)) stop("Arg omicName missed")
              
              SEObject <- getRflomicsSE(object, omicName)
              
              checkExpDesignCompleteness(SEObject, sampleList = sampleList)
          })



# ---- PLOTS ----
# ---- plotDataOverview ----
#' @title Dataset overview plot
#' @description This function plot an overview of the loaded datasets displaying per sample 
#' (n=number of entities (genes/metabolites/proteins); k=number of samples)
#' @param object An object of class \link{RflomicsMAE-class}
#' @param omicNames a vector with dataset names
#' @param realSize booleen value
#' @exportMethod plotDataOverview
#' @return plot
#' @rdname plotDataOverview
setMethod(f         = "plotDataOverview",
          signature = "RflomicsMAE",
          definition <- function(object, omicNames=NULL, realSize=FALSE){
              
              if(length(object) == 0) stop("object is NULL")
              
              object <- subRflomicsMAE(object, omicNames)
              
              if(is.null(object)) return(NULL)
              
              Groups <- getDesignMat(object)
              
              nb_entities <- lapply(names(object), function(SE){ dim(object[[SE]])[1] }) %>% unlist()
              names(nb_entities) <- names(object)
              
              data <- data.frame(nb_entities = nb_entities, assay = names(nb_entities)) %>%
                  full_join(data.frame(object@sampleMap), by="assay") %>%
                  mutate(y.axis = paste0(assay, "\n", "n=", nb_entities)) %>% arrange(primary)
              
              data$primary <- factor(data$primary, levels = levels(Groups$samples)) 
              
              nb_entities_ord <- select(data, y.axis, nb_entities) %>% unique() %>% arrange(desc(nb_entities))
              nb_entities_ord$nb_entities <- log(nb_entities_ord$nb_entities)
              tmp.vec <- c(0)
              breaks  <- vector()
              for(i in seq_len(length(nb_entities_ord$nb_entities))){ 
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
                              labs(x=paste0("Samples (k=", length(unique(object@sampleMap$primary)), ")"), y="") +
                              scale_y_continuous(breaks = (breaks), labels = nb_entities_ord$y.axis)
                          
                      },
                      "FALSE" = {
                          p <- ggplot(data, aes(x=primary, y=y.axis)) +
                              geom_tile(aes(fill = y.axis), colour = "grey50") +
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.ticks = element_blank(), legend.position="none",
                                    axis.text.x = element_text(angle = 90, hjust = 1)) +
                              labs(x=paste0("Samples (k=", length(unique(object@sampleMap$primary)), ")"), y="")
                      }
              )
              return(p)
          })



# ---- plotConditionsOverview ----
#' @title plotConditionsOverview
#' @description 
#'  A complete design and at least one biological and one batch factors are 
#' required for using RFLOMICS workflow.
#' @param An object of class \link{RflomicsMAE-class}
#' @param omicNames a vector with list of dataset names
#' @return a gg plot object
#' \itemize{
#'  \item{"plot:"}{ plot of count data.frame.}
#'  }
#'  
#' @exportMethod plotConditionsOverview
#' @noRd

setMethod(f         = "plotConditionsOverview",
          signature = "RflomicsMAE",
          definition <- function(object, omicNames = NULL){
              
              # check presence of bio factors
              if (!length(getBioFactors(object)) %in% seq_len(3)){ 
                  stop("No bio factor! or nbr of bio factors exceed 3!") }
              if (!length(getBatchFactors(object)) %in% c(1,2)){ 
                  stop("No replicates found!") }
              
              ####################
              
              BioFact <- getBioFactors(object)
              coldata <- getDesignMat(object) %>%
                  mutate(samples=rownames(.))
              #coldata <- tibble::as_tibble(coldata)
              coldata <- sampleMap(object) %>% as.data.frame() %>% 
                  left_join(coldata, by = c("primary" = "samples"))
              
              all_combin_cond <- lapply(BioFact, function(x){ 
                  df <- unique(coldata[x])
                  rownames(df) <- seq_len(nrow(df))
                  return(df) 
              }) %>% reduce(merge)
              
              counts <- coldata %>% select(assay, all_of(BioFact)) %>% 
                  unique() %>% 
                  group_by_at(BioFact) %>% count(name = "Count") %>% 
                  right_join(all_combin_cond, by = BioFact) %>% 
                  mutate_at(.vars = "Count", .funs = function(x) {
                      if_else(is.na(x), 0, x)})
              
              
              ####################                 
              
              counts <- counts %>% 
                  mutate(status = if_else(Count == length(object) , "all_data", if_else(Count == 0 , "no_data", "some_data")))
              
              #list of factor names
              factors <- names(counts)[seq_len(dim(counts)[2]-2)]
              
              col.panel <- c("all_data", "some_data", "no_data")
              names(col.panel) <- c("#00BA38", "orange", "red")
              
              col.panel.u <- col.panel[col.panel %in% unique(counts$status)]
              
              switch (length(factors),
                      "1" = { p <- ggplot(counts, aes_string(x = factors[1], y = 1)) + 
                          theme(axis.text.y = element_blank()) + ylab("") },
                      "2" = { p <- ggplot(counts, aes_string(x = factors[1], y = factors[2])) },
                      "3" = {
                          #get factor with min conditions -> to select for "facet_grid"
                          factors.l <- lapply(factors, function(x){ length(unique(counts[[x]])) }) %>% unlist()
                          names(factors.l) <- factors
                          factor.min <- names(factors.l[factors.l == min(factors.l)][1])
                          
                          factors <- factors[factors != factor.min]
                          
                          #add column to rename facet_grid
                          counts <- counts %>% 
                              mutate(grid = paste0(factor.min, "=",get(factor.min)))
                          
                          p <- ggplot(counts ,aes_string(x = factors[1], y = factors[2])) +
                              facet_grid(grid~.) })
              
              p <- p + geom_tile(aes(fill = status), 
                                 color = "white", size = 1,
                                 width = 1, height = 1)  + 
                  geom_text(aes(label = Count)) +
                  theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        axis.ticks = element_blank(), 
                        axis.text.x=element_text(angle=90, hjust=1),
                        legend.position = "none")
              return(p)
              
          })

# ---- plotExpDesignCompleteness ----
#' @title plotExpDesignCompleteness
#' @description This method checks that experimental design constraints are satisfied and plot 
#' a summary of the design.
#' A complete design (all combinations of factor modalities with at least 2 replicates for each have to be present) 
#' with at least one biological and one batch factors are required to use the RFLOMICS workflow.
#' @param An object of class \link{RflomicsSE-class}
#' @param sampleList list of samples to check.
#' @return plot
#' @exportMethod plotExpDesignCompleteness
#' @rdname plotExpDesignCompleteness

setMethod(f         = "plotExpDesignCompleteness",
          signature = "RflomicsSE",
          definition <- function(object, sampleList=NULL){
              
              # reduce object to sample list
              object <- runSampleFiltering(object, samples = sampleList)
              
              check <- checkExpDesignCompleteness(object)
              
              # Only works with bio and batch factors for the rest of the 
              # function
              ExpDesign <- getDesignMat(object)
              bio.fact <- getBioFactors(object)
              
              group_count <- .countSamplesPerCondition(ExpDesign, bio.fact)
              
              plot <- .plotExperimentalDesign(counts = group_count, 
                                              message= check$message)
              
              return(plot)
          })

#' @title plotExpDesignCompleteness
#' @param omicName a character string with the name of the dataset
#' @param sampleList list of samples to take into account
#' @exportMethod plotExpDesignCompleteness
#' @noRd
setMethod(f         = "plotExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){
              
              SEObject <- getRflomicsSE(object, omicName)
              
              plotExpDesignCompleteness(SEObject, sampleList = sampleList)
          })

