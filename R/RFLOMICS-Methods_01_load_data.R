### ============================================================================
### [01_Load_Data] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------


#' @importFrom dplyr full_join mutate arrange select group_by_at count left_join
#' right_join mutate_at if_else 
#' @importFrom ggplot2 ggplot aes element_blank element_text geom_col theme 
#' labs scale_y_continuous geom_tile scale_fill_manual ylab xlab labs
#' facet_grid 
#' @importFrom purrr reduce


#' @name RflomicsMAE-accessors
#' @title A group of functions to access and modify the metadata slot of an object of class 
#' \link{RflomicsMAE}
#'  
#' @aliases getProjectName getDesignMat getDatasetNames getOmicsTypes getFactorNames 
#' getFactorTypes getBioFactors getBatchFactors getMetaFactors getFactorModalities getRflomicsSE 
#' subRflomics getModelFormula getSelectedContrasts getValidContrasts getContrastMatrix 
#' getTransSettings getFilterSettings getFilteredFeatures getNormSettings getCoeffNorm 
#' getDiffSettings getCoexpSettings getClusterEntities getCoseqClusters
#'
#' @description A set of getters and setters generic functions to access and modify 
#' objects of the slot metadata of a \link{RflomicsMAE} object. 
#' 
#' Getter methods:
#' 
#' \itemize{
#'    \item getProjectName: vector. Name of the project.
#'    \item getDesignMat: matrix: Get the design matrix associated the omic analysis.
#'    \item getDatasetNames: Names of the loaded omics datasets.
#'    \item getOmicsTypes: Access to the type of omics datasets: "RNAseq", "Proteomics", "Metabolomics"
#'    \item getFactorNames: List the experimental factor names.
#'    \item getFactorTypes: List the experimental factor types ("bio", "batch" or "meta") containing in the design list of the metadata.
#'    \item getBioFactors: List the biological factors.
#'    \item getBatchFactors: List the batch factors.
#'    \item getMetaFactors: List the metadata factors.
#'    \item getFactorModalities: Access to the factor modalities 
#'    \item getRflomicsSE: Access to a particular RflomicsSE object 
#'    \item subRflomics: Extract a subset of a RflomicsMAE object.
#'    \item getModelFormula: Access to the model formula of the statistical analysis
#'    \item getSelectedContrasts: List the selected contrasts 
#'    \item getValidContrasts: List the valid contrasts
#'    \item getContrastMatrix: Get the contrasts matrix 
#'    \item getTransSettings: List: Transformation settings of a given omic dataset 
#'    \item getFilterSettings: Access to the filtering settings of a given omic dataset 
#'    \item getFilteredFeatures: Access to the filtered features of a given omic dataset 
#'    \item getNormSettings: Access to the normalization settings of a given omic dataset 
#'    \item getCoeffNorm: Access to the normalization coefficients of a given omic dataset 
#'    \item getDiffSettings: Access to the differential expression analysis settings of a given omic dataset
#'    \item getCoexpSettings : Access to the co-expression analysis settings of a given omic dataset
#' }
#' 
#' Setter methods:
#' 
#' \itemize{
#'     \item setModelFormula: Set the model formula stored in \code{metadata} slot.
#'     \item setSelectedContrasts: Set the selected contrasts stored in \code{metadata} slot.
#'     \item setValidContrasts: Set the valid contrasts stored in \code{metadata} slot.
#'     \item setTrans: Set the transformation settings stored in \code{metadata} slot.
#'     \item setNorm: Set the normalization settings stored in \code{metadata} slot.
#'     \item setCoeffNorm: Set the normalization settings stored in \code{metadata} slot.
#' }
#' @param object A \code{RflomicsMAE} object
#'
#' @return See the itemized list in the description section for details
#' @return Setters: A \code{RflomicsMAE}  object
NULL


#' @name RflomicsSE-accessors
#' @title A group of functions to access and modify the metadata slot of an object of class 
#' \link{RflomicsSE}
#' 
#' @aliases omicType design DataProcessing PCAlist DiffExpAnal CoExpAnal DiffExpEnrichAnal CoExpEnrichAnal
#' getOmicType getDesign getDesignMat getDatasetNames getOmicsTypes getFactorNames getFactorTypes 
#' getBioFactors getBatchFactors getMetaFactors getFactorModalities getRflomicsSE subRflomics 
#' getModelFormula getSelectedContrasts getValidContrasts getContrastMatrix getTransSettings
#' getFilterSettings getFilteredFeatures getNormSettings getCoeffNorm getDiffSettings getCoexpSettings
#' getCoseqClusters getClusterEntities
#' 
#' @description A set of getters and setters generic functions to modify 
#' a \code{RflomicsSE} object: the metadata slot. 
#'
#' Accessor methods
#' 
#' \itemize{
#'    \item getDesignMat: data.frame: Design matrix associated to the omic analysis.
#'    \item getDatasetNames: vector: Names of the loaded omics datasets.
#'    \item getOmicsTypes: vector. Access to the type of omics datasets: "RNAseq", "Proteomics", "Metabolomics"
#'    \item getFactorNames: vector: Experimental factor names.
#'    \item getFactorTypes: named vector: Experimental factor types ("bio", "batch" or "meta") containing in the design list of the metadata.
#'    \item getBioFactors: list: Biological factors.
#'    \item getBatchFactors: list. Bbatch factors.
#'    \item getMetaFactors: list. Metadata factors.
#'    \item getFactorModalities: Access to the factor modalities 
#'    \item getRflomicsSE: Access to a particular RflomicsSE object by giving its name in argument. 
#'    \item subRflomics: Extract a subset of a RflomicsMAE object.
#'    \item getModelFormula: Access to the model formula of the statistical analysis
#'    \item getSelectedContrasts: List the selected contrasts 
#'    \item getValidContrasts: List the valid contrasts
#'    \item getContrastMatrix: Get the contrasts matrix 
#'    \item getTransSettings: List: Transformation settings of a given omic dataset 
#'    \item getFilterSettings: List: Access to the filtering settings of a given omic dataset 
#'    \item getFilteredFeatures: Access to the filtered features of a given omic dataset 
#'    \item getNormSettings: Access to the normalization settings of a given omic dataset 
#'    \item getCoeffNorm: Access to the normalization coefficients of a given omic dataset 
#'    \item getDiffSettings: Access to the differential expression analysis settings of a given omic dataset
#'    \item getCoexpSettings : Access to the co-expression analysis settings of a given omic dataset
#' }
#' 
#' Setter methods:
#' 
#' \itemize{
#'     \item setTrans: Set the transformation settings stored in \code{metadata} slot.
#'     \item setNorm: Set the normalization settings stored in \code{metadata} slot.
#'     \item setCoeffNorm: Set the normalization settings stored in \code{metadata} slot.
#' }
#'
#' @param object A \code{RflomicsSE} object
#'
#' @return See the itemized list in the description section for details
#' @return Setters: A \code{RflomicsSE} object
NULL


# ---- ACCESSORS ----
# ---- getProjectName ----
# ---- getDesignMat:    get colData ----
#' @exportMethod getProjectName
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getProjectName",
          signature  = "RflomicsMAE",
          definition <- function(object){
              return(object@metadata$projectName)
          })


# ---- getDesignMat:    get colData ----
#' @exportMethod getDesignMat
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getDesignMat",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              return(as.data.frame(object@colData))
          })



#' @exportMethod getDesignMat
#' @rdname RflomicsSE-accessors
setMethod(f          = "getDesignMat",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(as.data.frame(object@colData))
          })

# ---- getDatasetNames: get experiment names from RflomicsMAE object ----

#' @exportMethod getDatasetNames
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getDatasetNames",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              ExperimentNames <- unlist(object@metadata$omicList)
              names(ExperimentNames) <- NULL
              
              return(ExperimentNames)
          })


#' @exportMethod getDatasetNames
#' @rdname RflomicsSE-accessors
setMethod(f          = "getDatasetNames",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(names(object@metadata$omicType))
          })

# ---- getOmicsTypes:   get omics types of experiments from RflomicsMAE object ----


#' @importFrom purrr reduce
#' @exportMethod getOmicsTypes
#' @rdname RflomicsMAE-accessors
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


#' @exportMethod getOmicsTypes
#' @rdname RflomicsSE-accessors
setMethod(f          = "getOmicsTypes",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(object@metadata$omicType)
          })


# ---- getFactorNames:  get Factor names ----
#' @exportMethod getFactorNames
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getFactorNames",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              return(names(object@metadata$design$Factors.Type))
          })


#' @exportMethod getFactorNames
#' @rdname RflomicsSE-accessors
setMethod(f          = "getFactorNames",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              return(names(object@metadata$design$factorType))
          })


# ---- getFactorTypes:  get Factor types ----
#' @exportMethod getFactorTypes
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getFactorTypes",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              return(object@metadata$design$Factors.Type)
          })



#' @exportMethod getFactorTypes
#' @rdname RflomicsSE-accessors
setMethod(f          = "getFactorTypes",
          signature  = "RflomicsSE",
          definition <- function(object){
              return(object@metadata$design$factorType)
          })


# ---- getBioFactors:   get bio factor ----
#' @exportMethod getBioFactors
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getBioFactors",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BIO"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

#' @exportMethod getBioFactors
#' @rdname RflomicsSE-accessors
setMethod(f          = "getBioFactors",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BIO"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

# ---- getBatchFactors: get batch factor names ----
#' @exportMethod getBatchFactors
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getBatchFactors",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BATCH"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })


#' @exportMethod getBatchFactors
#' @rdname RflomicsSE-accessors
setMethod(f          = "getBatchFactors",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "BATCH"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

# ---- getMetaFactors:  get meta Factor names ----
#' @exportMethod getMetaFactors
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getMetaFactors",
          signature  = "RflomicsMAE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "META"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })


#' @exportMethod getMetaFactors
#' @rdname RflomicsSE-accessors
setMethod(f          = "getMetaFactors",
          signature  = "RflomicsSE",
          definition <- function(object){
              
              factVect <- toupper(getFactorTypes(object))
              res <- names(factVect)[factVect == "META"]
              
              if(length(res) == 0) return(NULL)
              return(res)
          })

# ---- getRflomicsSE:   get RflomicsSE object of one omic dataset ----
#' @exportMethod getRflomicsSE
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getRflomicsSE",
          signature  = "RflomicsMAE",
          definition <- function(object, datasetName = NULL){
              
              if(is.null(datasetName)) return(NULL)
              
              return(object[[datasetName]])
          })

# ---- getFactorModalities: ----
#' @exportMethod getFactorModalities
#' @rdname RflomicsMAE-accessors
setMethod(f          = "getFactorModalities",
          signature  = "RflomicsMAE",
          definition <- function(object, factorName){
              
              if(is.null(factorName)) return(NULL)
              if(!factorName %in% getFactorNames(object)) return(NULL) 
              
              return(levels(getDesignMat(object)[[factorName]]))
          })


#' @exportMethod getFactorModalities
#' @rdname RflomicsSE-accessors
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
#' @exportMethod subRflomicsMAE
#' @rdname RflomicsMAE-accessors
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
#' @param object An object of class \link{RflomicsSE-class}
#' @param sampleList list of samples to check.
#' @return a string with message
#' @exportMethod checkExpDesignCompleteness
#' @rdname checkExpDesignCompleteness
#' @aliases checkExpDesignCompleteness
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
#' @param omicName the name of the data the normalization have to be applied to. 
#' @exportMethod checkExpDesignCompleteness
#' @aliases checkExpDesignCompleteness
#' @rdname checkExpDesignCompleteness
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
#' @aliases plotDataOverview
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
#' @param object An object of class \link{RflomicsMAE-class}
#' @param omicNames a vector with list of dataset names
#' @return a gg plot object
#' \itemize{
#'  \item{"plot:"}{ plot of count data.frame.}
#'  }
#'  
#' @exportMethod plotConditionsOverview
#' @aliases plotConditionsOverview

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
#' @param object An object of class \link{RflomicsSE-class}
#' @param sampleList list of samples to check.
#' @return plot
#' @exportMethod plotExpDesignCompleteness
#' @rdname plotExpDesignCompleteness
#' @aliases plotExpDesignCompleteness

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
#' @aliases plotExpDesignCompleteness
#' @rdname plotExpDesignCompleteness
setMethod(f         = "plotExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){
              
              SEObject <- getRflomicsSE(object, omicName)
              
              plotExpDesignCompleteness(SEObject, sampleList = sampleList)
          })

