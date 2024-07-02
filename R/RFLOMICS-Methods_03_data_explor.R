### ============================================================================
### [03_data_processing] accessors and methods for RflomicsMAE and RflomicsSE classes
### ----------------------------------------------------------------------------
# D. Charif, 
# N. Bessoltane, 
# A. Hulot

##==== DATA PROCESSING METHOD ====

###==== runDataProcessing ====

#' @title Data processing
#' @rdname runDataProcessing
#' @name runDataProcessing
#' @description
#' These functions applied a data processing (filtering, normalization 
#' and/or transformation, PCA) on RNAseq, proteomic, or metabolomic data.
#' 
#' runDataProcessing() calls the following functions:
#' @param object An object of class \link{RflomicsSE} or class \link{RflomicsSE} 
#' @param samples samples to keep.
#' @param lowCountFiltering_strategy strategy of RNAseq low count filtering. 
#' Mandatory for RNAseq data. Default value : "NbReplicates".
#' @param lowCountFiltering_CPM_Cutoff CPM cutoff for RNAseq low count filtering.
#'  Mandatory for RNAseq data. Default value : 1.
#' @param normMethod of normalisation. Mandatory for RNAseq data. 
#' Default value : RNAseq = TMM.
#' @param transformMethod method of transformation.
#' @return An object of class \link{RflomicsSE} or class \link{RflomicsSE} 
#' @exportMethod runDataProcessing
#' @seealso 
#' \link{checkExpDesignCompleteness} 
#' \link{getProcessedData} 
#' \link{getTransSettings} 
#' \link{getFilterSettings}
#' \link{getFiltredFeatures}
#' \link{getCoeffNorm}
#' \link{getNormSettings}
#' \link{plotLibrarySize}
#' \link{plotDataDistribution}
#' \link{plotOmicsPCA}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. 
#' DiCoExpress: a tool to process multifactorial
#' RNAseq experiments from quality controls to co-expression analysis through 
#' differential analysis based on contrasts inside GLM models. 
#' Plant Methods 16, 68 (2020).
#' 
#' ref2
#' @example inst/examples/runDataProcessing.R
setMethod(f          = "runDataProcessing",
          signature  = "RflomicsSE",
          definition = function(object, 
                                samples=NULL, 
                                lowCountFiltering_strategy = "NbReplicates", 
                                lowCountFiltering_CPM_Cutoff = 1, 
                                normMethod = "none", 
                                transformMethod = "none")
          {
            
            # keep selected samples
            message("[RFLOMICS] #    => select samples...")
            object <- runSampleFiltering(object, samples)
            
            if(nrow(getDesignMat(object)) == 0) stop("no samples in object!")
            
            # spported values:
            lowCountFiltering_strategy.sup     <- c("NbReplicates","NbConditions")
            transformation_method.sup          <- c("log2", "none")
            normalisation_method.abundance.sup <- c("median", "totalSum", "none")
            normalisation_method.count.sup     <- c("TMM")
            
            # apply data processing
            switch(object@metadata$omicType,
                   "RNAseq" = {
                     
                     # Filter low abundance
                     message("[RFLOMICS] #    => Low counts Filtering...")
                     if(is.null(lowCountFiltering_strategy)   || !lowCountFiltering_strategy %in% lowCountFiltering_strategy.sup) 
                       stop("the low count filtering strategy : ", lowCountFiltering_strategy, " isn't supported by rflomics package. Supported values : ",  paste(lowCountFiltering_strategy.sup, collapse = ", "))
                     
                     if(is.null(lowCountFiltering_CPM_Cutoff) || !is.numeric(lowCountFiltering_CPM_Cutoff)) 
                       stop(lowCountFiltering_CPM_Cutoff, " must be a integer value > 1")
                     
                     SE.processed <- filterLowAbundance(object = object, filterMethod= "CPM", filterStrategy = lowCountFiltering_strategy, cpmCutoff = lowCountFiltering_CPM_Cutoff)
                     
                     # Run Normalisation 
                     message("[RFLOMICS] #    => Counts normalization...")
                     if(is.null(normMethod) || normMethod != "TMM"){
                       normMethod <- "TMM"
                       warning("only ", normalisation_method.count.sup, 
                               " method is supported for ", 
                               object@metadata$omicType, " normalisation.")
                     }
                     
                     SE.processed <- runNormalization(SE.processed, normMethod = normMethod)
                   },
                   {
                     message("[RFLOMICS] #    => transformation data...")
                     if(is.null(transformMethod)) transformMethod <- "none"
                     if(! transformMethod %in% transformation_method.sup) 
                       stop("the transformation method ", transformMethod,
                            " is not support in rflomics package.Supported values : ", 
                            paste(transformation_method.sup, collapse = ", "))
                     SE.processed <- runTransformData(object, transformMethod = transformMethod)
                     
                     message("[RFLOMICS] #    => Run normalization...")
                     if(is.null(normMethod)) normMethod <- "none"
                     if(! normMethod %in% normalisation_method.abundance.sup) 
                       stop("the normalisation method ", normMethod,
                            " is not support in rflomics package. Supported values : ", 
                            paste(normalisation_method.abundance.sup, collapse = ", "))
                     SE.processed <- runNormalization(SE.processed, normMethod = normMethod)
                   }
            )
            
            # Run PCA for filtred & normalized data 
            message("[RFLOMICS] #    => Compute PCA ")
            SE.processed <- runOmicsPCA(SE.processed,ncomp = 5, raw = FALSE)  
            SE.processed@metadata$DataProcessing[["done"]] <- TRUE
            
            return(SE.processed)
          })


#' @rdname runDataProcessing
#' @param SE.name SE.name the name of the dataset if the input object 
#' is a \link{RflomicsMAE}
#' @exportMethod runDataProcessing
setMethod(f          = "runDataProcessing",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                samples=NULL, 
                                lowCountFiltering_strategy = "NbReplicates", 
                                lowCountFiltering_CPM_Cutoff = 1, 
                                normMethod= "none", 
                                transformMethod = "none", 
                                SE.name){
            
            # if paste0(SE.name, ".raw") exist
            if (!paste0(SE.name, ".raw") %in% names(object)){
              stop("SE name must be part of this list of names: ",
                   getDatasetNames(object))
            }
            
            # if paste0(SE.name, ".raw") exist
            if(is.null(object[[paste0(SE.name, ".raw")]])) {
              stop("SE name must be part of this list of names: ",
                   getDatasetNames(object))
            }
            
            SE.raw       <- object[[paste0(SE.name, ".raw")]]
            SE.processed <- 
              runDataProcessing(
                object = SE.raw,
                samples = samples, 
                lowCountFiltering_strategy = lowCountFiltering_strategy, 
                lowCountFiltering_CPM_Cutoff = lowCountFiltering_CPM_Cutoff,
                normMethod= normMethod, 
                transformMethod = transformMethod)
            
            
            # remove SE processed if exist
            if (SE.name %in% names(object)) object <- object[,, -which(names(object) == SE.name)]
            
            # add new SE with processed data
            object <- eval(parse(text = paste0('c( object ,', SE.name, ' = SE.processed )')))
            
            return(object)
          })

###==== filterLowAbundance ====

# METHOD to filter data

# Cette method est propre au RNASEQ => Est-ce que c'est vraiment ce que l'on souhaite ?
# Plutot qu'une fonction interface pour tous les omics ?
# Pourquoi ne pas avoir utilisée directement la fonction de edgeR ?
#' @rdname runDataProcessing
#' @name filterLowAbundance
#' @description 
#' \itemize{
#' \item{filterLowAbundance:}{ This function aims at removing transcript from the count data matrix of an omic of type "RNAseq".
#' by applying filtering criterion described in reference.}
#' }
#' @param filterMethod The filtering model ("CPM")
#' @param filterStrategy The filtering strategy ("NbConditions" or "NbReplicates")
#' @param cpmCutoff The CPM cutoff.
#' @details
#' filterLowAbundance(): By default, gene/transcript with 0 count are removed from the data. The function then
#' computes the count per million or read (CPM) for each gene in each sample and gives by
#' genes the number of sample(s) which are over the cpmCutoff (NbOfsample_over_cpm).
#' Then Two filtering strategies are proposed:
#' \itemize{
#' \item{NbConditions: }{keep gene if the NbOfsample_over_cpm >= NbConditions}
#' \item{NbReplicates: }{keep gene if the NbOfsample_over_cpm >= min(NbReplicat)}
#' \item{filterByExpr:} {the default filtering method implemented in the edgeR filterByExpr() function.}
#' }
#' @exportMethod filterLowAbundance
#' @importFrom edgeR DGEList filterByExpr cpm
#' @seealso edgeR::filterByExpr
setMethod(f         = "filterLowAbundance",
          signature = "RflomicsSE",
          definition = function(object, 
                                filterMethod= "CPM", 
                                filterStrategy = "NbConditions", 
                                cpmCutoff = 5){
            
            if (isFALSE(filterStrategy %in% c("NbReplicates","NbConditions"))) {
              stop("filterStrategy argument must be one of this tow options : NbReplicates or NbConditions")
            }
            
            objectFilt <- object
            assayFilt  <- assay(objectFilt)
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
                   "NbConditions" = { keep <- rowSums(cpm(assayFilt) >= cpmCutoff) >=  NbConditions },
                   "NbReplicates" = { keep <- rowSums(cpm(assayFilt) >= cpmCutoff) >=  min(NbReplicate) },
                   "filterByExpr" = { dge  <- DGEList(counts = assayFilt, genes = rownames(assayFilt))
                   keep <- filterByExpr(dge)
                   }
            )
            
            ## nbr of genes filtered
            genes_flt1  <- objectFilt[!keep]@NAMES
            
            object <- objectFilt[keep]
            
            Filtering <- list()
            Filtering[["setting"]][["method"]]           <- filterMethod
            Filtering[["setting"]][["filterStrategy"]]   <- filterStrategy
            Filtering[["setting"]][["cpmCutoff"]]        <- cpmCutoff
            Filtering[["results"]][["filteredFeatures"]] <- c(genes_flt0, genes_flt1)
            
            if (is.null(object@metadata$DataProcessing)) {
              object@metadata$DataProcessing <- list()
            }
            object@metadata$DataProcessing$Filtering <- Filtering
            
            return(object)
            
          })


#' @rdname runDataProcessing
#' @name filterLowAbundance
#' @exportMethod filterLowAbundance
setMethod(f          = "filterLowAbundance",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name, 
                                filterMethod= "CPM", 
                                filterStrategy = "NbConditions", 
                                cpmCutoff = 5){
            
            if (getOmicsTypes(object[[SE.name]]) == "RNAseq") {
              object[[SE.name]] <-  filterLowAbundance(object          = object[[SE.name]], 
                                                       filterStrategy = filterStrategy,
                                                       filterMethod   = filterMethod,
                                                       cpmCutoff      = cpmCutoff)
              return(object)
            }else{
              message("Can't apply this method to omics types other than RNAseq.")
            }
            
          })

###==== runSampleFiltering ====
# filtering per sample

#' @rdname runDataProcessing
#' @name runSampleFiltering
#' @description 
#' \itemize{
#' \item{runSampleFiltering:}
#' { This function applied sample filtering on an dataset.}
#' }
#' @exportMethod runSampleFiltering
#' @importFrom dplyr filter
setMethod(f          = "runSampleFiltering",
          signature  = "RflomicsSE",
          definition = function(object, 
                                samples=NULL) {
            
            # if no samples to filter
            if (is.null(samples)) return(object) 
            
            # check if samples overlap with
            samples <- intersect(colnames(object), samples)
            
            # if no iverlap with data colnames
            if (length(intersect(samples, colnames(object))) == 0) {
              warning("No overlap between samples and colnames of data")
              return(object)
            }
            
            # if 100% ovelap!
            if (length(samples) == length(colnames(object))) return(object)
            
            # keep selected samples
            # keep only samples in data matrix, and colData
            SE.new <- object[, object$samples %in% samples]
            
            # if we remove all samples :
            if (nrow(SE.new@colData) == 0) SE.new@colData <- SE.new@colData[c("samples", "groups")]
            
            for (factor in c(getBioFactors(object), getBatchFactors(object))){
              
              # if only one category remains after the filter, it's will be removed
              if (length(unique(SE.new@colData[[factor]])) <= 1 ) {
                SE.new@colData[[factor]] <- NULL
                factor.types <- getFactorTypes(SE.new)
                SE.new@metadata$design$factorType <- factor.types[which(names(factor.types) != factor)]
                # replace with setFactorTypes
              }
              else{
                F.levels <- levels(SE.new@colData[[factor]])
                SE.new@colData[[factor]] <- factor(SE.new@colData[[factor]], 
                                                   levels = intersect(F.levels, unique(SE.new@colData[[factor]])))
              }
            }
            colData <- as.data.frame(SE.new@colData)
            order_levels <- with(colData, do.call(order, colData[c(getBioFactors(SE.new), getBatchFactors(SE.new))]))
            SE.new$samples <- factor(SE.new$samples, levels = unique(SE.new$samples[order_levels]))
            SE.new$groups  <- factor(SE.new$groups,  levels = unique(SE.new$groups[order_levels]))
            
            # à retirer dès que je remplace Groups par colData
            #SE.new@metadata$Groups <- filter(SE.new@metadata$Groups, samples %in% SE.new$samples)
            
            SE.new@metadata$DataProcessing$filteredSamples <- setdiff(colnames(object), samples)
            
            selectedContrasts <- getSelectedContrasts(object)
            
            SE.new <- setSelectedContrasts(SE.new, selectedContrasts)
            
            return(SE.new)
          })


#' @rdname runDataProcessing
#' @exportMethod runSampleFiltering
setMethod(f          = "runSampleFiltering",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name,
                                samples=NULL) {
            
            object[[SE.name]] <- 
              runSampleFiltering(object[[SE.name]], samples = samples)
            
            return(object)
          })

###==== runTransformData ====
# METHOD to transform data

#' @rdname runDataProcessing
#' @name runTransformData
#' @description
#' \itemize{
#' \item{runTransformData:}
#' { This function applied a transformation on dataset. The transformation 
#' method is chosen according to dataset omic type 
#' (RNaseq: none, metabolomics/proteomics: log2)}
#' }
#' @param transformMethod The transformation to store in the metadata or to 
#' store and apply if modifyAssay is TRUE.
#' @param modifyAssay Boolean. Do the transformation need to be applied on 
#' the data? The raw data will be replaced by the transformed ones.
#' @exportMethod runTransformData
setMethod(f          = "runTransformData",
          signature  = "RflomicsSE",
          definition = function(object, 
                                transformMethod = NULL, 
                                modifyAssay = FALSE){
            
            transformation.setting <- getAnalysis(object,
                                                  name = "DataProcessing",
                                                  subName = "Transformation")
            
            if (is.null(transformMethod)) {
              if (!modifyAssay && getOmicsTypes(object) == "RNAseq") {
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
            
            transformation.setting[["setting"]][["method"]] <- transformMethod
            
            object <- setElementToMetadata(object,
                                           name = "DataProcessing",
                                           subName = "Transformation",
                                           content = transformation.setting)
            
            if (modifyAssay) {
              object <-  .applyTransformation(object)
            }
            
            return(object)
            
          })

#' @rdname runDataProcessing
#' @name runTransformData
#' @exportMethod runTransformData
setMethod(f          = "runTransformData",
          signature  = "RflomicsMAE",
          definition = function(object, 
                                SE.name, 
                                transformMethod = NULL, 
                                modifyAssay = FALSE){
            
            object[[SE.name]] <-  runTransformData(object[[SE.name]], 
                                                   transformMethod = transformMethod, 
                                                   modifyAssay = modifyAssay)
            
            return(object)
            
          })

###==== runNormalization ====
# METHOD to normalize data
# Function non generique pour les autres data

#' @rdname runDataProcessing
#' @name runNormalization
#' @description 
#' \itemize{
#' \item{runNormalization:}
#' { This function applied a normalization on a dataset. 
#' The normalisation method is chosen according to dataset omic type 
#' (RNAseq: TMM, metabolics/proteomics: median)}
#' }
#' @param normMethod Normalization method
#' @param modifyAssay Does the normalization have to be applied or just stored for later? Recommended it stays FALSE.
#' @return An object of class \link{RflomicsSE}
#' The applied normalization method and computed scaling factors (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "method" and "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{RflomicsSE} object.
#' @importFrom stats median
#' @exportMethod runNormalization
setMethod(f          = "runNormalization",
          signature  = "RflomicsSE",
          definition = function(object, normMethod = NULL, modifyAssay = FALSE){
            
            Groups     <- getDesignMat(object)
            
            if (.isNorm(object)) 
              warning("Data were already normalized before!")
            
            if (.isTransformed(object)) 
              warning("Data were transformed before!")
            
            if (is.null(normMethod)) {
              if (getOmicsTypes(object) == "RNAseq" && getTransSettings(object)$method == "none") {
                message("Using TMM normalization for RNAseq (counts) data")
                normMethod <- "TMM"
              } else {
                normMethod <- "none"
              }
            }
            
            if (getOmicsTypes(object) == "RNAseq" && normMethod != "TMM") {
              message("Forcing TMM normalization for RNAseq (counts) data")
              normMethod <- "TMM"
            }
            
            object2 <- object
            
            # Run normalization on transformed data (except for RNAseq data, transformation is expected to be "none" anyway)
            if (!.isTransformed(object) && getTransSettings(object)$method != "none")
              object2 <-  .applyTransformation(object2)
            
            switch(normMethod,
                   "TMM"        = {coefNorm  <- .tmmNormalization(assay(object2), Groups$groups) },
                   "median"     = {coefNorm  <- apply(assay(object2), 2, FUN = function(sample_vect) {median(sample_vect)}) },
                   "totalSum"   = {coefNorm  <- apply(assay(object2), 2, FUN = function(sample_vect) {sum(sample_vect^2)}) },
                   "none"       = {coefNorm  <- rep(1, ncol(assay(object2))) },
                   { message("Could not recognize the normalization method, applying 'none'. Please check your parameters.")
                     normMethod <- "none"
                     coefNorm  <-  rep(1, ncol(assay(object2)))
                   }
            )
            
            Normalization <- list()
            Normalization[["setting"]][["method"]]   <- normMethod
            Normalization[["results"]][["coefNorm"]] <- coefNorm
            Normalization[["normalized"]]            <- FALSE
            
            if (is.null(object@metadata$DataProcessing)) { 
              object@metadata$DataProcessing <- list()
            }
            object@metadata$DataProcessing[["Normalization"]] <- Normalization
            
            if (modifyAssay) object <- .applyNorm(object)
            
            return(object)
          })

#' @rdname runDataProcessing
#' @name runNormalization
#' @exportMethod runNormalization
setMethod(f          = "runNormalization",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, normMethod = NULL, modifyAssay = FALSE){
            
            object[[SE.name]] <-  runNormalization(object       = object[[SE.name]],
                                                   normMethod   = normMethod,
                                                   modifyAssay = modifyAssay)
            return(object)
          })

###==== runOmicsPCA ====

#' @title runOmicsPCA
#' @description 
#' \itemize{
#' \item{runOmicsPCA:}
#' {This function performs a principal component analysis on omic 
#' data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a 
#' "Normalization" slot is present in the metadata slot, then data are 
#' normalized before running the PCA according to the indicated transform 
#' method.}
#' }
#' This function performs a principal component analysis on omic 
#' data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a 
#' "Normalization" slot is present in the metadata slot, then data are 
#' normalized before running the PCA according to the indicated transform 
#' method.
#' @param object An object of class \link{RflomicsSE-class}.
#' @param ncomp Number of components to compute. Default is 5.
#' @param raw boolean. Does the pca have to be ran on raw data or transformed 
#' and normalized data? Default is FALSE, pca is ran on transformed and 
#' normalized data.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runOmicsPCA
#' @importFrom FactoMineR PCA
#' @rdname runOmicsPCA
#' @aliases runOmicsPCA
#'
setMethod(
  f          = "runOmicsPCA",
  signature  = "RflomicsSE",
  definition = function(object, ncomp = 5, raw = FALSE) {
    object2 <- .checkTransNorm(object, raw = raw)
    pseudo  <- assay(object2)
    if (raw)
      object@metadata[["PCAlist"]][["raw"]] <-
      PCA(t(pseudo), ncp = ncomp, graph = FALSE)
    else
      object@metadata[["PCAlist"]][["norm"]] <-
      PCA(t(pseudo), ncp = ncomp, graph = FALSE)
    return(object)
    
  }
)

#' @rdname runOmicsPCA
#' @aliases runOmicsPCA
#' @title runOmicsPCA
#' @param SE.name the name of the data the normalization have to be applied to.
#' @exportMethod runOmicsPCA
setMethod(
  f          = "runOmicsPCA",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        ncomp = 5,
                        raw = FALSE) {
    object[[SE.name]] <-  runOmicsPCA(object[[SE.name]],
                                      ncomp = ncomp,
                                      raw  = raw)
    
    return(object)
    
  }
)

## ---- checkExpDesignCompleteness ----

#' @name checkExpDesignCompleteness
#' @rdname Rflomics-accessors
#' @description
#' \itemize{
#'    \item checkExpDesignCompleteness: return a string with message.
#'    This method checks some experimental design characteristics.
#'    A complete design (all combinations of factor modalities with at 
#'    least 2 replicates for each have to be present) with
#'    at least one biological and one batch factors are required to use the 
#'    RFLOMICS workflow.}
#' @exportMethod checkExpDesignCompleteness
#' @param sampleList list of samples to check.
#' @exportMethod checkExpDesignCompleteness
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

#' @rdname Rflomics-accessors
#' @param omicName the name of the data the normalization have to be applied to. 
#' @exportMethod checkExpDesignCompleteness
setMethod(f         = "checkExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){
            
            if(is.null(omicName)) stop("Arg omicName missed")
            
            SEObject <- getRflomicsSE(object, omicName)
            
            checkExpDesignCompleteness(SEObject, sampleList = sampleList)
          })

##==== ACCESSORS ====

###==== getProcessedData ====

#' @rdname Rflomics-accessors
#' @name getProcessedData
#' @description 
#' Getters and Setters for data processing analysis:
#' \itemize{
#'    \item getProcessedData: return a processed data matrix  
#'    (filtering, normalization and/or transformation)}
#' @exportMethod getProcessedData
#' @examples
#' # See runDataProcessing for an example that includes getProcessedData
setMethod(f          = "getProcessedData",
          signature  = "RflomicsSE",
          definition = function(object)
          {
            object2 <- .checkTransNorm(object, raw = FALSE)
            pseudo <- assay(object2)
            
            return(data.frame(pseudo))
          })

#' @rdname Rflomics-accessors
#' @name getProcessedData
#' @exportMethod getProcessedData
setMethod(f          = "getProcessedData",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            
            if (!SE.name %in% getDatasetNames(object)){
              stop("SE name must be part of this list of names: ",
                   getDatasetNames(object))
            }
            
            pseudo <- getProcessedData(object[[SE.name]])
            
            return(pseudo)
          })

###==== getTransSettings ====

# Get transformation parameters
#' @rdname Rflomics-accessors
#' @name getTransSettings
#' @description 
#' \itemize{
#'    \item getTransSettings: return a list of transformation settings 
#'    of a given omic dataset}
#' @exportMethod getTransSettings
#' @examples
#' # See runDataProcessing for an example that includes getTransSettings
setMethod(f          = "getTransSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(object@metadata$DataProcessing$Transformation$setting)   
          })

#' @rdname Rflomics-accessors
#' @name getTransSettings
#' @exportMethod getTransSettings
setMethod(f          = "getTransSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getTransSettings(object = object[[SE.name]])
          })

###==== getFilterSettings ====

# Get filtering parameters
#' @rdname Rflomics-accessors
#' @name getFilterSettings
#' @description
#' \itemize{
#'    \item getFilterSettings: return a list the filtering settings of a given 
#'    omic dataset}
#' @exportMethod getFilterSettings
#' @examples
#' # See runDataProcessing for an example that includes getFilterSettings
setMethod(f          = "getFilterSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(object@metadata$DataProcessing$Filtering$setting)   
          })

#' @rdname Rflomics-accessors
#' @exportMethod getFilterSettings
setMethod(f          = "getFilterSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getFilterSettings(object = object[[SE.name]])
          })

###==== getFilteredFeatures ====

# Get filtred features
#' @rdname Rflomics-accessors
#' @name getFilteredFeatures
#' @description
#' \itemize{
#'    \item getFilteredFeatures: return a vector of filtered features of a given 
#'    omic dataset}
#' @exportMethod getFilteredFeatures 
#' @examples
#' # See runDataProcessing for an example that includes getFilteredFeatures
setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(object@metadata$DataProcessing$Filtering$results$filteredFeatures)   
          })

#' @rdname Rflomics-accessors
#' @exportMethod getFilteredFeatures 
setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getFilteredFeatures(object = object[[SE.name]])
          })

###==== getCoeffNorm ====

#' @rdname Rflomics-accessors
#' @name getCoeffNorm
#' @description
#' \itemize{
#'    \item getCoeffNorm: return a named vector with normalization coefficients 
#'    of a given omic dataset}
#' @exportMethod getCoeffNorm
#' @examples
#' # See runDataProcessing for an example that includes getCoeffNorm
setMethod(f          = "getCoeffNorm",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(metadata(object)[["DataProcessing"]][["Normalization"]][["results"]][["coefNorm"]])
          })

#' @rdname Rflomics-accessors
#' @exportMethod getCoeffNorm 

setMethod(f          = "getCoeffNorm",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getCoeffNorm(object = object[[SE.name]])
          })

###==== getNormSettings ====
# Get normalizationparameters

#' @rdname Rflomics-accessors
#' @name getNormSettings
#' @description
#' \itemize{
#'    \item getNormSettings: return a list of normalization settings 
#'    of a given omic dataset}
#' @exportMethod getNormSettings
#' @examples
#' # See runDataProcessing for an example that includes getNormSettings
setMethod(f          = "getNormSettings",
          signature  = "RflomicsSE",
          definition = function(object){
            return(object@metadata$DataProcessing$Normalization$setting)
          })

#' @rdname Rflomics-accessors
#' @exportMethod getNormSettings
setMethod(f          = "getNormSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getNormSettings(object = object[[SE.name]])
          })

##==== GRAPHICAL METHOD ====

###==== plotLibrarySize ====

#' @rdname Rflomics-plots
#' @name plotLibrarySize
#' @description
#' Graphical functions for data processing analysis:
#' \itemize{
#'    \item plotLibrarySize: return barplot of library size by sample.}
#' @param raw a boolean 
#' @exportMethod plotLibrarySize
#' @importFrom ggplot2 ggplot aes ggtitle element_text 
#' @importFrom ggplot2 theme labs ylab geom_bar
#' @importFrom dplyr full_join arrange
#' @examples
#' # See runDataProcessing for an example that includes plotLibrarySize
setMethod(f          = "plotLibrarySize",
          signature  = "RflomicsSE",
          definition = function(object, raw = FALSE){
            
            if (getOmicsTypes(object) != "RNAseq") {
              stop("Data are not RNAseq!")
            }
            
            abundances <- assay(object)
            Groups     <- getDesignMat(object)
            samples    <- colnames(abundances)
            
            if (raw) {
              pseudo <- assay(object) %>% colSums(., na.rm = TRUE)
              title  <- "Raw data"
            }else {
              # RNAseq, not expected to find any 
              # transformation method.
              if (!.isNorm(object)) {
                pseudo <- assay(.applyNorm(object))
              }
              
              pseudo <- pseudo %>% colSums(., na.rm = TRUE)
              title <- paste0("Filtered and normalized (",
                              getNormSettings(object)$method, 
                              ") data")
            }
            
            libSizeNorm <-  full_join(Groups, 
                                      data.frame("value" = pseudo , 
                                                 "samples" = names(pseudo)), 
                                      by = "samples") %>%
              arrange(groups)
            
            libSizeNorm$samples <- factor(libSizeNorm$samples, 
                                          levels = levels(Groups$samples))
            
            p <- ggplot(libSizeNorm,  aes(x = samples, y = value, fill = groups)) + 
              geom_bar(stat = "identity" ) + 
              ylab(ylab) + 
              theme(axis.text.x =  element_text(angle = 45, hjust = 1), legend.position  = "none") + 
              labs(x = "", y = "Total read count per sample") + 
              ggtitle(title)
            
            return(p)
            
          })

#' @rdname Rflomics-plots
#' @name plotLibrarySize
#' @exportMethod plotLibrarySize
setMethod(f          = "plotLibrarySize",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, 
                                raw = FALSE){
            
            if (getOmicsTypes(object[[SE.name]]) == "RNAseq") {
              return(plotLibrarySize(object[[SE.name]], 
                                     raw = raw))
            }else{
              stop("This function only applies to RNAseq data")
            }
          })

###==== plotDataDistribution ====

#' @rdname Rflomics-plots
#' @name plotDataDistribution
#' @description
#' \itemize{
#'    \item plotDataDistribution: return boxplot or density plot of expression
#'    or abundance distribution.}
#' @param plot plot type ("boxplot" or "density")
#' @param raw boolean. Plot the raw data or the transformed ones (raw = FALSE)
#' @exportMethod plotDataDistribution
#' @importFrom ggplot2 ggplot geom_boxplot geom_density aes
#' @importFrom ggplot2 xlab theme element_text ylab margin ggtitle
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join arrange
#' @examples
#' # See runDataProcessing for an example that includes plotDataDistribution
setMethod(
  f = "plotDataDistribution",
  signature = "RflomicsSE",
  definition = function(object, plot = "boxplot", raw = FALSE) {
    
    object2 <- .checkTransNorm(object, raw = raw)
    pseudo <- assay(object2)
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
      
      if (.isTransformed(object2)) {
        title <- paste0("Transformed (", 
                        getTransSettings(object2)$method, ") ", 
                        title)
      }
      
      if (.isNorm(object2)) {
        title <- paste0(title, " - normalization: ", 
                        getNormSettings(object2)$method)
      }
      
      if (omicsType == "RNAseq") {
        x_lab <- paste0("log2(", omicsType, " data)")
      }
    }
    
    pseudo.gg <- pseudo %>% reshape2::melt()
    colnames(pseudo.gg) <- c("features", "samples", "value")
    
    pseudo.gg <- pseudo.gg %>%
      full_join(Groups, by = "samples") %>%
      arrange(groups)
    
    pseudo.gg$samples <- factor(pseudo.gg$samples, 
                                levels = unique(pseudo.gg$samples))
    
    switch(plot,
           "density" = {
             p <- ggplot(pseudo.gg) +
               geom_density( aes(x = value, group = samples, color = groups), trim = FALSE) +
               xlab(x_lab) +
               theme(legend.position = "none") +
               ggtitle(title)
           },
           "boxplot" = {
             p <-  ggplot(pseudo.gg,  aes(x = samples, y = value, label = features)) +
               geom_boxplot( aes(fill = groups), outlier.size = 0.3) +
               theme(axis.text.x =  element_text(angle = 45, hjust = 1), legend.position = "none",
                     plot.margin= margin(0.5,0.5,0.5,1,"cm")) +
               xlab("") +
               ylab(x_lab) +
               ggtitle(title) #+
             # geom_point(alpha = 1/100,size=0)
           }
    )
    
    return(p)
  }
)

#' @rdname Rflomics-plots
#' @name plotDataDistribution
#' @exportMethod plotDataDistribution
setMethod(
  f = "plotDataDistribution",
  signature = "RflomicsMAE",
  definition = function(object, SE.name, 
                        plot = "boxplot", 
                        raw = FALSE) {
    plotDataDistribution(
      object = object[[SE.name]],
      plot = plot,
      raw = raw
    )
  }
)

### ---- plotOmicsPCA ----
#' @name plotOmicsPCA
#' @description 
#' \itemize{
#'    \item plotOmicsPCA: 
#' This function plot the factorial map from a PCA object stored
#' in a \link{RflomicsSE-class} object. By default, samples are
#' colored by groups (all combinations of level's factor)}
#' @param raw This argument indicates whether the scaled PCA has to be 
#' performed on raw [\sQuote{raw}] or normalized [\sQuote{norm}] data.
#' @param axes A vector giving the two axis that have to be drawn for the 
#' factorial map
#' @param groupColor All combination of level's factor
#' @importFrom dplyr mutate full_join select right_join
#' @importFrom FactoMineR coord.ellipse
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text aes xlab ylab 
#' geom_hline geom_vline geom_vline element_text ggtitle geom_polygon
#' @exportMethod plotOmicsPCA
#' @rdname Rflomics-plots
#' @examples
#' # See runDataProcessing for an example that includes plotOmicsPCA
setMethod(f = "plotOmicsPCA",
          signature = "RflomicsSE",
          definition = function(object,
                                raw = c("raw", "norm"),
                                axes = c(1, 2),
                                groupColor = "groups") {
            if (length(axes) != 2) {
              stop("PCA axes must be a vector of length 2")
            }
            
            ExpDesign <- getDesignMat(object)
            
            PC1 <- paste("Dim.", axes[1], sep = "")
            PC2 <- paste("Dim.", axes[2], sep = "")
            
            if (PC1 == PC2) PC2 <- PC1 + 1
            
            score <- metadata(object)$PCAlist[[raw]]$ind$coord[, axes] %>% 
              as.data.frame() %>%
              mutate(samples = row.names(.)) %>% 
              right_join(., ExpDesign, by = "samples")
            
            var1 <- round(metadata(object)$PCAlist[[raw]]$eig[axes, 2][1], digits = 3)
            var2 <- round(metadata(object)$PCAlist[[raw]]$eig[axes, 2][2], digits = 3)
            
            omicsType <- getOmicsTypes(object)
            
            switch(raw,
                   "raw"  = {
                     title <- paste0("Raw ", omicsType, " data")
                   },
                   "norm" = {
                     title <- switch (
                       omicsType,
                       "RNAseq" = {
                         paste0(
                           "Filtred and normalized ",
                           omicsType,
                           " data (",
                           getNormSettings(object)$method,
                           ")"
                         )
                       },
                       "proteomics" = {
                         paste0(
                           "Transformed and normalized ",
                           omicsType,
                           " data (",
                           getTransSettings(object)$method,
                           " - norm: ",
                           getNormSettings(object)$method,
                           ")"
                         )
                       },
                       "metabolomics" = {
                         paste0(
                           "Transformed and normalized ",
                           omicsType,
                           " data (",
                           getTransSettings(object)$method,
                           " - norm: ",
                           getNormSettings(object)$method,
                           ")"
                         )
                       }
                     )
                   })
            
            
            
            p <- ggplot(score,
                        aes_string(x = PC1, y = PC2, color = groupColor))  +
              geom_point(size = 2) +
              geom_text(
                aes(label = samples),
                size = 2,
                vjust = "inward",
                hjust = "inward"
              ) +
              xlab(paste(PC1, " (", var1, "%)", sep =
                           "")) +
              ylab(paste(PC2, " (", var2, "%)", sep =
                           "")) +
              geom_hline(yintercept = 0,
                         linetype = "dashed",
                         color = "red") +
              geom_vline(xintercept = 0,
                         linetype = "dashed",
                         color = "red") +
              theme(
                strip.text.x =  element_text(size = 8, 
                                             face = "bold.italic"),
                strip.text.y =  element_text(size = 8, 
                                             face = "bold.italic")
              ) +
              ggtitle(title)
            
            # ellipse corr
            aa <- select(score, all_of(groupColor), PC1, PC2)
            bb <- coord.ellipse(aa, bary = TRUE)
            p <- p + geom_polygon(
              data = bb$res,
              aes_string(x = PC1, y = PC2, fill = groupColor),
              show.legend = FALSE,
              alpha = 0.1
            )
            
            return(p)
            
          })

#' @rdname Rflomics-plots
#' @name plotOmicsPCA
#' @exportMethod plotOmicsPCA
setMethod(
  f          = "plotOmicsPCA",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        raw,
                        axes = c(1, 2),
                        groupColor = "groups") {
    plotOmicsPCA(object[[SE.name]], raw, axes, groupColor)
    
  }
)

### ---- plotExpDesignCompleteness ----
#' @name plotExpDesignCompleteness
#' @description 
#' \itemize{
#'    \item plotExpDesignCompleteness: 
#' This method checks that experimental design constraints are satisfied and 
#' plot a summary of the design.
#' A complete design (all combinations of factor modalities with at least 2 
#' replicates for each have to be present) 
#' with at least one biological and one batch factors are required to use the 
#' RFLOMICS workflow.}
#' @param sampleList list of samples to check.
#' @exportMethod plotExpDesignCompleteness
#' @rdname Rflomics-plots
#' @examples
#' # See runDataProcessing for an example that includes plotExpDesignCompleteness
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

#' @param omicName a character string with the name of the dataset
#' @exportMethod plotExpDesignCompleteness
#' @rdname Rflomics-plots
setMethod(f         = "plotExpDesignCompleteness",
          signature = "RflomicsMAE",
          definition <- function(object, omicName, sampleList=NULL){
            
            SEObject <- getRflomicsSE(object, omicName)
            
            plotExpDesignCompleteness(SEObject, sampleList = sampleList)
          })
