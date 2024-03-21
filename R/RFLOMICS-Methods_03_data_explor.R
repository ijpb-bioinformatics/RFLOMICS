

################################# EXPLORATION OF BIOLOGICAL AND TECHNICAL VARIABILITY ##################################


##### Statistical METHODS for exploring biological and technical variability

##### Graphical METHODS for exploring biological and technical variability


#' plotLibrarySize
#' 
#' @param object An object of class \link{RflomicsSE}. 
#' @param raw a boolean 
#' @return plot
#' @exportMethod plotLibrarySize
#' @importFrom ggplot2 ggplot aes ggtitle labs element_text theme 
#' ylab geom_bar
#' @importFrom dplyr full_join arrange
#' @rdname plotLibrarySize
#' @aliases plotLibrarySize
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"),
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")), 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest", "protetest"),
#'                          omicsTypes  = c("RNAseq","proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest", "protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE |>  runTransformData(SE.name = "protetest", transformMethod = "log2") |>
#'  filterLowAbundance(SE.name = "RNAtest")                           |>    
#'  runNormalization(SE.name = "RNAtest",   normMethod = "TMM")       |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")
#'  
#'  # plotLibrarySize
#'  plotLibrarySize(MAE@ExperimentList[[1]], raw=FALSE)
#'  plotLibrarySize(MAE@ExperimentList[[1]], raw=TRUE)
#'  
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

#' @rdname plotLibrarySize
#' @aliases plotLibrarySize
#' @title plotLibrarySize
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod plotLibrarySize
setMethod(f          = "plotLibrarySize",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, raw = FALSE){
              
              if (getOmicsTypes(object[[SE.name]]) == "RNAseq") {
                  return(plotLibrarySize(object[[SE.name]], 
                                         raw = raw))
              }else{
                  stop("This function only applies to RNAseq data")
              }
              
          })


#' @title plotDataDistribution
#'
#' @param object An object of class \link{RflomicsSE}
#' @param plot plot type ("boxplot" or "density")
#' @param raw boolean. Plot the raw data or the transformed ones (raw = FALSE)
#' @export
#' @exportMethod plotDataDistribution
#' @importFrom ggplot2 ggplot aes geom_density xlab theme ggtitle geom_boxplot 
#' element_text ylab margin
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join arrange
#' @rdname plotDataDistribution
#' @aliases plotDataDistribution
#' @return A ggplot object
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")), 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest", "protetest"),
#'                          omicsTypes  = c("RNAseq","proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest", "protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE,
#'                                      formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  filterLowAbundance(SE.name = "RNAtest")                           |>    
#'  runNormalization(SE.name = "RNAtest",   normMethod = "TMM")       |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")
#'  
#'  # plotLibrarySize
#'  plotDataDistribution(MAE@ExperimentList[[1]], raw=FALSE)
#'  plotDataDistribution(MAE@ExperimentList[[1]], raw=TRUE)
#'  

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

#' @rdname plotDataDistribution
#' @aliases plotDataDistribution
#' @title plotDataDistribution
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod plotDataDistribution
setMethod(
    f = "plotDataDistribution",
    signature = "RflomicsMAE",
    definition = function(object, SE.name, plot = "boxplot", raw = FALSE) {
        plotDataDistribution(
            object = object[[SE.name]],
            plot = plot,
            raw = raw
        )
    }
)

########################################## TRANSFORM DATA #################

#### METHOD to transform data

#' @title runTransformData
#'
#' @param object An object of class \link{RflomicsSE}
#' @param transformMethod The transformation to store in the metadata or to 
#' store and apply if modifyAssay is TRUE.
#' @param modifyAssay Boolean. Do the transformation need to be applied on 
#' the data? The raw data will be replaced by the transformed ones.
#'
#' @return An object of class \link{RflomicsSE}
#'
#' @exportMethod runTransformData
#' @rdname runTransformData
#' @aliases runTransformData
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"),omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2")
#'  
#'  # Get the transformation settings
#'  getTransSettings(MAE,SE.name = "protetest")
#'  
setMethod(f          = "runTransformData",
          signature  = "RflomicsSE",
          definition = function(object, 
                                transformMethod = NULL, 
                                modifyAssay = FALSE){
              
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
              
              object <- setTrans(object, transformMethod)
              # object@metadata[["DataProcessing"]][["Transformation"]]
              # [["setting"]][["method"]] <- transformMethod
              
              if (modifyAssay) {
                  object <-  .applyTransformation(object)
              }
              
              return(object)
              
          })

#' @rdname runTransformData
#' @title runTransformData
#' @aliases runTransformData
#' @param SE.name the name of the data the normalization have to be applied to. 
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


# ---- Get transformation parameters ----

#' @title get transformation  parameters
#'
#' @param object of class RflomicsSE
#' @exportMethod getTransSettings
#' @rdname getTransSettings
#' @aliases getTransSettings
#' @return List of transformation parametres.


setMethod(f          = "getTransSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
              return(object@metadata$DataProcessing$Transformation$setting)   
          })

#' @rdname getTransSettings
#' @title getTransSettings
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @aliases getTransSettings
#' @exportMethod getTransSettings

setMethod(f          = "getTransSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getTransSettings(object = object[[SE.name]])
          })


#' @title set transformation  method
#' @param object of class RflomicsSE
#' @param method the transformation method to set
#' @rdname setTrans
#' @aliases setTrans
#' @return a RflomicsSE or RflomicsMAE object.
#' @exportMethod setTrans
setMethod(f="setTrans", 
          signature=c("RflomicsSE"),
          definition=function(object, method = "none"){ 
              metadata(object)[["DataProcessing"]][["Transformation"]][["setting"]][["method"]] <- method
              return(object) 
          }) 

#' @title set transformation  method
#' @rdname setTrans
#' @title setTrans
#' @aliases setTrans
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @exportMethod setTrans

setMethod(f          = "setTrans",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, method = "none"){
              setTrans(object = object[[SE.name]], method=method ) 
          })

########################################## FILTER DATA #################

#### METHOD to filter data

# Cette method est propre au RNASEQ => Est-ce que c'est vraiment ce que l'on souhaite ?
# Plutot qu'une fonction interface pour tous les omics ?
# Pourquoi ne pas avoir utilisée directement la fonction de edgeR ?

#' @title filterLowAbundance
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
#' @exportMethod filterLowAbundance
#' @importFrom edgeR DGEList filterByExpr cpm
#' @seealso edgeR::filterByExpr
#' @rdname filterLowAbundance
#' @aliases filterLowAbundance
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest"), omicsTypes  = c("RNAseq"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' MAE |>  filterLowAbundance(SE.name = "RNAtest") 
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


#' @rdname filterLowAbundance
#' @title filterLowAbundance
#' @aliases filterLowAbundance
#' @param SE.name the name of the data the normalization have to be applied to. 
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



# ---- Get filtering parameters ----

#' @title get Filter setting parameters
#'
#' @param object of class RflomicsSE
#' @return List of differential analysis setting parametres.
#' @exportMethod getFilterSettings
#' @aliases getFilterSettings
#' @rdname getFilterSettings
#'

setMethod(f          = "getFilterSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
              return(object@metadata$DataProcessing$Filtering$setting)   
          })

#' @rdname getFilterSettings
#' @title getFilterSettings
#' @param SE.name the name of the data to fetch in the object if the object is 
#' a RflomicsMAE
#' @exportMethod getFilterSettings
#' @aliases getFilterSettings

setMethod(f          = "getFilterSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getFilterSettings(object = object[[SE.name]])
          })



# ---- Get filtred features ----

#' @title get Filtered Features 
#' @param object of class RflomicsSE
#' @return List of differential analysis setting parametres.
#' @exportMethod getFilteredFeatures 
#' @aliases getFilteredFeatures
#' @rdname getFilteredFeatures 
#'

setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsSE",
          
          definition = function(object){
              return(object@metadata$DataProcessing$Filtering$results$filteredFeatures)   
          })

#' @rdname getFilteredFeatures 
#' @title getFilteredFeatures 
#' @param SE.name the name of the data to fetch in the object if the 
#' object is a RflomicsMAE
#' @aliases getFilteredFeatures
#' @exportMethod getFilteredFeatures 

setMethod(f          = "getFilteredFeatures",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getFilteredFeatures(object = object[[SE.name]])
          })



######### NORMALIZATION #################

#### METHOD to normalize data

# Function non generique pour les autres data

#' @title runNormalization
#' @description This function applied a normalization method on an omic data sets stored in an object of
#' class \link{RflomicsSE} or \link{RFLomicsMAE}. In this case, the normalization method is applied to the data 
#' set and 
#' 
#' \itemize{
#' \item{For RNAseq data:}{the TMM function of edgeR is proposed by default, see the ref}
#' \item{For Proteomic data:}{}
#' \item{For Metabolomic data:}{}
#' }
#' @param object An object of class \link{RflomicsSE}
#' @param normMethod Normalization method
#' @param modifyAssay Does the normalization have to be applied or just stored for later? Recommended it stays FALSE.
#' @return An object of class \link{RflomicsSE}
#' The applied normalization method and computed scaling factors (by samples) are stored as a named list
#' ("normalization") of two elements (respectively "method" and "coefNorm") in the metadata slot of a
#' given data set, stored itself in the ExperimentList slot of a \link{RflomicsSE} object.
#' @importFrom stats median
#' @exportMethod runNormalization
#' @rdname runNormalization
#' @aliases runNormalization
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. 
#' DiCoExpress: a tool to process multifactorial
#' RNAseq experiments from quality controls to co-expression analysis through 
#' differential analysis based on contrasts inside GLM models. 
#' Plant Methods 16, 68 (2020).
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest"), omicsTypes  = c("RNAseq"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' MAE |>  filterLowAbundance(SE.name = "RNAtest") |>
#' runNormalization(SE.name = "RNAtest", normMethod = "TMM") 
#'
#' # Get the normalization method and the scaling factors
#' getCoeffNorm(MAE@ExperimentList[[1]])
#' 
#' # Get the normalization settings
#' getNormSettings(MAE@ExperimentList[[1]])
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

#' @rdname runNormalization
#' @title runNormalization
#' @aliases runNormalization
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @exportMethod runNormalization
setMethod(f          = "runNormalization",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, normMethod = NULL, modifyAssay = FALSE){
              
              object[[SE.name]] <-  runNormalization(object       = object[[SE.name]],
                                                     normMethod   = normMethod,
                                                     modifyAssay = modifyAssay)
              
              return(object)
              
          })


#' @title get Coeff Norm
#' @param object of class RflomicsSE
#' @return Coefficient of normalization
#' @exportMethod getCoeffNorm
#' @aliases getCoeffNorm
#' @rdname getCoeffNorm 
#'

setMethod(f          = "getCoeffNorm",
          signature  = "RflomicsSE",
          
          definition = function(object){
              return(metadata(object)[["DataProcessing"]][["Normalization"]][["results"]][["coefNorm"]])
          })

#' @rdname getCoeffNorm 
#' @aliases getCoeffNorm
#' @title getCoeffNorm 
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod getCoeffNorm 

setMethod(f          = "getCoeffNorm",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getCoeffNorm(object = object[[SE.name]])
          })

#' @title set Normalization coefficient
#' @param object of class RflomicsSE
#' @param coeff normalization coefficient
#' @rdname setCoeffNorm
#' @aliases setCoeffNorm
#' @exportMethod setCoeffNorm
#' @return object of class RflomicsSE
setMethod(f="setCoeffNorm", 
          signature=c("RflomicsSE"),
          definition=function(object, coeff = NULL){ 
              metadata(object)[["DataProcessing"]][["Normalization"]][["results"]][["coefNorm"]] <- coeff
              return(object) 
          }) 

#' @rdname setCoeffNorm
#' @title setCoeffNorm
#' @aliases setCoeffNorm
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod setCoeffNorm
setMethod(f          = "setCoeffNorm",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name = NULL, coeff = NULL){
              setCoeffNorm(object = object[[SE.name]], coeff = coeff) 
          })

# ---- Get normalizationparameters ----

#' @title Get normalization parameters
#' @param object of class RflomicsSE
#' @return List of differential analysis setting parametres.
#' @exportMethod getNormSettings
#' @aliases getNormSettings
#' @rdname getNormSettings
#'

setMethod(f          = "getNormSettings",
          signature  = "RflomicsSE",
          definition = function(object){
              return(object@metadata$DataProcessing$Normalization$setting)
          })

#' @rdname getNormSettings
#' @title getNormSettings
#' @aliases getNormSettings
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @exportMethod getNormSettings

setMethod(f          = "getNormSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
              getNormSettings(object = object[[SE.name]])
          })


#' @title set Normalization  method
#' @param object of class RflomicsSE
#' @param method the Normalization method to set
#' @rdname setNorm
#' @aliases setNorm
#' @exportMethod setNorm
#' @return object of class RflomicsSE
setMethod(f="setNorm", 
          signature=c("RflomicsSE"),
          definition=function(object, method = "none"){ 
              metadata(object)[["DataProcessing"]][["Normalization"]][["setting"]][["method"]] <- method
              return(object) 
          }) 

#' @rdname setNorm
#' @title setNorm
#' @aliases setNorm
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @exportMethod setNorm
#' @return object of class RflomicsMAE

setMethod(f          = "setNorm",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name = NULL, method = "none"){
              setNorm(object = object[[SE.name]], method=method ) 
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
#' @param lowCountFiltering_strategy strategy of RNAseq low count filtering. 
#' Mandatory for RNAseq data. Default value : "NbReplicates".
#' @param lowCountFiltering_CPM_Cutoff CPM cutoff for RNAseq low count filtering.
#'  Mandatory for RNAseq data. Default value : 1.
#' @param normMethod of normalisation. Mandatory for RNAseq data. 
#' Default value : RNAseq = TMM.
#' @param transformMethod method of transformation.
#' @return An object of class \link{RflomicsSE}
#' @aliases runDataProcessing
#' @exportMethod runDataProcessing
#' @seealso runSampleFiltering
#' @seealso  filterLowAbundance
#' @seealso runNormalization
#' @seealso runTransformData
#' @rdname runDataProcessing
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")), 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest.raw", "protetest.raw"),
#'                          omicsTypes  = c("RNAseq","proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest.raw", "protetest.raw")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' MAE <- runDataProcessing(MAE, SE.name = "RNAtest",  
#'                          lowCountFiltering_strategy = "NbReplicates", 
#'                          lowCountFiltering_CPM_Cutoff = 1, 
#'                          normMethod = "TMM", 
#'                          transformMethod = "log2") 
#' MAE <- runDataProcessing(MAE, SE.name = "protetest", 
#'                          normMethod = "median", 
#'                          transformMethod = "log2")


setMethod(f          = "runDataProcessing",
          signature  = "RflomicsSE",
          definition = function(object, samples=NULL, 
                                lowCountFiltering_strategy = "NbReplicates", 
                                lowCountFiltering_CPM_Cutoff = 1, 
                                normMethod = "none", transformMethod = "none")
          {
              
              # keep selected samples
              message("#    => select samples...")
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
                         message("#    => Low counts Filtering...")
                         if(is.null(lowCountFiltering_strategy)   || !lowCountFiltering_strategy %in% lowCountFiltering_strategy.sup) 
                             stop("the low count filtering strategy : ", lowCountFiltering_strategy, " isn't supported by rflomics package. Supported values : ",  paste(lowCountFiltering_strategy.sup, collapse = ", "))
                         
                         if(is.null(lowCountFiltering_CPM_Cutoff) || !is.numeric(lowCountFiltering_CPM_Cutoff)) 
                             stop(lowCountFiltering_CPM_Cutoff, " must be a integer value > 1")
                         
                         SE.processed <- filterLowAbundance(object = object, filterMethod= "CPM", filterStrategy = lowCountFiltering_strategy, cpmCutoff = lowCountFiltering_CPM_Cutoff)
                         
                         # Run Normalisation 
                         message("#    => Counts normalization...")
                         if(is.null(normMethod) || normMethod != "TMM"){
                             normMethod <- "TMM"
                             warning("only ", normalisation_method.count.sup, " method is supported for ", object@metadata$omicType, " normalisation.")
                         }
                         
                         SE.processed <- runNormalization(SE.processed, normMethod = normMethod)
                     },
                     {
                         message("#    => transformation data...")
                         if(is.null(transformMethod)) transformMethod <- "none"
                         if(! transformMethod %in% transformation_method.sup) 
                             stop("the transformation method ", transformMethod," is not support in rflomics package.Supported values : ", paste(transformation_method.sup, collapse = ", "))
                         SE.processed <- runTransformData(object, transformMethod = transformMethod)
                         
                         message("#    => Run normalization...")
                         if(is.null(normMethod)) normMethod <- "none"
                         if(! normMethod %in% normalisation_method.abundance.sup) 
                             stop("the normalisation method ", normMethod," is not support in rflomics package. Supported values : ", paste(normalisation_method.abundance.sup, collapse = ", "))
                         SE.processed <- runNormalization(SE.processed, normMethod = normMethod)
                     }
              )
              
              #### Run PCA for filtred & normalized data ####
              message("#    => Compute PCA ")
              
              SE.processed <- runOmicsPCA(SE.processed,ncomp = 5, raw = FALSE)  
              
              SE.processed@metadata$DataProcessing[["done"]] <- TRUE
              
              return(SE.processed)
          })


#' @rdname runDataProcessing
#' @title runDataProcessing
#' @param SE.name the name of the data the normalization have to be applied to. 
#' @aliases runDataProcessing
#' @exportMethod runDataProcessing
setMethod(f          = "runDataProcessing",
          signature  = "RflomicsMAE",
          definition = function(object, samples=NULL, lowCountFiltering_strategy = "NbReplicates", lowCountFiltering_CPM_Cutoff = 1, 
                                normMethod= "none", transformMethod = "none", SE.name){
              
              # if paste0(SE.name, ".raw") exist
              if (!paste0(SE.name, ".raw") %in% names(object)){
                  
                  if (!SE.name %in% names(object)) {
                      stop("no ", SE.name, " SE object.")
                  }
                  stop("raw SE must be tagged by .raw")
              }
              
              # if paste0(SE.name, ".raw") exist
              if(is.null(object[[paste0(SE.name, ".raw")]])) {
                  stop("raw SE must be tagged by .raw")
              }
              
              SE.raw       <- object[[paste0(SE.name, ".raw")]]
              SE.processed <- runDataProcessing(object = SE.raw,
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

# ------ filtering per sample -----

#' @title runSampleFiltering
#' @description This function applied sample filtering on an omic data sets stored in an object of
#' class \link{RflomicsSE}.
#' @param object An object of class \link{RflomicsSE}
#' @param samples samples to keep.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runSampleFiltering
#' @importFrom dplyr filter
#' @rdname runSampleFiltering
#' @aliases runSampleFiltering
#' @examples
#' #' Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0( datPath, "transcriptome_ecoseed.txt")), 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("RNAtest.raw", "protetest.raw"),
#'                          omicsTypes  = c("RNAseq","proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("RNAtest.raw", "protetest.raw")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#' SE.RNA <- runSampleFiltering(MAE[[1]],samples=colnames(MAE[[1]])[1:6])
#' SE.prot <- runSampleFiltering(MAE[[2]],samples=colnames(MAE[[2]])[1:6])

setMethod(f          = "runSampleFiltering",
          signature  = "RflomicsSE",
          definition = function(object, samples=NULL) {
              
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
              
              return(SE.new)
          })


