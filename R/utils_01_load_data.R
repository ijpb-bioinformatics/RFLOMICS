### ============================================================================
### Create Rflomics object / RflomicsSE and RflomicsMAE
### ----------------------------------------------------------------------------

###### RflomicsMAE CLASS Constructor for managing omics DATA and RESULTS

#' @title createRflomicsMAE is creator for the class \link{createRflomicsMAE-class}
#' @description This function initializes an object of class \link{createRflomicsMAE-class}
#' from a list of omics data.
#' @param projectName Project name
#' @param omicsData list of omics dataset.
#' @param omicsNames vector of dataset names
#' @param omicsTypes vector of dataset types
#' @param ExpDesign a data.frame. Row names give the name of each sample which has been to be construct
#' @param factorRef data.frame describing experimental factors.
#' \itemize{
#' \item{factorName:}{factor names}
#' \item{factorRef:}{factor references}
#' \item{factorType:}{factor type : "Bio", "batch", "Meta"}
#' \item{factorLevels:}{levels of each factor with "," separation.}
#' }
#' @return An object of class \link{createRflomicsMAE-class}
#' @examples
#' 
#' factorRef <- data.frame(factorName  = c("Repeat", "temperature" , "imbibition"),
#' factorRef   = c("rep1",   "Low",          "DS"),
#' factorType  = c("batch",  "Bio",          "Bio"),
#' factorLevels= c("rep1,rep2,rep3", "Low,Medium,Elevated", "DS,EI,LI"))
#' 
#' omicsData <- list(
#'   RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
#'   RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")),
#'   RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")))
#' 
#' MAE <- RFLOMICS::createRflomicsMAE(projectName = "Tests",
#'                                                omicsData   = omicsData,
#'                                                omicsNames  = c("RNAtest", "metatest", "protetest"),
#'                                                omicsTypes  = c("RNAseq","metabolomics","proteomics"),
#'                                                ExpDesign   = ExpDesign,
#'                                                factorRef   = factorRef)
#' 
#'
#' @name createRflomicsMAE
#' @rdname createRflomicsMAE
#' @export
#'
createRflomicsMAE <- function(projectName=NULL, omicsData=NULL, omicsNames=NULL, omicsTypes=NULL, ExpDesign=NULL, factorRef=NULL){
  
  #check arg
  ##projectName
  if(is.null(projectName)) stop("projectName is mandatory.")
  projectName <- stringr::str_replace_all(string = projectName, pattern = "[# /-]", replacement = "")
  
  ## omicsNames
  if(is.null(omicsNames)) stop("list of omicsNames is mandatory.")
  nb_omicsData <- length(omicsNames)
  omicsNames <- stringr::str_replace_all(string = omicsNames, pattern = "[# /-]", replacement = "")
  if (isTRUE(any(duplicated(omicsNames)))) stop("presence of duplicates in the omicsNames")
  
  ## omicsData
  if(!is.list(omicsData) || length(omicsData) == 0) stop("the omicsData list is mandatory.")
  if(nb_omicsData != length(omicsData)) stop("the number of omicsData matrix must match the number of omicsNames.")
  names(omicsData) <- omicsNames
  
  ## omicsTypes
  if(is.null(omicsTypes)) stop("the list of omicsTypes is mandatory.")
  if(nb_omicsData != length(omicsTypes)) stop("the number of omicsData matrix must match the number of omicsTypes")
  if(isTRUE(any(!unique(omicsTypes) %in% c("RNAseq","metabolomics","proteomics")))) stop("omicsTypes must be part of RNAseq, metabolomics, or proteomics.")
  names(omicsTypes) <- omicsNames
  
  ## ExpDesign
  if (is.null(ExpDesign)) stop("the ExpDesign is mandatory.")
  if (nrow(ExpDesign) == 0 || ncol(ExpDesign) == 0) stop("the ExpDesign is mandatory.")
  designRownames <- stringr::str_replace_all(string = rownames(ExpDesign), pattern = "[*# -/]", replacement = "")
  if (isTRUE(any(duplicated(designRownames)))) stop("presence of duplicates in the ExpDesign colnames")
  rownames(ExpDesign) <- designRownames
  
  ## factorRef
  if (is.null(factorRef)) stop("data.frame factorRef is mandatory.")
  if (is.null(factorRef$factorName)) stop("factorRef$factorName is mandatory")
  if (any(!factorRef$factorName %in% colnames(ExpDesign))) stop("factorRef$factorName don't match ExpDesign colnames")
  
  if (is.null(factorRef$factorType)) stop("factorRef$factorType is mandatory.")
  if (any(!unique(factorRef$factorType) %in% c("batch", "Bio", "Meta"))) stop("factorRef$factorType must be part of batch, Bio or Meta")
  
  factorBio   <- dplyr::filter(factorRef, factorType == "Bio")$factorName
  factorBatch <- dplyr::filter(factorRef, factorType == "batch")$factorName
  
  ## set ref and levels to ExpDesign
  for (i in 1:nrow(factorRef)){
    
    # set ref 
    if (!is.null(factorRef$factorRef)){
      
      if(!factorRef[i,]$factorRef %in% ExpDesign[[factorRef[i,]$factorName]]) stop(paste0("The factor ref : ", factorRef[i,]$factorRef, " don't exist"))
      ref <- factorRef[i,]$factorRef
    }
    else{
      ref <- sort(ExpDesign[[factorRef[i,]$factorName]])[1]
    }
    ExpDesign <- ExpDesign[order(row.names(ExpDesign)), ]
    ExpDesign[[factorRef[i,]$factorName]] <- relevel(as.factor(ExpDesign[[factorRef[i,]$factorName]]), ref=ref)
    
    # set level
    if (!is.null(factorRef$factorLevels)){
      
      levels <- stringr::str_split(factorRef[i,]$factorLevels, ",") |> unlist() %>% stringr::str_remove(" ")
      if(any(!levels %in% ExpDesign[[factorRef[i,]$factorName]])) stop(paste0("The factor levels : ", factorRef[i,]$factorLevels, " don't exist"))
      
      ExpDesign[[factorRef[i,]$factorName]] <- factor(ExpDesign[[factorRef[i,]$factorName]], levels = levels)
    }
  }
  
  ## consctuct ExpDesign object
  refList  <- factorRef$factorRef;  names(refList)  <- factorRef$factorName
  typeList <- factorRef$factorType; names(typeList) <- factorRef$factorName
  #Design   <- ExpDesign.constructor(ExpDesign = ExpDesign, refList = refList, typeList = typeList)
  
  # Create the List.Factors list with the choosen level of reference for each factor
  names(typeList) <- names(ExpDesign)
  
  
  Design <- list(Factors.Type  = typeList, 
                 Model.formula = vector(), 
                 Contrasts.Sel = data.frame())
  
  
  #
  ExpDesign   <- dplyr::mutate(ExpDesign, samples=row.names(ExpDesign)) |>
    tidyr::unite("groups", all_of(factorBio), sep = "_", remove = FALSE)
  
  order_levels      <- with(ExpDesign, do.call(order, ExpDesign[c(factorBio, factorBatch)]))
  ExpDesign$samples <- factor(ExpDesign$samples, levels = unique(ExpDesign$samples[order_levels]))
  ExpDesign$groups  <- factor(ExpDesign$groups,  levels = unique(ExpDesign$groups[order_levels]))
  
  ## create SE object of each dataset
  SummarizedExperimentList <- list()
  listmap  <- list()
  omicList <- list()
  k <- 0
  
  for(data in omicsNames){
    
    k <- k+1
    
    omicType <- omicsTypes[data]
    
    RflomicsSE <- createRflomicsSE(omicsData[[data]], omicType, ExpDesign, typeList)
    
    #### run PCA for raw count
    # SummarizedExperimentList[[data]] <- RunPCA(SE, raw = TRUE)
    SummarizedExperimentList[[data]] <- RflomicsSE
    
    
    # metadata for sampleMap for RflomicsMAE
    listmap[[data]] <- data.frame(primary = as.vector(SummarizedExperimentList[[data]]@colData$samples),
                                  colname = as.vector(SummarizedExperimentList[[data]]@colData$samples),
                                  stringsAsFactors = FALSE)
    
    # 
    colnames <- c(names(omicList[[omicType]]), k)
    omicList[[omicType]] <- c(omicList[[omicType]] ,data)
    names(omicList[[omicType]]) <- colnames
    
  }
  
  prepMAE <- MultiAssayExperiment::prepMultiAssay( ExperimentList = SummarizedExperimentList,
                                                   sampleMap      = MultiAssayExperiment::listToMap(listmap),
                                                   colData        = ExpDesign, outFile = stdout())
  
  
  MAE <- MultiAssayExperiment::MultiAssayExperiment(experiments = prepMAE$experiments,
                                                    colData     = prepMAE$colData,
                                                    sampleMap   = prepMAE$sampleMap,
                                                    metadata    = list(omicList = omicList, projectName = projectName, design = Design)) 
  
  rflomicsMAE <- new("RflomicsMAE")
  for(slot in slotNames(MAE)) {
    slot(rflomicsMAE, slot) <- slot(MAE, slot)
  }
  
  # tag as raw data (le temps de trouver une solution pour ne pas faire co-exister les raw et les process)
  names(rflomicsMAE) <- paste(names(rflomicsMAE), "raw", sep = ".")
  return(rflomicsMAE)
}



###### createRflomicsSE CLASS Constructor for managing omics DATA and RESULTS

#' @title createRflomicsSE is creator for the class \link{createRflomicsSE-class}
#' @description This function initializes an object of class \link{createRflomicsSE-class}
#' from a list of omics data.
#' @param omicsData omics dataset.
#' @return An object of class \link{createRflomicsSE-class}
#' @name createRflomicsSE
#' @rdname createRflomicsSE
#' @export
#'
createRflomicsSE <- function(omicData, omicType, ExpDesign, design){
  
  factorBio   <- names(design[design == "Bio"])
  factorBatch <- names(design[design == "batch"])
  
  # check overlap between design and data
  sample.intersect <- intersect(row.names(ExpDesign), colnames(omicData))
  if(length(sample.intersect) == 0) stop("samples in omics data should match the names in experimental design")
  
  # select abundance from design table and reorder
  omicData <- dplyr::select(omicData, tidyselect::all_of(sample.intersect))
  
  # remove row with sum == 0
  matrix <- as.matrix(omicData)
  # nbr of genes with 0 count
  genes_flt0  <- rownames(matrix[rowSums(matrix) <= 0, ])
  # remove 0 count
  matrix.filt  <- matrix[rowSums(matrix)  > 0, ]
  
  # create SE object
  colData   <- dplyr::mutate(ExpDesign, samples=row.names(ExpDesign)) |>
    dplyr::filter(samples %in% sample.intersect) |> 
    tidyr::unite("groups", all_of(factorBio), sep = "_", remove = FALSE)
  
  for (factor in c(factorBio, factorBatch)){
    
    F.levels <- levels(colData[[factor]])
    colData[[factor]] <- factor(colData[[factor]], levels = intersect(F.levels, unique(colData[[factor]])))
  }
  
  order_levels <- with(colData, do.call(order, colData[c(factorBio, factorBatch)]))
  colData$samples <- factor(colData$samples, levels = unique(colData$samples[order_levels]))
  colData$groups  <- factor(colData$groups,  levels = unique(colData$groups[order_levels]))
  
  metadata <- list(omicType = omicType, Groups = colData, 
                   design = list(factorType = design[intersect(names(design), names(colData))]), 
                   DataProcessing = list(rowSumsZero = genes_flt0,
                                         Filtering = NULL, 
                                         Normalization =  list(setting = list(method = "none"), results = NULL,  normalized = FALSE), 
                                         Transformation = list(setting = list(method = "none"), results = NULL,  transformed = FALSE)))
  
  SE <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(abundance = as.matrix(matrix.filt)), 
                                                   colData = DataFrame(colData), metadata = metadata)
  
  rflomicsSE <- new("RflomicsSE")
  for(slot in slotNames(SE)) {
    slot(rflomicsSE, slot) <- slot(SE, slot)
  }
  
  return(rflomicsSE)
}



