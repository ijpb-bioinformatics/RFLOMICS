##### Global import ####

#' @importFrom ggplot2 ggplot geom_col theme_classic aes
#' theme element_text element_blank ylab xlab ggtitle
#' scale_fill_gradientn geom_tile theme_bw guides scale_fill_gradient2
#' guide_colourbar labs

# @export
#' @importFrom magrittr "%>%" 
magrittr::`%>%`

### ============================================================================
### Create Rflomics object / RflomicsSE and RflomicsMAE
### ----------------------------------------------------------------------------

# ---- RflomicsMAE CLASS Constructor for managing omics DATA and RESULTS ----

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
    SummarizedExperimentList[[data]] <- runOmicsPCA(RflomicsSE, raw = TRUE)
    #SummarizedExperimentList[[data]] <- RflomicsSE
    
    
    # metadata for sampleMap for RflomicsMAE
    listmap[[data]] <- data.frame(primary = as.vector(SummarizedExperimentList[[data]]@colData$samples),
                                  colname = as.vector(SummarizedExperimentList[[data]]@colData$samples),
                                  stringsAsFactors = FALSE)
    
    # 
    colnames <- c(names(omicList[[omicType]]), k)
    omicList[[omicType]] <- c(omicList[[omicType]] ,data)
    names(omicList[[omicType]]) <- colnames
    
  }
  
  # prepMAE <- MultiAssayExperiment::prepMultiAssay( ExperimentList = SummarizedExperimentList,
  #                                                  sampleMap      = MultiAssayExperiment::listToMap(listmap),
  #                                                  colData        = ExpDesign, outFile = stdout())
  # 
  # 
  # MAE <- MultiAssayExperiment::MultiAssayExperiment(experiments = prepMAE$experiments,
  #                                                   colData     = prepMAE$colData,
  #                                                   sampleMap   = prepMAE$sampleMap,
  #                                                   metadata    = list(omicList = omicList, projectName = projectName, design = Design)) 
  # 
  # rflomicsMAE <- new("RflomicsMAE")
  # for(slot in slotNames(MAE)) {
  #   slot(rflomicsMAE, slot) <- slot(MAE, slot)
  # }
  
  RfMAE <- RflomicsMAE(experiments = SummarizedExperimentList,
                              colData     = ExpDesign,
                              sampleMap   = MultiAssayExperiment::listToMap(listmap),
                              metadata    = list(omicList = omicList, projectName = projectName, design = Design))
  
  # tag as raw data (le temps de trouver une solution pour ne pas faire co-exister les raw et les process)
  names(RfMAE) <- paste(names(RfMAE), "raw", sep = ".")
  return(RfMAE)
}



#' @title RflomicsMAE is RflomicsMAE object constructor.
#' @description
#' A short description...
#' 
#' @param name description
#' @seealso MultiAssayExperiment
#' @return An object of class \link{RflomicsMAE-class}
#' @name RflomicsMAE
#' @rdname RflomicsMAE
#' @export
#'
RflomicsMAE <- function(experiments = NULL, colData = NULL, sampleMap = NULL, metadata = NULL, ...){
  
  MAE <- MultiAssayExperiment::MultiAssayExperiment(experiments, colData, sampleMap, metadata, ...)
  
  rflomicsMAE <- new("RflomicsMAE")
  for(slot in slotNames(MAE)) {
    slot(rflomicsMAE, slot) <- slot(MAE, slot)
  }
  
  return(rflomicsMAE)
}



# ---- createRflomicsSE CLASS Constructor for managing omics DATA and RESULTS ----

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
                                         Transformation = list(setting = list(method = "none"), results = NULL,  transformed = FALSE)
                                         ))
  
  SE <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(abundance = as.matrix(matrix.filt)), 
                                                   colData = DataFrame(colData), metadata = metadata)

  rflomicsSE <- new("RflomicsSE")
  for(slot in slotNames(SE)) {
    slot(rflomicsSE, slot) <- slot(SE, slot)
  }
  
  return(rflomicsSE)
}







#' @title Read Experimental Design
#'
#' @param file path to experimental design file
#' @return data.frame
#' @importFrom dplyr  mutate across 
#' @importFrom tidyselect where
#' @importFrom stringr str_remove_all fixed
#' @importFrom purrr reduce
#' @importFrom vroom vroom
#' @export
#' 
read_exp_design <- function(file){
  
  if (missing(file)) {
    stop('Please provide a file path')
  }
  
  if(!file.exists(file))
  {
    stop(file, " don't exist!")
    return(NULL)
  }
  
  # read design and remove special characters
  # remove "_" from modality and factor names
  data <- vroom::vroom(file, delim = "\t", show_col_types = FALSE) %>%
    dplyr::mutate(dplyr::across(.cols = tidyselect::where(is.character), ~stringr::str_remove_all(.x, pattern = "[.,;:#@!?()§$€%&<>|=+-/]"))) %>%
    dplyr::mutate(dplyr::across(.cols = tidyselect::where(is.character), ~stringr::str_remove_all(.x, pattern = "[\\]\\[\'\"\ ]"))) %>%
    dplyr::mutate(dplyr::across(.cols = tidyselect::where(is.character), ~stringr::str_remove_all(.x, pattern = stringr::fixed("\\")))) %>% 
    dplyr::mutate(dplyr::across(.cols = c(-1), ~stringr::str_remove_all(.x, pattern = stringr::fixed("_")))) %>% 
    dplyr::mutate(dplyr::across(.cols = tidyselect::where(is.character), as.factor)) 
  
  names(data)  <- stringr::str_remove_all(string = names(data), pattern = "[.,;:#@!?()§$€%&<>|=+-/\\]\\[\'\"\ _]") %>%
    stringr::str_remove_all(., pattern = stringr::fixed("\\"))
  
  # check if there is duplication in sample names
  sample.dup <- as.vector(data[which(table(data[1]) > 1),1])[[1]]
  
  if (length(sample.dup) != 0) {
    
    stop("Duplicated sample names: ", paste0(sample.dup, collapse = ","))
  }
  
  # check if there is duplication in factor names
  # factor.dup <- as.vector(data[which(table(names(data[-1])) > 1),1])[[1]]
  factor.dup <- names(data[-1])[duplicated(names(data[-1]))]
  if (length(factor.dup) != 0) {
    stop("Duplicated factor name: ", paste0(factor.dup, collapse = ","))
  }
  
  # check if same name of moralities are used in diff factor
  mod.list <- sapply(names(data[-1]), function(x){ 
    unique(data[-1][[x]])
  }) %>% purrr::reduce(c)
  
  mod.dup <- mod.list[duplicated(mod.list)]
  if(length(mod.dup) != 0) {
    
    stop("Modality used in more than one factor: ", paste0(mod.dup[1:10], collapse = ", "))
  }
  
  # warning if number of factors exceed n = 10
  n <- 10
  if (dim(data)[2]-1 >= n){
    
    data <- data[, 1:n]
    warning("Large number of columns! only the first ", n," will be displayed")
  }
  
  # check nbr of modality of the 5th fist columns
  index <- sapply(names(data[-1]), function(x){ if(length(unique(data[-1][[x]]))>n){ FALSE }else{ TRUE } })
  F.mod <- names(data[-1])[index]
  
  ratio <- length(F.mod)/length(names(data[-1]))
  
  if(ratio != 1)
  {
    warning("The select input contains a large number of options")
  }
  
  data            <- data.frame(data) 
  row.names(data) <- data[,1]
  data            <- data[,-1]
  return(data)
}



#' @title Read omics data 
#'
#' @param file omics data matrix
#' @return data.frame
#' @importFrom vroom vroom
#' @importFrom stringr str_remove_all fixed
#' @export
#'

read_omics_data <- function(file){
  
  if(!file.exists(file))
  {
    stop(file, " don't exist!")
  }
  
  # read omics data and remove special characters
  data <- vroom::vroom(file, delim = "\t", show_col_types = FALSE)
  names(data)  <- stringr::str_remove_all(string = names(data), pattern = "[.,;:#@!?()§$€%&<>|=+-/\\]\\[\'\"\ ]") %>%
    stringr::str_remove_all(., pattern = stringr::fixed("\\"))
  
  # check if there is duplication in sample names
  sample.dup <- as.vector(data[which(table(names(data[-1])) > 1),1])[[1]]
  
  if (length(sample.dup) !=0){
    
    stop("Duplicated sample names: ", paste0(sample.dup, collapse = ","))
  }
  
  # check if there is duplication in factor names
  entity.dup <- as.vector(data[which(table(data[1]) > 1),1])[[1]]
  
  if (length(entity.dup) !=0){
    
    stop("Duplicated feature names: ", paste0(entity.dup, collapse = ","))
  }
  
  data            <- data.frame(data) 
  row.names(data) <- data[,1]
  data            <- data[,-1]
  return(data)
}



#' GetModelFormulae
#'
#' From a vector of character giving the name of the factors of an omics experiment,
#' and their type of effect: biological or batch, it returns all models formulae
#' that can be formulated in association with this factors. Batch effect factors do
#' not appear in interaction terms with biological factor. Model formulae stop in
#' second order interaction.
#'
#' @param FacBio a vector of character giving the name of the bio factors.
#' @param FacBatch a vector of character giving the name of the batch factors.
#' @param MAE a RflomicsMAE produced by RFLOMICS. Default is null. If not null, overwrite FacBio and FacBatch.
#'
#' @return a named list of object of class formula
#' @export
#' @noRd
#' @examples
#'
#' GetModelFormulae(FacBio=c("Genotype","Temperature","Environment"), FacBatch=c("Replicat"))
#' GetModelFormulae(FacBio=c("Genotype","Temperature"), FacBatch=c("Replicat"))
#' GetModelFormulae(FacBio=c("Genotype"), FacBatch=c("Replicat"))
#'
#' GetModelFormulae(FacBio=c("Genotype","Temperature"), FacBatch=c("Replicat", "laboratory"))
#' 
GetModelFormulae <- function(FacBio=NULL, FacBatch=NULL, MAE = NULL){
  
  MAE <- MAE
  
  if (!is.null(MAE)) { 
    facTypes <- getFactorTypes(MAE)
    FacBio   <- bioFactors(MAE) #names(facTypes)[facTypes == "Bio"]
    FacBatch <- batchFactors(MAE) #names(facTypes)[facTypes == "batch"]
  }
  
  # Initialize
  formulae <- list()
  
  # Verify that nbr of bio factors are between 1 and 3.
  if(!length(FacBio) %in% 1:3) stop(".... !")
  
  # Verify that nbr of batch factors are between 1 and 2.
  if(!length(FacBatch) %in% 1:2) stop(".... !")
  
  nFac <- length(FacBio)
  
  # get formulae without interation
  formulae[[1]] <- update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+",paste(FacBio,collapse="+"))),new=~.)
  
  # get formulae with interation if nbr of FacBio > 1
  if(nFac !=1)
    formulae[[2]] <- update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+","(",paste(FacBio,collapse="+"),")^2")),new=~.)
  
  formulae <- unlist(formulae)
  names(formulae) <- unlist(as.character(formulae))
  
  # Sort formulae
  
  formulae <- formulae[order(unlist(lapply(names(formulae),nchar)),decreasing=TRUE)]
  
  return(formulae)
}





#' Plot the balance of data in an experimental design
#'
#' This function provides easy visualization of the balance of data in a data set given a specified experimental design. This function is useful for identifying
#' missing data and other issues. The core of this function is from the function ezDesign in the package ez.
#'
#' @param counts : the number of data in each cell of the design
#' @param cell_border_size : Numeric value specifying the size of the border seperating cells (0 specifies no border)
#'
#' @return A printable/modifiable ggplot2 object.
#' @export
#' @importFrom ggplot2 aes_string theme facet_grid labs element_rect geom_rect scale_x_continuous scale_y_continuous facet_grid
#' @importFrom dplyr mutate if_else
#' @noRd
plotExperimentalDesign <- function(counts, cell_border_size = 10, message=""){
  if (names(counts)[ncol(counts)] != "Count"){
    stop("the last column of the input data frame must be labelled Count")
  }
  if(ncol(counts) < 2){
    stop("data frame with less than 2 columns")
  }
  
  # #add color column
  # # #00BA38
  
  counts <- counts %>% dplyr::mutate(status = dplyr::if_else(Count > 2 , "pass", dplyr::if_else(Count == 2 , "warning", "error")))
  
  #list of factor names
  factors <- names(counts)[1:(dim(counts)[2]-2)]
  
  col.panel <- c("pass", "warning", "error")
  names(col.panel) <- c("#00BA38", "orange", "red")
  
  col.panel.u <- col.panel[col.panel %in% unique(counts$status)]
  
  switch (length(factors),
          "1" = { p <- ggplot2::ggplot(counts ,ggplot2::aes_string(x = factors[1], y = 1)) + ggplot2::theme(axis.text.y = ggplot2::element_blank()) + ggplot2::ylab("") },
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
                   axis.ticks = ggplot2::element_blank(), axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
    ggplot2::ggtitle(message)
  
  return(p)
}



######## INTERNAL - CHECKS FUNCTIONS ###########

# .check_NA: checks if there are NA/nan in the RflomicsSE assay
#' @title .checkNA
#'
#' @param object An object of class \link{RflomicsSE}
#' @return boolean. if TRUE, NA/nan are detected in the SE::assay.
#' @keywords internal
#' @noRd
#'
.checkNA <- function(object) {
  NA_detect <- ifelse(any(is.na(assay(object))), TRUE, FALSE)
  return(NA_detect)
}


# ----- INTERNAL - Check if character vectors are tags Names : -----

#' @title Check if character vectors are tags Names
#'
#' @param object a MAE object or a SE object (produced by Flomics). If it's a RflomicsSE, expect to find
#'  a slot of differential analysis.
#' @param tagName vector of characters.
#' @return boolean. TRUE if all of tagName are indeed tags Names.
#' @noRd
#' @keywords internal
isTagName <- function(object, tagName) {
  df_contrasts <- getSelectedContrasts(object)
  
  search_match <- sapply(tagName, FUN = function(cn) {
    grep(cn, df_contrasts$tag, fixed = TRUE)
  })
  search_success <- sapply(search_match, identical, integer(0)) # if TRUE, not a success at all.
  
  if (!any(search_success)) {
    # Congratulations, it's a tag name!
    return(TRUE)
  } else {
    return(FALSE)
  }
}




# ---- INTERNAL - get variable name and type from omicstype ----

#' @title Omics Dictionary
#'
#' @param object a MAE object or a SE object (produced by Flomics). Expect to find a omicsType somewhere.
#' @param SE.name if object is a MAE, expect to find the experiment name from which the omics info has to be retrieved.
#' @return list of two elements: variableName and valueType.
#' @noRd
#' @keywords internal

omicsDic <- function(object, SE.name = NULL){
  
  if (!is(object, "RflomicsSE") && !is(object, "RflomicsMAE")) {
    stop("Object must be a RflomicsSE or a RflomicsMAE, not a ",
         class(object))
  }
  
  if (is(object, "RflomicsMAE")) {
    if (missing(SE.name)) {
      stop("Please provide an Experiment name (SE.name).")
    }
    
    object <- object[[SE.name]]
  }
  
  omicsType <- getOmicsTypes(object)
  
  valReturn <- switch(omicsType,
                      "RNAseq"       =  list("variableName" = "transcripts",
                                             "valueType" = "counts"),
                      "proteomics"   =  list("variableName" = "proteins",
                                             "valueType" = "XIC"),
                      "metabolomics" =  list("variableName" = "metabolites",
                                             "valueType" = "XIC")
  )
  
  return(valReturn)
  
}



