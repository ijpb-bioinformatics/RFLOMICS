
# ---- generateReport ----
#' @title Generate RFLOMICS rmarkdown report
#' @description
#' This function is used to generate a html report from a RFLOMICS
#' \link{RflomicsMAE}.
#' @param object a \link{RflomicsMAE} produced by RFLOMICS.
#' @param fileName Name of the html document.
#' @param archiveName name of result archive
#' @param export boolean value (default: FALSE)
#' @param tmpDir temporary dir (default: getwd()) 
#' @importFrom rmarkdown render
#' @exportMethod generateReport
#'
methods::setMethod(f          = "generateReport",
                   signature  = "RflomicsMAE",
                   definition = function(object, fileName = NULL, archiveName = NULL, 
                                         export = FALSE, tmpDir = getwd(), ...){
  
  # Copy the report file to a temporary directory before processing it, in
  # case we don't have write permissions to the current working dir (which
  # can happen when deployed).
  tempReport <-  file.path(path.package("RFLOMICS"), "/RFLOMICSapp/report.Rmd")
  
  # Check if the object is properly filled.
  ## function
  
  # project name
  projectName <- getProjectName(object)
  RDataName   <- paste0(projectName, ".MAE.RData")
  
  # tmp dir
  if(file.access(tmpDir, 2) != 0) stop("No writing access in ", tmpDir)
  tmpDir <- file.path(tmpDir, paste0(projectName, "_report"))
  dir.create(tmpDir, showWarnings=FALSE)
  
  # html name
  if(is.null(fileName)) fileName <- file.path(tmpDir, paste0(projectName, "_report.html"))
  
  # save FE rflomics.MAE in .Rdata and load it during report execution
  rflomics.MAE <- object
  save(rflomics.MAE, file = file.path(tmpDir, RDataName))
  
  # Set up parameters to pass to Rmd document
  param.list <- list(FEdata = file.path(tmpDir, RDataName),
                     title  = paste0(projectName, " project"),
                     outDir = tmpDir)
  
  # Knit the document, passing in the `params` list, and eval it in a
  # child of the global environment (this isolates the code in the document
  # from the code in this app).
  render(input             = tempReport, 
        output_file       = fileName,
        params            = param.list,
        knit_root_dir     = tmpDir,
        intermediates_dir = tmpDir,
        envir = new.env(parent = globalenv()), ...)
  
  #Export results
  if(isTRUE(export)){
    
    if(is.null(archiveName)) 
      archiveName <- file.path(dirname(tmpDir), 
                               paste0(projectName, "_archive.tar.gz"))
    
    # cp html in tmpDir
    file.copy(from = fileName, to = tmpDir)
    cmd <- paste0("tar -C ", dirname(tmpDir), " -czf ", archiveName, " ", basename(tmpDir))
    system(cmd)
    print(cmd)
    
  }else{
    file.copy(from = fileName, to = dirname(tmpDir))
  }
  unlink(tmpDir, recursive = TRUE)
  
})

# ---- runOmicsPCA ----
#' @title runOmicsPCA
#' @description This function performs a principal component analysis on omic data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a "Normalization" slot is present in the metadata slot, then data are normalized before running the PCA according to the indicated transform method.
#' @param object An object of class \link{RflomicsSE-class}.
#' @param ncomp Number of components to compute. Default is 5.
#' @param raw boolean. Does the pca have to be ran on raw data or transformed and normalized data? Default is FALSE, pca is ran on transformed and normalized data.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runOmicsPCA
#' @importFrom FactoMineR PCA
#' @rdname runOmicsPCA
#'
methods::setMethod(f          = "runOmicsPCA",
                   signature  = "RflomicsSE",
                   definition = function(object, ncomp = 5, raw = FALSE){
                     object2 <- RFLOMICS:::.checkTransNorm(object, raw = raw)
                     pseudo  <- SummarizedExperiment::assay(object2)
                     if (raw) object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = ncomp, graph = FALSE)
                     else    object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = ncomp, graph = FALSE)
                     return(object)

                   }
)

#' @rdname runOmicsPCA
#' @title runOmicsPCA
#' @param SE.name the name of the data the normalization have to be applied to.
#' @exportMethod runOmicsPCA
methods::setMethod(f          = "runOmicsPCA",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, ncomp = 5, raw = FALSE){

                     object[[SE.name]] <-  runOmicsPCA(object[[SE.name]],
                                                  ncomp = ncomp,
                                                  raw  = raw)

                     return(object)

                   })

# runOmicsPCA <- function(object, ncomp = 5, raw = FALSE){
#   object2 <- RFLOMICS:::.checkTransNorm(object, raw = raw)
#   pseudo  <- SummarizedExperiment::assay(object2)
#   if (raw) object@metadata[["PCAlist"]][["raw"]] <- FactoMineR::PCA(t(pseudo), ncp = ncomp, graph = FALSE)
#   else    object@metadata[["PCAlist"]][["norm"]] <- FactoMineR::PCA(t(pseudo), ncp = ncomp, graph = FALSE)
#   return(object)
# } 

# ---- plotPCA ----
#' @title plotPCA
#' @description This function plot the factorial map from a PCA object stored
#' in a \link{RflomicsSE-class} object. By default, samples are
#' colored by groups (all combinations of level's factor)
#' @param object An object of class \link{RflomicsSE-class}
#' @param raw This argument indicates whether the scaled PCA has to be performed on raw [\sQuote{raw}] or normalized [\sQuote{norm}] data.
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
                   definition <- function(object, raw=c("raw","norm"), axes=c(1,2), groupColor="groups"){
                     
                     if(length(axes) != 2) stop("PCA axes must be a vector of length 2")
                     
                     ExpDesign <- getDesignMat(object)
                     
                     #
                     PC1 <- paste("Dim.",axes[1], sep="")
                     PC2 <- paste("Dim.",axes[2], sep="")
                     
                     if(PC1 == PC2) PC2 <- PC1+1
                     
                     score     <- object@metadata$PCAlist[[raw]]$ind$coord[, axes] %>% as.data.frame() %>%
                       dplyr::mutate(samples=row.names(.)) %>% dplyr::right_join(., ExpDesign, by="samples")
                     
                     var1 <- round(object@metadata$PCAlist[[raw]]$eig[axes,2][1], digits=3)
                     var2 <- round(object@metadata$PCAlist[[raw]]$eig[axes,2][2], digits=3)
                     
                     omicsType <- getOmicsTypes(object)
                     
                     switch (raw,
                             "raw"  = {title <- paste0("Raw ", omicsType, " data")},
                             "norm" = {title <- switch (omicsType,
                                                        "RNAseq" = { paste0("Filtred and normalized ", omicsType, 
                                                                            " data (", getNormSettings(object)$method, ")")  },
                                                        "proteomics" = {paste0("Transformed and normalized ", omicsType, 
                                                                               " data (", getTransSettings(object)$method, 
                                                                               " - norm: ", getNormSettings(object)$method, ")")},
                                                        "metabolomics" = {paste0("Transformed and normalized ", omicsType, 
                                                                                 " data (", getTransSettings(object)$method, 
                                                                                 " - norm: ", getNormSettings(object)$method, ")")}
                             )}
                     )
                     
                     
                     
                     p <- ggplot2::ggplot(score, ggplot2::aes_string(x=PC1, y=PC2, color=groupColor))  +
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
                     aa <- dplyr::select(score, tidyselect::all_of(groupColor), PC1, PC2)
                     bb <- FactoMineR::coord.ellipse(aa, bary = TRUE)
                     p <- p + ggplot2::geom_polygon(data = bb$res, ggplot2::aes_string(x=PC1, y=PC2, fill = groupColor),
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
                   definition = function(object, SE.name, raw, axes=c(1,2), groupColor="groups"){
                     
                     plotPCA(object[[SE.name]], raw, axes, groupColor)
                     
                   })

# ---- resetFlomicsMultiAssay ----
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


