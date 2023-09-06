#' @title mvQCdesign
#' @description mvQCdesign is for multivariate quality check of design. For each design factor (one color for each),
#' and each PCA axis this function plot the coordinates of the sample in a PCA axis (y-axis) in an
#' increasing order along the x-axis. It allows to have a quick view of the variability associated to each factor.
#' @param object An object of class \link{MultiAssayExperiment}
#' @param data The name of the omic data for which the PCA plot has to be drawn
#' @param axis The number of PCA axis to keep
#' @param PCA This argument indicates which PCA results to plot: raw ("raw") or normalised ("normalised")
#' @param pngFile The name of the png file for saving the plot.
#' @exportMethod mvQCdesign
#'
#' @rdname mvQCdesign

methods::setMethod(f="mvQCdesign",
                   signature="MultiAssayExperiment",
                   definition <- function(object, data, PCA=c("raw","norm"), axis=5, pngFile=NULL){
                     
                     # Stop if the PCA object does not exist
                     resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
                     
                     if(is.null(resPCA)){
                       stop(paste0(PCA,"PCA does not exist for the ",data," data"))
                     }
                     
                     n_dFac <-  object@metadata$colDataStruc["n_dFac"]
                     
                     # corespondance between coordinate and factor modalities thanks to the sample's name
                     tab_tmp <- merge(as.data.frame(resPCA$ind$coord),as.data.frame(object@colData),by='row.names',all=TRUE)
                     
                     df <- list()
                     bigdf <- list()
                     
                     for(i in 1:n_dFac){
                       
                       Factor <-  object@colData[,i]
                       FactorName <- names(object@colData)[i]
                       nF <- length(Factor)
                       
                       for(j in 1:axis){
                         Dim=paste0("Dim.",j)
                         # select the coordinates and the factor columns
                         df[[j]] <- dplyr::select(tab_tmp, all_of(FactorName),all_of(Dim)) %>%
                           # rename
                           dplyr::rename(.,"Levels"=FactorName,"y"=starts_with("Dim.")) %>%
                           # sort by factor modalities then by coordinate
                           dplyr::arrange(Levels,y) %>%
                           # add column
                           dplyr::mutate(.,"x"=1:nF,
                                         "Axis"=rep(paste("PCA",j, "\n(",round(resPCA$eig[j,2],1),"%)",sep=""),nF),
                                         "FactorN"=rep(FactorName,nF))
                         
                       }
                       bigdf[[i]] <- dplyr::bind_rows(df)
                     }
                     big <- dplyr::bind_rows(bigdf)
                     
                     out <- by(data = big, INDICES = big$FactorN, FUN = function(m) {
                       
                       m <- droplevels(m)
                       
                       m <- ggplot(m,aes(y=y,x=x,colour=Levels))+
                         ggplot2::geom_bar(stat = "identity",position = ggplot2::position_dodge(),aes(fill=Levels))+
                         facet_grid(as.factor(FactorN)~Axis) +
                         labs(x = "Samples", y="Coordinates \n on the PCA axis")+
                         theme(axis.title.y = element_text(size = 10),
                               axis.title.x=element_text(size = 10),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
                     })
                     p <- do.call(gridExtra::grid.arrange, c(out,ncol=1))
                     print(p)
                     
                     
                     if(! is.null(pngFile)){
                       ggplot2::ggsave(filename = pngFile,  plot = p)
                     }
                   })

# copier coldataStruct dans metadata
# tester que PCA existe
#

#' @title mvQCdata
#' @description mvQCdata is for multivariate quality check of metadata.
#' This function helps to control if some experimental parameters given as metadata (numeric one) in input
#' explain much variability than expected in the data or if their effect could be confused with biological one.
#' This function correlates quantitative variable describing technical aspect for each sample with
#' their coordinate on the PCA axis.
#' \itemize{
#' \item{Technical parameters from sample preparation as the day of the RNAseq library preparation}
#' \item{Statistics results after the bioinformatics workflow as the percent of sequences with primers or \% of rrna in the library}
#'  }
#' @param object An object of class \link{MultiAssayExperiment-class}
#' @param data The name of the omic data for which the PCA plot has to be drawn
#' @param axis The number of PCA axis to keep
#' @param PCA This argument indicates which type of PCA results to take: on raw ("raw") or normalized ("norm") data.
#' @param pngFile The name of the png file for saving the plot.
#' @exportMethod mvQCdata
#' @rdname mvQCdata
methods::setMethod(f="mvQCdata",
                   signature="MultiAssayExperiment",
                   definition <- function(object, data, PCA=c("raw","norm"),axis=3, pngFile=NULL){
                     
                     resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
                     cc <- c(RColorBrewer::brewer.pal(9, "Set1"))
                     
                     n_dFac <- object@metadata$colDataStruc["n_dFac"]
                     n_qcFac <- dim(object[[data]]@colData[,-c(1:2)])[2]
                     var <- names(object[[data]]@colData[,-c(1:2)])
                     
                     corA=list()
                     VarAxis = list()
                     for(i in 1:axis){
                       corA[[i]] <- cor(resPCA$ind$coord[,i],
                                        as.data.frame(object[[data]]@colData[,-c(1:2)]),
                                        method="spearman")
                       VarAxis[[i]] <- paste("\n(Var=",round(resPCA$eig[i,2],1),")",sep="")
                     }
                     
                     df = data.frame("QCparam"=rep(var,axis),
                                     "Spearman"= unlist(corA),
                                     "Axis"=rep(paste(rep("Axis",axis),1:axis,VarAxis,sep=""), each=n_qcFac))
                     
                     p <- ggplot(df,aes(x=Axis, y=abs(Spearman),fill=QCparam))+
                       geom_bar(stat="identity",position=ggplot2::position_dodge(),width=0.7)+ggplot2::ylim(0,1)+
                       ggplot2::labs(x = "Axis number", y="Cor(Coord_dFactor_PCA,QCparam)")
                     
                     print(p)
                     
                     if (!is.null(pngFile)){
                       ggplot2::ggsave(filename = pngFile, plot = p)
                     }
                     
                   })


# ---- Old Stuff ----

#' Enrichment.plot
#'
#' @param object An object of class \link{SummarizedExperiment}
#' @param Over_Under "overrepresented" or "underrepresented" (default : overrepresented)
#' @param top level of enriched terms to display
#' @param listNames vector of DGEs or cluster names (default : all enriched lists)
#' @param from  "DiffExpAnal" or "CoExpAnal". default "DiffExpAnal"
#' @return plot
#' @export
#' @exportMethod Enrichment.plot
#' @importFrom dplyr desc
Enrichment.plot <- function(object, Over_Under = "overrepresented", top = 50,
                            domain = NULL, listNames = NULL, from = c("DiffExpEnrichAnal", "CoExpEnrichAnal")) {
  Decision <- Pvalue_over <- Pvalue_under <- Pvalue <- NULL
  Term <- Domain <- Trial_Success <- scale_size <- tail <- NULL
  
  # if Over_Under are not recognized we choose default value == overrepresented
  if (!Over_Under %in% c("overrepresented", "underrepresented")) {
    Over_Under <- "overrepresented"
  }
  
  if (!is.numeric(top)) {
    top <- "NA"
    Top.tag <- ""
  } else {
    Top.tag <- paste0("Top ", top)
  }
  
  # if listNames is null we take all results
  if (is.null(listNames)) {
    listNames <- names(object@metadata[[from]][["results"]])
  }
  
  p <- list()
  for (listname in listNames) {
    data <- object@metadata[[from]][["results"]][[listname]][["Over_Under_Results"]] %>% dplyr::filter(Domain %in% domain)
    
    data_ord <- switch(Over_Under,
                       "overrepresented" = {
                         dplyr::filter(data, Decision == Over_Under) %>%
                           dplyr::arrange(desc(Pvalue_over)) %>%
                           dplyr::mutate(Pvalue = Pvalue_over)
                       },
                       "underrepresented" = {
                         dplyr::filter(data, Decision == Over_Under) %>%
                           dplyr::arrange(desc(Pvalue_under)) %>%
                           dplyr::mutate(Pvalue = Pvalue_under)
                       }
    )
    
    data_ord$Term <- factor(data_ord$Term, levels = data_ord$Term)
    
    Urn_effective <- data$Urn_effective[1]
    Trial_effective <- data$Trial_effective[1]
    
    p[[listname]] <- ggplot2::ggplot(data = tail(data_ord, n = top), aes(x = sort(Trial_Success), y = Term, size = Urn_Success, color = Pvalue)) +
      geom_point(alpha = 0.5) +
      scale_size(range = c(0.1, 10)) +
      scale_color_gradient(low = "blue", high = "red") +
      ylab("") +
      xlab("Count") +
      ggtitle(paste0(listname, " :\n ", Over_Under, " ", Top.tag, "terms in ", domain, "\n", " (Urn effective = ", Urn_effective, "; Trial effective = ", Trial_effective, ")"))
  }
  
  return(p)
}

#' @title runAnnotationEnrichment
#' @description This function performs enrichment test from functional gene annotation data. This data could be
#' GO, KEGG or other... For instance, the hypergeometric test is applied. Parameters used are those
#' recommended in DiCoExpress workflow (see the paper in reference)
#' @param object An object of class \link{SummarizedExperiment}
#' @param CoExpListNames A list of clusters names.
#' @param from  "DiffExpAnal" or "CoExpAnal". default "DiffExpAnal"
#' @param alpha The pvalue cut-off
#' @param probaMethod The probabilistic method to use.
#' @param annotation The gene annotation file.
#' @return An object of class \link{SummarizedExperiment}
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runAnnotationEnrichment
#'
methods::setMethod(
  f = "runAnnotationEnrichment",
  signature = "SummarizedExperiment",
  definition = function(object,
                        annotation,
                        alpha = 0.01,
                        probaMethod = "hypergeometric",
                        ListNames = object@metadata$DiffExpAnal[["contrasts"]]$contrastName,
                        from = "DiffExpAnal") {
    EnrichAnal <- list()
    EnrichAnal[["list.names"]] <- ListNames
    EnrichAnal[["alpha"]] <- alpha
    EnrichAnal[["proba.test"]] <- probaMethod
    
    
    ## list of gene list to annotate
    geneLists <- list()
    if (from == "DiffExpAnal") {
      geneLists <- lapply(ListNames, function(listname) {
        row.names(object@metadata$DiffExpAnal[["TopDEF"]][[listname]])
      })
      names(geneLists) <- ListNames
    } else if (from == "CoExpAnal") {
      # geneLists.coseq <- lapply(CoExpListNames, function(listname){
      
      geneLists <- object@metadata[["CoExpAnal"]][["clusters"]][ListNames]
      # })
      # names(geneLists.coseq) <- CoExpListNames
    }
    
    
    Results <- list()
    count.na <- 0
    for (geneList in names(geneLists)) {
      if (length(intersect(geneLists[[geneList]], annotation$geneID)) != 0) {
        Results[[geneList]] <- switch(probaMethod,
                                      "hypergeometric" = EnrichmentHyperG(annotation, geneLists[[geneList]], alpha = 0.01)
        )
      } else {
        Results[[geneList]] <- NULL
        count.na <- count.na + 1
      }
    }
    
    if (count.na == length(names(geneLists))) {
      EnrichAnal[["results"]] <- NULL
    } else {
      EnrichAnal[["results"]] <- Results
    }
    
    if (from == "DiffExpAnal") {
      object@metadata[["DiffExpEnrichAnal"]] <- EnrichAnal
    } else if (from == "CoExpAnal") {
      object@metadata[["CoExpEnrichAnal"]] <- EnrichAnal
    }
    
    return(object)
  }
)