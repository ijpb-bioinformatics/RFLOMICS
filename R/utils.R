#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' GetDesignFromNames
#'
#' @param samples_name a vector of sample names giving the designs factor, each
#' separated by "_"
#'
#' @return a named list of two elements
#' * nb_dfac: the number of design factors
#' * tmpDesign: a data frame with the names of design factors in columns (ie: dFac1) and all
#' the modalities in row
#' @md
#' @export
#'
#' @examples
#'
#'
GetDesignFromNames <- function(samples_name){

  # Get the number of design factor and the factors from the names of the count matrix
  nb_dFac <- stringr::str_count(samples_name,pattern="_")+1
  # Test if the number of factor are the same for all sample names
  try(if(var(nb_dFac) != 0 ) stop("Column names do not have the same level of factor"))
  nb_dFac <- nb_dFac[1]
  names(nb_dFac) <- "n_dFac"

  #  Get all the factors from the names of the count matrix
  tmpDesign <- tibble::tibble(design=samples_name) %>%
    tidyr::separate(.,design,into=paste("dFac",1:nb_dFac["n_dFac"],sep="")) %>%
    dplyr::mutate_all(.,as.factor)

  return(list("nb_dFac"=nb_dFac,"tmpDesign"=tmpDesign))
}

#' GetModelFormulae
#'
#' @param Factors.Name
#' @param Factors.Type
#'
#' @return a named list of formulae
#' @export
#'
#' @examples
#'
GetModelFormulae <- function(Factors.Name,Factors.Type){


  formulae <- list()

  FacBio <- Factors.Name[which(Factors.Type == "Bio")]
  FacBatch <- Factors.Name[which(Factors.Type == "batch")]

  nFac <- length(FacBio)

  getF <- function(x,FacBatch){
    update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+","(",paste(x,collapse="+"),")^2")),new=~.)
  }
  getF2 <- function(x,FacBatch){
    update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+",paste(x,collapse="+"))),new=~.)
  }

  for(i in 1:nFac){

    formulae <- c(formulae, apply(combn(FacBio,i),2,getF, FacBatch=FacBatch))
    if(i !=1){
      formulae <- c(formulae, apply(combn(FacBio,i),2,getF2, FacBatch=FacBatch))
    }
    #formulae[[i]]  <- apply(combn(FacBio,i),2,getF, FacBatch=FacBatch)
  }

  formulae[[nFac+1]]<-apply(combn(FacBio,i),2,getF2,FacBatch=FacBatch)

  formulae <- unlist(formulae)
  names(formulae) <- unlist(as.character(formulae))

  #formulae[[nFac+1]]<- apply(combn(FacBio,i),2,getF2,FacBatch=FacBatch)

  #formulae <- unlist(formulae)
  #names(formulae) <- unlist(as.character(formulae))

  return(formulae)
}


#' GetContrasts
#'
#' @param An object of class design
#'
#' @return An object of class design
#' @export
#'
#' @examples
#'
GetContrasts <- function(design){
  # function qui est en cours d'ecriture par Christine
  data.frame()
}

#' @title TMM.Normalization
#' Interface to the calcNormFactors functionof the edgeR package  with the choosen TMM parameters as the normalization method
#' @param counts numeric matrix of read counts
#' @param groups vector or factor giving the experimental group/condition for each sample/library.
#' @return a data.frame with a row for each sample and columns group, lib.size and norm.factors containing the group labels, library sizes and normalization factors. Other columns can be optionally added to give more detailed sample information.
#' @export
#'
#' @examples
TMM.Normalization <- function(counts, groups){
  dge <- edgeR::DGEList(counts=counts, group=groups)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  nf  <- dge$samples
  return(nf)
}


#' @title edgeR.AnaDiff
#'
#' @param object an object of class [\code{\link{MultiAssayExperiment}]
#' @param data Omic data type
#' @param FDR The false discovery rate to apply
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @return A list of object of class [\code{\link{DGELRT}]
#' @export
#'
#' @examples
edgeR.AnaDiff <- function(object, data, FDR, clustermq){


  # retrieve the design matrix
  model_matrix <- model.matrix(as.formula(object@metadata$design@Model.formula),
                               data=as.data.frame(object@metadata$design@List.Factors))

  # Construct the DGE obect
  dge <- edgeR::DGEList(counts=assay(object@ExperimentList[[data]]),
                        group=object@ExperimentList[[data]]@metadata$Normalization$coefNorm$group,
                        lib.size =object@ExperimentList[[data]]@metadata$Normalization$coefNorm$lib.size,
                        norm.factors = object@ExperimentList[[data]]@metadata$Normalization$coefNorm$norm.factors)

  # Run the model
  dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)
  dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)
  dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)
  fit.f <- edgeR::glmFit(dge,design=model_matrix)

  # test clustermq

  if(clustermq == TRUE){

 # Fonction to run on contrast per job
 # y is the model, Contrasts are stored in a matrix, by columns

  fx <- function(x){
    edgeR::glmLRT(y, contrast = z[,x])
  }

  ListRes <- clustermq::Q(fx, x=1:length(object@metadata$design@Contrasts.Sel),
    export=list(y=fit.f,z=object@metadata$design@Contrasts.Coeff),
    n_jobs=length(object@metadata$design@Contrasts.Sel),pkgs="edgeR")
  }

  else{
   ListRes <-  lapply(object@metadata$design@Contrasts.Sel, function(x){
     resglm <- edgeR::glmLRT(fit.f, contrast = object@metadata$design@Contrasts.Coeff[,x])
     return(resglm)
   })
  }

  print(ListRes)
  names(ListRes) <- object@metadata$design@Contrasts.Sel
  return(ListRes)
}



#' @title colorPlot
#'
#' @param design
#' @param ColData
#' @param condition
#'
#' @return
#' @export
#'
#' @examples
colorPlot <- function(design, ColData, condition="samples"){

  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

    if(condition == "samples"){

      # combine only bio fact
      groups <- design@List.Factors[design@Factors.Type == "Bio"] %>%
        as.data.frame() %>% unite(col="groups", sep="_")
      list.cond <- factor(groups$groups)
    }
    else{

      list.cond <- design@List.Factors[[condition]]
    }

    len.cond  <- length(levels(list.cond))

    colors    <- getPalette(len.cond)

    #col <- colors[list.cond]
    col <- colors[levels(list.cond)]
    #names(col) <- row.names(ColData)
    names(col) <- row.names(levels(list.cond))

    return(col)
}


#' plotLibSize
#'
#' @param abundances
#' @param dataName
#' @param pngFile 
#' @return plot
#' @export
#'
#' @examples
plotLibSize <- function(abundances, dataName, pngFile=NULL){

  samples     <- colnames(abundances)
  libSizeNorm <- data.frame ( value = colSums(abundances, na.rm = TRUE) , samples=samples)

  libSizeNorm$samples <- factor(libSizeNorm$samples, levels = libSizeNorm$samples)

  p <- ggplot(libSizeNorm, aes(x=samples,y=value, fill=samples)) + geom_bar( stat="identity" ) +
    xlab(paste0(dataName, " samples")) + ylab("Library Size") + 
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          legend.position  = "none")
          #axis.text.x     = element_blank(),
          #axis.ticks      = element_blank())
          #legend.key.size = unit(0.3, "cm"))
          #legend.text     = element_text(size=5))
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }
}


#' plotDistr
#'
#' @param abundances matrix or dataframe of feature/gene abundances/counts
#' @param dataName name of dataset
#' @param pngFile png file name
#' @return plot
#' @export
#'
#' @examples plotDistr(assay(MAE), "dataset1", "tmp/countDist.png")
plotDistr <- function(abundances, dataName, pngFile=NULL){

  pseudo_counts <- log2(abundances+1) %>% reshape2::melt()
  colnames(pseudo_counts) <- c("features", "samples", "value")

  p <- ggplot(pseudo_counts) +
    geom_density(aes(value, color=samples) ) +
    xlab(paste0(dataName, " log2(feature abundances)")) + 
    theme(legend.position='none')
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }
}




#' pvalue.plot
#'
#' @param data 
#' @param contrast
#' @param pngFile 
#' @return plot
#' @export
#'
#' @examples
pvalue.plot <- function(data, contrast, pngFile=NULL){

  p <- ggplot(data=data) + geom_histogram(aes(x=PValue), bins = 200)
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }

}




