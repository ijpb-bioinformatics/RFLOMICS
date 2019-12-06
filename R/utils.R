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
  
  FacBio <- Factors.Name[which(Factors.Type == "Bio")]
  FacBatch <- Factors.Name[which(Factors.Type == "batch")]
  
  nFac <- length(FacBio)
  
  getF <- function(x,FacBatch){
    update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+","(",paste(x,collapse="+"),")^2")),new=~.)
  }
  getF2 <- function(x,FacBatch){
    update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+",paste(x,collapse="+"))),new=~.)
  }
  
  formulae <- list()
  for(i in 1:nFac){
    
    formulae <- c(formulae, apply(combn(FacBio,i),2,getF, FacBatch=FacBatch))
    if(i !=1){
      formulae <- c(formulae, apply(combn(FacBio,i),2,getF2, FacBatch=FacBatch))
    }
    #formulae[[i]]  <- apply(combn(FacBio,i),2,getF, FacBatch=FacBatch)
  }
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
#'
#' @param counts
#'
#' @return
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
#' @param counts
#'
#' @return
#' @export
#'
#' @examples
edgeR.AnaDiff <- function(object,data,FDR, clustermq=FALSE){

  
  # retrieve the matrix design
  model_matrix <- model.matrix(as.formula(object@metadata$design@Model.formula),
                               data=as.data.frame(object@metadata$design@List.Factors))
  
  # Construct the DGE obect
  dge <- edgeR::DGEList(counts=assay(object@ExperimentList[[data]]),
                        group=object@ExperimentList[[data]]@metadata$Normalization$coefNorm$group,
                        lib.size =object@ExperimentList[[data]]@metadata$Normalization$coefNorm$lib.size,
                        norm.factors = object@ExperimentList[[data]]@metadata$Normalization$coefNorm$norm.factors)
  
  # Run the model
  dge <- estimateGLMCommonDisp(dge, design=model_matrix)
  dge <- estimateGLMTrendedDisp(dge, design=model_matrix)
  dge <- estimateGLMTagwiseDisp(dge, design=model_matrix)
  fit.f<-glmFit(dge,design=model_matrix)
  
  ListRes <-  lapply(object@metadata$design@Contrasts.Sel, function(x){

    resglm <- glmLRT(fit.f, contrast = object@metadata$design@Contrasts.Coeff[,x])
    return(resglm)
  })

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
#' @param design
#' @param colData
#' @param pngDir 
#' @return
#' @export
#'
#' @examples
plotLibSize <- function(abundances, pngDir){

  samples     <- colnames(abundances)
  libSizeNorm <- data.frame ( value = colSums(abundances, na.rm = TRUE) , samples=samples)

  libSizeNorm$samples <- factor(libSizeNorm$samples, levels = libSizeNorm$samples)

  p <- ggplot(libSizeNorm, aes(x=samples,y=value, fill=samples)) + geom_bar( stat="identity" ) +
    xlab("") + ylab("Library Size") +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          legend.position  = "none")
          #axis.text.x     = element_blank(),
          #axis.ticks      = element_blank())
          #legend.key.size = unit(0.3, "cm"))
          #legend.text     = element_text(size=5))
  print(p)
  ggsave(filename = "LibSize.png", path = pngDir, plot = p)
}


#' plotDistr
#'
#' @param abundances
#' @param design
#' @param colData
#' @param pngDir omic data
#' @return
#' @export
#'
#' @examples
plotDistr <- function(abundances, pngDir){


  pseudo_counts <- log2(abundances+1) %>% reshape2::melt()
  colnames(pseudo_counts) <- c("features", "samples", "value")

  p <- ggplot(pseudo_counts) +
    geom_density(aes(value, color=samples) ) +
    xlab("log2(feature abundances)") + theme(legend.position='none')
  print(p)
  ggsave(filename = "CountDist.png", path = pngDir, plot = p)
}




#' pvalue.plot
#'
#' @param data 
#' @param tag 
#' @param pngDir 
#'
#' @return
#' @export
#'
#' @examples
pvalue.plot <- function(data, tag , pngDir){

  p <- ggplot(data=data) + geom_histogram(aes(x=PValue), bins = 200)
  print(p)
  ggsave(filename = paste0("PvalueDistribution_", tag, ".png" ), path = pngDir)

}
