##### Global import ####

#' @importFrom ggplot2 ggplot geom_col theme_classic aes
#' theme element_text element_blank ylab xlab ggtitle
#' scale_fill_gradientn geom_tile theme_bw guides scale_fill_gradient2
#' guide_colourbar labs

# @export
#' @importFrom magrittr "%>%" 
magrittr::`%>%`

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
#' @noRd
#' @importFrom tibble tibble
#' @importFrom tidyr separate
#' @importFrom dplyr mutate_all
#'
#'
GetDesignFromNames <- function(samples_name){
  
  # Get the number of design factor and the factors from the names of the count matrix
  nb_dFac <- stringr::str_count(samples_name, pattern = "_") + 1
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
#' From a vector of character giving the name of the factors of an omics experiment,
#' and their type of effect: biological or batch, it returns all models formulae
#' that can be formulated in association with this factors. Batch effect factors do
#' not appear in interaction terms with biological factor. Model formulae stop in
#' second order interaction.
#'
#' @param FacBio a vector of character giving the name of the bio factors.
#' @param FacBatch a vector of character giving the name of the batch factors.
#' @param MAE a MultiAssayExperiment produced by RFLOMICS. Default is null. If not null, overwrite FacBio and FacBatch.
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



#' @title TMM.Normalization
#' Interface to the calcNormFactors function of the edgeR package  with the choosen TMM parameters as the normalization method
#' @param counts numeric matrix of read counts
#' @param groups vector or factor giving the experimental group/condition for each sample/library.
#' @return a data.frame with a row for each sample and columns group, lib.size and norm.factors containing the group labels, library sizes and normalization factors. Other columns can be optionally added to give more detailed sample information.
#' @export
#' @importFrom edgeR DGEList calcNormFactors
#' @noRd

TMM.Normalization <- function(counts, groups){
  dge <- edgeR::DGEList(counts=counts, group=groups)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  nf  <- dge$samples
  return(nf)
}


#' @title edgeR.AnaDiff
#'
#' @param object an object of class \link{SummarizedExperiment}]
#' @param design an object of class \link{ExpDesign-class}
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @param parallel boolean. Compute parallel differential analyses (only when clustermq is FALSE)
#' @param nworkers integer. Number of core to use for the parallel operations. Only used when parallel is TRUE.
#' @return A list of object of class \link{DGELRT}
#' @export
#' @importFrom stats model.matrix as.formula 
#' @importFrom edgeR DGEList estimateGLMCommonDisp estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit glmLRT topTags
#' @importFrom clustermq Q 
#' @importFrom parallel mclapply
#' @importFrom dplyr filter rename
#' @noRd
#'

edgeR.AnaDiff <- function(count_matrix, model_matrix, group, lib.size, norm.factors, Contrasts.Sel, Contrasts.Coeff, FDR, clustermq = FALSE,
                          parallel = FALSE, nworkers = 1, cmd = FALSE){
  
  z <- y <- NULL
  
  ListRes <- list()
  
  # check clustermq and parallel
  if(clustermq && parallel) parallel <- FALSE 
  if(!parallel) nworkers <- 1
  
  # model_matrix <- model_matrix[match(colnames(count_matrix), rownames(model_matrix)),]
  # if (!identical(colnames(count_matrix), rownames(model_matrix), attrib.as.set = FALSE)) 
  #   stop("MisMatch in samples names orders.")
  
  # Construct the DGE obect
  dge <- edgeR::DGEList(counts       = count_matrix,
                        group        = group,
                        lib.size     = lib.size,
                        norm.factors = norm.factors)
  
  # Run the model
  if(cmd) print("[cmd] dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)")
  dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)
  if(cmd) print("[cmd] dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)")
  dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)
  if(cmd) print("[cmd] dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)")
  dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)
  if(cmd) print("[cmd] fit.f <- edgeR::glmFit(dge,design=model_matrix)")
  fit.f <- edgeR::glmFit(dge,design=model_matrix)
  
  
  # test clustermq
  if(clustermq == TRUE){
    
    # Fonction to run on contrast per job
    # y is the model, Contrasts are stored in a matrix, by columns
    
    fx <- function(x){
      
      try_rflomics <- function(expr) {
        warn <- err <- NULL
        value <- withCallingHandlers(
          tryCatch(expr,
                   error    =function(e){ err <- e
                   NULL
                   }),
          warning =function(w){ warn <- w
          invokeRestart("muffleWarning")}
        )
        list(value=value, warning=warn, error=err)
      }
      
      try_rflomics(edgeR::glmLRT(y, contrast = unlist(z[x,])))
    }
    
    ResGlm <- clustermq::Q(fx, x=1:length(Contrasts.Sel$contrast),
                           export=list(y=fit.f,z=Contrasts.Coeff),
                           n_jobs=length(Contrasts.Sel$contrast),pkgs="edgeR")
    
  }
  else{
    if(cmd) print("[cmd] apply model to each contrast")
    # ResGlm <-  BiocParallel::bplapply(Contrasts.Sel$contrast, function(x){
    #   #print(unlist(Contrasts.Coeff[x,]))
    #   try_rflomics(edgeR::glmLRT(fit.f, contrast = unlist(Contrasts.Coeff[x,])))
    #   
    # }, BPOPTIONS = bpoptions(workers = nworkers))
    ResGlm <-  parallel::mclapply(Contrasts.Sel$contrast, function(x){
      #print(unlist(Contrasts.Coeff[x,]))
      try_rflomics(edgeR::glmLRT(fit.f, contrast = unlist(Contrasts.Coeff[x,])))
      
    }, mc.cores = nworkers)
    
  }
  
  # Create a table of jobs summary
  error.list <- unlist(lapply(ResGlm, function(x){
    ifelse(is.null(x$error),"success",as.character(x$error))
  }))
  
  jobs.tab <- data.frame(H=Contrasts.Sel$contrast , error.message=as.factor(error.list))
  
  jobs.tab.error <- jobs.tab %>% dplyr::filter(., error.message != "success")
  
  # If no error
  if(dim(jobs.tab.error)[1]==0){
    
    ListRes[[1]] <- lapply(ResGlm,function(x){
      x$value
    })
    
    # Name the table of raw results
    names(ListRes[[1]]) <- Contrasts.Sel$contrastName
    
    # ListRes[[2]] => TOPDGE => TopDFE
    
    TopDGE <- lapply(ListRes[[1]], function(x){
      
      res <- edgeR::topTags(x, n = dim(x)[1])
      
      DEGs<- res$table[res$table$FDR <= FDR,]
      #DEGs<-res$table
      return(DEGs)
    })
    
    ListRes[[2]] <- TopDGE
    
    names(ListRes[[2]]) <- names(ListRes[[1]])
    
    # Mutate column name to render the anadiff results generic
    # Initial column Name:  gene_name  logFC      logCPM        LR        PValue           FDR
    ListRes[[2]] <- lapply(ListRes[[2]], function(x){
      #dplyr::rename(x,"Abundance"="logCPM","StatTest"="LR","pvalue"="PValue","Adj.pvalue"="FDR")
      dplyr::rename(x,"Abundance"="logCPM","pvalue"="PValue","Adj.pvalue"="FDR")
    })
    
    names(ListRes) <- c("RawDEFres","TopDEF")
  }
  else{
    ListRes[[1]] <- NULL
    ListRes[[2]] <- jobs.tab.error
    names(ListRes) <- c("RawDEFres","ErrorTab")
  }
  
  return(ListRes)
}


#' @title limma.AnaDiff
#'
#' @param object an object of class \link{SummarizedExperiment}
#' @param design an object of class \link{ExpDesign-class}
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @return A list
#' @export
#' @importFrom stats model.matrix as.formula
#' @importFrom dplyr filter rename
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @noRd
#'


limma.AnaDiff <- function(count_matrix, model_matrix, Contrasts.Sel, Contrasts.Coeff, 
                          Adj.pvalue.cutoff, Adj.pvalue.method,clustermq, cmd = FALSE){
  
  ListRes <- list()
  
  # Run the model
  fit <- limma::lmFit(count_matrix, model_matrix)
  
  # test clustermq
  if(clustermq == TRUE){
    
    fx <- function(x){
      
      try_rflomics <- function(expr) {
        warn <- err <- NULL
        value <- withCallingHandlers(
          tryCatch(expr,
                   error    =function(e){ err <- e
                   NULL
                   }),
          warning =function(w){ warn <- w
          invokeRestart("muffleWarning")}
        )
        list(value=value, warning=warn, error=err)
      }
      
      # Fonction to run on contrast per job
      # y is the model, Contrasts are stored in a matrix, by columns
      
      try_rflomics(limma::contrasts.fit(y, contrasts  = unlist(z[x,])))
    }
    
    ResGlm  <- clustermq::Q(fx, x=1:length(Contrasts.Sel$contrast),
                            export=list(y=fit,z=Contrasts.Coeff),
                            n_jobs=length(Contrasts.Sel$contrast),pkgs="edgeR")
    
  }
  else{
    if(cmd) print("[cmd] fit contrasts")
    ResGlm <-  lapply(Contrasts.Sel$contrast, function(x){
      #print(paste0(x," : ",as.vector(unlist(Contrasts.Coeff[x,]))))
      try_rflomics(limma::contrasts.fit(fit, contrasts  = as.vector(unlist(Contrasts.Coeff[x,]))))
    })
  }
  
  # Construct a table of jobs summary
  error.list <- unlist(lapply(ResGlm, function(x){
    ifelse(is.null(x$error),"success",as.character(x$error))
  }))
  
  jobs.tab <- data.frame(H=Contrasts.Sel$contrast , error.message=as.factor(error.list))
  
  jobs.tab.error <- jobs.tab %>% dplyr::filter(., error.message != "success")
  
  # If no error
  if(dim(jobs.tab.error)[1]==0){
    
    ListRes[[1]] <- lapply(ResGlm,function(x){
      x$value
    })
    
    # Name the table of raw results
    
    names(ListRes[[1]]) <- Contrasts.Sel$contrastName
    
    # ListRes[[2]] TopDPE with column names common to all AnaDiff function
    
    ListRes[[2]] <- lapply(ListRes[[1]], function(x){
      
      fit2 <- limma::eBayes(x, robust=TRUE)
      res <- limma::topTable(fit2, adjust.method = Adj.pvalue.method, number=Inf, sort.by="AveExpr")
      DEPs <- res[res$adj.P.Val <= Adj.pvalue.cutoff,]
      return(DEPs)
    })
    
    
    # Mutate column name to render the anadiff results generic
    # Initial column Name:  logFC  AveExpr         t      P.Value    adj.P.Val            B
    ListRes[[2]] <- lapply(ListRes[[2]], function(x){
      #dplyr::rename(x,"Abundance"="AveExpr","StatTest"="t","pvalue"="P.Value","Adj.pvalue"="adj.P.Val")
      dplyr::rename(x,"Abundance"="AveExpr","pvalue"="P.Value","Adj.pvalue"="adj.P.Val")
    })
    
    names(ListRes[[2]]) <- names(ListRes[[1]])
    
    names(ListRes) <- c("RawDEFres","TopDEF")
    
  }
  else{
    ListRes[[1]] <- NULL
    ListRes[[2]] <- jobs.tab.error
    names(ListRes) <- c("RawDEFres","ErrorTab")
  }
  
  
  return(ListRes)
}



#' @title colorPlot
#'
#' @param design design
#' @param ColData ColData from summarizedExperiment
#' @param condition Factor for the color of the plot. Default is samples, it takes all bio factors modalities to color the plot.
#'
#' @return a color palette
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom tidyr unite
#' @noRd

colorPlot <- function(design, ColData, condition="samples"){
  
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  
  if(condition == "samples"){
    
    # combine only bio fact
    groups <- design@List.Factors[design@Factors.Type == "Bio"] %>%
      as.data.frame() %>% tidyr::unite(col="groups", sep="_")
    list.cond <- factor(groups$groups)
  }
  else{
    
    list.cond <- design@List.Factors[[condition]]
  }
  
  len.cond  <- length(levels(list.cond))
  
  colors    <- getPalette(len.cond)
  
  col <- colors[levels(list.cond)]
  names(col) <- row.names(levels(list.cond))
  
  return(col)
}


#' plotDistr
#'
#' @param abundances matrix or dataframe of feature/gene abundances/counts
#' @export
#' @importFrom ggplot2 geom_density xlab
#' @importFrom reshape2 melt
#' @noRd
plotDistr <- function(abundances, dataType, transform_method){
  
  
  value <- samples <- NULL
  
  switch(dataType,
         "RNAseq" = {
           pseudo_counts <- log2(abundances+1) %>% reshape2::melt()
           colnames(pseudo_counts) <- c("features", "samples", "value")
           x_lab <-"log2(feature abundances)"
         },
         "proteomics"={
           pseudo_counts <- abundances %>% reshape2::melt()
           colnames(pseudo_counts) <- c("features", "samples", "value")
           x_lab <- switch (transform_method,
                            "log2"= "log2(feature abundances)",
                            "none"="feature abundances")
         },
         "metabolomics"={
           pseudo_counts <- abundances %>% reshape2::melt()
           colnames(pseudo_counts) <- c("features", "samples", "value")
           x_lab <- switch (transform_method,
                            "log2"= "log2(feature abundances)",
                            "none"="feature abundances")
         }
  )
  
  p <- ggplot2::ggplot(pseudo_counts) + geom_density(aes(value, color=samples) ) + xlab(x_lab) +
    theme(legend.position='none')
  print(p)
  
}


#' pvalue.plot
#'
#' @param data dataframe (ggplot2)
#' @param hypothesis the contrast, useful for plot title
#' @return plot
#' @export 
#' @importFrom ggplot2 geom_histogram theme_bw labs ggplot
#' @noRd
pvalue.plot <- function(data, hypothesis=hypothesis){
  
  PValue <- NULL
  
  p <- ggplot2::ggplot(data=data) +
    ggplot2::geom_histogram(aes(x=pvalue), bins = 200) +
    ggplot2::labs(x = expression(p - value), y = "count", title = hypothesis )+
    ggplot2::theme_bw(base_size = 10)
  
  return(p)
}

# TODO do we use this?
# utils::globalVariables(names(data))

#' MA.plot
#'
#' @param data dataframe (ggplot2)
#' @param Adj.pvalue.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff log2FC cutoff (absolute value)
#' @param hypothesis the contrast, useful for plot title
#' @return MA plot
#' @export
#' @importFrom ggplot2 aes geom_point scale_colour_manual ggsave theme_linedraw
#' @importFrom ggpubr ggmaplot
#' @importFrom dplyr select rename
#' @noRd
MA.plot <- function(data, Adj.pvalue.cutoff, logFC.cutoff, hypothesis=hypothesis){

  Abundance <- logFC <- Adj.pvalue <- NULL
  tmp <-dplyr::select(data,"Abundance","logFC","Adj.pvalue") %>% 
    dplyr::rename(., baseMeanLog2=Abundance, log2FoldChange=logFC, padj=Adj.pvalue)
  p <- ggpubr::ggmaplot(tmp, main = hypothesis,
                        fdr = Adj.pvalue.cutoff, fc = 2^logFC.cutoff, size = 0.4,
                        ylab = bquote(~Log[2] ~ "fold change"),
                        xlab = bquote(~Log[2] ~ "mean expression"),
                        palette = c("#B31B21", "#1465AC", "grey30"),
                        select.top.method=c("padj","fc"),
                        legend = "bottom", top = 20,
                        font.label = c("plain", 7),
                        font.legend = c(11, "plain", "black"),
                        font.main = c(11, "bold", "black"),
                        caption = paste("logFC cutoff=",logFC.cutoff, " and " ,"FDR cutoff=",Adj.pvalue.cutoff,sep=""),
                        ggtheme = ggplot2::theme_linedraw())
  
  
  return(p)
  
}




#' Title
#'
#' @param data dataframe (ggplot2)
#' @param Adj.pvalue.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff log2FC cutoff (absolute value)
#' @param hypothesis the contrast, useful for plot title
#'
#' @return a volcano plot, made with the \link{EnhancedVolcano} package.
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @export
#' @noRd
#'
Volcano.plot <- function(data, Adj.pvalue.cutoff, logFC.cutoff, hypothesis){
  
  # Modified 221123
  # Find pvalue corresponding to the FDR cutoff for the plot (mean between the last that passes the cutoff
  # and the first that is rejected to plot the line in the middle of the two points)
  # if pvalcutoff is 1 (no cutoff) no need to adjust
  
  
  if(Adj.pvalue.cutoff > 1){stop("Adj.pvalue.cutoff must be between 0 and 1")}
  
  pval1 <- data$pvalue[data$Adj.pvalue<Adj.pvalue.cutoff] %>% dplyr::last()
  pval2 <- data$pvalue[data$Adj.pvalue>Adj.pvalue.cutoff] %>% dplyr::first()
  pvalCutoff <- (pval1 + pval2)/2
  
  
  # If too low pvalues, unable to plot (error in if(d>0)...)
  # If drawconnectors is FALSE, it "works", with ylim being infinity, it doesn't look like anything. 
  # Modifiying the 0 pvalues to make sure it's working
  nz_pval <- data$pvalue[data$pvalue != 0][1] * 10^-1 # default replacement in EnhancedVolcanoPlot
  if (nz_pval == 0){
    data$pvalue[data$pvalue == 0] <- data$pvalue[data$pvalue != 0][1] 
    # message("10^-1 * current lowest non-zero p-value is still 0, all 0 pvalues are set to the lowest non-zero pvalue.")
  }
  
  Abundance <- logFC <- Adj.pvalue <- NULL
  p <- EnhancedVolcano::EnhancedVolcano(toptable = data,
                                        lab = rownames(data),
                                        x = 'logFC',
                                        y = 'pvalue',
                                        # pCutoff = Adj.pvalue.cutoff,
                                        pCutoff = pvalCutoff,
                                        FCcutoff = logFC.cutoff,
                                        axisLabSize=10,
                                        pointSize = 1.5,
                                        labSize = 2,
                                        title = hypothesis,
                                        titleLabSize=11,
                                        subtitle = "",
                                        subtitleLabSize = 10,
                                        caption = paste("logFC cutoff=", logFC.cutoff, " & " ,"FDR cutoff=", Adj.pvalue.cutoff, sep=""),
                                        legendPosition = "bottom",
                                        legendLabSize = 10,
                                        legendIconSize=1.5,
                                        captionLabSize=10,
                                        col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                                        colAlpha = 0.5,
                                        drawConnectors = TRUE,
                                        widthConnectors = 0.5)
  
  return(p)
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

######################################## get contrast ##################################################

# it is possible to define contrast combinations that are specifically suited to a particular experimental design and set of research questions

#Contrasts are used to test whether the levels of an effect are significantly different from one another. You
#can specify a contrast for each factor in the model. Contrasts represent linear combinations of the
#parameters.

# Comparison (or contrast) procedures are used to test more specific hypotheses about differences between means
# comparisons of the two drug groups to the control group (  a complex comparison)
# The contrasts are formed by applying a set of weights, called contrast coefficients, to the means
# A contrast is a test of the difference between the means of two groups from the ANOVA. There are two categories of contrasts among the groups tested by ANOVA, simple and complex.
# A simple contrast is a test of the difference between any two pairs, such as Experimental Group 1 and Control Group 2. A complex contrast is a test of the difference between combinations of groups.
# An example of a complex contrast is a test of the difference between a subgroup created by combining Experimental Groups 1, 2 and 4 combined, and a subgroup created by combining Control Groups 1 and 3

# pairwise comparisons between the groups
# biological_factors_list, interaction_in_model, model_matrix

## functions to define part of simple contrasts
#' define part of simple contrast data frame
#'
#' @param treatmentFactorsList list
#' @param i i
#' @param j j
#'
#' @return a dataframe of contrasts
#' @export
#' @importFrom utils combn
#' @importFrom data.table setDT setnames
#' @importFrom tidyr unite
#' @importFrom dplyr mutate
#' @importFrom tidyselect all_of
#' @noRd
#' @author Christine Paysant-Le Roux
define_partOfSimpleContrast_df <- function (treatmentFactorsList, i, j) {
  
  contrastPart <- fixFactor <- NULL
  
  comparisonPart <- treatmentFactorsList
  # combn(x,2) generate all combinations of the elements of x taken 2 at a time
  vectorFromCombn <- combn(treatmentFactorsList[[i]],2)[j,]
  comparisonPart[[i]] <- vectorFromCombn
  df_comparisonPart <- expand.grid(comparisonPart)
  data.table::setDT(df_comparisonPart)
  
  # paste all the column of the data table
  #df_comparisonPart[, contrastPart := do.call(paste, c(.SD, sep = "_")), .SDcols = names(df_comparisonPart)]
  df_comparisonPart <- df_comparisonPart %>% tidyr::unite(contrastPart, sep="_", remove=F)
  
  #colnameFactor_i <- names(df_comparisonPart)[i]
  colnameFactor_i <- names(treatmentFactorsList)[i]
  
  #df_comparisonPart[, comparisonPart := df_comparisonPart[[colnameFactor_i]]]
  df_comparisonPart <- df_comparisonPart %>% dplyr::mutate(comparisonPart = df_comparisonPart[[colnameFactor_i]])
  
  if(length(names(treatmentFactorsList)) != 1){
    colnamesToKeep <- setdiff(names(df_comparisonPart),c("contrastPart", "comparisonPart", colnameFactor_i))
    #df_comparisonPart[, fixFactor := do.call(paste, c(.SD, sep = "_")), .SDcols = colnamesToKeep]
    df_comparisonPart <- df_comparisonPart %>% tidyr::unite(fixFactor, all_of(colnamesToKeep), sep = "_", remove = F)
  }else{
    df_comparisonPart <- df_comparisonPart %>% dplyr::mutate(fixFactor= NA)
  }
  
  
  colnamesToDelete <- names(treatmentFactorsList)
  #df_comparisonPart[, (colnamesToDelete) := NULL]
  df_comparisonPart <- df_comparisonPart %>% dplyr::select(-all_of(colnamesToDelete))
  
  nameColumnContrast <- paste0("contrastPart", j)
  nameColumnComparison <- paste0("comparisonPart", j)
  data.table::setnames(df_comparisonPart, c("contrastPart", "comparisonPart"), c(nameColumnContrast, nameColumnComparison))
  return(df_comparisonPart)
}
#' compute a data table with all pairwise comparisons of one factor
#'
#' @param treatmentFactorsList list
#' @param i i
#' @noRd
#' @return a data frame with all simple contrasts
#' @export
#' @importFrom dplyr all_of
#' @author Christine Paysant-Le Roux
simpleContrastForOneFactor <- function (treatmentFactorsList, i){
  
  fixFactor <- groupComparison <- NULL
  
  contrastPart1 <- contrastPart2 <- contrastPart3 <- contrastPart4 <-  NULL
  comparisonPart1 <- comparisonPart2 <- NULL
  fixPart1 <- fixPart3 <- fixFactor1 <- fixFactor3 <- NULL
  comparisonPart3 <- comparisonPart4 <- NULL
  
  df_FirstComparisonPart <- define_partOfSimpleContrast_df(treatmentFactorsList,i,2)
  #df_FirstComparisonPart[,fixFactor := NULL]
  df_FirstComparisonPart <- df_FirstComparisonPart %>% dplyr::select(-fixFactor)
  
  df_SecondComparisonPart <- define_partOfSimpleContrast_df(treatmentFactorsList,i,1)
  df_simpleContrasts_factor <- cbind(df_FirstComparisonPart, df_SecondComparisonPart)
  
  #df_simpleContrasts_factor[, contrast := paste0("(", contrastPart2, " - ", contrastPart1, ")")]
  #df_simpleContrasts_factor[, groupComparison := paste0("(", comparisonPart2, " - ", comparisonPart1, ")")]
  
  df_simpleContrasts_factor <- df_simpleContrasts_factor %>%
    dplyr::mutate(contrast        = paste0("(", contrastPart2,   " - ", contrastPart1, ")"),
                  groupComparison = paste0("(", comparisonPart2, " - ", comparisonPart1, ")"))
  
  
  # case where fixFactor column is empty (NA inside)
  emptycolFixFactor <- unique(is.na(df_simpleContrasts_factor$fixFactor))
  if(emptycolFixFactor){
    #df_simpleContrasts_factor[, contrastName := groupComparison]
    df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::mutate(contrastName = groupComparison)
    
  } else {
    #df_simpleContrasts_factor[, contrastName := paste0(groupComparison, " in ", fixFactor )]
    df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::mutate(contrastName = paste0(groupComparison, " in ", fixFactor ))
  }
  #df_simpleContrasts_factor[, type := "simple"]
  df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::mutate(type = "simple")
  
  colnamesToDelete <- c("contrastPart2", "comparisonPart2", "contrastPart1", "comparisonPart1")
  
  #df_simpleContrasts_factor[, (colnamesToDelete) := NULL]
  #data.table::setcolorder(df_simpleContrasts_factor, c(names(df_simpleContrasts_factor)[2:length(names(df_simpleContrasts_factor))], names(df_simpleContrasts_factor)[1]))
  df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::select(-all_of(colnamesToDelete)) %>%
    dplyr::select("contrast", "groupComparison", "contrastName", "type", "fixFactor")
  
  
  
  return(df_simpleContrasts_factor)
}

#' define all simple contrasts
#'
#' @param treatmentFactorsList list of treatment factors
#'
#' @return a data frame with all simple contrasts
#' @export
#' @noRd
#' @author Christine Paysant-Le Roux
defineAllSimpleContrasts <- function(treatmentFactorsList){
  # create a data table with 5 columns
  # empty data table
  allSimpleContrast_df <- data.table::data.table(contrast = character(), groupComparison = factor(), contrastName = character(), type = character(), fixFactor = factor())
  # create each data frame and rbind it to the allSimpleContrast_df
  for(i in seq_along(treatmentFactorsList)){
    dataTableToCreate <- simpleContrastForOneFactor(treatmentFactorsList, i)
    allSimpleContrast_df <- rbind(allSimpleContrast_df, dataTableToCreate)
  }
  #allSimpleContrast_df[,contrastCoeff := sapply(contrast, function(x) defineCoefficient(x, colnamesGLMdesign))]
  return(allSimpleContrast_df[])
}

#' Define averaged contrasts
#'
#' @param allSimpleContrast_df : a data frame with all the simple contrasts comparisons (test of the difference between any two pairs of groups)
#'
#' @return allAveragedContrasts_df : a data frame with all the averaged contrasts
#' @export
#' @noRd
#' @importFrom dplyr add_tally group_by mutate select 
#' @importFrom data.table data.table
#' @author Christine Paysant-Le Roux
define_averaged_contrasts <- function(allSimpleContrast_df){
  
  groupComparison <- contrast <- n <- fixFactor <- contrastName <- NULL
  
  #allAveragedContrasts_df <- allSimpleContrast_df[, list(contrast = paste0(paste0(paste0("(", paste(contrast, collapse=" + ")),")/"),.N), meanIn = paste(fixFactor, collapse=" + ")), by = groupComparison]
  #allAveragedContrasts_df[, type :="mean"]
  allAveragedContrasts_df <- allSimpleContrast_df %>% dplyr::group_by(groupComparison) %>% dplyr::add_tally() %>%
    dplyr::mutate(contrast= paste0(paste0("(", paste(contrast, collapse=" + ")),")/", n),
                  meanIn  = paste(fixFactor, collapse=" + "),
                  type    = "mean") %>%
    dplyr::select(-contrastName, -fixFactor, -n) %>% unique() %>% data.table::data.table()
  
  #allAveragedContrasts_df[, contrastName := paste(groupComparison, "mean", sep = " in ")]
  allAveragedContrasts_df <- allAveragedContrasts_df %>% dplyr::mutate(contrastName = paste(groupComparison, "mean", sep = " in "))
  
  data.table::setcolorder(allAveragedContrasts_df, c("contrast", "groupComparison", "contrastName", "type", "meanIn"))
  
  #  allAveragedContrasts_df <- allSimpleContrast_df[, list(meanIn = paste(fixFactor, collapse=" + ")), by = groupComparison]
  return(allAveragedContrasts_df[])
}
# define contrasts for interactions

#' define a data frame with part of interaction contrast
#'
#' @param treatmentFactorsList list
#' @param i i 
#' @param j j
#' @param k k
#' @param row_i row_i
#' @param row_j row_j
#'
#' @return a dataframe with part of the interaction contrasts definition
#' @export
#' @importFrom data.table setDT
#' @noRd
#' @author Christine Paysant-Le Roux
define_partOfInteractionContrast_df <- function (treatmentFactorsList, i, j, k, row_i, row_j) {
  
  contrastPart_bis <- outsideGroup_bis <- fixFactor_bis <- NULL
  
  
  comparisonPart <- treatmentFactorsList
  # combn(x,2) generate all combinations of the elements of x taken 2 at a time
  comparisonPart [[i]] <- combn(treatmentFactorsList[[i]],2)[row_i,]
  comparisonPart [[j]] <- combn(treatmentFactorsList[[j]],2)[row_j,]
  df_comparisonPart <- expand.grid(comparisonPart)
  data.table::setDT(df_comparisonPart)
  # paste all the column of the data table
  #df_comparisonPart[, contrastPart := do.call(paste, c(.SD, sep = "_")), .SDcols = names(df_comparisonPart)]
  df_comparisonPart <- df_comparisonPart %>% tidyr::unite(contrastPart_bis, names(df_comparisonPart), sep="_", remove=F) %>%
    dplyr::mutate(contrastPart = contrastPart_bis) %>% dplyr::select(-contrastPart_bis)
  
  colnameFactor_i <- names(df_comparisonPart)[i]
  colnameFactor_j <- names(df_comparisonPart)[j]
  
  #data.table::setDT(df_comparisonPart)
  #df_comparisonPart[, comparisonPart := df_comparisonPart[[colnameFactor_i]]]
  #df_comparisonPart[, fixPart := df_comparisonPart[[colnameFactor_j]]]
  df_comparisonPart <- df_comparisonPart %>%
    dplyr::mutate(comparisonPart = df_comparisonPart[[colnameFactor_i]]) %>%
    dplyr::mutate(fixPart        = df_comparisonPart[[colnameFactor_j]])
  
  
  colnamesToKeep <- setdiff(names(treatmentFactorsList),c(colnameFactor_i, colnameFactor_j))
  #df_comparisonPart[, outsideGroup := do.call(paste, c(.SD, sep = "_")), .SDcols = colnamesToKeep]
  if(length(colnamesToKeep)){
    
    df_comparisonPart <- df_comparisonPart %>% tidyr::unite(outsideGroup_bis, all_of(colnamesToKeep), sep="_", remove=F) %>%
      dplyr::mutate(outsideGroup = outsideGroup_bis) %>% dplyr::select(-outsideGroup_bis)
  }else{
    
    df_comparisonPart <- df_comparisonPart %>% dplyr::mutate(outsideGroup = NA)
  }
  
  colnamesToKeep <- setdiff(names(df_comparisonPart),c("contrastPart", "comparisonPart", "fixPart", "outsideGroup", colnameFactor_i))
  #df_comparisonPart[, fixFactor := do.call(paste, c(.SD, sep = "_")), .SDcols = colnamesToKeep]
  df_comparisonPart <- df_comparisonPart %>% tidyr::unite(fixFactor_bis, all_of(colnamesToKeep), sep="_", remove=F) %>%
    dplyr::mutate(fixFactor = fixFactor_bis) %>% dplyr::select(-fixFactor_bis)
  
  colnamesToDelete <- names(treatmentFactorsList)
  #df_comparisonPart[, (colnamesToDelete) := NULL]
  df_comparisonPart <- df_comparisonPart %>% dplyr::select(-all_of(colnamesToDelete))
  
  
  nameColumnContrast <- paste0("contrastPart", k)
  nameColumnComparison <- paste0("comparisonPart", k)
  nameFixFactor <- paste0("fixFactor", k)
  namePartFixFactor <- paste0("fixPart", k)
  nameOutsideGroup <- paste0("outsideGroup", k)
  data.table::setnames(df_comparisonPart, c("contrastPart", "comparisonPart", "fixFactor", "fixPart", "outsideGroup"),
                       c(nameColumnContrast, nameColumnComparison, nameFixFactor, namePartFixFactor, nameOutsideGroup))
  return(df_comparisonPart)
}

#' define interaction constrast for pairs of biological factors
#'
#' @param treatmentFactorsList list
#' @param i i
#' @param j j
#'
#' @return a dataframe
#' @export
#' @noRd
#' @author Christine Paysant-Le Roux
defineInteractionConstrastForPairsOfFactors <- function(treatmentFactorsList, i, j){
  
  contrastPart1 <- contrastPart2 <- contrastPart3 <- contrastPart4 <- NULL
  comparisonPart1 <- comparisonPart2 <- comparisonPart3 <- comparisonPart4 <- NULL
  fixPart1 <- fixPart3 <- fixFactor1 <- fixFactor3 <- NULL
  
  df_part1 <-  define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 1, 2, 2)
  df_part2 <-  define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 2, 1, 2)
  df_part3 <-  define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 3, 2, 1)
  df_part4 <-  define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 4, 1, 1)
  df_interactionContrasts <- cbind(df_part1, df_part2, df_part3, df_part4)
  #df_interactionContrasts[, contrast := paste0("(", "(", contrastPart1, " - ", contrastPart2, ")"," - ", "(", contrastPart3, " - ", contrastPart4, ")", ")")]
  #df_interactionContrasts[, groupComparison := paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " vs ", "(", fixPart1, " - ", fixPart3, ")")]
  #df_interactionContrasts[, contrastName := paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " in ", fixFactor1, " - ", "(", comparisonPart3, " - ", comparisonPart4, ")", " in ", fixFactor3 )]
  #df_interactionContrasts[, type := "interaction"]
  df_interactionContrasts <- df_interactionContrasts %>% dplyr::mutate(contrast = paste0("(", "(", contrastPart1, " - ", contrastPart2, ")"," - ",
                                                                                         "(", contrastPart3, " - ", contrastPart4, ")", ")"),
                                                                       groupComparison = paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " vs ",
                                                                                                "(", fixPart1, " - ", fixPart3, ")"),
                                                                       contrastName  = paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " in ", fixFactor1, " - ",
                                                                                              "(", comparisonPart3, " - ", comparisonPart4, ")", " in ", fixFactor3 ),
                                                                       type = "interaction")
  
  colnamesToDelete <- c("contrastPart1",  "comparisonPart1", "fixFactor1", "fixPart1", "outsideGroup1",
                        "contrastPart2", "comparisonPart2", "fixFactor2", "fixPart2", "outsideGroup2",
                        "contrastPart3", "comparisonPart3", "fixFactor3", "fixPart3", "outsideGroup3",
                        "contrastPart4", "comparisonPart4", "fixFactor4", "fixPart4")
  #df_interactionContrasts[, (colnamesToDelete) := NULL]
  df_interactionContrasts <- df_interactionContrasts %>% dplyr::select(-tidyselect::all_of(colnamesToDelete))
  
  
  data.table::setnames(df_interactionContrasts, "outsideGroup4", "outsideGroup")
  
  #df_interactionContrasts[,groupInteraction := paste0(names(treatmentFactorsList)[i], " vs ", names(treatmentFactorsList)[j])]
  df_interactionContrasts <- df_interactionContrasts %>% dplyr::mutate(groupInteraction = paste0(names(treatmentFactorsList)[i], " vs ", names(treatmentFactorsList)[j]))
  
  # if 3 factors bio / outsideGroup exist
  if(!is.null(df_interactionContrasts$outsideGroup)){
    
    ## add factor name of outsideGroup modality
    if(length(names(treatmentFactorsList)[-c(i,j)]) != 0){
      
      df_interactionContrasts <- df_interactionContrasts %>% dplyr::mutate(outsideGroup = names(treatmentFactorsList)[-c(i,j)]) %>% 
        dplyr::group_by(outsideGroup, groupComparison) %>% dplyr::add_tally() %>% 
        dplyr::mutate(contrast= paste0(paste0("(", paste(contrast, collapse=" + ")),")/", n), 
                      contrastName=paste0(groupComparison, " in ", outsideGroup)) %>% 
        dplyr::select(-n) %>% unique()
    }
  }
  
  return(df_interactionContrasts)
}

#' define all interaction contrasts
#'
#' @param treatmentFactorsList list
#' @param groupInteractionToKeep tokeep
#'
#' @return a dataframe with all the interaction contrasts
#' @export
#' @noRd
#' @author Christine Paysant-Le Roux
defineAllInteractionContrasts <- function(treatmentFactorsList, groupInteractionToKeep = NULL){
  
  groupInteraction <- NULL
  
  allInteractionsContrasts_df <- data.table::data.table(contrast = character(), groupComparison = factor(), groupInteraction = character(),
                                                        outsideGroup = character(),contrastName = character(), type = character())
  # combn(names(treatmentFactorsList),2)
  # cat(paste("\ntreatment factors names:\n"))
  # print(as.character(names(treatmentFactorsList)))
  vecFori <- combn(length(treatmentFactorsList),2)[1,]
  vecForj <- combn(length(treatmentFactorsList),2)[2,]
  for(k in seq_along(vecFori)){
    i <- vecFori[k]
    j <- vecForj[k]
    #print(paste("i:",i))
    #print(paste("j:",j))
    dataTableToCreate <-  defineInteractionConstrastForPairsOfFactors(treatmentFactorsList, i, j)
    allInteractionsContrasts_df <- rbind(allInteractionsContrasts_df, dataTableToCreate)
  }
  if(!missing(groupInteractionToKeep)){
    allInteractionsContrasts_df <- subset(allInteractionsContrasts_df, (groupInteraction %in% groupInteractionToKeep))
  }
  return(allInteractionsContrasts_df)
}


################## function for getContrastMatrix ExpDesign method  ########################


#' compute group binary vector
#'
#' Compute a binary vector (a single vector of 0 and 1) returning the matched group string(s) from a grepl match on the design model matrix colnames.
#'
#' @param biologicalGroups: factor giving group membership (treatment condition associated to one sample)
#' @param colnamesMatrixDesign : vector giving the column names of the model design matrix
#' @param interactionPresent: logical. If TRUE interaction is include in the design model matrix
#'
#' @return a binary vector (a single vector of 0 and 1) returning the matched group string(s) from a grepl match on the design model matrix colnames
#' @export
#' @noRd
#' @author Christine Paysant-Le Roux
computeGroupVector <- function(treatmentGroups, colnamesMatrixDesign, interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction) {
  vectorLength <- length(colnamesMatrixDesign)
  groupVector <- rep(0, vectorLength)
  if(interactionPresent){
    simples <- unique(unlist(strsplit(treatmentGroups, "_")))
    order2interaction <- combn(simples,2, FUN=paste, collapse=':')
    toMatchList <- c(simples, order2interaction)
    if(isThreeOrderInteraction){
      order3interaction <- combn(simples,3, FUN=paste, collapse=':')
      toMatchList <- c(toMatchList, order3interaction)
    }
    toMatchList <- paste("^", toMatchList, "$", sep = "")
    #grepl return TRUE for each matched pattern
    pos <- grepl(paste(toMatchList, collapse = "|"), x= colnamesMatrixDesign)
    groupVector <- as.numeric(pos)
    # for intercept
    groupVector[1] <- 1
  } else {
    simples <- unlist(strsplit(treatmentGroups, "_"))
    toMatchList <- simples
    toMatchList <- paste("^", toMatchList, "$", sep = "")
    #toMatch <- paste("^", treatmentGroups, "$", sep = "")
    pos <- grepl(paste(toMatchList, collapse = "|"), x= colnamesMatrixDesign)
    groupVector <- as.numeric(pos)
    # for intercept
    groupVector[1] <- 1
  }
  return(groupVector)
}
#' Assign binary vector to groups
#'
#' @param treatmentFactorsList: list of treatment factors levels
#' @param interaction_in_model: logical. If TRUE interaction is include in the design model matrix
#' @param modelMatrix: numeric matrix giving the design matrix of the GLM.
#' @param treatmentCondenv: the environment to use
#'
#' @return binary vector
#' @export
#' @noRd
#' @author Christine Paysant-Le Roux
assignVectorToGroups <- function(treatmentFactorsList = treatmentFactorsList, modelMatrix = modelMatrix,  interactionPresent = interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction, treatmentCondenv = treatmentCondenv){
  # treatment conditions (group) compatible with colnames of the design model matrix
  treatmentGroups <- do.call(paste, c(expand.grid(treatmentFactorsList), sep = "_"))
  # assign binary vector to each group
  modelMatrixColnames <- colnames(modelMatrix)
  # treatmentCondenv <- new.env()
  binaryVectorsList <- lapply(treatmentGroups, function(x)  computeGroupVector(x, modelMatrixColnames, interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction))
  names(binaryVectorsList) <- treatmentGroups
  groupDF <- as.data.frame(binaryVectorsList, row.names = modelMatrixColnames)
  mapply(function(x, value) assign(x, value, pos = treatmentCondenv), treatmentGroups, binaryVectorsList)
  # http://r.789695.n4.nabble.com/Using-assign-with-mapply-td4681789.html
  # for (i in 1:n) assign(levels[i], indicatorForModelWithInteraction(levels[i], colnamesGLMdesign = colnamesGLMdesign), pos = levelsenv)
}

#' return contrast coefficients
#'
#' @param contrast contrast considered
#' @param colnamesGLMdesign colnames
#' @param treatmentCondenv: the environment to use
#' @noRd
#' @return the contrast vector
#' @export
#' @author Christine Paysant-Le Roux
returnContrastCoefficients <- function(contrast, colnamesGLMdesign, treatmentCondenv){
  expression <- NULL
  if (!is.null(contrast)) {
    expression <- as.character(contrast)
    contrastVector <-rep(0,length(colnamesGLMdesign))
    ej <- parse(text = expression[1])
    # eval evaluates the expression argument in the environment specified by envir and returns the computed value
    contrastVector <- eval(ej, envir = treatmentCondenv)
  }
  return(contrastVector)
}

####################################### CO-EXPRESSION ##############################

#' return gene list
#'
#' @param matrix matrix of DE results
#' @param colnames colnames
#' @param mergeType either union or intersection.
#' @return list of genes
#' @export
#' @noRd
#'
getDEGlist_for_coseqAnalysis <- function(matrix, colnames = colnames(matrix)[-1], mergeType="union"){
  
  if (length(colnames) == 0 ){ return(NULL) }
  
  matrix_sum <- matrix %>% dplyr::mutate(sum = dplyr::select(., tidyselect::all_of(colnames)) %>% rowSums(.))
  
  DEG_list <- switch(mergeType,
                     
                     "union"={        dplyr::filter(matrix_sum, sum != 0) },
                     "intersection"={ dplyr::filter(matrix_sum, sum == length(colnames)) }
  )
  
  if (length(DEG_list$DEF) == 0 ){ return(NULL) }
  
  return(DEG_list$DEF)
}



#' @title try_rflomics
#' @details
#' This function come from https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#  The author indicated that he merged Martins solution (https://stackoverflow.com/a/4952908/2161065) and
#  the one from the R-help mailing list you get with demo(error.catching).
#' @param expr The expression that has been to be evaluated.
#' @return a named list
#' \itemize{
#' \item{\code{value:} }{The results of the expr evaluation or NULL if an error occured }
#' \item{\code{warning:} }{warning message or NULL}
#' \item{\code{error:} }{error message or NULL}
#' }
#' @export
#' @keywords internal
#' @noRd

try_rflomics <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr,
             error    = function(e){ err <- e
             NULL
             }),
    warning = function(w){ warn <- w
    invokeRestart("muffleWarning")}
  )

  return(list(value = value,
              warning = warn,
              error = err)
  )

}


#' @title coseq.error.manage
#' @param coseq.res.list list of coseq object
#' @param K number of group
#' @param replicates number of replication to run
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @export
#' @importFrom dplyr n group_by summarise mutate filter
#' @importFrom coseq ICL
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @noRd
#'
coseq.error.manage <- function(coseq.res.list, K, replicates, cmd = FALSE){
  
  # Create a table of jobs summary
  error.list <- unlist(lapply(coseq.res.list, function(x){
    ifelse(is.null(x$error), "success", as.character(x$error))
  }))
  
  # status of jobs
  nK_success.job <- table(error.list)["success"]
  
  if (is.na(nK_success.job)) { nK_success.job <- 0 }
  
  # if at least one failed job
  # => generate table with error summary
  K.list <- rep(paste0("K",min(K), "-", max(K)), each = replicates)
  
  jobs.tab <- data.frame(K = K.list, error.message = as.factor(error.list))
  
  jobs.tab.sum1 <- jobs.tab %>% 
    dplyr::group_by(K, error.message) %>%
    dplyr::summarise(n = dplyr::n()) %>%  
    dplyr::mutate(prop.failed = round((n/replicates)*100)) %>%
    dplyr::filter(error.message != "success")
  
  jobs.tab.sum <- jobs.tab.sum1
  
  if (nK_success.job != 0) {
    
    # Generate the list of results
    coseq.res.list[["value"]] <- list()
    
    for (x in names(coseq.res.list)) {
      
      if (!is.null(coseq.res.list[[x]]$value)) {
        coseq.res.list[["value"]][[x]] <- coseq.res.list[[x]]$value
      }
    }
    
    if (cmd) print("#     => error management: level 2 ")
    ICL.vec <- unlist(lapply(1:nK_success.job, function(x){ (coseq::ICL(coseq.res.list[["value"]][[x]])) })) %>%
      lapply(., function(x){ ifelse(is.na(x), "failed", "success") }) %>% unlist()
    
    nK_success <- table(ICL.vec)["success"]
    
    replicates <- nK_success.job
    
    # expected list of cases
    K.list.ex <- rep(K, each = replicates)
    
    # observed list of cases
    K.list.ob <- stringr::str_replace(string = names(ICL.vec), pattern = "K=", replacement = "") %>% 
      as.numeric() %>% 
      sort()
    
    # missed cases
    if (length(K.list.ob) != length(K.list.ex)) {
      
      missed.K.vec <- names(table(K.list.ob)[table(K.list.ob) < nK_success.job])
      
      ICL.vec.bis <- rep("failed", length(missed.K.vec))
      names(ICL.vec.bis) <- paste0("K=", missed.K.vec)
      
      ICL.vec <- c(ICL.vec, ICL.vec.bis)
    }
    
    jobs.tab <- data.frame(K = names(ICL.vec), error.message = as.factor(ICL.vec))
    
    jobs.tab.sum2 <- jobs.tab %>% 
      dplyr::group_by(K, error.message) %>%
      dplyr::summarise(n = dplyr::n()) %>%  
      dplyr::mutate(prop.failed = round((n/replicates)*100)) %>%
      dplyr::filter(error.message != "success")
    
    # if (dim(jobs.tab.sum1)[1] == 0){ jobs.tab.sum <- jobs.tab.sum2 }
    # else if(dim(jobs.tab.sum2)[1] == 0){ jobs.tab.sum <- jobs.tab.sum1 }
    # else{ jobs.tab.sum <- rbind(jobs.tab.sum1, jobs.tab.sum2) }
    
    jobs.tab.sum <- data.table::rbindlist(list(jobs.tab.sum1, jobs.tab.sum2), use.names = TRUE) %>% 
      tibble::tibble()
    
  }
  else{
    nK_success <- 0
  }
  
  return(list(jobs.tab.sum = jobs.tab.sum, 
              nK_success = nK_success, 
              coseq.res.list.values = coseq.res.list[["value"]]))
}



#' @title coseq.results.process
#' @param coseqObjectList list of coseq object
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @importFrom coseq ICL likelihood clusters
#' @export
#' @noRd
#'
coseq.results.process <- function(coseqObjectList, K, conds){
  
  # ICL plot
  ICL.list <- list()
  
  ICL.vec <- lapply(1:length(coseqObjectList), function(x){ coseq::ICL(coseqObjectList[[x]]) }) %>% unlist()
  ICL.list[["ICL.vec"]] <- ICL.vec
  
  ICL.min.per.cluster <- lapply(1:length(coseqObjectList), function(x){
    ICL.vec.min <- coseq::ICL(coseqObjectList[[x]])
    ICL.vec.min[which(ICL.vec.min == min(ICL.vec.min))]
    })  %>% unlist()
  
  # TODO tests 231117
  # ICL.min.per.cluster <- lapply(K, function(x){
  #   # x = 2
  #   ICL.vec.min <- ICL.vec[names(ICL.vec)==paste0("K=", x)]
  #   # ICL.vec.min[]
  #   which(ICL.vec.min == min(ICL.vec.min))
  #   })  %>% unlist()
  
  ICL.tab <- data.frame(K = stringr::str_replace(names(ICL.vec), "K=", ""), ICL = ICL.vec) %>%
    dplyr::mutate(K = as.numeric(K))
  ICL.list[["ICL.tab"]] <- ICL.tab
  
  ICL.n <- ICL.tab  %>% 
    dplyr::group_by(.,K) %>% 
    dplyr::filter(!is.na(ICL)) %>%
    dplyr::summarise(median = median(ICL, na.rm = TRUE), n = dplyr::n()) %>%
    dplyr::mutate(K = as.numeric(K))
  ICL.list[["ICL.n"]] <- ICL.n

  # min ICL
  K.ICL.median.min <- ICL.n[which.min(ICL.n$median),]$K
  index  <- which(names(ICL.min.per.cluster) == paste0("K=", K.ICL.median.min))
  index2 <- which(ICL.min.per.cluster[index] == min(ICL.min.per.cluster[index]))
  
  # coseq object with the min ICL
    coseq.res <- coseqObjectList[index][index2][[1]]
  # coseq.res <- coseqObjectList[[index2]]@allResults[[index]]
  
  # K.ICL.min <- min(ICL.vec[names(ICL.vec) == paste0("K=", K.ICL.median.min)], na.rm = TRUE)
  # 
  # 
  # index <- sapply(names(coseqObjectList), function(x){
  #   
  #   any(coseq::ICL(coseqObjectList[[x]]) == K.ICL.min)
  #   })
  # coseq.res <- coseqObjectList[index][[1]]
  # coseq.res <- coseqObjectList[[which.min(ICL.vec)]]
  
  # logLike plot
  logLike.list <- list()
  
  logLike.vec <- lapply(1:length(coseqObjectList), function(x){ coseq::likelihood(coseqObjectList[[x]]) }) %>% unlist()
  logLike.list[["logLike.vec"]] <- logLike.vec
  
  logLike.tab <- data.frame(K = stringr::str_replace(names(logLike.vec), "K=", ""), logLike = logLike.vec) %>% 
    dplyr::mutate(K = as.numeric(K))
  logLike.list[["logLike.tab"]] <- logLike.tab
  
  logLike.n <- logLike.tab %>% 
    dplyr::group_by(.,K) %>% 
    dplyr::filter(!is.na(logLike)) %>%
    dplyr::summarise(median = median(logLike), n = dplyr::n()) %>%
    dplyr::mutate(K = as.numeric(K))
  logLike.list[["logLike.n"]] <- logLike.n
  
  # list of genes per cluster
  clusters <- lapply(1:length(table(coseq::clusters(coseq.res))), function(i){
    names(coseq::clusters(coseq.res)[coseq::clusters(coseq.res) == i])
  })
  names(clusters) <- paste("cluster", 1:length(table(coseq::clusters(coseq.res))), sep = ".")
  
  # nbr of cluster
  # Gestion des NA dans les ICLs
  ICLv <- na.omit(coseq.res@metadata$ICL)
  nb_cluster <- na.omit(coseq.res@metadata$nbCluster)[min(ICLv) == ICLv]
  
  #output
  CoExpAnal <- list()
  CoExpAnal[["results"]]      <- TRUE
  CoExpAnal[["coseqResults"]] <- coseq.res
  CoExpAnal[["clusters"]]     <- clusters
  CoExpAnal[["cluster.nb"]]   <- nb_cluster
  CoExpAnal[["plots"]]        <- list("ICL" = ICL.list, "logLike" = logLike.list)
  
  #### Plots
  
  #### plot ICL
  # ICL.p   <- ggplot2::ggplot(data = ICL.tab) + 
  #   ggplot2::geom_boxplot(ggplot2::aes(x = as.factor(K), y = ICL, group = K)) +
  #   ggplot2::geom_text(data = ICL.n, ggplot2::aes(x = 1:length(K), y = max(ICL.vec, na.rm = TRUE), 
  #                                                 label = paste0("n=", n)), col = 'red', size = 4) +
  #   ggplot2::ylim(min(ICL.vec, na.rm = TRUE), max(ICL.vec, na.rm = TRUE)) +
  #   ggplot2::xlab("K")
  
  #### plot logLike
  # logLike.p   <- ggplot2::ggplot(data = logLike.tab) + 
  #   ggplot2::geom_boxplot(ggplot2::aes(x = as.factor(K), y = logLike, group = K)) + 
  #   ggplot2::xlab("K") +
  #   ggplot2::geom_text(data = logLike.n, ggplot2::aes(x = 1:length(K), y = max(logLike.vec, na.rm = TRUE), 
  #                                                     label = paste0("n=", n)), col = 'red', size = 4)
  
  #### coseq plots
  # plot.coseq.res <- coseq::plot(coseq.res, conds = conds, collapse_reps = "average",
  #                               graphs = c("profiles", "boxplots", "probapost_boxplots",
  #                                          "probapost_barplots", "probapost_histogram")) 
  
  # CoExpAnal[["plots"]] <- plot.coseq.res
  # CoExpAnal[["plots"]][["ICL"]]     <- ICL.p
  # CoExpAnal[["plots"]][["logLike"]] <- logLike.p
  
  return(CoExpAnal)
}


#' @title run Coseq for co-expression analysis on cluster
#' @param counts matrix
#' @param K number of groups
#' @param replicates number of replication to run
#' @param param.list list of coseq parameters
#' @return coseqResults
#' @importFrom coseq coseq
#' @importFrom clustermq Q_rows 
#' @export
#' @noRd
#'
runCoseq_clustermq <- function(counts, conds, K=2:20, replicates = 5, param.list, silent = TRUE, cmd = FALSE){

  # iter <-  rep(K, each = replicates)
  # seed_arg = rep(1:replicates, max(K) - 1)  
  
  df_args <- data.frame(x = rep(K, each = replicates),
                        seed_arg = rep(1:replicates, max(K) - 1) )
  # colnames(df_args) = c("x", "seed_arg")
  
  # nbr_iter <- length(iter)
  nbr_iter <- nrow(df_args)
  coseq.res.list <- list()
  set.seed(12345)
  
  # setting to run coseq on clustermq
  param.list[["object"]] <- counts
  param.list[["K"]] <- K
  
  fx <- function(x, seed_arg){
    
    try_rflomics <- function(expr) {
      warn <- err <- NULL
      value <- withCallingHandlers(
        tryCatch(expr,
                 error    = function(e){ err <- e
                 NULL
                 }),
        warning = function(w){ warn <- w
        invokeRestart("muffleWarning")}
      )
      list(value = value, warning = warn, error = err)
    }
    if (silent) { 
      co <- suppressMessages(capture.output(
        res <- try_rflomics(coseq::coseq(object = param.list[["object"]], 
                                         K = x,
                                         model = param.list$model,
                                         transformation = param.list$transformation,
                                         GaussianModel = param.list$GaussianModel,
                                         normFactors = param.list$normFactors,
                                         meanFilterCutoff = param.list$meanFilterCutoff,
                                         seed = seed_arg, 
                                         verbose = FALSE))
      ))
    }else{
      res <- try_rflomics(coseq::coseq(object = param.list[["object"]], 
                                       K = x,
                                       model = param.list$model,
                                       transformation = param.list$transformation,
                                       GaussianModel = param.list$GaussianModel,
                                       normFactors = param.list$normFactors,
                                       meanFilterCutoff = param.list$meanFilterCutoff,
                                       seed = seed_arg))
    }
    return(res)
  }
  
  # coseq.res.list <- clustermq::Q(fx,
  #                                x = iter, seed_arg = rep(1:replicates, replicates), 
  #                                export = param.list, n_jobs = nbr_iter, pkgs = "coseq")
  # 
  coseq.res.list <- clustermq::Q_rows(fun = fx,
                                      df = df_args, 
                                      export = param.list, n_jobs = nbr_iter, pkgs = "coseq")

  names(coseq.res.list) <- c(1:nbr_iter)
  
  CoExpAnal <- list()
  
  if (cmd) print("#     => error management ")
  
  # Create a table of jobs summary
  error.list <- unlist(lapply(coseq.res.list, function(x){
    ifelse(is.null(x$error), "success", as.character(x$error))
  }))
  
  nK_success <- table(error.list)["success"]
  if (cmd) print(paste0("#     => nbr of success jobs: ", nK_success))
  
  K.list <- rep(paste0("K=", K), each = replicates)
  
  jobs.tab <- data.frame(K = K.list, error.message = as.factor(error.list))
  
  jobs.tab.sum <- jobs.tab %>% dplyr::group_by(K, error.message) %>%
    dplyr::summarise(n=dplyr::n()) %>%  dplyr::mutate(prop.failed=round((n/replicates)*100)) %>%
    dplyr::filter(error.message != "success")
  
  
  # If they are at least the half of K which succeed, valid results
  if(nK_success !=0 ){
    
    if (cmd) print("#     => process results ")
    # Generate the list of results
    #coseq.res.list[["value"]] <- lapply(coseq.res.list,function(x){x$value})
    coseq.res.list[["value"]] <- list()
    for(x in names(coseq.res.list)){
      
      if(!is.null(coseq.res.list[[x]]$value)){
        coseq.res.list[["value"]][[x]] <- coseq.res.list[[x]]$value
      }
    }
    
    CoExpAnal <- coseq.results.process(coseq.res.list[["value"]], K = K, conds = conds)
    CoExpAnal[["warning"]] <- coseq.res.list$warning
    
    if(nK_success/length(iter) < 0.8){
      
      CoExpAnal[["error"]] <- TRUE
    }
    
  }
  # Réinitialisation de l'objet CoExpAnal
  else{
    CoExpAnal[["results"]] <- FALSE
    CoExpAnal[["error"]] <- TRUE
    
  }
  
  CoExpAnal[["stats"]] <- jobs.tab.sum
  
  return(CoExpAnal)
}




#' @title run Coseq for co-expression analysis on cluster
#' @param counts matrix
#' @param K number of groups
#' @param replicates number of replication to run
#' @param param.list list of coseq parameters
#' @return coseqResults
#' @importFrom coseq coseq clusters
#' @export
#' @noRd
#'
runCoseq_local <- function(counts, conds, K=2:20, replicates = 5, param.list, silent = TRUE, cmd = FALSE){
  
  iter <- rep(K, replicates)
  coseq.res.list <- list()
  
  # set.seed(12345)
  if (silent) {
    
    coseq.res.list <- lapply(1:replicates, function(x){
      
      co <- capture.output(suppressMessages(
        res <- try_rflomics(
          coseq::coseq(counts, K = K, parallel = TRUE,
                       model            = param.list[["model"]],
                       transformation   = param.list[["transformation"]],
                       meanFilterCutoff = param.list[["meanFilterCutoff"]],
                       normFactors      = param.list[["normFactors"]],
                       GaussianModel    = param.list[["GaussianModel"]],
                       seed = x, 
                       verbose = FALSE)
        )))
      return(res)
    })  
    
  }else{
    coseq.res.list <- lapply(1:replicates, function(x){
      
      try_rflomics(coseq::coseq(counts, K = K, parallel = TRUE,
                                model            = param.list[["model"]],
                                transformation   = param.list[["transformation"]],
                                meanFilterCutoff = param.list[["meanFilterCutoff"]],
                                normFactors      = param.list[["normFactors"]],
                                GaussianModel    = param.list[["GaussianModel"]],
                                seed = x))
    })  
  }
  
  names(coseq.res.list) <- c(1:replicates)
  
  CoExpAnal <- list()
  
  # error managment
  if (cmd) print("#     => error management: level 1 ")
  coseq.error.management <- coseq.error.manage(coseq.res.list = coseq.res.list, 
                                               K = K, 
                                               replicates = replicates,
                                               cmd = cmd)
  
  nK_success   <- coseq.error.management$nK_success
  
  # If they are at least the half of jobs succeed, valid results
  if (nK_success != 0) {
    
    CoExpAnal <-  coseq.results.process(coseqObjectList = coseq.error.management$coseq.res.list.values, 
                                        K = K,
                                        conds = conds)
    CoExpAnal[["warning"]] <- coseq.res.list$warning
    
    if (nK_success/length(iter) < 0.8) {
      CoExpAnal[["error"]] <- TRUE
    }
    
    if (cmd) { 
      print(paste0("#     => Number of clusters: ", 
                   max(unique(coseq::clusters(CoExpAnal$coseqResults)))))
    }
    
  }else{
    
    CoExpAnal[["results"]] <- FALSE
    CoExpAnal[["error"]] <- TRUE
  }
  
  CoExpAnal[["stats"]] <- coseq.error.management$jobs.tab.sum
  
  return(CoExpAnal)
}



#' @title profile boxplot per cluster
#' @param coseq.res coseq object
#' @param selectedCluster cluster num
#' @param conds condition matrix
#' @return boxplot profiles.
#' @export
#' @importFrom ggplot2 geom_boxplot facet_wrap theme element_blank
#' @importFrom dplyr arrange
#' @importFrom purrr reduce
#' @noRd
#'
coseq.y_profile.one.plot <- function(coseq.res, selectedCluster, conds){
  
  samples <- variable <- value <- cluster <- NULL
  
  nb_cluster <- coseq.res@metadata$nbCluster[min(coseq.res@metadata$ICL) == coseq.res@metadata$ICL]
  groups <- conds %>% dplyr::arrange(factor(samples, levels = names(coseq.res@y_profiles)))
  y_profiles <- list()
  
  for (i in 1:nb_cluster) {
    y_profiles[[i]] <- coseq.res@y_profiles[coseq.res@allResults[[paste0("K=", nb_cluster)]][,i] > 0.8,] %>%
      data.frame() %>% 
      reshape2::melt() %>%  
      dplyr::rename(samples = variable) %>%
      dplyr::full_join(conds , by = "samples") %>% 
      dplyr::mutate(cluster = i)
  }
  
  y_profiles.gg <-  y_profiles %>% purrr::reduce(rbind)
  y_profiles.gg$groups <- factor(y_profiles.gg$groups, levels = unique(conds$groups))
  y_profiles.gg$samples <- factor(y_profiles.gg$samples, levels = unique(conds$samples))
  
  
  p <- ggplot2::ggplot(data = dplyr::filter(y_profiles.gg, cluster == selectedCluster)) +
    ggplot2::geom_boxplot(ggplot2::aes(x = samples, y = value, fill = groups), outlier.size = 0.3) + 
    ggplot2::facet_wrap(~cluster) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  print(p)
}





## contrastName to name of contrast directory
# Exemple:
# contrastName
# "(temperatureMedium - temperatureElevated) in imbibitionEI - (temperatureMedium - temperatureElevated) in imbibitionDS"
# contrastDir
# "temperatureMedium-temperatureElevated_in_imbibitionEI_vs_temperatureMedium-temperatureElevated_in_imbibitionDS"

#' Title
#'
#' @param a string: contrastName
#'
#' @return a string: contrastDir
#' @importFrom stringr str_replace_all str_remove_all
#' @export
#' @noRd
contrastName2contrastDir <- function(contrastName){
  # remplacement des comparaisons centrales
  tmp <- stringr::str_replace_all(contrastName,"[:blank:]-[:blank:]\\(","_vs_")
  # remplacement des comparaisons dans les parenthèses
  tmp <- stringr::str_replace_all(tmp,"[:blank:]-[:blank:]","-")
  # suppression des parenthèses
  tmp <- stringr::str_remove_all(tmp,c("\\(|\\)"))
  # remplcement des espaces par des _
  tmp <- stringr::str_replace_all(tmp,"[:blank:]","_")
  return(tmp)
}

