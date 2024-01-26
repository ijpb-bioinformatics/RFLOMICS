##### Global import ####

#' @importFrom ggplot2 ggplot geom_col theme_classic aes
#' theme element_text element_blank ylab xlab ggtitle
#' scale_fill_gradientn geom_tile theme_bw guides scale_fill_gradient2
#' guide_colourbar labs

# @export
#' @importFrom magrittr "%>%" 
magrittr::`%>%`




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
#' @param object an object of class \link{RflomicsSE}]
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
#' @param object an object of class \link{RflomicsSE}
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
  
  # if there is at least one failed job
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
    
    # if ICL == NA | loglike == 0
    
    # if (cmd) print("#     => error management: level 2 ")
    like.vec <- unlist(lapply(1:nK_success.job, function(x){ 
      coseq::likelihood(coseq.res.list[["value"]][[x]])
    })) %>%
      lapply(., function(x){ 
        ifelse(is.na(x) | (x==0), "failed", "success") }) %>% unlist()
    
    nK_success <- table(like.vec)["success"]
    
    replicates <- nK_success.job
    
    # expected list of cases
    K.list.ex <- rep(K, each = replicates)
    
    # observed list of cases
    K.list.ob <- stringr::str_replace(string = names(like.vec), pattern = "K=", replacement = "") %>% 
      as.numeric() %>% 
      sort()
    
    # missed cases
    if (length(K.list.ob) != length(K.list.ex)) {
      
      missed.K.vec <- names(table(K.list.ob)[table(K.list.ob) < nK_success.job])
      
      like.vec.bis <- rep("failed", length(missed.K.vec))
      names(like.vec.bis) <- paste0("K=", missed.K.vec)
      
      like.vec <- c(like.vec, like.vec.bis)
    }
    
    jobs.tab <- data.frame(K = names(like.vec), error.message = as.factor(like.vec))
    
    jobs.tab.sum2 <- jobs.tab %>% 
      dplyr::group_by(K, error.message) %>%
      dplyr::summarise(n = dplyr::n()) %>%  
      dplyr::mutate(prop.failed = round((n/replicates)*100)) %>%
      dplyr::filter(error.message != "success")
    
    jobs.tab.sum <- data.table::rbindlist(list(jobs.tab.sum1, jobs.tab.sum2), use.names = TRUE) %>% 
      tibble::tibble()
    
  } else{
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
  
  # get ICL as one vector for all replicates
  
  ICL.vec <- lapply(1:length(coseqObjectList), function(x){ coseq::ICL(coseqObjectList[[x]]) }) %>% unlist()
  ICL.list[["ICL.vec"]] <- ICL.vec
  
  # Find the ICL min by replicates
  
  ICL.min.per.rep <- lapply(1:length(coseqObjectList), function(x){
    ICL.vec.min <- coseq::ICL(coseqObjectList[[x]])
    ICL.vec.min[which(ICL.vec.min == min(ICL.vec.min))]
  })  %>% unlist()
  
  # Construct a table of results with ICL and K:
  
  ICL.tab <- data.frame(K = stringr::str_replace(names(ICL.vec), "K=", ""), ICL = ICL.vec) %>%
    dplyr::mutate(K = as.numeric(K))
  ICL.list[["ICL.tab"]] <- ICL.tab
  
  # Summarize the table: by K, compute the median of the replicate's ICL.
  
  ICL.n <- ICL.tab  %>% 
    dplyr::group_by(.,K) %>% 
    dplyr::filter(!is.na(ICL)) %>%
    dplyr::summarise(median = median(ICL, na.rm = TRUE), n = dplyr::n()) %>%
    dplyr::mutate(K = as.numeric(K))
  
  ICL.list[["ICL.n"]] <- ICL.n
  
  # Search for a replicate with a ICL min corresponding to the K with the min median 
  
  K.ICL.median.min <- ICL.n[which.min(ICL.n$median),]$K
  
  index  <- which(names(ICL.min.per.rep) == paste0("K=", K.ICL.median.min))
  
  # Case where the median.min has a K.min.rep
  if(length(index)>0){
    # Case where there is several rep with a min.rep, the min of them is taken
    index2 <- which(ICL.min.per.rep[index] == min(ICL.min.per.rep[index]))
    coseq.res <- coseqObjectList[index][index2][[1]]
    
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
  } else{
    # Pb of convergence: if there is no K.min.rep which correspond to the median.min, return an error
    CoExpAnal <- list()
    CoExpAnal[["results"]]      <- FALSE
    CoExpAnal[["error"]]  <- "No min.median correspond to min.ICL.rep"
  }
  
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
  
  # If they are more than 60 % of succeeded jobs, valid results, find min.median.ICL
  
  if (nK_success/length(iter) >=0.8) {
    
    CoExpAnal <-  coseq.results.process(coseqObjectList = coseq.error.management$coseq.res.list.values, 
                                        K = K,
                                        conds = conds)
    # If ICL.median has been found
    
    if(CoExpAnal[["results"]]==TRUE){
      
      CoExpAnal[["results"]] <- TRUE
      CoExpAnal[["warning"]] <- coseq.res.list$warning
      
      if (cmd) { 
        print(paste0("#     => Number of clusters: ", 
                     max(unique(coseq::clusters(CoExpAnal$coseqResults)))))
      }
    }
  } else{
    CoExpAnal[["results"]] <- FALSE
    CoExpAnal[["error"]] <- "nK_success/length(iter) < 0.8"
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



