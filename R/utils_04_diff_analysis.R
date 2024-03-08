
######### Stat functions for differential analysis #########

#' @title .edgeRAnaDiff
#'
#' @param object an object of class \link{RflomicsSE}]
#' @param clustermq A boolean indicating if the constrasts have to be computed 
#' in local or in a distant machine
#' @param parallel boolean. Compute parallel differential analyses 
#' (only when clustermq is FALSE)
#' @param nworkers integer. Number of core to use for the parallel operations.
#'  Only used when parallel is TRUE.
#' @return A list of object of class \link{DGELRT}
#' @keywords internal
#' @importFrom stats model.matrix as.formula 
#' @importFrom edgeR DGEList estimateGLMCommonDisp estimateGLMTrendedDisp 
#' estimateGLMTagwiseDisp glmFit glmLRT topTags
#' @importFrom clustermq Q 
#' @importFrom parallel mclapply
#' @importFrom dplyr filter rename
#' @noRd
#'

.edgeRAnaDiff <- function(count_matrix, model_matrix, group, 
                          lib.size, norm.factors, Contrasts.Sel, 
                          Contrasts.Coeff, FDR, clustermq = FALSE,
                          parallel = FALSE, nworkers = 1, cmd = FALSE){
    
    z <- y <- NULL
    
    ListRes <- list()
    
    # check clustermq and parallel
    if (clustermq && parallel) parallel <- FALSE 
    if (!parallel) nworkers <- 1
    
    
    # Construct the DGE obect
    dge <- DGEList(counts       = count_matrix,
                   group        = group,
                   lib.size     = lib.size,
                   norm.factors = norm.factors)
    
    # Run the model
    if (cmd) message("[cmd] dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)")
    dge <- estimateGLMCommonDisp(dge, design=model_matrix)
    if (cmd) message("[cmd] dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)")
    dge <- estimateGLMTrendedDisp(dge, design=model_matrix)
    if (cmd) message("[cmd] dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)")
    dge <- estimateGLMTagwiseDisp(dge, design=model_matrix)
    if (cmd) message("[cmd] fit.f <- edgeR::glmFit(dge,design=model_matrix)")
    fit.f <- glmFit(dge,design=model_matrix)
    
    
    # test clustermq
    if (clustermq == TRUE){
        
        # Fonction to run on contrast per job
        # y is the model, Contrasts are stored in a matrix, by columns
        
        fx <- function(x){
            
            .tryRflomics <- function(expr) {
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
            
            .tryRflomics(glmLRT(y, contrast = unlist(z[x,])))
        }
        
        ResGlm <- Q(fx, x= seq_len(length(Contrasts.Sel$contrast)),
                    export=list(y=fit.f,z=Contrasts.Coeff),
                    n_jobs=length(Contrasts.Sel$contrast),pkgs="edgeR")
        
    }
    else{
        if(cmd) message("[cmd] apply model to each contrast")
        ResGlm <-  mclapply(Contrasts.Sel$contrast, function(x){
            .tryRflomics(glmLRT(fit.f, contrast = unlist(Contrasts.Coeff[x,])))
        }, mc.cores = nworkers)
        
    }
    
    # Create a table of jobs summary
    error.list <- unlist(lapply(ResGlm, function(x){
        ifelse(is.null(x$error),"success",as.character(x$error))
    }))
    
    jobs.tab <- data.frame(H=Contrasts.Sel$contrast, 
                           error.message=as.factor(error.list))
    
    jobs.tab.error <- jobs.tab %>% 
        filter(., error.message != "success")
    
    # If no error
    if(dim(jobs.tab.error)[1]==0){
        
        ListRes[[1]] <- lapply(ResGlm,function(x){
            x$value
        })
        
        # Name the table of raw results
        names(ListRes[[1]]) <- Contrasts.Sel$contrastName
        
        # ListRes[[2]] => TOPDGE => TopDFE
        
        TopDGE <- lapply(ListRes[[1]], function(x){
            
            res <-  topTags(x, n = dim(x)[1])
            
            DEGs <- res$table[res$table$FDR <= FDR,]
            #DEGs<-res$table
            return(DEGs)
        })
        
        ListRes[[2]] <- TopDGE
        
        names(ListRes[[2]]) <- names(ListRes[[1]])
        
        # Mutate column name to render the anadiff results generic
        # Initial column Name:  
        # gene_name  logFC      logCPM        LR        PValue           FDR
        ListRes[[2]] <- lapply(ListRes[[2]], function(x){
            rename(x,"Abundance"="logCPM",
                   "pvalue"="PValue",
                   "Adj.pvalue"="FDR")
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


#' @title .limmaAnaDiff
#'
#' @param object an object of class \link{RflomicsSE}
#' @param clustermq A boolean indicating if the constrasts have to be 
#' computed in local or in a distant machine
#' @return A list
#' @keywords internal
#' @importFrom stats model.matrix as.formula
#' @importFrom dplyr filter rename
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @noRd
#'


.limmaAnaDiff <- function(count_matrix, model_matrix, 
                          Contrasts.Sel, Contrasts.Coeff, 
                          p.adj.cutoff, p.adj.method,clustermq,
                          cmd = FALSE){
    
    ListRes <- list()
    
    # Run the model
    fit <- lmFit(count_matrix, model_matrix)
    
    # test clustermq
    if(clustermq == TRUE){
        
        fx <- function(x){
            
            try_rflomics <- function(expr) {
                warn <- err <- NULL
                value <- withCallingHandlers(
                    tryCatch(expr,
                             error    = function(e){ err <- e
                             NULL
                             }),
                    warning =function(w){ warn <- w
                    invokeRestart("muffleWarning")}
                )
                list(value=value, warning=warn, error=err)
            }
            
            # Fonction to run on contrast per job
            # y is the model, Contrasts are stored in a matrix, by columns
            
            .tryRflomics(contrasts.fit(y, contrasts  = unlist(z[x,])))
        }
        
        ResGlm  <- Q(fx,
                     x = seq_len(length(Contrasts.Sel$contrast)),
                     export = list(y = fit, z = Contrasts.Coeff),
                     n_jobs = length(Contrasts.Sel$contrast),
                     pkgs = "edgeR")
        
    }
    else{
        if(cmd) message("[cmd] fit contrasts")
        ResGlm <-  lapply(Contrasts.Sel$contrast, function(x){
            .tryRflomics(contrasts.fit(fit, contrasts  = as.vector(unlist(Contrasts.Coeff[x,]))))
        })
    }
    
    # Construct a table of jobs summary
    error.list <- unlist(lapply(ResGlm, function(x){
        ifelse(is.null(x$error),"success",as.character(x$error))
    }))
    
    jobs.tab <- data.frame(H=Contrasts.Sel$contrast , 
                           error.message=as.factor(error.list))
    
    jobs.tab.error <- jobs.tab %>% filter(., error.message != "success")
    
    # If no error
    if(dim(jobs.tab.error)[1]==0){
        
        ListRes[[1]] <- lapply(ResGlm,function(x){
            x$value
        })
        
        # Name the table of raw results
        
        names(ListRes[[1]]) <- Contrasts.Sel$contrastName
        
        # ListRes[[2]] TopDPE with column names common to all AnaDiff function
        
        ListRes[[2]] <- lapply(ListRes[[1]], function(x){
            
            fit2 <- eBayes(x, robust=TRUE)
            res <- topTable(fit2, adjust.method = p.adj.method, 
                            number=Inf, sort.by="AveExpr")
            DEPs <- res[res$adj.P.Val <= p.adj.cutoff,]
            return(DEPs)
        })
        
        
        # Mutate column name to render the anadiff results generic
        # Initial column Name:  
        # logFC  AveExpr         t      P.Value    adj.P.Val            B
        ListRes[[2]] <- lapply(ListRes[[2]], function(x){
            rename(x,"Abundance"="AveExpr",
                   "pvalue"="P.Value",
                   "Adj.pvalue"="adj.P.Val")
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


# ----- Get Summary for diffExpAnalysis : -----

#' @title Get summary table from diffExpAnalysis analysis
#'
#' @param object a SE object (produced by Flomics) or a MAE
#' @return a table
#' @importFrom dplyr filter
#' @export
sumDiffExp <- function(object, SE.name = NULL) {
    if (is(object, "RflomicsMAE")) {
        if (!is.null(SE.name)) {
            object <- object[[SE.name]]
        }
    }
    
    pcut <- getDiffSettings(object)$p.adj.cutoff
    lcut <- getDiffSettings(object)$abs.logFC.cutoff
    
    df_sim <- lapply(object@metadata$DiffExpAnal$DEF,
                     FUN = function(tab) {
                         tab <- tab %>%
                             filter(Adj.pvalue < pcut) %>%
                             filter(abs(logFC) > lcut)
                         
                         return(c("All" = nrow(tab),
                                  "Up" = nrow(tab %>% filter(logFC > 0)),
                                  "Down" = nrow(tab %>% filter(logFC < 0))
                         ))
                     }
    )
    return(do.call("rbind", df_sim))
}

######### Plot functions for differential analysis #########


#'.plotPValue
#'
#' @param data dataframe (ggplot2)
#' @param contrastName the contrast, useful for plot title
#' @return plot
#' @keywords internal
#' @importFrom ggplot2 geom_histogram theme_bw labs ggplot aes
#' @noRd
.plotPValue <- function(data, contrastName = contrastName) {
    PValue <- NULL
    
    p <- ggplot(data = data) +
        geom_histogram(aes(x = pvalue), bins = 100) +
        labs(x = expression(p - value),
             y = "count",
             title = contrastName) +
        theme_bw(base_size = 10)
    
    return(p)
}



#' MA.plot
#'
#' @param data dataframe (ggplot2)
#' @param p.adj.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff |log2FC| cutoff (absolute value)
#' @param contrastName the contrast, useful for plot title
#' @return MA plot
#' @keywords internal
#' @importFrom ggplot2 aes geom_point scale_colour_manual ggsave theme_linedraw
#' @importFrom ggpubr ggmaplot
#' @importFrom dplyr select rename
#' @noRd

.plotMA <- function(data,
                    p.adj.cutoff,
                    logFC.cutoff,
                    contrastName = contrastName) {
    Abundance <- logFC <- Adj.pvalue <- NULL
    
    tmp <- select(data, "Abundance", "logFC", "Adj.pvalue") %>%
        rename(.,
               baseMeanLog2 = Abundance,
               log2FoldChange = logFC,
               padj = Adj.pvalue
        )
    
    p <- ggmaplot(tmp,
                  main = contrastName,
                  fdr = p.adj.cutoff,
                  fc = 2 ^ logFC.cutoff,
                  size = 0.4,
                  ylab = bquote( ~ Log[2] ~ "fold change"),
                  xlab = bquote( ~ Log[2] ~ "mean expression"),
                  palette = c("#B31B21", "#1465AC", "grey30"),
                  select.top.method = c("padj", "fc"),
                  legend = "bottom",
                  top = 20,
                  font.label = c("plain", 7),
                  label.rectangle = TRUE,
                  font.legend = c(11, "plain", "black"),
                  font.main = c(11, "bold", "black"),
                  caption = paste(
                      "logFC cutoff=",
                      logFC.cutoff,
                      " and " ,
                      "FDR cutoff=",
                      p.adj.cutoff,
                      sep = ""
                  ),
                  ggtheme = theme_linedraw()
    )
    
    
    return(p)
    
}




#' Title
#'
#' @param data dataframe (ggplot2)
#' @param p.adj.cutoff adjusted pvalue cutoff
#' @param logFC.cutoff log2FC cutoff (absolute value)
#' @param contrastName the contrast, useful for plot title
#' @return a volcano plot, made with the \link{EnhancedVolcano} package.
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom dplyr first last 
#' @keywords internal
#' @noRd
#'
.plotVolcanoPlot <- function(data,
                             p.adj.cutoff,
                             logFC.cutoff,
                             contrastName) {
    # Find pvalue corresponding to the FDR cutoff for the plot
    # (mean between the last that passes the cutoff
    # and the first that is rejected to plot the line in the middle of
    #  the two points)
    # if pvalcutoff is 1 (no cutoff) no need to adjust
    
    if (p.adj.cutoff > 1) {
        stop("p.adj.cutoff must be between 0 and 1")
    }
    
    pval1 <- data$pvalue[data$Adj.pvalue < p.adj.cutoff] %>% last()
    pval2 <- data$pvalue[data$Adj.pvalue > p.adj.cutoff] %>% first()
    pvalCutoff <- (pval1 + pval2) / 2
    
    # If too low pvalues, unable to plot (error in if(d>0)...)
    # If drawconnectors is FALSE, it "works", with ylim being infinity,
    # it doesn't look like anything.
    # Modifiying the 0 pvalues to make sure it's working
    # default replacement in EnhancedVolcanoPlot
    nz_pval <- data$pvalue[data$pvalue != 0][1] * 10 ^ -1
    if (nz_pval == 0) {
        data$pvalue[data$pvalue == 0] <- data$pvalue[data$pvalue != 0][1]
    }
    
    Abundance <- logFC <- Adj.pvalue <- NULL
    p <- EnhancedVolcano(toptable = data,
                         lab = rownames(data),
                         x = 'logFC',
                         y = 'pvalue',
                         # pCutoff = p.adj.cutoff,
                         pCutoff = pvalCutoff,
                         FCcutoff = logFC.cutoff,
                         axisLabSize = 10,
                         pointSize = 1.5,
                         labSize = 2.5,
                         title = contrastName,
                         titleLabSize = 11,
                         subtitle = "",
                         subtitleLabSize = 10,
                         caption = paste("logFC cutoff=",
                                         logFC.cutoff,
                                         " & " , "FDR cutoff=",
                                         p.adj.cutoff,
                                         sep = ""),
                         legendPosition = "bottom",
                         legendLabSize = 10,
                         legendIconSize = 1.5,
                         captionLabSize = 10,
                         col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                         colAlpha = 0.5,
                         drawConnectors = FALSE,
                         widthConnectors = 0.5,
                         max.overlaps = 15
    )
    
    return(p)
}



