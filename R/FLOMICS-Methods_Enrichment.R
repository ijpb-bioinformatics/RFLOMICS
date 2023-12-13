######################## ANNOTATION USING CLUSTERPROFILER ######################

#' @title runAnnotationEnrichment
#' @description This function performs overrepresentation analysis (ORA) using
#' clusterprofiler functions. It can be used with custom annotation file
#' (via enricher), GO (enrichGO) or KEGG (enrichKEGG) annotations.
#' @param object An object of class \link{SummarizedExperiment} or
#' \link{MultiAssayExperiment}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results.
#' @param SE.name name of the experiment to consider if object is a
#' MultiAssayExperiment.
#' @param nameList name of contrasts (tags or names) from which to extract DE
#' genes if from is DiffExpAnal.
#' @param list_args list of arguments to pass to the enrichment function.
#' These arguments must match the ones from the clusterprofiler package.
#' E.g: universe, keytype, pvalueCutoff, qvalueCutoff, etc.
#' @param from indicates if ListNames are from differential analysis results
#' (DiffExpAnal) or from the co-expression analysis results (CoExpAnal)
#' @param dom.select: is it a custom annotation, GO or KEGG annotations
#' @return A list of results from clusterprofiler.
#' @export
#' @importFrom dplyr filter select
#' @importFrom clusterProfiler enrichKEGG enrichGO enricher
#' @exportMethod runAnnotationEnrichment
#' @rdname runAnnotationEnrichment
methods::setMethod(
  f = "runAnnotationEnrichment",
  signature = "SummarizedExperiment",
  definition = function(object,
                        nameList = NULL,
                        list_args = list(
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 1,
                          # minGSSize = 10,
                          minGSSize = 3,
                          maxGSSize = 500,
                          universe = names(object)
                        ),
                        from = "DiffExp",
                        dom.select = "custom",
                        Domain = "no-domain",
                        col_term = "term",
                        col_gene = "gene",
                        col_name = "name",
                        col_domain = NULL,
                        annot = NULL) {
    
    EnrichAnal <- list()
    
    if (is.null(object@metadata$DiffExpAnal)) {
      stop("There is no differential analysis. Please run a differential 
           analysis before running enrichment")
    }
    
    searchFrom <- as.character(c(1,2)[c(grepl("DIFFEXP", toupper(from)), 
                                        grepl("COEXP", toupper(from)))])
    if (length(searchFrom) < 1) searchFrom <- 3
    
    switch(searchFrom, 
           "1" = { 
             from <- "DiffExp"
             contrasts <- NULL
             
             if (is.null(getValidContrasts(object))) {
               contrasts <- getSelectedContrasts(object)$contrastName
             } else contrasts <- getValidContrasts(object)$contrastName                 
             
             if (!is.null(nameList)) {
               if (isTagName(object, nameList)) nameList <- convertTagToContrast(object, nameList)
               
               contrasts <- intersect(contrasts, nameList)
             }
             
             geneLists <- lapply(contrasts, function(contrastName) {
               row.names(object@metadata$DiffExpAnal[["TopDEF"]][[contrastName]])
             })
             names(geneLists) <- contrasts
             
           },
           "2" = { 
             from <- "CoExp"   
             namesClust <- names(object@metadata[["CoExpAnal"]][["clusters"]])
             if (!is.null(nameList)) namesClust <- intersect(namesClust, nameList)
             
             geneLists <- lapply(namesClust, function(namClust) {
               object@metadata[["CoExpAnal"]][["clusters"]][[namClust]]
             })
             names(geneLists) <- namesClust
           },
           {
             message("Argument from is detected to be neither DiffExp nor CoExp, 
                     taking DiffExp results.")
             from <- "DiffExp"
           })
    
    if (is.null(Domain)) Domain <- "no-domain"
    
    if (dom.select == "custom") {
      if (is.null(annot)) {
        stop("You need an annotation file for a custom enrichment")
      }else if (length(intersect(c(col_term, col_gene), colnames(annot))) != 2) {
        stop("The name of columns for gene and term names don't match the ones of the annotation files")
      }
      if (!is.null(col_domain) && is.null(annot[[col_domain]])) {
        stop("The column you indicated for the domain in your annotation file doesn't seem to exist.")
      } else if (!is.null(col_domain)) {
        Domain <- unique(annot[[col_domain]])
        Domain <- Domain[!is.na(Domain)]
      }
      
      annotation <- annot
      
    } 
    
    # for each list
    results_list <- lapply(names(geneLists), FUN = function(listname) {
      list_args$gene <- geneLists[[listname]]
      results_ont <- lapply(Domain, FUN = function(ont) {
        switch(dom.select,
               "GO" = {
                 func_to_use <- "enrichGO"
                 list_args$ont <- ont
               },
               "KEGG" = {
                 func_to_use <- "enrichKEGG"
               },
               "custom" = {
                 func_to_use <- "enricher"
                 
                 list_args$TERM2NAME <- NA
                 annotation2 <- annotation
                 
                 if (ont != "no-domain") {
                   annotation2 <- filter(annotation, get(col_domain) == ont)
                 }
                 
                 list_args$TERM2GENE <- list("term" = annotation2[[col_term]], "gene" = annotation2[[col_gene]])
                 
                 if (!is.null(annotation2[[col_name]])) {
                   tmp <- dplyr::select(annotation2, tidyselect::all_of(c(col_term, col_name))) %>% unique()
                   list_args$TERM2NAME <- list("term" = tmp[[col_term]], "name" = tmp[[col_name]])
                 }
               }
        )
        do.call(getFromNamespace(func_to_use, ns = "clusterProfiler"), list_args)
      })
      names(results_ont) <- Domain
      return(results_ont)
    })
    names(results_list) <- names(geneLists)
    
    overview_list <- list()
    term.list <- list()
    for (listname in names(results_list)) {
      # TODO why is it a for?!
      
      for (ont in names(results_list[[listname]])) {
        if (!is.null(results_list[[listname]][[ont]])) {
          res <- results_list[[listname]][[ont]]@result
          res.n <- nrow(res[res$p.adjust < list_args$pvalueCutoff, ])
          overview_list[[listname]][[ont]] <- res.n
          
          if (res.n == 0) {
            results_list[[listname]][[ont]] <- NULL
          } else {
            term.list[[ont]][[listname]] <- data.frame(
              term = row.names(res[res$p.adjust < list_args$pvalueCutoff, ]),
              bin = rep(1, res.n)
            )
            names(term.list[[ont]][[listname]]) <- c("term", listname)
          }
        } else {
          results_list[[listname]][[ont]] <- NULL
        }
      }
    }
    
    if (length(overview_list) == 0) {
      EnrichAnal[["summary"]] <- NULL
    } else {
      dt_res <- as.data.frame(do.call("rbind", overview_list))
      if (!any(is.na(dt_res[, -1]))) {
        dt_res <- dt_res %>%
          dplyr::mutate(Contrast = rownames(.)) %>%
          dplyr::relocate(Contrast)
        EnrichAnal[["summary"]] <- dt_res
      }
    }
    
    EnrichAnal[["list_args"]] <- list_args[names(list_args) %in% c("universe", "keyType", "pvalueCutoff", "qvalueCutoff", "OrgDb", "organism")]
    EnrichAnal[["list_args"]] <- c(EnrichAnal[["list_args"]], list("Domain"=Domain))
    EnrichAnal[["enrichResult"]] <- results_list
    
    switch(from, 
           "DiffExp" = {
             object@metadata[["DiffExpEnrichAnal"]][[dom.select]] <- EnrichAnal
           },
           "CoExp" = {
             object@metadata[["CoExpEnrichAnal"]][[dom.select]] <- EnrichAnal
           }
    )
    
    return(object)
  }
)

#' @rdname runAnnotationEnrichment
#' @title runAnnotationEnrichment
#' @exportMethod runAnnotationEnrichment
methods::setMethod(
  f = "runAnnotationEnrichment",
  signature = "MultiAssayExperiment",
  definition = function(object,
                        SE.name,
                        nameList = NULL,
                        list_args = list(
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 1,
                          # minGSSize = 10,
                          minGSSize = 3,
                          maxGSSize = 500,
                          universe = names(object)
                        ),
                        from = "DiffExp",
                        dom.select = "custom",
                        Domain = "no-domain",
                        col_term = "term",
                        col_gene = "gene",
                        col_name = "name",
                        col_domain = NULL,
                        annot = NULL) {
    
    object[[SE.name]] <- runAnnotationEnrichment(object = object[[SE.name]],
                                                 nameList = nameList,
                                                 list_args = list_args,
                                                 from = from,
                                                 dom.select = dom.select,
                                                 Domain = Domain, 
                                                 col_term = col_term,
                                                 col_gene = col_gene,
                                                 col_domain = col_domain,
                                                 annot = annot)
    
    return(object)
    
  })

#' @title plotCPRKEGG
#' @description TODO
#' @param object An object of class \link{SummarizedExperiment}. It is expected the SE object is produced by rflomics previous analyses, as it relies on their results.
#' @param contrast the name of the contrast to consider. For Co expression analysis, it is expected to be one of "cluster.1", "cluster.2", etc.
#' @param ont the ontology (GO, KEGG or custom)
#' @param Domain if the ontology is GO, expect one of BP, MF or CC. Default is NULL.
#' @param from what type of analysis to consider? One of 'DiffExpAnal' or 'CoExpAnal'
#' @param type type of plot. Define the function used inside. One of dotplot, heatplot or cnetplot.
#' @param showCategory max number of terms to show.
#' @param searchExpr expression to search in the showCategory terms.
#' @param node_label same as in enrichplot::cnetplot function, defines the labels on the graph. One of "all", "category" or "gene". Default is 'all'.
#' @param pvalueCutoff pvalueCutoff to define the enrichment threshold. Default is the one find in the clusterprofiler results in object.
#' @param ... Not in use at the moment
#'
#' @return A plot.
#' @export
#' @exportMethod plotCPRKEGG
methods::setMethod(
  f = "plotCPRKEGG",
  signature = "SummarizedExperiment",
  definition = function(object,
                        contrast,
                        pathway_id = NULL,
                        species = "ath",
                        gene_idtype = "kegg",
                        from = "DiffExp",
                        pvalueCutoff = metadata(object)[["DiffExpEnrichAnal"]][["KEGG"]]$list_args$pvalueCutoff,
                        ...) {
    
    
    searchFrom <- as.character(c(1,2)[c(grepl("DIFFEXP", toupper(from)), 
                                        grepl("COEXP", toupper(from)))])
    if (length(searchFrom) < 1) searchFrom <- 3
    
    log2FC_vect <- NULL
    if (searchFrom == 1) {
      log2FC_vect <- object@metadata$DiffExpAnal[["TopDEF"]][[contrast]][["logFC"]]
      names(log2FC_vect) <- rownames(object@metadata$DiffExpAnal[["TopDEF"]][[contrast]])
    }
    
    if (isTagName(object, contrast)) contrast <- convertTagToContrast(object, contrast)
    
    see_pathview(
      gene.data = log2FC_vect,
      pathway.id = pathway_id,
      species = species,
      gene.idtype = gene_idtype,
      map.symbol = FALSE,
      same.layer = FALSE,
      low = list(gene = "blue"),
      mid = list(gene = "gray"),
      high = list(gene = "red"),
      na.col = "transparent"
      # cex = 1 # too much
    )
    
    return()
  }
)

#' @title plotCPR
#' @description TODO
#' @param object An object of class \link{SummarizedExperiment}. It is expected the SE object is produced by rflomics previous analyses, as it relies on their results.
#' @param contrast the name of the contrast to consider. For Co expression analysis, it is expected to be one of "cluster.1", "cluster.2", etc.
#' @param ont the ontology (GO, KEGG or custom)
#' @param Domain if the ontology is GO, expect one of BP, MF or CC. Default is NULL.
#' @param from what type of analysis to consider? One of 'DiffExpAnal' or 'CoExpAnal'
#' @param type type of plot. Define the function used inside. One of dotplot, heatplot or cnetplot.
#' @param showCategory max number of terms to show.
#' @param searchExpr expression to search in the showCategory terms.
#' @param node_label same as in enrichplot::cnetplot function, defines the labels on the graph. One of "all", "category" or "gene". Default is 'all'.
#' @param pvalueCutoff pvalueCutoff to define the enrichment threshold. Default is the one find in the clusterprofiler results in object.
#' @param ... additionnal parameters for cnetplot, heatplot or enrichplot functions.
#'
#' @return A plot.
#' @export
#' @importFrom enrichplot cnetplot heatplot dotplot
#' @importFrom ggrepel geom_label_repel
#' @exportMethod plotCPR
methods::setMethod(
  f = "plotCPR",
  signature = "SummarizedExperiment",
  definition = function(object,
                        contrast,
                        ont,
                        Domain = NULL,
                        from = "DiffExp",
                        type = "dotplot",
                        showCategory = 15,
                        searchExpr = "",
                        node_label = "all",
                        pvalueCutoff = object@metadata$DiffExpEnrichAnal[[ont]]$list_args$pvalueCutoff,
                        ...) {
    
    # if (isTagName(contrast)) contrast <- convertTagToContrast(object, contrast)
    
    searchFrom <- as.character(c(1,2)[c(grepl("DIFFEXP", toupper(from)), 
                                        grepl("COEXP", toupper(from)))])
    if (length(searchFrom) < 1) searchFrom <- 3
    
    log2FC_vect <- NULL
    switch(searchFrom, 
           "1" = { 
             from <- "DiffExp"
             dataPlot <- object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult[[contrast]]
             log2FC_vect <- object@metadata$DiffExpAnal[["TopDEF"]][[contrast]][["logFC"]]
             names(log2FC_vect) <- rownames(object@metadata$DiffExpAnal[["TopDEF"]][[contrast]])
           },
           "2" = { 
             from <- "CoExp"   
             dataPlot <- object@metadata$CoExpEnrichAnal[[ont]]$enrichResult[[contrast]]
           },
           {
             message("Argument from is detected to be neither DiffExp nor CoExp, 
                     taking DiffExp results.")
             from <- "DiffExp"
             dataPlot <- object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult[[contrast]]  
             log2FC_vect <- object@metadata$DiffExpAnal[["TopDEF"]][[contrast]][["logFC"]]  
             names(log2FC_vect) <- rownames(object@metadata$DiffExpAnal[["TopDEF"]][[contrast]])
           })
    
    if (ont == "GO") {
      if (is.null(Domain)) {
        stop("Ontology is GO, non null Domain (BP, CC or MF) is expected.")
      } else {
        dataPlot <- dataPlot[[Domain]]
      }
    } else if (ont == "custom" || ont == "KEGG") {
      if (is.null(Domain)) {
        Domain <- "no-domain"
      }
      
      if (!Domain %in% names(dataPlot)) {
        stop("Domain is expected to be one of ", paste(names(dataPlot), collapse = ","))
      } else {
        dataPlot <- dataPlot[[Domain]]
      }
    }
    
    # Select categories to show
    dataTab <- dataPlot@result[dataPlot@result$p.adjust < pvalueCutoff, ]
    Categories <- dataTab$Description
    if (searchExpr != "") Categories <- Categories[grep(toupper(searchExpr), toupper(Categories))]
    NbtoPlot <- min(length(Categories), showCategory)
    if (NbtoPlot == 0) stop("There is no terms to show")
    
    Categories <- Categories[1:NbtoPlot]
    
    # Create the plot
    type <- tolower(type)
    returnplot <- NULL
    
    returnplot <-  switch(type, 
                          "cnetplot" = {
                            cnetplot(dataPlot, showCategory = Categories,
                                     color.params = list(foldChange = log2FC_vect), 
                                     node_label = node_label, ...) +
                              guides(colour = guide_colourbar(title = "log2FC")) +
                              scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 
                          },
                          "heatplot" = { 
                            heatplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, ...) +
                              labs(fill = "log2FC") +
                              scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                              theme(axis.text.y = element_text(size = 10))},
                          {
                            tryCatch(dotplot(dataPlot, showCategory = Categories, ...),
                                     error = function(e) e,
                                     warnings = function(w) w)
                          })
    
    return(returnplot)
  }
)



#' @title plotEnrichComp
#' @description TODO
#' @param object An object of class \link{SummarizedExperiment}. It is expected the SE object is produced by rflomics previous analyses, as it relies on their results.
#' @param from indicates if the enrichment results are taken from differential analysis results (DiffExpAnal) or from the co-expression analysis results (CoExpAnal)
#' @param ont is it a custom annotation, GO or KEGG annotations
#' @param domain domain from the ontology (eg GO has three domains, BP, CC and MF)
#' @param matrixType Heatmap matrix to plot, one of GeneRatio, p.adjust or presence.
#' @param nClust number of separate cluster to plot on the heatmap, based on the clustering. 
#' @param decorate one of stars or GeneRatio. Decoration of the heatmap. Default is NULL, no decoration.
#' @param ... more arguments for ComplexHeatmap::Heatmap. 
#' @return A ggplot object
#' @export
#' @importFrom reshape2 recast
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap ht_opt
#' @importFrom stringr str_wrap
#' @exportMethod plotEnrichComp
#' @rdname plotEnrichComp
methods::setMethod(
  f = "plotEnrichComp",
  signature = "SummarizedExperiment",
  definition = function(object, 
                        from = "DiffExp", 
                        ont = NULL, 
                        domain = "no-domain",
                        matrixType = "GeneRatio",
                        nClust = NULL,
                        ...){
    
    allData <- switch(toupper(from), 
                      "DIFFEXP"           = { object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult },
                      "DIFFEXPANAL"       = { object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult },
                      "DIFFEXPENRICHANAL" = { object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult },
                      "COEXP"             = { object@metadata$CoExpEnrichAnal[[ont]]$enrichResult   },
                      "COEXPANAL"         = { object@metadata$CoExpEnrichAnal[[ont]]$enrichResult   },
                      "COEXPENRICHANAL"   = { object@metadata$CoExpEnrichAnal[[ont]]$enrichResult   },
                      {
                        message("Argument from is detected to be neither DiffExp nor CoExp, 
                     taking DiffExp results.")
                        allData <- object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult
                      })
    
    if (length(allData) == 0) {
      stop("The selected ontology ", ont, " does not seem to exist in the object.")
    }
    
    allData <- allData[lengths(allData) > 0]
    if (length(allData) == 0) {
      stop("It seems there is no results here.")
    }
    
    pvalThresh <- allData[[1]][[1]]@pvalueCutoff
    
    domainPoss <- unique(unlist(lapply(allData, names)))
    if (missing(domain) || is.null(domain)) domain <- domainPoss
    if (any(!domain %in% domainPoss)) {
      stop("Trying to select a domain that does not exist in the object.")
    }
    
    extract <- do.call("rbind", lapply(1:length(allData), FUN = function(i){
      if (domain %in% names(allData[[i]])) {
        cprRes <- allData[[i]][[domain]]
        cprRes <- cprRes@result[cprRes@result$p.adjust < cprRes@pvalueCutoff,]
        cprRes$contrast <- names(allData)[i]
        cprRes$domain <- domain
        return(cprRes)
      }
    }))
    
    toKeep <- names(which(table(extract$ID) > 1)) 
    if (length(toKeep) == 0) stop("There is no common terms to show.")
    extract <- extract[extract$ID %in% toKeep,]
    
    # handling description and ID
    extract$ID <- switch(ont, 
                         "GO" = { 
                           if (!identical(extract$ID, extract$Description)) {
                             paste0("(", extract$ID, ")", "\n", extract$Description)
                           } else {extract$ID}
                         },
                         "KEGG" = {
                           gsub(" - .*", "", extract$Description)
                         },
                         "custom" = {  
                           if (!identical(extract$ID, extract$Description)) {
                             paste0("(", extract$ID, ")", "\n", extract$Description)
                           } else {extract$ID}
                         })
    
    extract$ID <- str_wrap(extract$ID, width = 20)
    extract$contrast <- str_wrap(extract$contrast, width = 30)
    
    extract$GeneRatio <- as.numeric(vapply(extract$GeneRatio, 
                                           FUN = function(x) eval(parse(text = x)),
                                           FUN.VALUE = 1))
    extract$BgRatio <- as.numeric(vapply(extract$BgRatio, 
                                         FUN = function(x) eval(parse(text = x)),
                                         FUN.VALUE = 1))
    extract$FC <- extract$GeneRatio/extract$BgRatio
    
    dat <- switch(matrixType, 
                  "GeneRatio" = {
                    inter <- recast(extract[, c("ID", "contrast", "GeneRatio")], 
                                    ID ~ contrast, 
                                    measure.var = "GeneRatio")
                    rownames(inter) <- inter$ID
                    inter <- inter[, colnames(inter) != "ID"]
                    inter[is.na(inter)] <- 0
                    inter
                  }, 
                  "p.adjust" = {
                    inter <- recast(extract[, c("ID", "contrast", "p.adjust")], 
                                    ID ~ contrast, 
                                    measure.var = "p.adjust")
                    rownames(inter) <- inter$ID
                    inter <- inter[, colnames(inter) != "ID"]
                    inter[is.na(inter)] <- 1
                    inter
                  },
                  "presence" = {
                    inter <- recast(extract[, c("ID", "contrast", "p.adjust")], 
                                    ID ~ contrast, 
                                    measure.var = "p.adjust")
                    rownames(inter) <- inter$ID
                    inter <- inter <- inter[, colnames(inter) != "ID"]
                    inter[!is.na(inter)] <- 1
                    inter[is.na(inter)]  <- 0
                    inter
                  },
                  "FC" = {
                    inter <- recast(extract[, c("ID", "contrast", "FC")], 
                                    ID ~ contrast, 
                                    measure.var = "FC")
                    rownames(inter) <- inter$ID
                    inter <- inter[, colnames(inter) != "ID"]
                    inter[is.na(inter)] <- 0
                    inter
                  },
                  "log2FC" = {
                    inter <- recast(extract[, c("ID", "contrast", "FC")], 
                                    ID ~ contrast, 
                                    measure.var = "FC")
                    rownames(inter) <- inter$ID
                    inter <- inter[, colnames(inter) != "ID"]
                    inter <- log2(inter)
                    inter[is.infinite(as.matrix(inter))] <- 0
                    # means FC is 0, shouldn't happen much...
                    inter[is.na(inter)] <- 0
                    # means it's not significant and not in the matrix. 
                    inter
                  })
    
    if (nrow(dat) > 1) {
      switch(matrixType, 
             "presence" = { 
               hcPlot <- hclust(dist(dat, method = "binary"), method = "complete") 
               hcCol <- hclust(dist(t(dat), method = "binary"), method = "complete")
             },
             { 
               hcPlot <- hclust(dist(inter, method = "euclidean"), method = "complete")
               hcCol <- hclust(dist(t(inter), method = "euclidean"), method = "complete")
             })
    } else {
      hcPlot <- FALSE
    }
    
    colors <- switch(matrixType,
                     "presence"   = {structure(c("white", "firebrick"), names = c("0", "1"))},
                     "GeneRatio"  = {colorRamp2(c(0, max(dat)), c("white", "firebrick"))},
                     "p.adjust"   = {colorRamp2(c(0, pvalThresh, 1), c("firebrick", "white", "white"))},
                     "FC"         = {colorRamp2(c(0, max(dat)), c("white", "firebrick"))},
                     "log2FC"     = {colorRamp2(c(-max(abs(dat)), 0, max(abs(dat))), c("blue", "white", "firebrick"))})
    
    ht_opt(DENDROGRAM_PADDING = unit(0.1, "cm"))
    
    suppressWarnings(
      Heatmap(dat, 
              col = colors,
              name = matrixType,
              cluster_columns = hcCol, 
              show_column_dend = FALSE,
              cluster_rows = hcPlot, 
              row_names_side = "left", 
              column_names_rot = 30,
              # column_names_rot = 0, 
              # column_names_centered = TRUE, 
              rect_gp = gpar(col = "gray50", lwd = 0.5),
              width =  ncol(dat)*5,
              height = nrow(dat)*5,
              heatmap_legend_param = list(direction = "horizontal"),
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              row_names_gp = gpar(fontsize = 10)))
    
  }
)
