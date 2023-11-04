######################## ANNOTATION USING CLUSTERPROFILER ########################

#' @title runAnnotationEnrichment
#' @description This function performs overrepresentation analysis (ORA) using clusterprofiler functions. It can be used with custom annotation file (via enricher), GO (enrichGO) or KEGG (enrichKEGG) annotations.
#' @param object An object of class \link{SummarizedExperiment} or \link{MultiAssayExperiment}. It is expected the SE object is produced by rflomics previous analyses, as it relies on their results.
#' @param SE.name name of the experiment to consider if object is a MultiAssayExperiment.
#' @param nameList name of contrasts (tags or names) from which to extract DE genes if from is DiffExpAnal. 
#' @param list_args list of arguments to pass to the enrichment function. These arguments must match the ones from the clusterprofiler package. E.g: universe, keytype, pvalueCutoff, qvalueCutoff, etc.
#' @param from indicates if ListNames are from differential analysis results (DiffExpAnal) or from the co-expression analysis results (CoExpAnal)
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
                        list_args = list(),
                        from = "DiffExpAnal",
                        dom.select = "custom",
                        Domain = "no-domain",
                        col_term = "term",
                        col_gene = "gene",
                        col_name = "name",
                        col_domain = NULL,
                        annot = NULL) {
  
    EnrichAnal <- list()

    if (is.null(object@metadata$DiffExpAnal)) {
      stop("There is no differential analysis. Please run a differential analysis before running enrichment")
    }
    
    # "Retrieving the lists of DE entities")
    switch(from,
           "DiffExpAnal" = {
             contrasts <- NULL
             
             if (is.null(getValidContrasts(object))) contrasts <- getSelectedContrasts(object)$contrastName
             else contrasts <- getValidContrasts(object)$contrastName                 
             
             if (!is.null(nameList)) {
               if (isTagName(object, nameList)) nameList <- convertTagToContrast(object, nameList)
               
               contrasts <- intersect(contrasts, nameList)
             }  
             
             geneLists <- lapply(contrasts, function(contrastName) {
               row.names(object@metadata$DiffExpAnal[["TopDEF"]][[contrastName]])
             })
             names(geneLists) <- contrasts
             
           },
           "CoExpAnal" = {
             
             namesClust <- names(object@metadata[["CoExpAnal"]][["clusters"]])
             if (!is.null(nameList)) namesClust <- intersect(namesClust, nameList)
             
             geneLists <- lapply(namesClust, function(namClust) {
               object@metadata[["CoExpAnal"]][["clusters"]][[namClust]]
             })
             names(geneLists) <- namesClust
           }
    )

    # Checks arguments
    if (dom.select == "custom") {
      if (is.null(annot)) {
        stop("You need an annotation file for a custom enrichment")
      }
      if (nrow(annot) < 1) {
        stop("Your annotation file seems to have 0 lines")
      }
      if (length(intersect(c(col_term, col_gene), colnames(annot))) != 2) {
        stop("The name of columns for gene and term names don't match the ones of the annotation files")
      }
    }

    # Change Domain if needed
    if (is.null(Domain)) {
      Domain <- "no-domain"
    }
    if (!is.null(col_domain) && dom.select == "custom") {
      if (is.null(annot[[col_domain]])) {
        stop("The column you indicated for the domain in your annotation file doesn't seem to exist.")
      } else {
        Domain <- unique(annot[[col_domain]])
        Domain <- Domain[!is.na(Domain)]
      }
    }

    # common parameters (is this useful ?)
    if (is.null(list_args$pvalueCutoff)) list_args$pvalueCutoff <- 0.05 # default in clusterprofiler
    if (is.null(list_args$qvalueCutoff)) list_args$qvalueCutoff <- 1 # no threshold on qvalue (default 0.2)
    # if(is.null(list_args$minGSSize))    list_args$minGSSize    <- 10 # default in clusterprofiler
    if (is.null(list_args$minGSSize)) list_args$minGSSize <- 3 # tried for SBML
    if (is.null(list_args$maxGSSize)) list_args$maxGSSize <- 500 # default in clusterprofiler
    if (is.null(list_args$universe)) list_args$universe <- names(object)

    annotation <- annot

    # for each list
    results_list <- lapply(names(geneLists), FUN = function(listname) {
      # listname <- names(geneLists)[1]
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

    ## summary

    # after filter
    #          BP  MF  CC
    # case1    10  22  4
    # case2    0   20  5
    # case3    0   0   0
    # case4    NA  10  3
    # case5    NA  NA  NA


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
    
    if (from == "DiffExpAnal") {
      object@metadata[["DiffExpEnrichAnal"]][[dom.select]] <- EnrichAnal
    } else if (from == "CoExpAnal") {
      object@metadata[["CoExpEnrichAnal"]][[dom.select]] <- EnrichAnal
    }

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
                        list_args = list(),
                        from = "DiffExpAnal",
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
                        from = "DiffExpEnrichAnal",
                        pvalueCutoff = NULL,
                        ...) {

    
    if (toupper(from) %in% toupper(c("DiffExpAnal", "DiffExpEnrichAnal"))) {
      from <- "DiffExpEnrichAnal"
    } else {
      from <- "CoExpEnrichAnal"
    }
    
    if (isTagName(object, contrast)) contrast <- convertTagToContrast(object, contrast)
    
    if (is.null(pvalueCutoff)) pvalueCutoff <- metadata(object)[[from]][["KEGG"]]$list_args$pvalueCutoff

    log2FC_vect <- NULL
    # Get the log2FC if appropriate
    if (from == "DiffExpEnrichAnal") {
      log2FC_vect <- object@metadata$DiffExpAnal[["TopDEF"]][[contrast]][["logFC"]]
      names(log2FC_vect) <- rownames(object@metadata$DiffExpAnal[["TopDEF"]][[contrast]])
    }

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
#' @exportMethod plotCPR
methods::setMethod(
  f = "plotCPR",
  signature = "SummarizedExperiment",
  definition = function(object,
                        contrast,
                        ont,
                        Domain = NULL,
                        from = "DiffExpAnal",
                        type = "dotplot",
                        showCategory = 15,
                        searchExpr = "",
                        node_label = "all",
                        pvalueCutoff = object@metadata$DiffExpEnrichAnal[[ont]]$list_args$pvalueCutoff,
                        ...) {
    # if from diffExpAnal, then takes the log2FC by default.
    # -> what if the user want something else printed ?! can modify it through scales ?
    # if from coexp, then no log2FC

    # dataPlot the enrichment results for correct ontology and contrast.
    
    # if (isTagName(contrast)) contrast <- convertTagToContrast(object, contrast)
    
    if (from == "DiffExpAnal") {
      dataPlot <- object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult[[contrast]]
    } else {
      dataPlot <- object@metadata$CoExpEnrichAnal[[ont]]$enrichResult[[contrast]]
    }

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


    log2FC_vect <- NULL
    # Get the log2FC if appropriate
    if (from == "DiffExpAnal") {
      log2FC_vect <- object@metadata$DiffExpAnal[["TopDEF"]][[contrast]][["logFC"]]
      names(log2FC_vect) <- rownames(object@metadata$DiffExpAnal[["TopDEF"]][[contrast]])
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
    if (type == "cnetplot") {
      suppressMessages( # delete warnings for scale fill replacement
        returnplot <- cnetplot(dataPlot, showCategory = Categories, color.params = list(foldChange = log2FC_vect), node_label = node_label, ...) +
          guides(colour = guide_colourbar(title = "log2FC")) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
      )
    } else if (type == "heatplot") {
      suppressMessages( # delete warnings for scale fill replacement
        returnplot <- heatplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, ...) +
          labs(fill = "log2FC") +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
          theme(axis.text.y = element_text(size = 10))
      )
    } else if (type == "dotplot") {
      returnplot <- dotplot(dataPlot, showCategory = Categories, ...)
    }

    return(returnplot)
  }
)




