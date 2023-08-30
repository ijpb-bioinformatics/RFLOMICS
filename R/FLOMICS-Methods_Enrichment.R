######################## ANNOTATION USING CLUSTERPROFILER ########################

#' @title runAnnotationEnrichment_CPR
#' @description This function performs overrepresentation analysis (ORA) using clusterprofiler functions. It can be used with custom annotation file (via enricher), GO (enrichGO) or KEGG (enrichKEGG) annotations.
#' @param object An object of class \link{SummarizedExperiment}. It is expected the SE object is produced by rflomics previous analyses, as it relies on their results.
#' @param list_args list of arguments to pass to the enrichment function. These arguments must match the ones from the clusterprofiler package. E.g: universe, keytype, pvalueCutoff, qvalueCutoff, etc.
#' @param from indicates if ListNames are from differential analysis results (DiffExpAnal) or from the co-expression analysis results (CoExpAnal)
#' @param dom.select: is it a custom annotation, GO or KEGG annotations
#' @return A list of results from clusterprofiler.
#' @export
#' @exportMethod runAnnotationEnrichment_CPR
methods::setMethod(
  f = "runAnnotationEnrichment_CPR",
  signature = "SummarizedExperiment",
  definition = function(object,
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

    # "Retrieving the lists of DE entities")
    switch(from,
      "DiffExpAnal" = {
        geneLists <- lapply(object@metadata$DiffExpAnal$Validcontrasts$contrastName, function(contrastName) {
          row.names(object@metadata$DiffExpAnal[["TopDEF"]][[contrastName]])
        })
        names(geneLists) <- object@metadata$DiffExpAnal$Validcontrasts$contrastName
      },
      "CoExpAnal" = {
        geneLists <- lapply(object@metadata[["CoExpAnal"]][["clusters"]], function(clusters) {
          clusters
        })
        names(geneLists) <- names(object@metadata[["CoExpAnal"]][["clusters"]])
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
              annotation2 <- dplyr::filter(annotation, get(col_domain) == ont)
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

    EnrichAnal[["list_args"]] <- list_args[names(list_args) %in% c("universe", "keytype", "pvalueCutoff", "qvalueCutoff", "OrgDb")]
    EnrichAnal[["enrichResult"]] <- results_list

    if (from == "DiffExpAnal") {
      object@metadata[["DiffExpEnrichAnal"]][[dom.select]] <- EnrichAnal
    } else if (from == "CoExpAnal") {
      object@metadata[["CoExpEnrichAnal"]][[dom.select]] <- EnrichAnal
    }

    return(object)
  }
)

#' @title plot.CPRKEGG_Results
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
#' @exportMethod plot.CPRKEGG_Results
methods::setMethod(
  f = "plot.CPRKEGG_Results",
  signature = "SummarizedExperiment",
  definition = function(object,
                        contrast,
                        pathway_id = NULL,
                        species = "ath",
                        gene_idtype = "kegg",
                        from = "DiffExpAnal",
                        pvalueCutoff = object@metadata$DiffExpEnrichAnal[["KEGG"]]$list_args$pvalueCutoff,
                        ...) {
    if (from == "DiffExpAnal") {
      dataPlot <- object@metadata$DiffExpEnrichAnal[["KEGG"]]$enrichResult[[contrast]]
    } else {
      dataPlot <- object@metadata$CoExpAnal[["KEGG"]]$enrichResult[[contrast]]
    }

    log2FC_vect <- NULL
    # Get the log2FC if appropriate
    if (from == "DiffExpAnal") {
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

#' @title plot.CPR_Results
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
#' @exportMethod plot.CPR_Results
methods::setMethod(
  f = "plot.CPR_Results",
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
    if (from == "DiffExpAnal") {
      dataPlot <- object@metadata$DiffExpEnrichAnal[[ont]]$enrichResult[[contrast]]
    } else {
      dataPlot <- object@metadata$CoExpAnal[[ont]]$enrichResult[[contrast]]
    }

    if (ont == "GO") {
      if (is.null(Domain)) {
        stop("Ontology is GO, non null Domain (BP, CC or MF) is expected.")
      } else {
        dataPlot <- dataPlot[[Domain]]
      }
    } else if (ont == "custom") {
      if (is.null(Domain)) {
        Domain <- "no-domain"
      }

      if (!Domain %in% names(dataPlot)) {
        stop(paste0("Domain is expected to be one of ", paste(names(dataPlot), collapse = ",")))
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

    # Select cateogries to show
    dataTab <- dataPlot@result[dataPlot@result$p.adjust < pvalueCutoff, ]
    NbtoPlot <- min(nrow(dataTab), showCategory)
    Categories <- dataTab$Description[1:NbtoPlot]
    if (searchExpr != "") Categories <- Categories[grep(toupper(searchExpr), toupper(Categories))]


    # Create the plot
    returnplot <- NULL
    if (type == "cnetplot") {
      suppressMessages( # delete warnings for scale fill replacement
        returnplot <- enrichplot::cnetplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, node_label = node_label, ...) +
          ggplot2::guides(colour = ggplot2::guide_colourbar(title = "log2FC")) +
          ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
      )
    } else if (type == "heatplot") {
      suppressMessages( # delete warnings for scale fill replacement
        returnplot <- enrichplot::heatplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, ...) +
          ggplot2::labs(fill = "log2FC") +
          ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
      )
    } else if (type == "dotplot") {
      returnplot <- enrichplot::dotplot(dataPlot, showCategory = Categories, ...)
    }

    return(returnplot)
  }
)




