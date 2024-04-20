######################## ANNOTATION USING CLUSTERPROFILER ######################

#' @title runAnnotationEnrichment
#' @description This function performs overrepresentation analysis (ORA) using
#' clusterprofiler functions. It can be used with custom annotation file
#' (via enricher), GO (enrichGO) or KEGG (enrichKEGG) annotations.
#' @param object An object of class \link{RflomicsSE} or
#' \link{RflomicsMAE}. It is expected the SE object is produced by
#' rflomics previous analyses, as it relies on their results.
#' @param SE.name name of the experiment to consider if object is a
#' RflomicsMAE.
#' @param nameList name of contrasts (tags or names) from which to extract DE
#' genes if from is DiffExpAnal.
#' @param list_args list of arguments to pass to the enrichment function.
#' These arguments must match the ones from the clusterprofiler package.
#' E.g: universe, keytype, pvalueCutoff, qvalueCutoff, etc.
#' @param from indicates if ListNames are from differential analysis results
#' (DiffExpAnal) or from the co-expression analysis results (CoExpAnal)
#' @param database is it a custom annotation, GO or KEGG annotations
#' @param domain subcatgory for the database (eg BP for GO)
#' @param col_term for custom annotation, column name of the file containing
#' term names or id.
#' @param col_gene for custom annotation, column name of the file containing
#' entities names.
#' @param col_name (optional) for custom annotation, column name of the 
#' file containing term names if col_term is used for ids.
#' @param col_domain for custom annotation, column name of the file containing
#' the domains.
#' @param annot for custom annotation, path of the annotation file.
#' @return A RflomicsMAE or a RflomicsS, depending on the class of object
#' parameter. The enrichment results are added to the metadata slot, either
#' in DiffExpEnrichAnal or CoExpEnrichAnal.
#' @importFrom dplyr filter select mutate relocate
#' @importFrom tidyselect all_of
#' @importFrom clusterProfiler enrichKEGG enrichGO enricher
#' @importFrom utils getFromNamespace
#' @exportMethod runAnnotationEnrichment
#' @rdname runAnnotationEnrichment
#' @aliases runAnnotationEnrichment
#' @examples
#' # Generate RflomicsMAE for example
#'  datPath <- paste0(system.file(package = "RFLOMICS"),
#'                    "/ExamplesFiles/ecoseed/")
#'
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#'facRef <- data.frame(
#'                 factorName   = c("Repeat", "temperature" , "imbibition"),
#'                 factorRef    = c("rep1",   "Low",          "DS"),
#'                 factorType   = c("batch",  "Bio",          "Bio"),
#'                 factorLevels = c("rep1,rep2,rep3",
#'                                  "Low,Medium,Elevated",
#'                                  "DS,EI,LI"))
#'                        
#' omicsData <- list(
#'    readOmicsData(file = paste0(datPath, "transcriptome_ecoseed.txt")),
#'    readOmicsData(file = paste0(datPath, "metabolome_ecoseed.txt")))
#'    
#' MAEtest <- createRflomicsMAE(projectName = "Tests",
#'                              omicsData   = omicsData,
#'                              omicsNames  = c("RNAtest",  "protetest"),
#'                              omicsTypes  = c("RNAseq", "proteomics"),
#'                              ExpDesign   = ExpDesign,
#'                              factorRef   = facRef
#'                              )
#' names(MAEtest) <- c("RNAtest", "protetest")
#' 
#' # Run GO annotation (enrichGO)
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(OrgDb = "org.At.tair.db",
#'                                          keyType = "TAIR",
#'                                          pvalueCutoff = 0.05),
#'                         from = "DiffExp", database = "GO",
#'                         domain = "CC")
#' getEnrichRes(MAEtest[["protetest"]])
#' 
#' # Run KEGG annotation (enrichKEGG)
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(organism = "ath",
#'                                          keyType = "kegg",
#'                                          pvalueCutoff = 0.05),
#'                         from = "DiffExp", database = "KEGG")
#' sumORA(MAEtest[["protetest"]], from = "DiffExp", database = "KEGG")
setMethod(
    f = "runAnnotationEnrichment",
    signature = "RflomicsSE",
    definition = function(object,
                          nameList = NULL,
                          list_args = list(
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 1,
                              minGSSize = 10,
                              maxGSSize = 500,
                              universe = names(object)
                          ),
                          from = "DiffExp",
                          database = "custom",
                          domain = "no-domain",
                          col_term = "term",
                          col_gene = "gene",
                          col_name = "name",
                          col_domain = NULL,
                          annot = NULL) {
        EnrichAnal <- list()
        
        if (is.null(object@metadata$DiffExpAnal)) {
            stop(
                "There is no differential analysis. Please run a differential
           analysis before running enrichment"
            )
        }
        
        from <- .determineFrom(from)
        
        switch(from,
               "DiffExp" = {
                   contrasts <- NULL
                   
                   if (is.null(getValidContrasts(object))) {
                       contrasts <- getSelectedContrasts(object)$contrastName
                   } else {
                       contrasts <- getValidContrasts(object)$contrastName
                   }
                   
                   if (!is.null(nameList)) {
                       if (.isTagName(object, nameList))
                           nameList <-
                               .convertTagToContrast(object, nameList)
                       
                       contrasts <- intersect(contrasts, nameList)
                   }
                   
                   geneLists <-
                       lapply(contrasts, function(contrastName) {
                           row.names(object@metadata$DiffExpAnal[["TopDEF"]][[contrastName]])
                       })
                   names(geneLists) <- contrasts
                   
               },
               "CoExp" = {
                   namesClust <- names(object@metadata[["CoExpAnal"]][["clusters"]])
                   if (!is.null(nameList))
                       namesClust <- intersect(namesClust, nameList)
                   
                   geneLists <- lapply(namesClust, function(namClust) {
                       object@metadata[["CoExpAnal"]][["clusters"]][[namClust]]
                   })
                   names(geneLists) <- namesClust
               },
               {
                   message(
                       "Argument from is detected to be neither DiffExp nor CoExp,
                     taking DiffExp results."
                   )
                   from <- "DiffExp"
               })
        
        # Domain
        if (is.null(domain))
            domain <- "no-domain"
        if (database == "GO" && "ALL" %in% toupper(domain)) {
            domain <- c("MF", "BP", "CC")
        }
        
        # If custom, make sure the annotation file matches the requirements
        if (database == "custom") {
            if (is.null(annot)) {
                stop("You need an annotation file for a custom enrichment")
            } else if (length(intersect(c(col_term, col_gene), colnames(annot))) != 2)
            {
                stop(
                    "The name of columns for gene and term names don't
             match the ones of the annotation files"
                )
            }
            if (!is.null(col_domain) &&
                is.null(annot[[col_domain]])) {
                stop(
                    "The column you indicated for the domain in your annotation
             file doesn't seem to exist."
                )
            } else if (!is.null(col_domain)) {
                domain <- unique(annot[[col_domain]])
                domain <- domain[!is.na(domain)]
            }
            annotation <- annot
        }
        
        geneLists <- geneLists[lengths(geneLists) > 0]
        
        # for each list
        results_list <-
            lapply(
                names(geneLists),
                FUN = function(listname) {
                    list_args$gene <- geneLists[[listname]]
                    results_ont <- lapply(
                        domain,
                        FUN = function(dom) {
                            switch(
                                database,
                                "GO" = {
                                    func_to_use <- "enrichGO"
                                    list_args$ont <- dom
                                },
                                "KEGG" = {
                                    func_to_use <- "enrichKEGG"
                                },
                                "custom" = {
                                    func_to_use <- "enricher"
                                    
                                    list_args$TERM2NAME <- NA
                                    annotation2 <- annotation
                                    
                                    if (dom != "no-domain") {
                                        annotation2 <- filter(annotation, get(col_domain) == dom)
                                    }
                                    
                                    list_args$TERM2GENE <-
                                        list("term" = annotation2[[col_term]],
                                             "gene" = annotation2[[col_gene]])
                                    
                                    if (!is.null(annotation2[[col_name]])) {
                                        tmp <- select(annotation2, all_of(c(col_term, col_name))) %>%
                                            unique()
                                        list_args$TERM2NAME <-
                                            list("term" = tmp[[col_term]],
                                                 "name" = tmp[[col_name]])
                                    }
                                }
                            )
                           res1 <- do.call(getFromNamespace(func_to_use, 
                                                            ns = "clusterProfiler"),
                                    list_args)
                           # delete heavy slots
                           slot(res1, name = "geneSets", 
                                check = FALSE) <- list("removed")
                           slot(res1, name = "universe", 
                                check = FALSE) <- c("removed")
                           return(res1)
                        }
                    )
                    names(results_ont) <- domain
                    return(results_ont)
                }
            )
        names(results_list) <- names(geneLists)
        
        overview_list <- list()
        term.list <- list()
        for (listname in names(results_list)) {
            # TODO why is it a for?!
            
            for (dom in names(results_list[[listname]])) {
                if (!is.null(results_list[[listname]][[dom]])) {
                    res <- results_list[[listname]][[dom]]@result
                    res.n <-
                        nrow(res[res$p.adjust < list_args$pvalueCutoff,])
                    overview_list[[listname]][[dom]] <- res.n
                    
                    if (res.n == 0) {
                        results_list[[listname]][[dom]] <- NULL
                    } else {
                        term.list[[dom]][[listname]] <- data.frame(term = row.names(res[res$p.adjust < list_args$pvalueCutoff,]),
                                                                   bin = rep(1, res.n))
                        names(term.list[[dom]][[listname]]) <-
                            c("term", listname)
                    }
                } else {
                    results_list[[listname]][[dom]] <- NULL
                }
            }
        }
        
        if (length(overview_list) == 0) {
            EnrichAnal[["summary"]] <- NULL
        } else {
            dt_res <- as.data.frame(do.call("rbind", overview_list))
            if (!any(is.na(dt_res[, -1]))) {
                dt_res <- dt_res %>%
                    mutate(Contrast = rownames(.)) %>%
                    relocate(Contrast)
                EnrichAnal[["summary"]] <- dt_res
            }
        }
        
        EnrichAnal[["list_args"]] <-
            list_args[names(list_args) %in%
                          c("universe",
                            "keyType",
                            "pvalueCutoff",
                            "qvalueCutoff",
                            "OrgDb",
                            "organism")]
        EnrichAnal[["list_args"]] <- c(EnrichAnal[["list_args"]],
                                       list("domain" = domain))
        EnrichAnal[["enrichResult"]] <- results_list
        
        switch(from,
               "DiffExp" = {
                   object@metadata[["DiffExpEnrichAnal"]][[database]] <- EnrichAnal
               },
               "CoExp" = {
                   object@metadata[["CoExpEnrichAnal"]][[database]] <- EnrichAnal
               })
        
        return(object)
    }
)

#' @rdname runAnnotationEnrichment
#' @aliases runAnnotationEnrichment
#' @title runAnnotationEnrichment
#' @exportMethod runAnnotationEnrichment
setMethod(
    f = "runAnnotationEnrichment",
    signature = "RflomicsMAE",
    definition = function(object,
                          SE.name,
                          nameList = NULL,
                          list_args = list(
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 1,
                              minGSSize = 10,
                              maxGSSize = 500,
                              universe = names(object)
                          ),
                          from = "DiffExp",
                          database = "custom",
                          domain = "no-domain",
                          col_term = "term",
                          col_gene = "gene",
                          col_name = "name",
                          col_domain = NULL,
                          annot = NULL) {
        object[[SE.name]] <-
            runAnnotationEnrichment(
                object = object[[SE.name]],
                nameList = nameList,
                list_args = list_args,
                from = from,
                database = database,
                domain = domain,
                col_term = col_term,
                col_gene = col_gene,
                col_domain = col_domain,
                annot = annot
            )
        
        return(object)
        
    }
)


#' @title plotKEGG
#' @description Plot the KEGG pathway using the pathview package. It overlays
#' the result of the differential analysis or the clusters entity on the
#' image.
#' @param object An object of class \link{RflomicsSE}.
#' It is expected the SE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param contrastName the name of the contrast to consider.
#' For Co expression analysis,
#' it is expected to be one of "cluster.1", "cluster.2", etc.
#' @param pathway_id the KEGG id pathway to plot.
#' @param species Kegg code (eg hsa or ath)
#' @param gene_idtype idtype (kegg, uniprot,...)
#' @param from differential analysis (diffExp) or coexpression analysis (coexp)
#' Used to link contrastName to the right metadata.
#' @param ... Not in use at the moment
#'
#' @return Only displays the KEGG pathway, it does not return any object.
#' @exportMethod plotKEGG
#' @rdname plotKEGG
#' @aliases plotKEGG
#' @examples
#' # Generate RflomicsMAE for example
#'  datPath <- paste0(system.file(package = "RFLOMICS"),
#'                    "/ExamplesFiles/ecoseed/")
#'
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#'facRef <- data.frame(
#'                 factorName   = c("Repeat", "temperature" , "imbibition"),
#'                 factorRef    = c("rep1",   "Low",          "DS"),
#'                 factorType   = c("batch",  "Bio",          "Bio"),
#'                 factorLevels = c("rep1,rep2,rep3",
#'                                  "Low,Medium,Elevated",
#'                                  "DS,EI,LI"))
#'                        
#' omicsData <- list(
#'    readOmicsData(file = paste0(datPath, "transcriptome_ecoseed.txt")),
#'    readOmicsData(file = paste0(datPath, "metabolome_ecoseed.txt")))
#'    
#' MAEtest <- createRflomicsMAE(projectName = "Tests",
#'                              omicsData   = omicsData,
#'                              omicsNames  = c("RNAtest",  "protetest"),
#'                              omicsTypes  = c("RNAseq", "proteomics"),
#'                              ExpDesign   = ExpDesign,
#'                              factorRef   = facRef
#'                              )
#' names(MAEtest) <- c("RNAtest", "protetest")
#' 
#' # Run KEGG annotation (enrichKEGG)
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                        list_args = list(organism = "ath",
#'                        keyType = "kegg",
#'                        pvalueCutoff = 0.05),
#'                        from = "CoExp", database = "KEGG")
#'                        
#' plotKEGG(MAEtest[["protetest"]], pathway_id = "ath00710", species = "ath",
#'           contrastName = "cluster.4", from = "Coexp")
#'
setMethod(
    f = "plotKEGG",
    signature = "RflomicsSE",
    definition = function(object,
                          contrastName,
                          pathway_id = NULL,
                          species = "ath",
                          gene_idtype = "kegg",
                          from = "DiffExp",
                          ...) {
        from <- .determineFromEnrich(from)
        
        log2FC_vect <- NULL
        switch(from,
               "DiffExpEnrichAnal" = {
                   obj <- metadata(object)$DiffExpAnal[["TopDEF"]][[contrastName]]
                   log2FC_vect <- obj[["logFC"]]
                   names(log2FC_vect) <- rownames(obj)
                   
                   .see_pathview(
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
                   )
               },
               "CoExpEnrichAnal" = {
                   entities <- getClusterEntities(object, contrastName)
                   log2FC_vect <- rep(1, length(entities))
                   names(log2FC_vect) <- entities
                   
                   .see_pathview(
                       gene.data = log2FC_vect,
                       pathway.id = pathway_id,
                       species = species,
                       gene.idtype = gene_idtype,
                       plot.col.key = FALSE,
                       map.symbol = FALSE,
                       same.layer = FALSE,
                       high = list(gene = "antiquewhite3"),
                       na.col = "transparent"
                   )
               }) # end switch
        
        return()
    }
)

#' @title plotClusterProfiler
#' @description Plot a dotplot, a cnetplot or an heatplot, using enrichplot
#' package. It is a wrapper method destined for the RflomicsSE class.
#' @param object An object of class \link{RflomicsSE}.
#' It is expected the SE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param contrastName the name of the contrast to consider.
#' For Co expression analysis, it is expected to be one of
#' "cluster.1", "cluster.2", etc.
#' @param database the database (GO, KEGG or custom)
#' @param domain if the database is GO, expect one of BP, MF or CC.
#' Default is NULL.
#' @param from what type of analysis to consider?
#' One of 'DiffExp' or 'CoExp'
#' @param plotType type of plot. Define the function used inside.
#' One of dotplot, heatplot or cnetplot.
#' @param showCategory max number of terms to show.
#' @param searchExpr expression to search in the showCategory terms.
#' @param nodeLabel same as in enrichplot::cnetplot function, defines
# the labels on the graph. One of "all", "category" or "gene". Default is 'all'.
#' @param p.adj.cutoff pvalueCutoff to define the enrichment threshold.
# Default is the one find in the clusterprofiler results in object.
#' @param ... additionnal parameters for cnetplot, heatplot or
#' enrichplot functions.
#'
#' @return A ggplot object.
#' @importFrom enrichplot cnetplot heatplot dotplot set_enrichplot_color
#' @importFrom ggplot2 scale_fill_gradient2 guide_colourbar
#' @importFrom ggrepel geom_label_repel
#' @exportMethod plotClusterProfiler
#' @rdname plotClusterProfiler
#' @aliases plotClusterProfiler
#' @examples
#' # Generate RflomicsMAE for example
#' MAEtest <- generateExample(annotation = FALSE, integration = FALSE)
#' # Run KEGG annotation on differential analysis results
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(organism = "ath",
#'                                          keyType = "kegg",
#'                                          pvalueCutoff = 0.05),
#'                         from = "DiffExp", database = "KEGG")
#'
#' # From differential analysis proteins lists:
#' plotClusterProfiler(MAEtest[["protetest"]],
#'         contrastName = "(temperatureElevated - temperatureMedium) in mean",
#'         database = "KEGG", from = "DiffExp",
#'         plotType = "heatplot", p.adj.cutoff = 0.05,
#'         domain = "no-domain")
#'
setMethod(
    f = "plotClusterProfiler",
    signature = "RflomicsSE",
    definition = function(object,
                          contrastName,
                          database,
                          domain = "no-domain",
                          from = "DiffExp",
                          plotType = "dotplot",
                          showCategory = 15,
                          searchExpr = "",
                          nodeLabel = "all",
                          p.adj.cutoff = NULL,
                          ...) {
        
        from <- .determineFromEnrich(from)
        
        log2FC_vect <- NULL
        switch(from,
               "DiffExpEnrichAnal" = {
                   dataPlot <- getEnrichRes(
                       object,
                       contrastName = contrastName,
                       database = database,
                       from = from
                   )
                   inter <-
                       metadata(object)$DiffExpAnal[["TopDEF"]][[contrastName]]
                   log2FC_vect <- inter[["logFC"]]
                   names(log2FC_vect) <- rownames(inter)
                   
                   if (is.null(p.adj.cutoff)) {
                       p.adj.cutoff <- getEnrichPvalue(object, from = from,
                                                       database = database)
                   }
               },
               "CoExpEnrichAnal" = {
                   dataPlot <- getEnrichRes(
                       object,
                       contrastName = contrastName,
                       database = database,
                       from = from
                   )
                   if (is.null(p.adj.cutoff)) {
                       p.adj.cutoff <- getEnrichPvalue(object, from = from,
                                                       database = database)
                   }
               },
               {
                   message(
                       "Argument from is detected to be neither DiffExp
                     nor CoExp, taking DiffExp results."
                   )
                   from <- "DiffExpEnrichAnal"
                   dataPlot <-
                       getEnrichRes(
                           object,
                           contrastName = contrastName,
                           database = database,
                           from = from
                       )
                   inter <-
                       metadata(object)$DiffExpAnal[["TopDEF"]][[contrastName]]
                   log2FC_vect <- inter[["logFC"]]
                   names(log2FC_vect) <- rownames(inter)
                   
                   if (is.null(p.adj.cutoff)) {
                       p.adj.cutoff <- getEnrichPvalue(object, from = from,
                                                       database = database)
                   }
                   
               })
        
        if (database == "GO") {
            if (is.null(domain)) {
                stop("database is GO, non null domain (BP, CC or MF) is expected.")
            } else {
                dataPlot <- dataPlot[[domain]]
            }
        } else if (database == "custom" || database == "KEGG") {
            if (is.null(domain)) {
                domain <- "no-domain"
            }
            
            if (!domain %in% names(dataPlot)) {
                stop("domain is expected to be one of ",
                     paste(names(dataPlot), collapse = ","))
            } else {
                dataPlot <- dataPlot[[domain]]
            }
        }
        
        # Select categories to show
        dataTab <-
            dataPlot@result[which(dataPlot@result$p.adjust < p.adj.cutoff),]
        Categories <- dataTab$Description
        if (searchExpr != "") {
            Categories <-
                Categories[grep(toupper(searchExpr), toupper(Categories))]
        }
        NbtoPlot <- min(length(Categories), showCategory)
        if (NbtoPlot == 0)
            stop("There is no terms to show")
        
        Categories <- Categories[seq_len(NbtoPlot)]
        
        # Create the plot
        plotType <- tolower(plotType)
        returnplot <- NULL
        
        returnplot <-  switch(
            plotType,
            "cnetplot" = {
                cnetplot(
                    dataPlot,
                    showCategory = Categories,
                    color.params = list(foldChange = log2FC_vect),
                    node_label = nodeLabel,
                    ...) +
                    guides(colour = guide_colourbar(title = "log2FC")) 
            },
            "heatplot" = {
                outgg <- heatplot(dataPlot,
                                  showCategory = Categories,
                                  foldChange = log2FC_vect,
                                  ...) 
                outgg$scales$scales <- list()
                outgg + labs(fill = "log2FC") +
                    scale_fill_gradient2(
                        low = "blue",
                        mid = "white",
                        high = "red",
                        midpoint = 0
                    ) +
                    theme(axis.text.y = element_text(size = 10))
            },
            {
                tryCatch(
                    dotplot(dataPlot,
                            showCategory = Categories,
                            ...),
                    error = function(e)
                        e,
                    warnings = function(w)
                        w
                )
            })
        return(returnplot)
    }
)



#' @title plotEnrichComp
#' @description Plot an heatmap of all the enriched term found for a given
#' database and a given source (differential analysis or coexpression clusters).
#' Allow for the comparison of several enrichment results.
#' @param object An object of class \link{RflomicsSE}.
#' It is expected the SE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param from indicates if the enrichment results are taken from differential
#' analysis results (DiffExpAnal) or from the co-expression analysis results
#'  (CoExpAnal)
#' @param database is it a custom annotation, GO or KEGG annotations
#' @param domain domain from the database (eg GO has three domains,
#' BP, CC and MF)
#' @param matrixType Heatmap matrix to plot, one of GeneRatio, p.adjust or
#' presence.
#' @param nClust number of separate cluster to plot on the heatmap, based o
#' n the clustering.
#' @param ... more arguments for ComplexHeatmap::Heatmap.
#' @return A complexHeatmap object.
#'
#' @importFrom reshape2 recast
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap ht_opt
#' @importFrom stats hclust dist
#' @importFrom stringr str_wrap
#' @exportMethod plotEnrichComp
#' @rdname plotEnrichComp
#' @aliases plotEnrichComp
#' @examples
#' # Generate RflomicsMAE for example
#' MAEtest <- generateExample(annotation = FALSE, integration = FALSE)
#' # Run KEGG annotation on differential analysis results
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(organism = "ath",
#'                                          keyType = "kegg",
#'                                          pvalueCutoff = 0.1),
#'                         from = "DiffExp", database = "KEGG")
#'
#' plotEnrichComp(MAEtest[["protetest"]], from = "DiffExp",
#'                database = "KEGG", matrixType = "FC")
#'
setMethod(
    f = "plotEnrichComp",
    signature = "RflomicsSE",
    definition = function(object,
                          from = "DiffExp",
                          database = NULL,
                          domain = "no-domain",
                          matrixType = "FC",
                          nClust = NULL,
                          ...) {
        from <- .determineFromEnrich(from)
        
        allData <- switch(
            toupper(from),
            "DIFFEXPENRICHANAL" = {
                object@metadata$DiffExpEnrichAnal[[database]]$enrichResult
            },
            "COEXPENRICHANAL"   = {
                object@metadata$CoExpEnrichAnal[[database]]$enrichResult
            },
            {
                message(
                    "Argument from is detected to be neither DiffExp nor CoExp,
                     taking DiffExp results."
                )
                allData <-
                    object@metadata$DiffExpEnrichAnal[[database]]$enrichResult
            }
        )
        
        if (length(allData) == 0) {
            stop("The selected database ",
                 database,
                 " does not seem to exist in the object.")
        }
        
        allData <- allData[lengths(allData) > 0]
        if (length(allData) == 0) {
            stop("It seems there is no results here.")
        }
        
        pvalThresh <- allData[[1]][[1]]@pvalueCutoff
        
        domainPoss <- unique(unlist(lapply(allData, names)))
        if (missing(domain) || is.null(domain))
            domain <- domainPoss
        if (any(!domain %in% domainPoss)) {
            stop("Trying to select a domain that does not exist in the object.")
        }
        
        extract <-
            do.call("rbind", lapply(
                seq_len(length(allData)),
                FUN = function(i) {
                    if (domain %in% names(allData[[i]])) {
                        cprRes <- allData[[i]][[domain]]
                        cprRes <-
                            cprRes@result[cprRes@result$p.adjust < cprRes@pvalueCutoff, ]
                        cprRes$contrastName <- names(allData)[i]
                        cprRes$domain <- domain
                        return(cprRes)
                    }
                }
            ))
        
        toKeep <- names(which(table(extract$ID) > 1))
        if (length(toKeep) == 0)
            stop("There is no common terms to show.")
        extract <- extract[extract$ID %in% toKeep, ]
        
        # handling description and ID
        extract$Description <- switch(database,
                                      "KEGG" = {
                                          gsub(" - .*", "", extract$Description)
                                      },
                                      {
                                          # if there is duplication in description
                                          # not necessarily duplicated in the ID
                                          # should not be grouped in the heatmap
                                          if (sum(duplicated(extract$Description)) > 0) {
                                              if (!identical(extract$Description,
                                                             extract$ID)) {
                                                  posDup <-
                                                      unique(which(
                                                          duplicated(extract$Description),
                                                          duplicated(extract$Description, fromLast = TRUE)
                                                      ))
                                                  extract$Description[posDup] <-
                                                      paste0("(",
                                                             extract$ID[posDup],
                                                             ")",
                                                             "\n",
                                                             extract$Description[posDup])
                                                  extract$Description
                                              } else{
                                                  extract$Description
                                              }
                                              
                                          } else {
                                              extract$Description
                                          }
                                      })
        
        extract$Description <-
            str_wrap(extract$Description, width = 20)
        extract$contrastName <-
            str_wrap(extract$contrastName, width = 30)
        
        extract$GeneRatio <- as.numeric(vapply(
            extract$GeneRatio,
            FUN = function(x)
                eval(parse(text = x)),
            FUN.VALUE = 1
        ))
        extract$BgRatio <- as.numeric(vapply(
            extract$BgRatio,
            FUN = function(x)
                eval(parse(text = x)),
            FUN.VALUE = 1
        ))
        extract$FC <- extract$GeneRatio / extract$BgRatio
        
        dat <- switch(
            matrixType,
            "GeneRatio" = {
                inter <- recast(extract[, c("Description", "contrastName", "GeneRatio")],
                                Description ~ contrastName,
                                measure.var = "GeneRatio")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "ID"]
                inter[is.na(inter)] <- 0
                inter
            },
            "p.adjust" = {
                inter <- recast(extract[, c("Description", "contrastName", "p.adjust")],
                                Description ~ contrastName,
                                measure.var = "p.adjust")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "Description"]
                inter[is.na(inter)] <- 1
                inter
            },
            "presence" = {
                inter <- recast(extract[, c("Description", "contrastName", "p.adjust")],
                                Description ~ contrastName,
                                measure.var = "p.adjust")
                rownames(inter) <- inter$Description
                inter <-
                    inter <- inter[, colnames(inter) != "Description"]
                inter[!is.na(inter)] <- 1
                inter[is.na(inter)]  <- 0
                inter
            },
            "FC" = {
                inter <- recast(extract[, c("Description", "contrastName", "FC")],
                                Description ~ contrastName,
                                measure.var = "FC")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "Description"]
                inter[is.na(inter)] <- 0
                inter
            },
            "log2FC" = {
                inter <- recast(extract[, c("Description", "contrastName", "FC")],
                                Description ~ contrastName,
                                measure.var = "FC")
                rownames(inter) <- inter$Description
                inter <- inter[, colnames(inter) != "Description"]
                inter <- log2(inter)
                inter[is.infinite(as.matrix(inter))] <- 0
                # means FC is 0, shouldn't happen much...
                inter[is.na(inter)] <- 0
                # means it's not significant and not in the matrix.
                inter
            }
        )
        
        if (nrow(dat) > 1) {
            switch(matrixType,
                   "presence" = {
                       hcPlot <- hclust(dist(dat, method = "binary"),
                                        method = "complete")
                       hcCol <-
                           hclust(dist(t(dat), method = "binary"),
                                  method = "complete")
                   },
                   {
                       hcPlot <- hclust(dist(inter, method = "euclidean"),
                                        method = "complete")
                       hcCol <-
                           hclust(dist(t(inter), method = "euclidean"),
                                  method = "complete")
                   })
        } else {
            hcPlot <- FALSE
        }
        
        colors <- switch(
            matrixType,
            "presence"   = {
                structure(c("white", "firebrick"), names = c("0", "1"))
            },
            "GeneRatio"  = {
                colorRamp2(c(0, max(dat)), c("white", "firebrick"))
            },
            "p.adjust"   = {
                colorRamp2(c(0, pvalThresh, 1),
                           c("firebrick", "white", "white"))
            },
            "FC"         = {
                colorRamp2(c(0, max(dat)), c("white", "firebrick"))
            },
            "log2FC"     = {
                colorRamp2(c(-max(abs(dat)), 0, max(abs(dat))),
                           c("blue", "white", "firebrick"))
            }
        )
        
        ht_opt(DENDROGRAM_PADDING = unit(0.1, "cm"))
        
        suppressWarnings(
            Heatmap(
                dat,
                col = colors,
                name = matrixType,
                cluster_columns = hcCol,
                show_column_dend = FALSE,
                cluster_rows = hcPlot,
                row_names_side = "left",
                column_names_rot = 30,
                rect_gp = gpar(col = "gray50", lwd = 0.5),
                width =  ncol(dat) * 5,
                height = nrow(dat) * 5,
                heatmap_legend_param = list(direction = "horizontal"),
                border = TRUE,
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = 10)
            )
        )
        
    }
)


# ---- Get a particular enrichment result ----
#
#' @title Get a particular enrichment result
#'
#' @param object a RflomicsSE object or a RflomicsMAE object. Should have at
#' least results from an enrichment, on differential analysis or co expression
#' clusters.
#' @param contrastName the contrast or cluster name on which the enrichment
#' was perform.
#' @param experiment if the object is a RflomicsMAE, then experiment is the
#' name of the RflomicsSE to look for.
#' @param from either diffexp or coexp, indicates where to search for the
#' results
#' @param database the database used for the enrichment (GO, KEGG or custom)
#' @param domain the subonology or subdomain for the database (eg CC, MF or
#' BP for GO.)
#' @return enrichment results given in the form of lists of clusterprofiler
#' results.
#' @exportMethod getEnrichRes
#' @rdname getEnrichRes
#' @aliases getEnrichRes
#' @examples
#' # Generate RflomicsMAE for example
#' MAEtest <- generateExample(annotation = FALSE, integration = FALSE)
#' # Run KEGG annotation on differential analysis results
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(organism = "ath",
#'                                          keyType = "kegg",
#'                                          pvalueCutoff = 0.1),
#'                         from = "DiffExp", database = "KEGG")
#' # Get all results from KEGG on differential expression lists:
#' getEnrichRes(MAEtest[["protetest"]],
#'               from = "diffexp", database = "KEGG")
setMethod(
    f = "getEnrichRes",
    signature = "RflomicsSE",
    definition = function(object,
                          contrastName = NULL,
                          from = "DiffExpEnrichAnal",
                          database = "GO",
                          domain = NULL) {
        from <- .determineFromEnrich(from)
        
        res_return <- .getEnrichResIntSE(
            object,
            contrastName = contrastName,
            from = from,
            database = database
        )
        
        if (!is.null(domain) && !is.null(contrastName)) {
            return(res_return[[domain]])
        } else {
            return(res_return)
        }
    }
)

#' @rdname getEnrichRes
#' @aliases getEnrichRes
#' @exportMethod getEnrichRes
setMethod(
    f = "getEnrichRes",
    signature = "RflomicsMAE",
    definition = function(object,
                          experiment,
                          contrastName = NULL,
                          from = "DiffExpEnrichAnal",
                          database = "GO",
                          domain = NULL) {
        if (missing(experiment)) {
            stop("Please indicate from which data you want to extract
           the enrichment results.")
        }
        
        from <- .determineFromEnrich(from)
        
        res_return <- .getEnrichResIntSE(
            object[[experiment]],
            contrastName = contrastName,
            from = from,
            database = database
        )
        
        if (!is.null(domain) && !is.null(contrastName)) {
            return(res_return[[domain]])
        } else {
            return(res_return)
        }
    }
)

# ---- Get summary from ORA : ----

#' @title Get summary tables from ORA analyses -
#' once an enrichment has been conducted.
#'
#' @param object a SE object (produced by Flomics)
#' @param database either NULL, GO, KEGG or custom.
#' if NULL, all tables are returned in a list.
#' @param from either DiffExpEnrichAnal or CoExpAnal.
#' @param contrastName the contrastName or clusterName to retrieve 
#' the results from. If NULL, all results are returned.
#' @return a list of tables or a table
#' @exportMethod sumORA
#' @rdname sumORA
#' @aliases sumORA
#' @examples
#' # Generate RflomicsMAE for example
#' MAEtest <- generateExample(annotation = FALSE, integration = FALSE)
#' # Run KEGG annotation on differential analysis results
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(organism = "ath",
#'                                          keyType = "kegg",
#'                                          pvalueCutoff = 0.1),
#'                         from = "DiffExp", database = "KEGG")
#' # Search for the pvalue cutoff:
#' sumORA(MAEtest[["protetest"]], from = "diffexp", database = "KEGG")
#'
setMethod(
    f = "sumORA",
    signature = "RflomicsSE",
    definition = function(object,
                          from = "DiffExpEnrichAnal",
                          database = NULL,
                          contrastName = NULL) {
        from <- .determineFromEnrich(from)
        
        # cat("|From: ", from, "\n")
        
        if (!is.null(contrastName)) {
            if (.isTagName(object, contrastName)) {
                contrastName <- .convertTagToContrast(object, contrastName)
            }
        }
        
        if (!is.null(database)) {
            toReturn <- object@metadata[[from]][[database]]$summary
            if (!is.null(contrastName)) {
                toReturn <- toReturn[which(toReturn$contrastName == contrastName),]
            }
            if (!is.null(toReturn) && from == "CoExpEnrichAnal") {
                colnames(toReturn)[1] <- "Cluster"
            }
            return(toReturn)
        } else {
            list_res <- lapply(
                names(object@metadata[[from]]),
                FUN = function(ontres) {
                    interRes <- object@metadata[[from]][[ontres]]$summary
                    if (!is.null(contrastName)) {
                        interRes <- interRes[which(interRes$contrastName == contrastName),]
                    }
                    interRes
                }
            )
            names(list_res) <- names(object@metadata[[from]])
            return(list_res)
        }
    }
)



# ---- Get a pvalue threshold used in enrichment analysis ----
#
#' @title Get the pvalue threshold used in enrichment analysis
#'
#' @param object a RflomicsSE object
#' @param from where to search for the results (either coexp or diffExp)
#' @param database which database (GO, KEGG, custom...)
#' @return the pvalue cutoff used for the analysis.
#' @exportMethod getEnrichPvalue
#' @rdname getEnrichPvalue
#' @aliases getEnrichPvalue
#' @examples
#' # Generate RflomicsMAE for example
#' MAEtest <- generateExample(annotation = FALSE, integration = FALSE)
#' # Run KEGG annotation on differential analysis results
#' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#'                         list_args = list(organism = "ath",
#'                                          keyType = "kegg",
#'                                          pvalueCutoff = 0.1),
#'                         from = "DiffExp", database = "KEGG")
#' # Search for the pvalue cutoff:
#' getEnrichPvalue(MAEtest[["protetest"]], from = "diffexp", database = "KEGG")
#'
setMethod(
    f = "getEnrichPvalue",
    signature = "RflomicsSE",
    definition = function(object,
                          from = "DiffExpEnrichAnal",
                          database = "GO") {
        from <- .determineFromEnrich(from)
        
        if (!database  %in% c("GO", "KEGG", "custom")) {
            stop(database,
                 " is not a valid value.
           Choose one of GO, KEGG or custom.")
        }
        pvalCutoff <-
            metadata(object)[[from]][[database]]$list_args$pvalueCutoff
        if (is.null(pvalCutoff)) {
            stop("P-value not found (returns NULL)")
        }
        
        return(pvalCutoff)
        
    }
)


#' 
#' # ---- Reset enrichment results ----
#' #
#' #' @title Set enrichment results to null
#' #' Users should not have to use this method. It is mainly planned for the
#' #' interface code.
#' #' @param object a RflomicsSE object
#' #' @param from where to search for the results (either coexp or diffExp)
#' #' @param database which database (GO, KEGG, custom...)
#' #' @return an RflomicsSE object with NULL value in the
#' #' enrichment results.
#' #' @exportMethod setEnrichNull
#' #' @rdname setEnrichNull
#' #' @aliases setEnrichNull
#' #' @examples
#' #' # Generate RflomicsMAE for example
#' #' MAEtest <- generateExample(annotation = FALSE, integration = FALSE)
#' #' # Run KEGG annotation on differential analysis results
#' #' MAEtest[["protetest"]] <- runAnnotationEnrichment(MAEtest[["protetest"]],
#' #'                         list_args = list(organism = "ath",
#' #'                                          keyType = "kegg",
#' #'                                          pvalueCutoff = 0.1),
#' #'                         from = "DiffExp", database = "KEGG")
#' #' # Search for the pvalue cutoff:
#' #' setEnrichNull(MAEtest[["protetest"]], from = "diffexp", database = "KEGG")
#' #'
#' setMethod(
#'     f = "setEnrichNull",
#'     signature = "RflomicsSE",
#'     definition = function(object,
#'                           from = "DiffExp",
#'                           database = NULL) {
#'         from <- .determineFromEnrich(from)
#'         
#'         if (is.null(database)) {
#'             metadata(object)[[from]] <- NULL
#'         } else {
#'             metadata(object)[[from]][[database]] <- NULL
#'         }
#'         
#'         return(object)
#'         
#'     }
#' )