#' @import methods

# ---- generateReport ----
#' @title Generate RFLOMICS rmarkdown report
#' @description
#' This function is used to generate a html report from a
#' \link{RflomicsMAE} object or archive with results.
#' @param object a \link{RflomicsMAE} produced by RFLOMICS.
#' @param fileName Name of the html report (default: date()_projectName.html).
#' @param archiveName name of archive with all analysis results 
#' (default: date()_projectName.tar.gz).
#' @param export boolean value to create archive (default: FALSE)
#' @param tmpDir temporary directory (default: working directory)
#' @param ... other arguments to pass into the render function.
#' @return An html report or archive (tar.gz)
#' @importFrom rmarkdown render
#' @exportMethod generateReport
#' @rdname generateReport
#' @aliases generateReport
#' @example inst/examples/generateReport.R
#'
setMethod(
  f          = "generateReport",
  signature  = "RflomicsMAE",
  definition = function(object,
                        fileName = NULL,
                        archiveName = NULL,
                        export = FALSE,
                        tmpDir = getwd(),
                        ...) {
    # Copy the report file to a temporary directory before processing it, in
    # case we don't have write permissions to the current working dir (which
    # can happen when deployed).
    tempReport <-
      file.path(path.package("RFLOMICS"), "/RFLOMICSapp/report.Rmd")
    
    # Check if the object is properly filled.
    ## function
    
    # project name
    projectName <- getProjectName(object)
    RDataName   <- paste0(projectName, ".MAE.RData")
    
    # tmp dir
    if (file.access(tmpDir, 2) != 0)
      stop("No writing access in ", tmpDir)
    
    tmpDir <-
      file.path(tmpDir, 
                paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName))
    # file.path(tmpDir, 
    #           paste0(projectName, "_report"))
    
    
    dir.create(tmpDir, showWarnings = FALSE)
    
    # html name
    if (is.null(fileName))
      fileName <- file.path(tmpDir,
                            paste0(format(Sys.time(), "%Y_%m_%d"), "_", 
                                   projectName, ".html"))
    
    # save FE rflomics.MAE in .Rdata and load it during report execution
    rflomics.MAE <- object
    save(rflomics.MAE, file = file.path(tmpDir, RDataName))
    
    # Set up parameters to pass to Rmd document
    param.list <-
      list(
        FEdata = file.path(tmpDir, RDataName),
        title  = paste0(projectName, " project"),
        outDir = tmpDir
      )
    
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the 
    # document from the code in this app).
    render(
      input             = tempReport,
      output_file       = fileName,
      params            = param.list,
      knit_root_dir     = tmpDir,
      intermediates_dir = tmpDir,
      envir = new.env(parent = globalenv()),
      ...
    )
    
    #Export results
    if (isTRUE(export)) {
      if (is.null(archiveName))
        archiveName <- file.path(dirname(tmpDir),
                                 paste0(format(Sys.time(),"%Y_%m_%d"),
                                        "_", projectName, 
                                        ".tar.gz"))
      
      
      
      # cp html in tmpDir
      file.copy(from = fileName, to = tmpDir)
      cmd <- paste0("tar -C ", dirname(tmpDir),
                    " -czf ", archiveName,  " ", basename(tmpDir))
      system(cmd)
      #message(cmd)
      
    } else{
      file.copy(from = fileName, to = dirname(tmpDir))
    }
    unlink(tmpDir, recursive = TRUE)
    
  }
)

# ---- runOmicsPCA ----
#' @title runOmicsPCA
#' @description This function performs a principal component analysis on omic 
#' data stored in an object of class \link{RflomicsSE-class}
#' Results are stored in the metadata slot of the same object. If a 
#' "Normalization" slot is present in the metadata slot, then data are 
#' normalized before running the PCA according to the indicated transform 
#' method.
#' @param object An object of class \link{RflomicsSE-class}.
#' @param ncomp Number of components to compute. Default is 5.
#' @param raw boolean. Does the pca have to be ran on raw data or transformed 
#' and normalized data? Default is FALSE, pca is ran on transformed and 
#' normalized data.
#' @return An object of class \link{RflomicsSE}
#' @exportMethod runOmicsPCA
#' @importFrom FactoMineR PCA
#' @rdname runOmicsPCA
#' @aliases runOmicsPCA
#'
setMethod(
  f          = "runOmicsPCA",
  signature  = "RflomicsSE",
  definition = function(object, ncomp = 5, raw = FALSE) {
    object2 <- .checkTransNorm(object, raw = raw)
    pseudo  <- assay(object2)
    if (raw)
      object@metadata[["PCAlist"]][["raw"]] <-
      PCA(t(pseudo), ncp = ncomp, graph = FALSE)
    else
      object@metadata[["PCAlist"]][["norm"]] <-
      PCA(t(pseudo), ncp = ncomp, graph = FALSE)
    return(object)
    
  }
)

#' @rdname runOmicsPCA
#' @aliases runOmicsPCA
#' @title runOmicsPCA
#' @param SE.name the name of the data the normalization have to be applied to.
#' @exportMethod runOmicsPCA
setMethod(
  f          = "runOmicsPCA",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        ncomp = 5,
                        raw = FALSE) {
    object[[SE.name]] <-  runOmicsPCA(object[[SE.name]],
                                      ncomp = ncomp,
                                      raw  = raw)
    
    return(object)
    
  }
)

# ---- plotOmicsPCA ----
#' @title plotOmicsPCA
#' @description This function plot the factorial map from a PCA object stored
#' in a \link{RflomicsSE-class} object. By default, samples are
#' colored by groups (all combinations of level's factor)
#' @param object An object of class \link{RflomicsSE-class}
#' @param raw This argument indicates whether the scaled PCA has to be 
#' performed on raw [\sQuote{raw}] or normalized [\sQuote{norm}] data.
#' @param axes A vector giving the two axis that have to be drawn for the 
#' factorial map
#' @param groupColor All combination of level's factor
#' @return PCA plot (ggplot2 object)
#' @importFrom dplyr mutate full_join select right_join
#' @importFrom FactoMineR coord.ellipse
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text aes xlab ylab 
#' geom_hline geom_vline geom_vline element_text ggtitle geom_polygon
#' @exportMethod plotOmicsPCA
#' @rdname plotOmicsPCA
#' @aliases plotOmicsPCA
#'
setMethod(f = "plotOmicsPCA",
          signature = "RflomicsSE",
          definition = function(object,
                                raw = c("raw", "norm"),
                                axes = c(1, 2),
                                groupColor = "groups") {
            if (length(axes) != 2) {
              stop("PCA axes must be a vector of length 2")
            }
            
            ExpDesign <- getDesignMat(object)
            
            PC1 <- paste("Dim.", axes[1], sep = "")
            PC2 <- paste("Dim.", axes[2], sep = "")
            
            if (PC1 == PC2) PC2 <- PC1 + 1
            
            score <- metadata(object)$PCAlist[[raw]]$ind$coord[, axes] %>% 
              as.data.frame() %>%
              mutate(samples = row.names(.)) %>% 
              right_join(., ExpDesign, by = "samples")
            
            var1 <- round(metadata(object)$PCAlist[[raw]]$eig[axes, 2][1], digits = 3)
            var2 <- round(metadata(object)$PCAlist[[raw]]$eig[axes, 2][2], digits = 3)
            
            omicsType <- getOmicsTypes(object)
            
            switch(raw,
                   "raw"  = {
                     title <- paste0("Raw ", omicsType, " data")
                   },
                   "norm" = {
                     title <- switch (
                       omicsType,
                       "RNAseq" = {
                         paste0(
                           "Filtred and normalized ",
                           omicsType,
                           " data (",
                           getNormSettings(object)$method,
                           ")"
                         )
                       },
                       "proteomics" = {
                         paste0(
                           "Transformed and normalized ",
                           omicsType,
                           " data (",
                           getTransSettings(object)$method,
                           " - norm: ",
                           getNormSettings(object)$method,
                           ")"
                         )
                       },
                       "metabolomics" = {
                         paste0(
                           "Transformed and normalized ",
                           omicsType,
                           " data (",
                           getTransSettings(object)$method,
                           " - norm: ",
                           getNormSettings(object)$method,
                           ")"
                         )
                       }
                     )
                   })
            
            
            
            p <- ggplot(score,
                        aes_string(x = PC1, y = PC2, color = groupColor))  +
              geom_point(size = 2) +
              geom_text(
                aes(label = samples),
                size = 2,
                vjust = "inward",
                hjust = "inward"
              ) +
              xlab(paste(PC1, " (", var1, "%)", sep =
                           "")) +
              ylab(paste(PC2, " (", var2, "%)", sep =
                           "")) +
              geom_hline(yintercept = 0,
                         linetype = "dashed",
                         color = "red") +
              geom_vline(xintercept = 0,
                         linetype = "dashed",
                         color = "red") +
              theme(
                strip.text.x =  element_text(size = 8, 
                                             face = "bold.italic"),
                strip.text.y =  element_text(size = 8, 
                                             face = "bold.italic")
              ) +
              ggtitle(title)
            
            # ellipse corr
            aa <- select(score, all_of(groupColor), PC1, PC2)
            bb <- coord.ellipse(aa, bary = TRUE)
            p <- p + geom_polygon(
              data = bb$res,
              aes_string(x = PC1, y = PC2, fill = groupColor),
              show.legend = FALSE,
              alpha = 0.1
            )
            
            return(p)
            
          })

#' @rdname plotOmicsPCA
#' @aliases plotOmicsPCA
#' @title plotOmicsPCA
#' @param SE.name the name of the data the normalization have to be applied to.
#' @exportMethod plotOmicsPCA
setMethod(
  f          = "plotOmicsPCA",
  signature  = "RflomicsMAE",
  definition = function(object,
                        SE.name,
                        raw,
                        axes = c(1, 2),
                        groupColor = "groups") {
    plotOmicsPCA(object[[SE.name]], raw, axes, groupColor)
    
  }
)

# ---- resetRflomicsMAE ----
#' resetRflomicsMAE allows for initializing the object or initializing a selection of results.
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param singleAnalyses vector of single omics analysis results names 
#' (c("DataProcessing", "PCAlist", "DiffExpAnal", 
#' "DiffExpEnrichAnal", "CoExpAnal", "CoExpEnrichAnal"))
#' @param multiAnalyses vector of multi omics analysis results names 
#' (c("IntegrationAnalysis"))
#' @param datasetNames dataset name. 
#' If dataset == NULL, all datasets will be reset
#' @return An object of class \link{RflomicsMAE}
#' @noRd
#' @keywords internal
#'
setMethod(f = "resetRflomicsMAE", 
          signature = "RflomicsMAE",
          definition = function(object,
                                singleAnalyses = NULL, 
                                multiAnalyses  = NULL,
                                datasetNames   = NULL) {
            
            all.datasets <- getDatasetNames(object)
            
            if(!is.null(datasetNames)) {
              datasetNames <- intersect(datasetNames, all.datasets)
              if(length(datasetNames) == 0) stop("")
            }
            
            if(is.null(datasetNames) && is.null(singleAnalyses) && is.null(multiAnalyses)) stop("")
            
            # if analyses is null and datasetNames not null -> remove 
            if(is.null(singleAnalyses)){
              
              if(!is.null(datasetNames)){
                
                index <- names(object) %in% datasetNames
                object <- object[,, !index]
                
              }
            }
            else{
              
              # if dataset is null we take all datasets 
              # present in RflomicsMAE object
              if (is.null(datasetNames)) {
                datasetNames <- getDatasetNames(object)
              }
              
              for (data in datasetNames) {
                
                for (analysis in singleAnalyses) {
                  
                  if(!is.null(object[[data]]@metadata[[analysis]])){
                    object[[data]]@metadata[[analysis]] <- list()
                  }
                }
              }
            }
            
            if(!is.null(multiAnalyses)){
              
              for (analysis in multiAnalyses) {
                
                if(!is.null(object@metadata[[analysis]])){
                  object@metadata[[analysis]] <- list()
                }
              }
            }
            
            return(object)
          })



# ---- getAnalyzedDatasetNames ----
#' @title getAnalyzedDatasetNames
#' @description
#' A short description...
#'
#' @param object An object of class \link{RflomicsMAE-class}
#' @param analyses vector of list of analysis name
#' @exportMethod getAnalyzedDatasetNames
#' @return a list of dataset nmaes per analysis
#' @keywords internal
#' @noRd
#' 
setMethod(
  f          = "getAnalyzedDatasetNames",
  signature  = "RflomicsMAE",
  definition = function(object, analyses = NULL) {
    
    all.analyses <- c("DataProcessing",
                      "DiffExpAnal", "DiffExpEnrichAnal", 
                      "CoExpAnal", "CoExpEnrichAnal")
    
    if(is.null(analyses)) analyses <- all.analyses
    
    df.list <- list()
    for (dataset in getDatasetNames(object)) {
      
      if(is.null(object[[dataset]])) next
      
      for(analysis in analyses){
        
        if(length(object[[dataset]]@metadata[[analysis]]) == 0)
          next
        
        switch (analysis,
                "DataProcessing" = {
                  if(isTRUE(object[[dataset]]@metadata[[analysis]]$done))
                    df.list[[analysis]] <- c(df.list[[analysis]], dataset)
                },
                "DiffExpAnal" = {
                  if(!is.null(object[[dataset]]@metadata[[analysis]]$Validcontrasts))
                    df.list[[analysis]] <- c(df.list[[analysis]], dataset)
                },
                "CoExpAnal" = {
                  if(isTRUE(object[[dataset]]@metadata[[analysis]]$results))
                    df.list[[analysis]] <- c(df.list[[analysis]], dataset)
                },
                {
                  for(db in names(object[[dataset]]@metadata[[analysis]])){
                    df.list[[analysis]][[db]] <- c(df.list[[analysis]][[db]], dataset)
                  }
                }
        )
      }
    }
    
    if(length(df.list) == 0) return(NULL)
    if(length(analyses) == 1) return(df.list[[1]])
    return(df.list)
  })


# ---- getDiffAnalysesSummary ----
#' @title getDiffAnalysesSummary
#' @description
#' A short description...
#'
#' @param object An object of class \link{RflomicsMAE-class}
#' @param plot FALSE or TRUE
#' @importFrom dplyr mutate
#' @importFrom purrr reduce
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_col geom_text
#' facet_grid scale_x_continuous labs position_stack
#' @exportMethod getDiffAnalysesSummary
#' @return a data.frame with differential analyses summary
#' @rdname getDiffAnalysesSummary
#' @aliases getDiffAnalysesSummary
setMethod(
  f          = "getDiffAnalysesSummary",
  signature  = "RflomicsMAE",
  definition = function(object, plot = FALSE, 
                        ylabelLength = 30) {
    # DataProcessing
    
    omicNames <- getAnalyzedDatasetNames(object, "DiffExpAnal")
    
    df.list <- list()
    for (dataset in omicNames) {
      
      Validcontrasts <- getValidContrasts(object[[dataset]])$contrastName
      
      if (is.null(Validcontrasts) || length(Validcontrasts) == 0){
        
        Validcontrasts <- getSelectedContrasts(object[[dataset]])$contrastName
        if (is.null(Validcontrasts) || length(Validcontrasts) == 0)
            next
      }
      
      p.adj <- getDiffSettings(object[[dataset]])$p.adj.cutoff
      logFC <- getDiffSettings(object[[dataset]])$abs.logFC.cutoff
      
      df.list[[dataset]] <-
        as.data.frame(
          object[[dataset]]@metadata$DiffExpAnal$stats)[Validcontrasts,] %>%
        mutate(dataset = dataset, contrasts = rownames(.), 
               settings = paste0("(p.adj: ", p.adj, " & logFC: ", logFC, ")"))
    }
    
    if (length(df.list) == 0)
      return(NULL)
    
    df <- reduce(df.list, rbind) |>
      melt(id = c("dataset", "settings", "contrasts", "All"), value.name = "Up_Down") |>
      filter(All != 0, !is.na(All)) |>
      mutate(percent = Up_Down / All * 100)
    
    if (isFALSE(plot))
      return(df)
    
    df$tabel <- 
      lapply(df$contrasts, function(x){
        vec <- seq(1, stringr::str_length(x), ylabelLength)
        stringr::str_sub(x, vec, vec+ylabelLength-1) %>% 
          paste(., collapse = "\n")
        }) %>% unlist()
    
    p <- ggplot(data = df, aes(y = tabel, x = percent, fill = variable)) +
      geom_col() +
      geom_text(aes(label = Up_Down), position = position_stack(vjust = 0.5)) +
      facet_wrap(.~ dataset+settings, ncol = 4) +
      scale_x_continuous(breaks = seq(0, 100, 25),
                         labels = paste0(seq(0, 100, 25), "%")) +
      labs(fill = NULL, x = "", y="")
    
    return(p)
  }
)

# ---- getAnnotAnalysesSummary ----
#' @title getAnnotAnalysesSummary
#' @description TODO
#' @param object An object of class \link{RflomicsMAE}. It is expected the SE 
#' object is produced by rflomics previous analyses, as it relies on their 
#' results.
#' @param from indicates if the enrichment results are taken from differential 
#' analysis results (DiffExpEnrichAnal) or from the co-expression analysis 
#' results (CoExpEnrichAnal)
#' @param matrixType Heatmap matrix to plot, one of GeneRatio, p.adjust 
#' or presence.
#' @param ... more arguments for ComplexHeatmap::Heatmap.
#' @return A list of heatmaps, one for each ontology/domain.
#' @importFrom reshape2 recast
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap ht_opt
#' @importFrom stringr str_wrap
#' @exportMethod getAnnotAnalysesSummary
#' @rdname getAnnotAnalysesSummary
#' @aliases getAnnotAnalysesSummary
setMethod(
  f = "getAnnotAnalysesSummary",
  signature = "RflomicsMAE",
  definition = function(object,
                        from = "DiffExpEnrichAnal",
                        matrixType = "presence",
                        ...) {
    extract.list <- list()
    
    for (data in getDatasetNames(object)) {
      if (is.null(object[[data]]))
        next
      if (length(object[[data]]@metadata[[from]]) == 0)
        next
      
      # for each database
      databases <- names(object[[data]]@metadata[[from]])
      for (database in databases) {
        if (is.null(object[[data]]@metadata[[from]][[database]]$enrichResult))
          next
        
        clusterNames <-
          names(object[[data]]@metadata[[from]][[database]]$enrichResult)
        pvalThresh   <-
          object[[data]]@metadata[[from]][[database]]$list_args$pvalueCutoff
        
        for (name in clusterNames) {
          domains <-
            names(object[[data]]@metadata[[from]][[database]]$enrichResult[[name]])
          
          for (dom in domains) {
            cprRes              <-
              object[[data]]@metadata[[from]][[database]]$enrichResult[[name]][[dom]]
            cprRes              <-
              cprRes@result[cprRes@result$p.adjust < cprRes@pvalueCutoff, ]
            cprRes$contrastName <- name
            cprRes$dataset      <- data
            
            extract.list[[database]][[dom]] <-
              rbind(extract.list[[database]][[dom]], cprRes)
            
          }
        }
      }
    }
    
    
    p.list <- list()
    
    for (database in names(extract.list)) {
      for (dom in names(extract.list[[database]])) {
        extract <- extract.list[[database]][[dom]]
        
        extract$Description <-
          str_wrap(extract$Description, width = 30)
        extract$contrastName <-
          str_wrap(extract$contrastName, width = 30)
        
        extract$GeneRatio <-
          as.numeric(vapply(
            extract$GeneRatio,
            FUN = function(x)
              eval(parse(text = x)),
            FUN.VALUE = 1
          ))
        extract$BgRatio <-
          as.numeric(vapply(
            extract$BgRatio,
            FUN = function(x)
              eval(parse(text = x)),
            FUN.VALUE = 1
          ))
        extract$FC <- extract$GeneRatio / extract$BgRatio
        
        extract$contrastNameLabel <- extract$contrastName
        extract$contrastName <-
          paste(extract$contrastName, extract$dataset, sep = "\n")
        
        split.df <- unique(extract[c("dataset", "contrastName", "contrastNameLabel")])
        split <- split.df$dataset
        names(split) <- split.df$contrastName
        
        if (nrow(extract) == 0)
          next
        
        dat <- switch(
          matrixType,
          "GeneRatio" = {
            inter <- recast(extract[, c("Description", "contrastName", "GeneRatio")],
                            Description ~ contrastName,
                            measure.var = "GeneRatio")
            rownames(inter) <-
              inter$Description
            #inter <- inter[, colnames(inter) != "Description"]
            inter <-
              select(inter, -"Description")
            inter[is.na(inter)] <- 0
            inter
          },
          "p.adjust" = {
            inter <-
              recast(extract[, c("Description", "contrastName", "p.adjust")],
                     Description ~ contrastName,
                     measure.var = "p.adjust")
            rownames(inter) <-
              inter$Description
            #inter <- inter[, colnames(inter) != "Description"]
            inter <-
              select(inter, -"Description")
            inter[is.na(inter)] <- 1
            inter
          },
          "presence" = {
            inter <- recast(extract[, c("Description", "contrastName", "p.adjust")],
                            Description ~ contrastName,
                            measure.var = "p.adjust")
            rownames(inter) <-
              inter$Description
            #inter <- inter[, colnames(inter) != "Description"]
            inter <-
              select(inter, -"Description")
            inter[!is.na(inter)] <- 1
            inter[is.na(inter)]  <- 0
            inter
          },
          "FC" = {
            inter <- recast(extract[, c("Description", "contrastName", "FC")],
                            Description ~ contrastName,
                            measure.var = "FC")
            rownames(inter) <-
              inter$Description
            #inter <- inter[, colnames(inter) != "Description"]
            inter <-
              select(inter, -"Description")
            inter[is.na(inter)] <- 0
            inter
          },
          "log2FC" = {
            inter <- recast(extract[, c("Description", "contrastName", "FC")],
                            Description ~ contrastName,
                            measure.var = "FC")
            rownames(inter) <-
              inter$Description
            #inter <- inter[, colnames(inter) != "Description"]
            inter <-
              select(inter, -"Description")
            inter <- log2(inter)
            inter[is.infinite(as.matrix(inter))] <-
              0
            # means FC is 0, shouldn't happen much...
            inter[is.na(inter)] <- 0
            # means it's not significant and not in the matrix.
            inter
          }
        )
        
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
            colorRamp2(c(-max(abs(
              dat
            )), 0, max(abs(
              dat
            ))),
            c("blue", "white", "firebrick"))
          }
        )
        
        ht_opt(DENDROGRAM_PADDING = unit(0.1, "cm"))
        
        p.list[[database]][[dom]] <- suppressWarnings(
          Heatmap(
            t(dat),
            col = colors,
            name = matrixType,
            row_split = split[names(dat)],
            #cluster_columns = hcCol,
            #cluster_rows = hcPlot,
            show_column_dend = FALSE,
            show_row_dend = FALSE,
            row_names_side = "left",
            #column_names_rot = 20,
            row_labels = split.df[which(split.df$contrastName %in% names(dat)), ]$contrastNameLabel,
            # column_names_rot = 0,
            # column_names_centered = TRUE,
            rect_gp = gpar(col = "gray80", lwd = 0.1),
            width =  ncol(dat) * 5,
            height = nrow(dat) * 5,
            heatmap_legend_param = list(direction = "horizontal"),
            border = TRUE,
            column_names_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 10)
          )
        )
        
        
      }
    }
    if (length(p.list) == 0)
      return(NULL)
    return(p.list)
  }
)

# ---- getCoExpAnalysesSummary ----
#' @title getCoExpAnalysesSummary
#' @description TODO
#' @param object An object of class \link{RflomicsMAE}. 
#' It is expected the SE object is produced by rflomics previous analyses, 
#' as it relies on their results.
#' @param omicNames the name of the experiment to summarize.
#' @param ... Not in use at the moment
#' @return A ggplot object
#' @importFrom dplyr filter mutate rename full_join group_by
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme element_text 
#' facet_grid labs vars
#' @exportMethod getCoExpAnalysesSummary
#' @rdname getCoExpAnalysesSummary
#' @aliases getCoExpAnalysesSummary
setMethod(
  f = "getCoExpAnalysesSummary",
  signature = "RflomicsMAE",
  definition = function(object, omicNames = NULL,
                        ...) {
    mean.y_profiles.list <- list()
    
    if (is.null(omicNames))
      omicNames <- getAnalyzedDatasetNames(object, "CoExpAnal")
    
    for (data in omicNames) {
      
      Groups     <- getDesignMat(object[[data]])
      cluster.nb <-
        object[[data]]@metadata$CoExpAnal$cluster.nb
      coseq.res  <-
        object[[data]]@metadata$CoExpAnal[["coseqResults"]]
      
      mean.y_profiles.list[[data]] <-
        lapply(seq_len(cluster.nb), function(cluster) {
          assays.data <- filter(as.data.frame(coseq.res@assays@data[[1]]),
                                get(paste0("Cluster_", cluster)) > 0.8)
          
          y_profiles.gg <-
            coseq.res@y_profiles[rownames(assays.data), ] %>%
            data.frame() %>%
            mutate(observations = rownames(.)) %>%
            melt(id = "observations", value.name = "y_profiles") %>%
            rename(samples = variable) %>%
            full_join(Groups , by = "samples")
          
          y_profiles.gg %>% group_by(groups) %>%
            summarise(mean = mean(y_profiles)) %>%
            mutate(cluster = paste0("cluster.", cluster))
          
        }) %>% reduce(rbind) %>% mutate(dataset = data)
      
    } %>% reduce(rbind)
    
    mean.y_profiles.gg <- reduce(mean.y_profiles.list, rbind)
    
    if (nrow(mean.y_profiles.gg) == 0)
      return(NULL)
    
    p <- ggplot(data = mean.y_profiles.gg, aes(x = groups, y = mean, group = 1)) +
      geom_line(aes(color = as.factor(cluster))) + 
      geom_point(aes(color = as.factor(cluster))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      facet_grid(cols = vars(cluster), rows = vars(dataset)) +
      labs(x = "Conditions", y = "Expression profiles mean")
    
    return(p)
    
  }
)
