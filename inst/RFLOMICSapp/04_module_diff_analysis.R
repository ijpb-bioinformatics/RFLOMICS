

DiffExpAnalysisUI <- function(id){

  #name space for id
  ns <- NS(id)

  tagList(
    fluidRow(
      uiOutput(ns("instruction"))),

    ### parametres for Diff Analysis
    fluidRow(
      column(3,uiOutput(ns("DiffParamUI")),
               uiOutput(ns("FilterPvalueUI"))),
      column(9,uiOutput(ns("ContrastsResults")),
               tags$br(),
               uiOutput(ns("ResultsMerge"))))

  )
}



DiffExpAnalysis <- function(input, output, session, dataset, rea.values){

  local.rea.values <- reactiveValues(dataset.SE = NULL)

  # list of tools for diff analysis
  MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")

  method <- switch (rea.values[[dataset]]$omicsType,
                    "RNAseq"       = MethodList[1],
                    "proteomics"   = MethodList[2],
                    "metabolomics" = MethodList[2])

  output$instruction <- renderUI({
    box(title = span(tagList(icon("cogs"), "  ",  a(names(method), href="https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf"), "    (Scroll down for instructions)"  )),
        solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
        p("Differential expression analysis is conducted for each hypothesis. There is just two options to set (the ajusted-pvalue cut-off and the |FC| cut-off).
          The results will appear in blocks (one per hypothesis) with 3 outputs:"),
        p("- the distribution of pvalue's : which has to be validated", a("(some help to identify the good shapes)", href="Pvalue_distrib.pdf"),""),
        p("- the MA plot (DE genes in red will varie with the p-value cutoff)"),
        p("- the table of statistics per gene/protein/metabolite (Number of stats displayed will varie with the p-value cutoff)")
    )
  })

  output$DiffParamUI <- renderUI({

    #we must run process before
    validate(
      need(rea.values[[dataset]]$process != FALSE, "Please run data processing")
    )

    validate(
      need(!is.null(rea.values$Contrasts.Sel), "Please run data processing ")
    )

    # #we must select list of contrast to test
    # validate(
    #   need(rea.values$analysis != FALSE, "Please select contrast")
    # )

    #design must be complete
    validate(
      need(rea.values[[dataset]]$compCheck != FALSE, session$userData$FlomicsMultiAssay@metadata$completeCheck[["error"]])
    )

    box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",

        fluidRow(column(12,
               ## list of contrasts to test
               pickerInput(
                   inputId  = session$ns("contrastList"),
                   label    = "Selected contrast :",
                   choices  = session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                   multiple = TRUE, selected = session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName),

               # method for Diff analysis
               selectInput(inputId  = session$ns("AnaDiffMethod"), label = "Method :",
                           choices  = method,
                           selected = method),

               # use of cluster. need setting step
               materialSwitch(inputId = session$ns("clustermq"),
                              label   =  popify(actionLink("infoCluster",paste0("Cluster: (?)")),"",
                                                "If there is a huge number of contrasts, the calculation can be send to the cluster to be run in parrallel",
                                                options=list(container="body"))
                              , value = FALSE, status = "success"),

               actionButton(session$ns("runAnaDiff"),"Run"))
    ))

  })

  # Run the differential analysis for each contrast set
  # Filter
  #   -> return a dynamic user interface with a collapsible box for each contrast
  #         - Pvalue graph
  #         - MAplot
  #         - Table of the DE genes
  #   -> combine data : union or intersection
  observeEvent(input$runAnaDiff, {

    # check list of genes
    if(length(input$contrastList) == 0){

      showModal(modalDialog( title = "Error message", "Please select at least 1 hypothesis"))
    }
    validate({
      need(length(input$contrastList) != 0, message="Please select at least 1 hypothesis")
    })

    local.rea.values$dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]

    print(paste("# 4- Diff Analysis...", dataset))

    rea.values[[dataset]]$diffAnal   <- FALSE
    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE

    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]] <- NULL


    if (dataset %in% rea.values$datasetDiff){
      rea.values$datasetDiff <- rea.values$datasetDiff[-which(rea.values$datasetDiff == dataset)]
      }

    local.rea.values$dataset.SE@metadata$DiffExpAnal <- list()
    local.rea.values$dataset.SE@metadata$CoExpAnal   <- list()
    local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal  <- list()
    local.rea.values$dataset.SE@metadata$CoExpEnrichAnal  <- list()

    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#



    # run diff analysis with selected method
    local.rea.values$dataset.SE <- RunDiffAnalysis(object =             local.rea.values$dataset.SE,
                                                   design =             session$userData$FlomicsMultiAssay@metadata$design,
                                                   contrastList =       session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                                                   Adj.pvalue.method =  "BH",
                                                   DiffAnalysisMethod = input$AnaDiffMethod,
                                                   clustermq =          input$clustermq)

    # error management
    if(isFALSE(local.rea.values$dataset.SE@metadata$DiffExpAnal[["results"]])){
      showModal(modalDialog( title = "Error message",
                             if(! is.null(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorStats"]])){
                               renderDataTable(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorStats"]],rownames = FALSE)
                             }
                             else{
                               as.character(local.rea.values$dataset.SE@metadata$DiffExpAnal[["Error"]])
                             }
      ))
    }

    if(is.null(local.rea.values$dataset.SE@metadata$DiffExpAnal[["RawDEFres"]])){

      showModal(modalDialog( title = "Error message",
                             if(! is.null(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorTab"]])){
                               renderDataTable(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorTab"]],rownames = FALSE)
                             }
                             else{
                               as.character(local.rea.values$dataset.SE@metadata$DiffExpAnal[["error"]])
                             }
      ))
    }

    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]] <- local.rea.values$dataset.SE

    rea.values[[dataset]]$diffAnal <- TRUE

    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#

  }, ignoreInit = TRUE)

  # filterin param
  output$FilterPvalueUI <- renderUI({

    if (rea.values[[dataset]]$diffAnal == FALSE) return()

    box(title = NULL,status = "warning", width = 14,
        numericInput(inputId = session$ns("Adj.pvalue.cutoff"),
                     label="Adjusted pvalue cutoff:",
                     value=0.05, min=0, max=1, 0.01),
        numericInput(inputId = session$ns("abs.FC.cutoff"),
                     label="|FC| cutoff:",
                     value=0, min=0, max=100, 0.1),
        actionButton(session$ns("validContrast"),"Validate"))
  })

  # display results per contrast
  output$ContrastsResults <- renderUI({

    #if (rea.values[[dataset]]$diffAnal == FALSE) return()
    if (rea.values[[dataset]]$diffAnal == FALSE ||
        is.null(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal)) return()

    ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
    local.rea.values$dataset.SE <- FilterDiffAnalysis(object = local.rea.values$dataset.SE,
                                                      Adj.pvalue.cutoff = input$Adj.pvalue.cutoff,
                                                      FC.cutoff = input$abs.FC.cutoff)
    #Contrasts.Sel <- local.rea.values$dataset.SE@metadata$DiffExpAnal$contrasts

    list(
      lapply(1:length(rea.values$Contrasts.Sel$contrast), function(i) {

      dataset.SE <- local.rea.values$dataset.SE
      vect     <- unlist(rea.values$Contrasts.Sel[i,])
      res      <- dataset.SE@metadata$DiffExpAnal[["RawDEFres"]][[vect["contrastName"]]]
      stats    <- dataset.SE@metadata$DiffExpAnal[["stats"]][[vect["contrastName"]]]

      diff.plots <- DiffAnal.plot(dataset.SE, hypothesis=vect["contrastName"],
                                  Adj.pvalue.cutoff = input$Adj.pvalue.cutoff, FC.cutoff = input$abs.FC.cutoff)

      fluidRow(
        column(10,
               box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning",
                   title = tags$h5(paste0(vect["tag"], " : ", vect["contrastName"],"  [#DE: ", stats$gDE," (up: ", stats$pgDEup,"%, ", "down: ", stats$pgDEdown,"%)]")),
                   
                   tabsetPanel(
                     
                     ### pvalue plot ###
                     tabPanel("Pvalue's distribution", renderPlot({ diff.plots$Pvalue.hist })),
                     
                     ### MAplot
                     tabPanel("MA plot", renderPlot({ suppressMessages(diff.plots$MA.plot) })),
                     
                     ### MAplot
                     tabPanel("Volcano plot", renderPlot({ suppressMessages(diff.plots$Volcano.plot) }, height = 600)),
                     
                     ### DEF result table ###
                     tabPanel("Table",
                              ### DEF result table ###
                              DT::renderDataTable({
                                resTable <- dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]]
                                resTable %>% DT::datatable(extensions = 'Buttons',
                                                           options = list(dom = 'lfrtipB',
                                                                          rownames = FALSE,
                                                                          pageLength = 10,
                                                                          buttons = c('csv', 'excel'),
                                                                          lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All")))) %>%
                                  formatStyle('logFC',
                                              backgroundColor = styleInterval(c(0, 0.01), c('royalblue', 'white', 'red2')),
                                              fontWeight = 'bold') %>% formatSignif(columns = 1:dim(resTable)[2], digits = 3)
                              }, server = FALSE)),
                     ### Heatmap ###
                     tabPanel("Heatmap",
                              renderUI({
                                renderPlot({
                                  resTable <- dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]]
                                  m.def <- assays(local.rea.values$dataset.SE)[[1]][,session$userData$FlomicsMultiAssay@metadata$design@Groups$samples]
                                  
                                  # Normalize counts (added 221123)
                                  if(dataset.SE@metadata$Normalization$methode == "TMM"){
                                    m.def <- log2(scale(m.def+1, center = FALSE,
                                                        scale = dataset.SE@metadata$Normalization$coefNorm$lib.size*dataset.SE@metadata$Normalization$coefNorm$norm.factors))
                                  }
                                  
                                  # filter by DE
                                  m.def.filter <- subset(m.def, rownames(m.def) %in% row.names(resTable))
                       
                                  # normalize count
                                  
                                  # Center
                                  # m.def.filter.center <- scale(m.def.filter,center=TRUE,scale=FALSE)
                                  m.def.filter.center <- t(scale(t(m.def.filter),center = TRUE, scale = FALSE)) # Modified 221123 : centered by genes and not by samples
                                  column_split.value <- if(input[[paste0(vect["contrastName"],"-","condColorSelect")]] != "none"){
                                    session$userData$FlomicsMultiAssay@metadata$design@Groups[,input[[paste0(vect["contrastName"],"-","condColorSelect")]]]
                                  }
                                  else{NULL}
                             
                                  # Color annotations
                                  df_annotation <- session$userData$FlomicsMultiAssay@metadata$design@Groups %>% dplyr::select(!samples & !groups)
                                  
                                  set.seed(10000) ; selectPal <- sample(rownames(RColorBrewer::brewer.pal.info),  size = ncol(df_annotation), replace = FALSE)
                                  
                                  color_list <- lapply(1:ncol(df_annotation), FUN = function(i){
                                    annot_vect <- unique(df_annotation[,i])
                                    
                                    col_vect <- RColorBrewer::brewer.pal(n = length(annot_vect), name = selectPal[i]) 
                                    names(col_vect) <- annot_vect 
                                    col_vect[!is.na(names(col_vect))] # RcolorBrewer::brewer.pal n is minimum 3, remove NA names if only 2 levels
                                  })
                                  names(color_list) <- colnames(df_annotation)
                                  
                                  column_ha <- ComplexHeatmap::HeatmapAnnotation(df = df_annotation, col = color_list)
                                  
                                  # Drawing heatmap
                                  ha <- ComplexHeatmap::Heatmap(m.def.filter.center, name = "normalized counts\nor XIC",
                                                                show_row_names= ifelse( dim(m.def.filter.center)[1] > 50, FALSE, TRUE),
                                                                row_names_gp = grid::gpar(fontsize = 8),
                                                                column_names_gp = grid::gpar(fontsize = 12),
                                                                row_title_rot = 0 ,
                                                                clustering_method_columns = "ward.D2",
                                                                cluster_column_slice=FALSE,
                                                                column_split = column_split.value,
                                                                top_annotation = column_ha)
                                  
                                  ComplexHeatmap::draw(ha, merge_legend = TRUE)
                                })
                              })
                              ,
                              renderText("Clustering method=ward.D2, center=TRUE, scale=FALSE")
                              ,
                              ## select cluster to plot
                              radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-","condColorSelect")),
                                           label = 'Levels :',
                                           choices = c("none", names(session$userData$FlomicsMultiAssay@colData)),
                                           selected = "none", inline = TRUE)

                     ),
                     ### PCA ###
                     tabPanel("PCA on DE",
                              
                              fluidRow(
                                renderPlot({
                                  
                                  resTable <- dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]]
                                  newDataset.SE <- dataset.SE[rownames(dataset.SE) %in% row.names(resTable)]
                                  newDataset.SE <- RFLOMICS::RunPCA(newDataset.SE)  

                                  RFLOMICS::plotPCA(newDataset.SE, 
                                                    PCA = "norm", 
                                                    PCs = c(as.numeric(input[["DEG_PCA_Firstaxis"]]), as.numeric(input[["DEG_PCA_Secondaxis"]])), 
                                                    condition = input$DEG_PCA_condColorSelect)
                                  
                                })
                              ),
                              fluidRow(
                                # 12/2022 : don't know how to make module common work...
                                column(width = 3, radioButtons(inputId = session$ns("DEG_PCA_condColorSelect"),
                                             label = 'Levels:',
                                             choices = c("groups", names(session$userData$FlomicsMultiAssay@colData)),
                                             selected = "groups")),
                                column(width = 3,   radioButtons(inputId  = session$ns("DEG_PCA_Firstaxis"),
                                                                 label    = "Choice of PCs:",
                                                                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                                                                 selected = 1,
                                                                 inline = TRUE)),
                                column(width = 3, 
                                       
                                       radioButtons(inputId  = session$ns("DEG_PCA_Secondaxis"),
                                                                 label    = "Choice of PCs:",
                                                                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                                                                 selected = 2,
                                                                 inline = TRUE))

                              ),
                              
                     )
                   ),
               )
        ),
        column(2,
               checkboxInput(session$ns(paste0("checkbox_", vect[["tag"]])), "OK", value = TRUE))
      )
    })
    )

  })

  # merge results on upset plot
  output$ResultsMerge <- renderUI({

    if (rea.values[[dataset]]$diffAnal == FALSE) return()

    # dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]

    index <- sapply(rea.values$Contrasts.Sel$tag, function(x){(input[[paste0("checkbox_",x)]])}) %>% unlist()

    H_selected <- rea.values$Contrasts.Sel$tag[index]

    DEF_mat <- as.data.frame(local.rea.values$dataset.SE@metadata$DiffExpAnal[["mergeDEF"]])

    rea.values[[dataset]]$DiffValidContrast <- dplyr::filter(rea.values$Contrasts.Sel, tag %in% H_selected)

    if (length(H_selected) > 1){

       box(width=12,  status = "warning",

           renderPlot({ UpSetR::upset(DEF_mat, sets = H_selected) })
       )
    }


  })

  # validate contrasts
  observeEvent(input$validContrast, {

    print(paste("# 9bis- Filter Diff Analysis...", dataset))

    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE

    # filter DEG according pvalue adj cut-off

    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal <-
      local.rea.values$dataset.SE@metadata$DiffExpAnal

    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]] <-
      rea.values[[dataset]]$DiffValidContrast

    rea.values[[dataset]]$diffValid <- TRUE
    rea.values$datasetDiff <- unique(c(rea.values$datasetDiff , dataset))


  }, ignoreInit = TRUE)


  return(input)
}



