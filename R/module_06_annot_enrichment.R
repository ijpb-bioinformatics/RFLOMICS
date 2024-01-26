# Step 6 : ANNOTAION

# ---- Module enrichment : main -----

#' @keywords internal
.modEnrichmentUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  htmltools::tagList(  
    fluidRow(
      box(title = span(tagList(icon("chart-pie"), " ", 
                               a("ClusterProfiler/", href = "https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html"), 
                               a("Pathview",         href = "https://bioconductor.org/packages/devel/bioc/vignettes/pathview/inst/doc/pathview.pdf"), 
                               "   " , 
                               tags$small("(Scroll down for instructions)"))),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
          div(  
            p("Analyses in this module are conducted using the clusterprofiler R-package. If you have more questions or interest in this package,
              please check the associated paper or the online vignette at https://yulab-smu.top/biomedical-knowledge-mining-book/index.html."),
            p(""),
            h4(tags$span("Parameters set up:", style = "color:orange")),
            p("Choose the lists of genes you want to run the enrichment for. Default option selects all the available lists (contrasts or co-expression clusters)."),
            p("Then choose the ontology you want to refer to for the analysis. Multiple choices are not allowed. 
            If you select custom, you'll have to enter an annotation file with at least two columns : 
              the names of the entity (same as the rownames of your dataset) and an id for an ontology term (eg: GO:0030198). 
              It can also contains a column for a more explicit name for the term (eg: extracellular matrix organization) 
              and a column for specifying the domain (eg: MF, BP or CC). 
              An enrichment analysis will be ran on each specified domain."),
            p("You will have to specify the names of the columns after validating the annotation file."),
            p("If you choose GO, the three GO:BP, GO:MF and GO:CC will be analyzed. You can chose to only analyze one of them by selecting the wanted ontology domain).
              It requires to indicate an R-package for the database, in the form of org.*db. (eg: org.At.tair.db for Arabidopsis thaliana), 
              and to specify which type of identifier is used in the data (eg: TAIR for arabidopsis)."),
            p("For KEGG analysis, only four identifiers are possible. Please check if the rownames correspond to one of them (eg: TAIR is also kegg identifiers)."),
            p("KEGG analysis uses an API to search for pathways online, it requires to have access to an internet connection."),
            p("Set the adjusted pvalue threshold. Only results below this threshold will be displayed."),
            
            h4(tags$span("Outputs:", style = "color:orange")),
            p("For each list, either contrast results or co-expression cluster, multiple table and plots are displayed."),
            p("- Overview: shows all results for all ontology domain (usefull when multiple domains) for all contrasts or all clusters. The blue line indicates the current results"),
            p("- Results table: for the current list, all terms that passed the adjusted pvalue threshold, ordered by increasing adjusted pvalue. You can change the domain at the bottom of the table.
              You can also order the table with the columns or search for a particular expression"),
            p("- Dotplot: for each domain, shows the 15 (default) first terms by adj. pvalue. You can also change the domain or search for a particular expression. When searching for the expression, 
              you might want to increase the number of terms to consider (if there is enough that passed the threshold)"),
            p("- Heatplot: for each domain, shows the 15 (default) first terms and the genes that are both part of the list and the pathway in the form of a heatmap. 
            For the contrasts lists, colors indicate the log2FC of the gene as found in the differential analysis. You can change the domain and search for a particular expression 
              (adjusting the number of terms to consider if you want to check further on the list of terms)"),
            p("- cnetplot: for each domain, shows the 15 (default) first terms and the genes that are both part of the list and the pathway in the form of a network. 
            As for the heatplot, only the contrasts list have colors, according to the log2FC of each genes.
              Default only shows the terms labels, you can turn on the genes names as well (it can be unreadable). You can also search for a particular expression."),
          )
      )
    ),
    ## parameters + input data
    fluidRow(uiOutput(ns("tabsetPanel_DB_UI")) )
  )
}

#' @importFrom data.table fread
#' @importFrom AnnotationDbi keytypes
#' @importFrom DT renderDataTable datatable
#' @keywords internal
.modEnrichment <- function(input, output, session, dataset, rea.values){
  
  ns <- session$ns
  
  local.rea.values <- reactiveValues()
  
  # UI
  output$tabsetPanel_DB_UI <- renderUI({
    
    # validate(
    #   need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis and validate your choices...")
    # )
    
    omicsType <- getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]])
    
    tabPanel.list <- list(tabPanel(title = "Custom annotation", .modEnrichmentDBUI(id = ns("custom"))))
    
    if (omicsType %in% c("RNAseq", "proteomics")) {
      
      tabPanel.list <- c(list(
        tabPanel(title = "Gene Onlology database", .modEnrichmentDBUI(id = ns("GO"))),
        tabPanel(title = "KEGG database", .modEnrichmentDBUI(id = ns("KEGG")))),
        tabPanel.list )
    }
    do.call(what = tabsetPanel, args = tabPanel.list)
  })
  
  # SERVEUR
  shiny::callModule(module = .modEnrichmentDB, id = "GO", dataset = dataset, database = "GO", rea.values)
  shiny::callModule(module = .modEnrichmentDB, id = "KEGG", dataset = dataset, database = "KEGG", rea.values)
  shiny::callModule(module = .modEnrichmentDB, id = "custom", dataset = dataset, database = "custom", rea.values)
  
}

# ---- Module enrichment : per database ----

#' @title .modEnrichmentDBUI
#' @keywords internal
#' @noRd
.modEnrichmentDBUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  htmltools::tagList(  
    br(),
    box(title = span(tagList(icon("sliders"), "  ", "Settings")),
        width = 3, status = "warning",
        uiOutput(ns("ParamDB_UI"))),
    box(title = NULL, headerBorder = FALSE,
        width = 9, status = "warning",
        uiOutput(ns("tabsetRes_UI")))
  )
  
}


#' @title .modEnrichmentDB
#' @importFrom RCurl url.exists
#' @keywords internal
#' @importFrom DT renderDataTable datatable
#' @noRd
.modEnrichmentDB <- function(input, output, session, dataset, database, rea.values){
  
  ns <- session$ns
  
  # UI
  
  ## setting
  output$ParamDB_UI <- renderUI({
    
    switch (database,
            "GO"     = { uiOutput(ns("AnnotParamGO_UI")) },
            "KEGG"   = { uiOutput(ns("AnnotParamKEGG_UI")) },
            "custom" = { uiOutput(ns("AnnotParamCustom_UI")) }
    )
    
  })
  
  ## display results of enrichment analysis for selected database
  output$tabsetRes_UI <- renderUI({
    
    ### for diff analysis results
    tabPanel.list <- list(tabPanel(title = "Differential expression lists", 
                                   br(),
                                   .modRunEnrichmentUI(id = ns("DiffExpEnrichAnal"))))
    
    ### for coexp analysis results if exist
    if(!isFALSE(rea.values[[dataset]]$coExpAnal)) 
      tabPanel.list <- c(tabPanel.list,
                         list(tabPanel(title = "Co-expression clusters",
                                       br(),
                                       .modRunEnrichmentUI(id = ns("CoExpEnrichAnal")))))
    
    do.call(what = tabsetPanel, args = tabPanel.list )
  })
  
  ## GO setting
  output$AnnotParamGO_UI <-   renderUI({
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    #if (omicsType %in% c("RNAseq", "proteomics")) {
    
    # box(title = span(tagList(icon("sliders"), "  ", "Settings")),
    #     width = 12, status = "warning",
    
    list(
      # select ontology
      pickerInput(
        inputId = ns("ont.select"), label = "Select Domain:",
        choices = c("ALL", "BP", "MF", "CC"),
        selected = "ALL",
        options = list(`actions-box` = TRUE, size = 10, 
                       `selected-text-format` = "count > 3")
      ),
      
      # select database
      pickerInput(
        inputId = ns("db.select"), label = "Select Data Base:",
        choices = c("", dir(.libPaths())[grep("org.*db", dir(.libPaths()))]),
        selected = "",
        options = list(`actions-box` = TRUE, size = 10, 
                       `selected-text-format` = "count > 3")
      ),
      
      # select keytype
      pickerInput(
        inputId = ns("keytype.go"), label = "Select id type:",
        choices = c("", "ENTREZID", "SYMBOL",  "TAIR", 
                    "ENSEMBL", "PMID", "REFSEQ", "GENENAME"),
        selected = "",
        options = list(`actions-box` = TRUE, size = 10, 
                       `selected-text-format` = "count > 3")
      ),
      #),
      
      ## alpha threshold
      numericInput(inputId = ns("pValue_GO"), 
                   label = "Adjusted p-value threshold:", 
                   value = 0.01 , min = 0, max = 1, step = 0.01)
    )
    
    #}
  })
  
  ## KEGG setting
  output$AnnotParamKEGG_UI <-   renderUI({
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    #if (omicsType %in% c("RNAseq", "proteomics")) {
    list(
        pickerInput(
          inputId = ns("keyTypeKEGG"), label = "Select id type:",
          choices = c("", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
          selected = "",
          options = list(`actions-box` = TRUE, size = 10, 
                         `selected-text-format` = "count > 3")),
        
        # type organism 
        textInput(inputId = ns("organism"), 
                  label = "Organism: (eg. ath)", value = "", 
                  width = NULL, placeholder = NULL),
        
        ## alpha threshold
        numericInput(inputId = ns("pValueKEGG"), 
                     label = "Adjusted p-value threshold:",
                     value = 0.01 , min = 0, max = 1, step = 0.01)
    )
    
    #}
  })
  
  ## custom setting
  output$AnnotParamCustom_UI <-   renderUI({
    
    if (rea.values[[dataset]]$diffValid == FALSE) return()
    
    list(
        
        fileInput(inputId = ns("annotationFileCPR"), label = "Annotation file:",
                  accept  = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),  
        
        #uiOutput(ns("useExampleFile_UI")),
        
        uiOutput(ns("selectColumnsFromCustomAnnot")),
        
        ## alpha threshold
        numericInput(inputId = ns("pValue_custom"), 
                     label = "Adjusted p-value threshold:", 
                     value = 0.01 , min = 0, max = 1, step = 0.01)
    )
  })
  observeEvent(input$annotationFileCPR, {
    local.rea.values$dataPathAnnot <- NULL
    local.rea.values$dataPathAnnot <- input$annotationFileCPR$datapath
    
    output$selectColumnsFromCustomAnnot <- renderUI({
      
      annotation <- data.table::fread(file = local.rea.values$dataPathAnnot, 
                                      sep = "\t", header = TRUE)
      
      column(width = 12,
             
             # Select the right columns for the analysis
             pickerInput(inputId = ns("col_geneName"), label = "Genes ID: *",
                         choices = c("", colnames(annotation)),
                         selected = ""),
             pickerInput(inputId = ns("col_termID"), label = "Terms IDs: *",
                         choices = c("", colnames(annotation)),
                         selected = ""),
             pickerInput(inputId = ns("col_termName"), label = "Term Names:",
                         choices = c("", colnames(annotation)),
                         selected = ""),
             pickerInput(inputId = ns("col_domain"), label = "Domain:",
                         choices = c("", colnames(annotation)),
                         selected = ""),
      )
    })
  })
  
  # SERVEUR
  shiny::callModule(module  = .modRunEnrichment, 
                    id = "DiffExpEnrichAnal", dataset = dataset,
                    database = database, 
                    listSource = "DiffExpEnrichAnal",
                    rea.values = rea.values, 
                    input2 = input)
  
  shiny::callModule(module  = .modRunEnrichment, 
                    id = "CoExpEnrichAnal", dataset = dataset,
                    database = database, 
                    listSource = "CoExpEnrichAnal",
                    rea.values = rea.values, 
                    input2 = input)
  
}

# ---- module enrichment : run per analysis type ----

#' @title .modRunEnrichmentUI
#' @keywords internal
#' @noRd
.modRunEnrichmentUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  htmltools::tagList(
    uiOutput(ns("setting")),
    uiOutput(ns("summary")),
    uiOutput(ns("AnnotResults"))
  )
}

#' @title .modRunEnrichment
#' @importFrom RCurl url.exists
#' @keywords internal
#' @importFrom DT renderDataTable datatable
#' @noRd
.modRunEnrichment <- function(input, output, session, dataset, database, 
                              listSource, rea.values, input2){
  ns <- session$ns
  
  switch(listSource,
         "DiffExpEnrichAnal" = {
           title <- "Summary - Enrichment from differential expression analysis"
           fromAnnot <- "diffAnnot"
           from <- "DiffExpAnal"
         },
         "CoExpEnrichAnal" = {
           title <- "Summary - Enrichment from  co-expression analysis"
           fromAnnot <- "coExpAnnot"
           from <- "CoExpAnal"
         }
  )
  
  # UI
  
  ## setting
  output$setting <- renderUI({
    
    ### lists
    data.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    ListNames <- switch (listSource,
                         "DiffExpEnrichAnal" = getValidContrasts(data.SE)$contrastName,
                         "CoExpEnrichAnal"   = names(.getCoseqClusters(data.SE)))
    
    tagList(
      column(width = 6,
             pickerInput(
               inputId  = ns("listToAnnot"),
               label    = "Select lists:",
               choices  = ListNames,
               multiple = TRUE,
               selected = ListNames,
               options  = list(`actions-box` = TRUE, size = 10,
                               `selected-text-format` = "count > 3")
             )
      ),
      #tag$br(),
      column(width = 6, 
             actionButton(inputId = ns("run"), label = "Run ORA"))
    )
  })
  
  ## display results
  output$summary <- renderUI({

    # if (rea.values[[dataset]][[fromAnnot]] == FALSE ||
    #     is.null(sumORA(session$userData$FlomicsMultiAssay[[dataset]],
    #                    from = listSource,
    #                    ont = database)))
    if (rea.values[[dataset]][[fromAnnot]] == FALSE || 
        is.null(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]]))
      return()

    results <- getEnrichRes(session$userData$FlomicsMultiAssay[[dataset]], listSource, database)
    
    # session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]]
    if (is.null(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]][["summary"]])) {

      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,
            collapsed = TRUE, status = "warning", title = title,
            "There is no results for enrichment analysis! Check geneID"))
    } else {
      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,
            collapsed = TRUE, status = "primary", title = title,
            renderDataTable({datatable(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]][["summary"]],
                                       rownames = FALSE,
                                       options = list(pageLength = 6, dom = 'tip'))
            })))
    }
  })
  output$AnnotResults <- renderUI({

    dataSE <- session$userData$FlomicsMultiAssay[[dataset]]

    if (rea.values[[dataset]][[fromAnnot]] == FALSE ||
        is.null(sumORA(dataSE, from = listSource, ont = database))) return()

    #foreach genes list selected (contrast)
    lapply(names(getEnrichRes(dataSE, from = listSource, ont = database)), function(listname) {
      if (length(getEnrichRes(dataSE, from = listSource, ont = database, contrast = listname)) != 0) {
        if (sum(unlist(sumORA(dataSE, from = listSource, ont = database)[-1]), na.rm = TRUE) == 0) {

          fluidRow(
            box(width = 12, title = paste0(listname, ": 0 enriched terms found"), status = "danger")
          )
        }
        else{

          choices <- names(getEnrichRes(dataSE, from = listSource,
                                        contrast = listname,
                                        ont = database))

          tabPanel.list <- list(
            # ---- Tab Panel : dotPlot : ----
            tabPanel("DotPlot",
                     br(),
                     renderUI({
                       outdot <- .doNotSpeak(plotCPR(dataSE,
                                                     contrast = listname,
                                                     from = from,
                                                     plotType = "dotplot",
                                                     ontology = database,
                                                     domain = input[[paste0(listname, "-domain")]],
                                                     showCategory = input[[paste0(listname, "-top.over")]],
                                                     searchExpr = input[[paste0(listname, "-grep")]])
                       )

                       if (is(outdot, "gg")) renderPlot({
                         warnOpt <- getOption("warn")
                         options(warn = -1) # avoid ggrepel warnings, or at least trying to
                         suppressMessages(suppressWarnings(print(outdot)))
                         options(warn = warnOpt)
                       })
                       else renderText({outdot$message})
                     }),
            ),
            # ---- Tab Panel : heatplot : ----
            tabPanel("Heatplot",
                     br(),
                     renderUI({
                       outheat <- .doNotSpeak(plotCPR(dataSE,
                                                      contrast = listname,
                                                      from = from,
                                                      plotType = "heatplot",
                                                      ontology = database,
                                                      domain = input[[paste0(listname, "-domain")]],
                                                      showCategory = input[[paste0(listname, "-top.over")]],
                                                      searchExpr = input[[paste0(listname, "-grep")]]))


                       if (is(outheat, "gg")) renderPlot({
                         warnOpt <- getOption("warn")
                         options(warn = -1) # avoid ggrepel warnings, or at least trying to
                         suppressMessages(suppressWarnings(print(outheat)))
                         options(warn = warnOpt)
                       })
                       else renderText({outheat$message})
                     })
            ),
            # ---- Tab Panel : cnetplot : ----
            tabPanel("Cnetplot",
                     br(),
                     verticalLayout(
                       renderUI({
                         nodeLabelArg <- "none"
                         if (input[[paste0(listname, "-genesLabels_cnet")]] &&
                             input[[paste0(listname, "-termsLabels_cnet")]]) {
                           nodeLabelArg <- "all"
                         }else if (input[[paste0(listname, "-genesLabels_cnet")]] &&
                                   !input[[paste0(listname, "-termsLabels_cnet")]]) {
                           nodeLabelArg <- "gene"
                         }else if (input[[paste0(listname, "-termsLabels_cnet")]] &&
                                   !input[[paste0(listname, "-genesLabels_cnet")]]) {
                           nodeLabelArg <- "category"
                         }

                         outcnet <- .doNotSpeak(plotCPR(dataSE,
                                                        contrast = listname,
                                                        from = from,
                                                        plotType = "cnetplot",
                                                        ontology = database,
                                                        domain = input[[paste0(listname, "-domain")]],
                                                        showCategory = input[[paste0(listname, "-top.over")]],
                                                        searchExpr = input[[paste0(listname, "-grep")]],
                                                        nodeLabel = nodeLabelArg))

                         if (is(outcnet, "gg")) renderPlot({
                           warnOpt <- getOption("warn")
                           options(warn = -1) # avoid ggrepel warnings, or at least trying to
                           suppressMessages(suppressWarnings(print(outcnet)))
                           options(warn = warnOpt)
                         })
                         else renderText({outcnet$message})

                       }),
                       fluidRow(

                         column(2, checkboxInput(inputId = ns(paste0(listname, "-genesLabels_cnet")),
                                                 label = "Genes Labels", value = FALSE)),
                         column(2, checkboxInput(inputId = ns(paste0(listname, "-termsLabels_cnet")),
                                                 label = "Terms Labels", value = TRUE))

                       )
                     )
            ),
            # ---- Tab Panel : Results table : ----
            tabPanel("Result Table",
                     br(),
                     verticalLayout(

                       tags$br(),
                       DT::renderDataTable({

                         dataPlot <- getEnrichRes(object = dataSE,
                                                  from = listSource,
                                                  contrast = listname,
                                                  ont = database,
                                                  domain = input[[paste0(listname, "-domain")]])

                         pvalue <- getEnrichPvalue(dataSE,
                                                   from = listSource, dom = database)
                         DT::datatable(dataPlot@result[dataPlot@result$p.adjust <  pvalue,],
                                       rownames = FALSE,
                                       options = list(pageLength = 5,
                                                      lengthMenu = c(5, 10, 15, 20),
                                                      scrollX = TRUE))


                       })
                     )
            )

          )# TabsetPanel

          # ---- Tab Panel : only for KEGG, pathview : ----
          if (database == "KEGG") {
            data <-   getEnrichRes(object = dataSE,
                                   contrast = listname, from = listSource,
                                   ont = "KEGG")[["no-domain"]]@result

            pvalue <- getEnrichPvalue(dataSE,
                                      from = listSource,
                                      dom = database)
            mapChoices <- sort(data$ID[data$p.adjust < pvalue])

            tabPanel.list <- c(tabPanel.list,
                               list(
                                 tabPanel("Pathview results",
                                          br(),
                                          fluidRow(
                                            column(3,

                                                   selectInput(
                                                     inputId = ns(paste0(listname, "-MAP.sel")), label = "Select map:",
                                                     choices = mapChoices, multiple = FALSE, selectize = FALSE,
                                                     size = 5,
                                                     selected = mapChoices[1]
                                                   ),
                                            ),
                                            column(9,
                                                   h5("Link to interactive map online"),
                                                   renderPrint({
                                                     link_to_map <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
                                                                           input[[paste0(listname, "-MAP.sel")]], "/",
                                                                           data[input[[paste0(listname, "-MAP.sel")]], "geneID"])
                                                     paste0(link_to_map)
                                                   }),
                                            )
                                          ),
                                          fluidRow(
                                            column(12,
                                                   renderUI({
                                                     link_to_map <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
                                                                           input[[paste0(listname, "-MAP.sel")]], "/",
                                                                           data[input[[paste0(listname, "-MAP.sel")]], "geneID"])

                                                     # test validity of URL
                                                     if (RCurl::url.exists(link_to_map)) {

                                                       renderPlot({
                                                         plotCPRKEGG(object = dataSE,
                                                                     contrast = listname,
                                                                     from = listSource,
                                                                     pathway_id = input[[paste0(listname, "-MAP.sel")]],
                                                                     species = input2$organism,
                                                                     gene_idtype = input2$keyTypeKEGG
                                                         )}, res = 300, width = 1000, height = 1000)
                                                     } else {
                                                       renderText("Please check your connection. It seems the URL does not exist, or you're not connected.")
                                                     }

                                                   })
                                            )
                                          )
                                 )
                               )
            )

          }
          # display results
          fluidRow(

            box(width = 12, solidHeader = TRUE, collapsible = TRUE,
                collapsed = TRUE, status = "success", title = listname,

                fluidRow(
                  column(4,
                         radioButtons(inputId = ns(paste0(listname, "-domain")),
                                      label = "Domain",
                                      choices = choices,
                                      selected = choices[1])),
                  column(4,
                         numericInput(inputId = ns(paste0(listname, "-top.over")),
                                      label = "Top terms:" , value=15 ,
                                      min = 1, max = 10000, step = 5)), # max code en dur : pas bien
                  column(4,
                         textInput(inputId = ns(paste0(listname, "-grep")),
                                   label = "Search Expression"))
                ),
                fluidRow(
                  column(width = 12,
                         do.call(what = tabsetPanel, args = tabPanel.list)
                  )
                )

            ) # box
          ) # fluidrow
        } # else
      } # if
    })# lapply
  })

  # SERVER

  ## run enrichment
  # ---- run Annotation  ----
  observeEvent(input$run, {

    # check list of genes
    if (length(input$listToAnnot) == 0) {
      showModal(modalDialog( title = "Error message", "Please select at least 1 variable list."))
    }
    validate({
      need(length(input$listToAnnot) != 0,
           message = "Please select at least 1 variable list")
    })

    if(database == "GO"){

      # check parameters
      if (input2$db.select == "") {
        showModal(modalDialog(title = "Error message", "Select databases..."))
      }
      validate({
        need(input2$db.select != "", "Select databases...")
      })

      library(input2$db.select, character.only = TRUE)
      accepted_keytypes <- AnnotationDbi::keytypes(get(input2$db.select))
      if (!input2$keytype.go %in% accepted_keytypes) {
        showModal(modalDialog(title = paste("Keytype must be one of: ",
                                            paste(accepted_keytypes, collapse = ", "), sep = " ")))
      }
      validate({
        need(input2$keytype.go %in% accepted_keytypes,
             message = paste("Keytype must be one of: ",
                             paste(accepted_keytypes, collapse = ", "), sep = " "))
      })

      # set list args
      database <- "GO"
      list_args <- list()
      list_args[["OrgDb"]]        <- input2$db.select
      list_args[["keyType"]]      <- input2$keytype.go
      list_args[["pvalueCutoff"]] <- input2$pValue_GO

      domain <- input2$ont.select
      if (input2$ont.select == "ALL") domain <- c("MF", "BP", "CC")

      paramList <- list(#method       = "ORA",
                        ontology     = database,
                        list_args    = list_args,
                        #OrgDb        = input2$db.select,
                        #keyType      = input2$keytype.go,
                        domain       = domain#,
                        #pvalueCutoff = input2$pValue_GO
      )
    }
    if(database == "KEGG"){

      # check param
      ## ID key
      if (!input2$keyTypeKEGG %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")) {
        showModal(modalDialog( title = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot"))
      }
      validate({
        need(input2$keyTypeKEGG %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
             message = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot")
      })

      ## code organism KEGG
      if (input2$organism == "") {
        showModal(modalDialog( title = "Error message", "set organism kegg code."))
      }
      validate({
        need(input2$organism != "", message = "Set organism kegg code.")
      })


      # set list args
      database <- "KEGG"
      list_args <- list()
      list_args[["organism"]] <- input2$organism
      list_args[["keyType"]]  <- input2$keyTypeKEGG
      list_args[["pvalueCutoff"]] <- input2$pValueKEGG

      # prevent multiple execution
      paramList <- list(#method       = "ORA",
                        list_args    = list_args,
                        ontology     = database,
                        domain       = "no-domain"#,
                        #organism     = input2$organism,
                        #keyType      = input2$keyTypeKEGG,
                        #pvalueCutoff = input2$pValueKEGG
      )

      # input2[["organism"]]    <- input2$organism
      # input2[["keyTypeKEGG"]] <- input2$keyTypeKEGG

    }
    if(database == "custom"){

      # check param
      ## if custom annotation file
      if (is.null(input2$dataPathAnnot)) {
        showModal(modalDialog(title = "Error message", "load custom annotation file."))
      }
      validate({
        need(!is.null(input2$dataPathAnnot),
             message = "load custom annotation file.")
      })
      # Check some custom elements and annotation file
      if (input2$col_geneName == "" || input2$col_termID == "") {
        showModal(modalDialog( title = "Error message",
                               "Please choose columns names for the gene names/ID and the ontology terms ID!"))
      }
      validate({
        need(input2$col_geneName != "" && input2$col_termID != "",
             message = "Please choose columns names for the gene names/ID and the ontology terms ID!")
      })

      # set param
      domain <- NULL
      annotation2 <- NULL
      col_domain_arg <- NULL
      list_args <- list()
      list_args[["pvalueCutoff"]] <- input2$pValue_custom

      annotation <- data.table::fread(file = input2$dataPathAnnot,
                                      sep = "\t", header = TRUE)
      listCharRm <- c(".", " ", "")
      annotation2 <- list()

      # if column with domain exist
      if (input2$col_domain != "") {

        # remove column with domain == c(".", " ", "")
        domainRm <- intersect(annotation[[input2$col_domain]], listCharRm)
        if (length(domainRm) != 0) {
          annotation <- dplyr::filter(annotation,
                                      !get(input2$col_domain) %in% listCharRm)
        }

        annotation2[["domain"]] <- annotation[[input2$col_domain]]
        domain <- unique(annotation2$domain)
        col_domain_arg <- "domain"
      }

      annotation2[["gene"]] <- annotation[[input2$col_geneName]]
      annotation2[["term"]] <- annotation[[input2$col_termID]]
      if (input2$col_termName != "") {
        annotation2[["name"]] <- annotation[[input2$col_termName]]
      }
      annotation2 <- data.frame(annotation2)

      ##
      paramList <- list(method       = "ORA",
                        list_args    = list_args,
                        ontology     = "custom",
                        domain       = domain,
                        annot        = annotation2,
                        col_domain   = col_domain_arg)
    }

    # prevent multiple execution
    #if (.checkRunORAExecution(session$userData$FlomicsMultiAssay[[dataset]], database, paramList) == FALSE) return()

    #input2[[database]] <- FALSE

    # ---- Annotation on diff results: ----
    # run analysis
    print(paste("# 11- ", database, " Enrichment Analysis of ", listSource ," lists...", dataset))

    session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]] <- NULL
    rea.values[[dataset]][[fromAnnot]] <- FALSE

    # run annotation diff analysis
    runRes <- tryCatch({
      paramList <- c(paramList, list(object = session$userData$FlomicsMultiAssay[[dataset]], 
                                     from = from, nameList = input$listToAnnot))
      
      paramList <<- paramList
      
      do.call(runAnnotationEnrichment, paramList)
      
      # Something went wrong: unused arguments (method = "ORA", organism = "ath", keyType = "kegg", pvalueCutoff = 0, listToAnnot = c("(temperatureLow - temperatureElevated) in mean", "(temperatureMedium - temperatureElevated) in mean", "(temperatureMedium - temperatureLow) in mean", "(imbibitionDS - imbibitionLI) in mean", "(imbibitionEI - imbibitionLI) in mean", "(imbibitionEI - imbibitionDS) in mean")
      
      # runAnnotationEnrichment(object = session$userData$FlomicsMultiAssay[[dataset]],
      #                         nameList = input$listToAnnot,
      #                         list_args = list_args,
      #                         from = from,
      #                         annot = annotation2,
      #                         ontology = database,
      #                         domain = domain,
      #                         col_domain = col_domain_arg)
    },
    warning = function(war) return(war),
    error   = function(err) return(err)
    )

    if (!is(runRes, "SummarizedExperiment")) {
      showModal(modalDialog(title = paste("Something went wrong: ", runRes$message)))
    }
    validate({need(is(runRes, "SummarizedExperiment"),
                   message = paste0("Something went wrong: ", runRes$message))})

    session$userData$FlomicsMultiAssay[[dataset]] <-  runRes
    rea.values[[dataset]][[fromAnnot]] <- TRUE
    runRes <<- runRes
  }, ignoreInit = TRUE)
}

# ---- Comparison between enrichment results ----

#' @title .modCompEnrichmentUI
#' @keywords internal
#' @noRd
.modCompEnrichmentUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  htmltools::tagList(  
    fluidRow(
      column(width = 12,
             uiOutput(ns("CompResults"))
      )
    ))
}

#' @title .modCompEnrichment
#' @importFrom RCurl url.exists
#' @importFrom grid unit
#' @keywords internal
#' @noRd
.modCompEnrichment <- function(input, output, session, dataset, database, 
                               listSource, rea.values, local.rea.values){
  
  ns <- session$ns
  
  switch(listSource,
         "DiffExpEnrichAnal" = {
           title <- "Enrichment from differential expression analysis"
           fromAnnot <- "diffAnnot"
           from <- "DiffExp"
         },
         "CoExpEnrichAnal" = {
           title <- "Enrichment from  co-expression analysis"
           fromAnnot <- "coExpAnnot"
           from <- "CoExp"
         }
  )
  
  output$CompResults <- renderUI({
    
    if (rea.values[[dataset]][[fromAnnot]] == FALSE || 
        is.null(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]][["summary"]])) return()
    
    # Possible domains of ontology:
    possDomain <-  unique(unlist(
      lapply(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]][["enrichResult"]], 
             FUN = function(x) names(x)))) 
    
    # display results
    verticalLayout(
      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE, 
            collapsed = TRUE, status = "success", title = database,
            
            fluidRow(
              column(4,
                     radioButtons(inputId = ns(paste0(database, "-compDomain")), 
                                  label = "Domain",
                                  choices = possDomain, selected = possDomain[1], 
                                  inline = TRUE)),
              column(4,
                     radioButtons(inputId = ns(paste0(database, "-compType")), 
                                  label = "Matrix Type",
                                  choices = c("presence", "FC", "log2FC"), 
                                  selected = "FC", inline = TRUE)),
            ),
            
            fluidRow(
              column(12,
                     renderUI({
                       
                       outHeatmap <- tryCatch({
                         outHeatmap <-  plotEnrichComp(session$userData$FlomicsMultiAssay[[dataset]],
                                                       from = from, ontology = database, 
                                                       domain = input[[paste0(database, "-compDomain")]],
                                                       matrixType = input[[paste0(database, "-compType")]])
                       },
                       warning = function(warn) warn,
                       error = function(err) err
                       )
                       
                       if (is(outHeatmap, "Heatmap")) renderPlot({ComplexHeatmap::draw(outHeatmap, 
                                                                                       heatmap_legend_side = "top", 
                                                                                       padding = unit(5, "mm"),
                                                                                       gap = unit(2, "mm"))}, 
                                                                 width = "auto", 
                                                                 height = min(400 + nrow(outHeatmap@matrix)*50, 1000))
                       else renderText({outHeatmap$message})
                     })
              ) # column
            ),  # fluidrow
        ), # box
      ) # fluidrow
    )
  })
}



# ---- check run enrichement analysis execution ----
#' @title .checkRunORAExecution
#' @noRd
#' @keywords internal
#' @importFrom dplyr setequal
.checkRunORAExecution <- function(object.SE, dbType = NULL, paramList = NULL){
  
  
  if (!is.null(paramList$diffList)) {
    
    if (is.null(object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]])) return(TRUE)
    if (isFALSE(dplyr::setequal(names(object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$enrichResult), paramList$diffList)) ) return(TRUE)
    
    # check param
    if (isFALSE(dplyr::setequal(paramList$domain, object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$domain))) return(TRUE)
    if (paramList$pvalueCutoff != object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$pvalueCutoff) return(TRUE)
    
    switch(dbType,
           "GO" = {
             if (paramList$OrgDb    != object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$OrgDb)    return(TRUE)
             if (paramList$keyType  != object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$keyType)  return(TRUE)
           },
           "KEGG" = {
             if (paramList$keyType  != object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$keyType)  return(TRUE)
             if (paramList$organism != object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$organism) return(TRUE)
           }
    )
    
  }
  if (!is.null(paramList$CoexpList)) {
    
    if (is.null(object.SE@metadata[["CoExpEnrichAnal"]][[dbType]])) return(TRUE)
    if (isFALSE(dplyr::setequal(names(object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$enrichResult), paramList$CoexpList)) ) return(TRUE)
    
    # check param
    if (isFALSE(dplyr::setequal(paramList$domain, object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$domain))) return(TRUE)
    if (paramList$pvalueCutoff != object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$pvalueCutoff) return(TRUE)
    
    switch(dbType,
           "GO" = {
             if (paramList$OrgDb    != object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$OrgDb)   return(TRUE)
             if (paramList$keyType  != object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$keyType) return(TRUE)
           },
           "KEGG" = {
             if (paramList$keyType  != object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$keyType)  return(TRUE)
             if (paramList$organism != object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$organism) return(TRUE)
           }
    )
    
  }
  
  return(FALSE)
}




