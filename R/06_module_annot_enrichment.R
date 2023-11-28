##########################################
# Part6 : ANNOTAION
##########################################

### subModule

module_runEnrichment_UI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  htmltools::tagList(  
    fluidRow(
      column(width = 12,
             uiOutput(ns("summary")),
             uiOutput(ns("AnnotResults"))
      )
    ))
}

module_runEnrichment <- function(input, output, session, dataset, dom.select, list.source, rea.values, local.rea.values){
  
  ns <- session$ns
  
  switch (list.source,
          "DiffExpEnrichAnal" = {
            title <- "Summary - Enrichment from differential expression analysis"
            diffAnnot_coExpAnal <- "diffAnnot"
            from <- "DiffExpAnal"
          },
          "CoExpEnrichAnal" = {
            title <- "Summary - Enrichment from  co-expression analysis"
            diffAnnot_coExpAnal <- "coExpAnnot"
            from <- "CoExpAnal"
          }
  )
  
  output$summary <- renderUI({
    
    if (rea.values[[dataset]][[diffAnnot_coExpAnal]] == FALSE || 
        is.null(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[list.source]][[dom.select]])) return()
    
    results <- session$userData$FlomicsMultiAssay[[dataset]]@metadata[[list.source]][[dom.select]]
    
    if (is.null(results[["summary"]])) {
      
      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = title,
            "There is no results for enrichment analysis! Check geneID"))
    }
    else{
      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", title = title,
            DT::renderDataTable({ DT::datatable(results[["summary"]], rownames = FALSE, options = list(pageLength = 6, dom = 'tip')) })))
    }
  })
  
  output$AnnotResults<- renderUI({
    
    if (rea.values[[dataset]][[diffAnnot_coExpAnal]] == FALSE || 
        is.null(session$userData$FlomicsMultiAssay[[dataset]]@metadata[[list.source]][[dom.select]][["summary"]])) return()
    
    #foreach genes list selected (contrast)
    lapply(names(getEnrichRes(session$userData$FlomicsMultiAssay[[dataset]], from = list.source, ont = dom.select)), function(listname) {
      
      if (length(getEnrichRes(session$userData$FlomicsMultiAssay[[dataset]], from = list.source, ont = dom.select, contrast = listname)) != 0) {
        if (sum(unlist(sumORA(session$userData$FlomicsMultiAssay[[dataset]], from = list.source, ont = dom.select)[-1]), na.rm = TRUE) == 0) {
          
          fluidRow(
            box(width = 12, title = paste0(listname, ": 0 enriched terms found"), status = "danger")
          )
        }
        else{
          
          choices <- names(getEnrichRes(session$userData$FlomicsMultiAssay[[dataset]], from = list.source, contrast = listname, ont = dom.select))
          
          tabPanel.list <- list(
            # ---- Tab Panel : dotPlot : ----
            tabPanel("DotPlot",
                     br(),
                     renderUI({
                       outdot <- .doNotSpeak({plotCPR(session$userData$FlomicsMultiAssay[[dataset]],
                                                      contrast = listname, 
                                                      from = from,
                                                      type = "dotplot",
                                                      ont = dom.select,
                                                      Domain = input[[paste0(listname, "-domain")]],
                                                      showCategory = input[[paste0(listname, "-top.over")]],
                                                      searchExpr = input[[paste0(listname, "-grep")]])
                       })
                       
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
                       outheat <- .doNotSpeak({  plotCPR(session$userData$FlomicsMultiAssay[[dataset]],
                                                         contrast = listname, 
                                                         from = from,
                                                         type = "heatplot",
                                                         ont = dom.select,
                                                         Domain = input[[paste0(listname, "-domain")]],
                                                         showCategory = input[[paste0(listname, "-top.over")]],
                                                         searchExpr = input[[paste0(listname, "-grep")]])})
                       
                       
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
                         node_label_arg <- "none"
                         if (input[[paste0(listname, "-genesLabels_cnet")]] && input[[paste0(listname, "-termsLabels_cnet")]]) {
                           node_label_arg <- "all"
                         }else if (input[[paste0(listname, "-genesLabels_cnet")]] && !input[[paste0(listname, "-termsLabels_cnet")]]) {
                           node_label_arg <- "gene"
                         }else if (input[[paste0(listname, "-termsLabels_cnet")]] && !input[[paste0(listname, "-genesLabels_cnet")]]) {
                           node_label_arg <- "category"
                         }
                         
                         outcnet <- .doNotSpeak({plotCPR(session$userData$FlomicsMultiAssay[[dataset]],
                                                         contrast = listname, 
                                                         from = from,
                                                         type = "cnetplot",
                                                         ont = dom.select,
                                                         Domain = input[[paste0(listname, "-domain")]],
                                                         showCategory = input[[paste0(listname, "-top.over")]],
                                                         searchExpr = input[[paste0(listname, "-grep")]],
                                                         node_label = node_label_arg)})
                         
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
                         
                         dataPlot <-  getEnrichRes(object = session$userData$FlomicsMultiAssay[[dataset]],
                                                   from = list.source,
                                                   contrast = listname,
                                                   ont = dom.select,
                                                   domain = input[[paste0(listname, "-domain")]])
                         
                         pvalue <- getEnrichPvalue(session$userData$FlomicsMultiAssay[[dataset]], from = list.source, dom = dom.select)
                         DT::datatable(dataPlot@result[dataPlot@result$p.adjust <  pvalue,], # sometimes it prints non significant results...
                                       rownames = FALSE,
                                       options = list(pageLength = 5, lengthMenu = c(5, 10, 15, 20), scrollX = TRUE))
                         
                         
                       })
                     )
            )
            
          )# TabsetPanel
          
          # ---- Tab Panel : only for KEGG, pathview : ----
          if (dom.select == "KEGG") {
            data <-   getEnrichRes(object = session$userData$FlomicsMultiAssay[[dataset]],
                                   contrast = listname, from = list.source,
                                   ont = "KEGG")[["no-domain"]]@result
            
            pvalue <- getEnrichPvalue(session$userData$FlomicsMultiAssay[[dataset]], from = list.source, dom = dom.select)
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
                                                   
                                                   # renderUI({
                                                   #   # Crée le lien URL
                                                   #   #link_to_map <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", input[[paste0(listname, "-MAP.sel")]], "/", data[input[[paste0(listname, "-MAP.sel")]], "geneID"])
                                                   #   #a("Link to interactive map online", href = link_to_map, target = "_blank")
                                                   # })
                                            )
                                          ),
                                          fluidRow(
                                            column(12,
                                                   renderPlot({
                                                     
                                                     plotCPRKEGG(object = session$userData$FlomicsMultiAssay[[dataset]],
                                                                 contrast = listname,
                                                                 from = list.source,
                                                                 pathway_id = input[[paste0(listname, "-MAP.sel")]],
                                                                 species = local.rea.values$KEGG_org,
                                                                 gene_idtype = local.rea.values$keytype.kegg
                                                     )
                                                   }, res = 300, width = 1000, height = 1000),
                                            )
                                          )
                                 )
                               )
            )
            
          }
          # display results
          fluidRow(
            
            box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success", title = listname,
                
                fluidRow(
                  column(4,
                         radioButtons(inputId = ns(paste0(listname, "-domain")), label = "Domain",
                                      choices = choices, selected = choices[1])),
                  column(4,
                         numericInput(inputId = ns(paste0(listname, "-top.over")), label = "Top terms:" , value=15 ,
                                      min = 1, max = 10000, step = 5)), # max code en dur : pas bien
                  column(4,
                         textInput(inputId = ns(paste0(listname, "-grep")), label = "Search Expression"))
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
  
}


### main 

AnnotationEnrichmentClusterProfUI <- function(id){
  
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
    fluidRow( uiOutput(ns("tabsetPanel_DB_UI")) )
  )
}

#' @importFrom data.table fread
#' @importFrom AnnotationDbi keytypes
#' @importFrom DT renderDataTable datatable
AnnotationEnrichmentClusterProf <- function(input, output, session, dataset, rea.values){
  
  ns <- session$ns
  
  local.rea.values <- reactiveValues()
  
  # ----- Panel for results -----
  output$tabsetPanel_DB_UI <- renderUI({
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis and validate your choices...")
    )
    
    omicsType <- getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]])
    
    tabPanel_custom.list <- list()
    tabPanel_genes.list  <- list()
    
    tabPanel_custom.list <- list(
      tabPanel(title = "Custom annotation",
               br(),
               column(3, uiOutput(ns("AnnotParamCustom_UI"))),
               column(9, 
                      tabsetPanel(
                        tabPanel(title = "Differential expression lists",
                                 br(),
                                 module_runEnrichment_UI(id = ns("custom_DiffExpEnrichAnal"))
                        ),
                        tabPanel(title = "Co-expression clusters",
                                 br(),
                                 module_runEnrichment_UI(id = ns("custom_CoExpEnrichAnal"))
                        )
                      )
               )
      )
    )
    
    tabPanel_genes.list <- list(
      tabPanel(title = "Gene Onlology database",
               br(),
               column(3, uiOutput(ns("AnnotParamGO_UI"))),
               column(9, 
                      tabsetPanel(
                        tabPanel(title = "Differential expression lists",
                                 br(),
                                 module_runEnrichment_UI(id = ns("GO_DiffExpEnrichAnal"))
                        ),
                        tabPanel(title = "Co-expression clusters",
                                 br(),
                                 module_runEnrichment_UI(id = ns("GO_CoExpEnrichAnal"))
                        )
                      )
               )
      ),
      tabPanel(title = "KEGG database",
               br(),
               column(3, uiOutput(ns("AnnotParamKEGG_UI"))),
               column(9, 
                      tabsetPanel(
                        tabPanel(title = "Differential expression lists",
                                 br(),
                                 module_runEnrichment_UI(id = ns("KEGG_DiffExpEnrichAnal"))
                        ),
                        tabPanel(title = "Co-expression clusters",
                                 br(),
                                 module_runEnrichment_UI(id = ns("KEGG_CoExpEnrichAnal"))
                        )
                      )
               )
      )
    )
    
    
    
    if (omicsType %in% c("RNAseq", "proteomics")) {
      
      do.call(what = tabsetPanel, args = c(tabPanel_genes.list, tabPanel_custom.list) )
    }
    else{
      do.call(what = tabsetPanel, args = c(tabPanel_custom.list) )
    }
    
    
    
    
  })
  
  # ---- settings: ----
  output$AnnotParamGO_UI <-   renderUI({
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    ## gene lists
    data.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    ListNames.diff <- getValidContrasts(data.SE)$contrastName
    omicsType <- getOmicsTypes(data.SE)
    
    if (omicsType %in% c("RNAseq", "proteomics")) {
      
      box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
          
          # select DE list
          pickerInput(
            inputId = ns("GeneList.diff_GO"), label = "Select DE lists:", 
            choices = ListNames.diff, multiple = TRUE, selected = ListNames.diff,
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          
          # select clusters
          uiOutput(ns("GeneList.coseqCPRUI_GO")),
          
          hr(),
          
          # select ontology
          pickerInput(
            inputId = ns("ont.select"), label = "Select Domain:",
            choices = c("ALL", "BP", "MF", "CC"),
            selected = "ALL",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          
          # select database
          pickerInput(
            inputId = ns("db.select"), label = "Select Data Base:",
            choices = c("", dir(.libPaths())[grep("org.*db", dir(.libPaths()))]),
            selected = "",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          
          # select keytype
          pickerInput(
            inputId = ns("keytype.go"), label = "Select id type:",
            choices = c("", "ENTREZID", "SYMBOL",  "TAIR", "ENSEMBL", "PMID", "REFSEQ", "GENENAME"),
            selected = "",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          #),
          
          ## alpha threshold
          numericInput(inputId = ns("pValue_GO"), label = "Adjusted p-value threshold:", value = 0.01 , min = 0, max = 1, step = 0.01),
          
          # run enrichment
          actionButton(inputId = ns("runEnrich_CPR_GO"), label = "Run")
      )
      
    }
  })
  
  output$AnnotParamKEGG_UI <-   renderUI({
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    ## gene lists
    data.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    ListNames.diff <- getValidContrasts(data.SE)$contrastName
    omicsType <- getOmicsTypes(data.SE)
    
    if (omicsType %in% c("RNAseq", "proteomics")) {
      
      box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
          
          # select DE list
          pickerInput(
            inputId = ns("GeneList.diff_KEGG"), label = "Select DE lists:", 
            choices = ListNames.diff, multiple = TRUE, selected = ListNames.diff,
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          
          # select clusters
          uiOutput(ns("GeneList.coseqCPRUI_KEGG")),
          
          hr(),
          
          pickerInput(
            inputId = ns("keytype.kegg"), label = "Select id type:",
            choices = c("", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
            selected = "",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")),
          
          # type organism 
          textInput(inputId = ns("KEGG_org"), 
                    label = "Organism: (eg. ath)", value = "", 
                    width = NULL, placeholder = NULL),
          
          ## alpha threshold
          numericInput(inputId = ns("pValue_KEGG"), label = "Adjusted p-value threshold:", value = 0.01 , min = 0, max = 1, step = 0.01),
          
          # run enrichment
          actionButton(inputId = ns("runEnrich_CPR_KEGG"), label = "Run")
      )
      
    }
  })
  
  output$AnnotParamCustom_UI <-   renderUI({
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    ## gene lists
    data.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    ListNames.diff <- getValidContrasts(data.SE)$contrastName
    omicsType <- getOmicsTypes(data.SE)
    
    box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
        
        # select DE list
        pickerInput(
          inputId = ns("GeneList.diff_custom"), label = "Select DE lists:", 
          choices = ListNames.diff, multiple = TRUE, selected = ListNames.diff,
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        # select clusters
        uiOutput(ns("GeneList.coseqCPRUI_custom")),
        
        hr(),
        
        fileInput(inputId = ns("annotationFileCPR"), label = "Annotation file:",
                  accept  = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        
        uiOutput(ns("selectColumnsFromCustomAnnot")),
        
        ## alpha threshold
        numericInput(inputId = ns("pValue_custom"), label = "Adjusted p-value threshold:", value = 0.01 , min = 0, max = 1, step = 0.01),
        
        # run enrichment
        actionButton(inputId = ns("runEnrich_CPR_Custom"), label = "Run")
    )
  })
  
  # ---- if custom annoattion ----
  observeEvent(input$annotationFileCPR$datapath, {
    
    output$selectColumnsFromCustomAnnot <- renderUI({
      
      annotation <- data.table::fread(file = input$annotationFileCPR$datapath, sep = "\t", header = TRUE)
      
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
  
  # ---- Coseq UI if there are results ----
  output$GeneList.coseqCPRUI_GO <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnal == FALSE) return()
    
    ListNames.coseq <- names(session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpAnal[["clusters"]])
    
    pickerInput(
      inputId = ns("GeneList.coseq_GO"), label = "Select Clusters:",
      choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
    )
  })
  
  output$GeneList.coseqCPRUI_KEGG <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnal == FALSE) return()
    
    ListNames.coseq <- names(session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpAnal[["clusters"]])
    
    pickerInput(
      inputId = ns("GeneList.coseq_KEGG"), label = "Select Clusters:",
      choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
    )
  })
  
  output$GeneList.coseqCPRUI_custom <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnal == FALSE) return()
    
    ListNames.coseq <- names(session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpAnal[["clusters"]])
    
    pickerInput(
      inputId = ns("GeneList.coseq_custom"), label = "Select Clusters:",
      choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
    )
  })
  
  ##================================= RUN ====================================##
  
  # ---- run Annotation  ----
  observeEvent(input$runEnrich_CPR_GO, {
    
    # check list of genes
    if (length(c(input$GeneList.diff_GO, input$GeneList.coseq_GO )) == 0) {
      showModal(modalDialog( title = "Error message", "Please select at least 1 variable list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff_GO, input$GeneList.coseq_GO)) != 0, message = "Please select at least 1 variable list") 
    })
    
    # check parameters
    if (input$db.select == "") {
      showModal(modalDialog(title = "Error message", "Select databases..."))
    }
    validate({ 
      need(input$db.select != "", "Select databases...") 
    })
    
    library(input$db.select, character.only = TRUE)
    accepted_keytypes <- AnnotationDbi::keytypes(get(input$db.select))
    if (!input$keytype.go %in% accepted_keytypes) {
      showModal(modalDialog(title = paste("Keytype must be one of: ", paste(accepted_keytypes, collapse = ", "), sep = " ")))
    }
    validate({ 
      need(input$keytype.go %in% accepted_keytypes, message = paste("Keytype must be one of: ", paste(accepted_keytypes, collapse = ", "), sep = " ")) 
    })
    
    # set list args
    dom.select <- "GO"
    list_args <- list()
    list_args[["OrgDb"]]        <- input$db.select
    list_args[["keyType"]]      <- input$keytype.go
    list_args[["pvalueCutoff"]] <- input$pValue_GO
    list_args[["minGSSize"]]    <- 10
    Domain <- input$ont.select; if (input$ont.select == "ALL") Domain <- c("MF", "BP", "CC")
    annotation2    <- NULL
    Domain.        <- NULL
    col_domain_arg <- NULL
    
    param.list <- list(method       = "ORA",
                       diffList     = input$GeneList.diff_GO,
                       CoexpList    = input$GeneList.coseq_GO,
                       OrgDb        = input$db.select,
                       keyType      = input$keytype.go,
                       Domain       = Domain,
                       pvalueCutoff = input$pValue_GO
    )
    
    # prevent multiple execution
    if (.checkRunORAExecution(session$userData$FlomicsMultiAssay[[dataset]], "GO", param.list) == FALSE) return()
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Annot", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    #local.rea.values[[dom.select]] <- FALSE
    
    # ---- Annotation on diff results: ----  
    if (length(input$GeneList.diff_GO) != 0) {
      
      # run analysis
      print(paste("# 11- Enrichment Analysis of DE lists...", dataset))
      
      session$userData$FlomicsMultiAssay[[dataset]]@metadata[["DiffExpEnrichAnal"]][[dom.select]] <- NULL
      
      # run annotation diff analysis
      runRes <- tryCatch({
        runAnnotationEnrichment(object = session$userData$FlomicsMultiAssay[[dataset]], 
                                nameList = input$GeneList.diff_GO,
                                list_args = list_args, 
                                from = "DiffExpAnal",
                                annot = annotation2,
                                dom.select = dom.select, 
                                Domain = Domain,
                                col_domain = col_domain_arg)
      },
      warning = function(war) return(war),
      error   = function(err) return(err)
      )
      
      if (!is(runRes, "SummarizedExperiment")) {
        showModal(modalDialog(title = paste("Something went wrong: ", runRes$message)))
      }
      validate({need(is(runRes, "SummarizedExperiment"), message = paste0("Something went wrong: ", runRes$message))})
      
      session$userData$FlomicsMultiAssay[[dataset]] <-  runRes
      rea.values[[dataset]]$diffAnnot <- TRUE
      shiny::callModule(module  = module_runEnrichment, id = "GO_DiffExpEnrichAnal", dataset = dataset,
                        dom.select = "GO", 
                        list.source = "DiffExpEnrichAnal",
                        rea.values = rea.values, local.rea.values = local.rea.values)
      
    }
    
    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep = ""))
    #----------------------#
    
    # ---- Annotation on COEXP results: ----
    if (length(input$GeneList.coseq_GO) != 0) {
      
      # run analysis
      print(paste("# 11- Enrichment Analysis of clusters...", dataset))
      
      session$userData$FlomicsMultiAssay[[dataset]]@metadata[["CoExpEnrichAnal"]][[dom.select]]   <- NULL
      
      runResCo <- tryCatch({
        session$userData$FlomicsMultiAssay[[dataset]] <-  
          runAnnotationEnrichment(session$userData$FlomicsMultiAssay[[dataset]], 
                                  nameList = input$GeneList.coseq_GO,
                                  list_args = list_args, 
                                  from = "CoExpAnal",
                                  annot = annotation2,
                                  dom.select = dom.select, 
                                  Domain = Domain,
                                  col_domain = col_domain_arg)
        
      },
      warning = function(war) return(war),
      error   = function(err) return(err)
      )
      
      if (!is(runResCo, "SummarizedExperiment")) {
        showModal(modalDialog(title = paste("Something went wrong: ", runResCo$message)))
      }
      validate({need(is(runResCo, "SummarizedExperiment"), message = paste0("Something went wrong: ", runResCo$message))})
      
      session$userData$FlomicsMultiAssay[[dataset]] <-  runResCo
      rea.values[[dataset]]$coExpAnnot <- TRUE
      shiny::callModule(module  = module_runEnrichment, id = "GO_CoExpEnrichAnal", dataset = dataset, 
                        dom.select = "GO", list.source = "CoExpEnrichAnal",
                        rea.values = rea.values, local.rea.values = local.rea.values)
    }
    
    #local.rea.values[[dom.select]] <- TRUE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep = ""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
  # ---- run Annotation  KEGG ----
  observeEvent(input$runEnrich_CPR_KEGG, {
    
    # check list of genes
    if (length(c(input$GeneList.diff_KEGG, input$GeneList.coseq_KEGG)) == 0) {
      showModal(modalDialog( title = "Error message", "Please select at least 1 variable list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff_KEGG, input$GeneList.coseq_KEGG)) != 0, message = "Please select at least 1 variable list") 
    })
    
    # check param
    ## ID key
    if (!input$keytype.kegg %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")) {
      showModal(modalDialog( title = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot"))
    }
    validate({ 
      need(input$keytype.kegg %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"), 
           message = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot") 
    })
    
    ## code organism KEGG_org
    if (input$KEGG_org == "") {
      showModal(modalDialog( title = "Error message", "set organism kegg code."))
    }
    validate({ 
      need(input$KEGG_org != "", message = "Set organism kegg code.") 
    })
    
    # set list args
    dom.select <- "KEGG"
    list_args <- list()
    Domain <- NULL
    annotation2 <- NULL
    col_domain_arg <- NULL
    list_args[["organism"]] <- input$KEGG_org
    list_args[["keyType"]]  <- input$keytype.kegg
    list_args[["pvalueCutoff"]] <- input$pValue_KEGG
    list_args[["minGSSize"]] <- 10
    
    # prevent multiple execution   
    param.list <- list(method       = "ORA",
                       diffList     = input$GeneList.diff_KEGG,
                       CoexpList    = input$GeneList.coseq_KEGG,
                       organism     = input$KEGG_org,
                       keyType      = input$keytype.kegg,
                       Domain       = "no-domain",
                       pvalueCutoff = input$pValue_KEGG
    )
    
    param.list <<- param.list
    
    if (.checkRunORAExecution(session$userData$FlomicsMultiAssay[[dataset]], "KEGG", param.list) == FALSE) return()

    local.rea.values[["KEGG_org"]]     <- input$KEGG_org
    local.rea.values[["keytype.kegg"]] <- input$keytype.kegg
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Annot", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    # ---- Annotation on diff results: ----  
    if (length(input$GeneList.diff_KEGG) != 0) {
      
      print(paste("# 11- Enrichment Analysis of DE lists...", dataset))
      
      session$userData$FlomicsMultiAssay[[dataset]]@metadata[["DiffExpEnrichAnal"]][[dom.select]] <- NULL
      
      # run annotation
      runRes <- tryCatch({
      session$userData$FlomicsMultiAssay[[dataset]] <-  
        runAnnotationEnrichment(object = session$userData$FlomicsMultiAssay[[dataset]], 
                                nameList = input$GeneList.diff_KEGG,
                                list_args = list_args, 
                                from = "DiffExpAnal",
                                annot = annotation2,
                                dom.select = dom.select, 
                                Domain = Domain,
                                col_domain = col_domain_arg)
      },
      warning = function(war) return(war),
      error   = function(err) return(err)
      )
      
      if (!is(runRes, "SummarizedExperiment")) {
        showModal(modalDialog(title = paste("Something went wrong: ", runRes$message)))
      }
      validate({need(is(runRes, "SummarizedExperiment"), message = paste0("Something went wrong: ", runRes$message))})
      
      session$userData$FlomicsMultiAssay[[dataset]] <-  runRes
      shiny::callModule(module  = module_runEnrichment, id = "KEGG_DiffExpEnrichAnal", 
                        dataset = dataset, dom.select = "KEGG", 
                        list.source = "DiffExpEnrichAnal",
                        rea.values = rea.values, local.rea.values = local.rea.values)
      
      rea.values[[dataset]]$diffAnnot <- TRUE
    }
    
    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep = ""))
    #----------------------#
    
    # ---- Annotation on COEXP results: ----
    if (length(input$GeneList.coseq_KEGG) != 0) {
      
      print(paste("# 11- Enrichment Analysis of clusters...", dataset))
      
      session$userData$FlomicsMultiAssay[[dataset]]@metadata[["CoExpEnrichAnal"]][[dom.select]]   <- NULL
      
      runResCo <- tryCatch({
      session$userData$FlomicsMultiAssay[[dataset]] <-  
        runAnnotationEnrichment(session$userData$FlomicsMultiAssay[[dataset]], 
                                nameList = input$GeneList.coseq_KEGG,
                                list_args = list_args, 
                                from = "CoExpAnal",
                                annot = annotation2,
                                dom.select = dom.select, 
                                Domain = Domain,
                                col_domain = col_domain_arg) 
      },
      warning = function(war) return(war),
      error   = function(err) return(err)
      )
      
      if (!is(runResCo, "SummarizedExperiment")) {
        showModal(modalDialog(title = paste("Something went wrong: ", runResCo$message)))
      }
      validate({need(is(runResCo, "SummarizedExperiment"), message = paste0("Something went wrong: ", runResCo$message))})
      
      session$userData$FlomicsMultiAssay[[dataset]] <-  runResCo
      shiny::callModule(module  = module_runEnrichment, id = "KEGG_CoExpEnrichAnal", 
                        dataset = dataset, dom.select = "KEGG", 
                        list.source = "CoExpEnrichAnal",
                        rea.values = rea.values, local.rea.values = local.rea.values)
      
      rea.values[[dataset]]$coExpAnnot <- TRUE
    }
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep = ""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
  # ---- run Annotation on custom file ----
  observeEvent(input$runEnrich_CPR_Custom, {
    
    # ---- Checks: ----
    
    # check list of genes
    if (length(c(input$GeneList.diff_custom, input$GeneList.coseq_custom )) == 0) {
      showModal(modalDialog( title = "Error message", "Please select at least 1 variable list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff_custom, input$GeneList.coseq_custom)) != 0, message = "Please select at least 1 variable list") 
    })
    
    # check param
    ## if custom annotation file 
    if (is.null(input$annotationFileCPR$datapath)) {
      showModal(modalDialog(title = "Error message", "load custom annotation file."))
    }
    validate({ 
      need(!is.null(input$annotationFileCPR$datapath), 
           message = "load custom annotation file.")
    })
    # Check some custom elements and annotation file
    if (input$col_geneName == "" || input$col_termID == "") {
      showModal(modalDialog( title = "Error message", 
                             "Please choose columns names for the gene names/ID and the ontology terms ID!"))
    }
    validate({ 
      need(input$col_geneName != "" && input$col_termID != "", 
           message = "Please choose columns names for the gene names/ID and the ontology terms ID!")
    })
    
    
    # set param
    dom.select <- "custom"
    Domain <- NULL
    annotation2 <- NULL
    col_domain_arg <- NULL
    list_args <- list()
    list_args[["pvalueCutoff"]] <- input$pValue_custom
    list_args[["minGSSize"]] <- 10
    
    annotation <- data.table::fread(file = input$annotationFileCPR$datapath, sep = "\t", header = TRUE)
    listCharRm <- c(".", " ", "")
    annotation2 <- list()
    
    # if column with domain exist
    if (input$col_domain != "") {
      
      # remove column with domain == c(".", " ", "")
      domainRm <- intersect(annotation[[input$col_domain]], listCharRm)
      if (length(domainRm) != 0) {
        annotation <- dplyr::filter(annotation, !get(input$col_domain) %in% listCharRm)
      }
      
      annotation2[["domain"]] <- annotation[[input$col_domain]]
      Domain <- unique(annotation2$domain)
      col_domain_arg <- "domain"
    }
    
    annotation2[["gene"]] <- annotation[[input$col_geneName]]
    annotation2[["term"]] <- annotation[[input$col_termID]]
    if (input$col_termName != "") annotation2[["name"]] <- annotation[[input$col_termName]]
    annotation2 <- data.frame(annotation2)
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Annot", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    # ---- Annotation on diff results: ----  
    if (length(input$GeneList.diff_custom) != 0) {
      
      print(paste("# 11- Enrichment Analysis of DE lists...", dataset))
      
      session$userData$FlomicsMultiAssay[[dataset]]@metadata[["DiffExpEnrichAnal"]][[dom.select]] <- NULL
      
      # run annotation
      runRes <- tryCatch({
      session$userData$FlomicsMultiAssay[[dataset]] <-  
        runAnnotationEnrichment(object = session$userData$FlomicsMultiAssay[[dataset]], 
                                nameList = input$GeneList.diff_custom,
                                list_args = list_args, 
                                from = "DiffExpAnal",
                                annot = annotation2,
                                dom.select = dom.select, 
                                Domain = Domain,
                                col_domain = col_domain_arg)  
      },
      warning = function(war) return(war),
      error   = function(err) return(err)
      )
      
      if (!is(runRes, "SummarizedExperiment")) {
        showModal(modalDialog(title = paste("Something went wrong: ", runRes$message)))
      }
      validate({need(is(runRes, "SummarizedExperiment"), message = paste0("Something went wrong: ", runRes$message))})
      
      session$userData$FlomicsMultiAssay[[dataset]] <-  runRes
      shiny::callModule(module  = module_runEnrichment, id = "custom_DiffExpEnrichAnal", 
                        dataset = dataset, dom.select = dom.select, 
                        list.source = "DiffExpEnrichAnal",
                        rea.values = rea.values, local.rea.values = local.rea.values)
      
      rea.values[[dataset]]$diffAnnot <- TRUE
    }
    
    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep = ""))
    #----------------------#
    
    # ---- Annotation on COEXP results: ----
    if (length(input$GeneList.coseq_custom) != 0) {
      
      print(paste("# 11- Enrichment Analysis of clusters...", dataset))
      
      session$userData$FlomicsMultiAssay[[dataset]]@metadata[["CoExpEnrichAnal"]][[dom.select]]   <- NULL
      
      runResCo <- tryCatch({
      session$userData$FlomicsMultiAssay[[dataset]] <-  
        runAnnotationEnrichment(session$userData$FlomicsMultiAssay[[dataset]], 
                                nameList = input$GeneList.coseq_custom,
                                list_args = list_args, 
                                from = "CoExpAnal",
                                annot = annotation2,
                                dom.select = dom.select, 
                                Domain = Domain,
                                col_domain = col_domain_arg) 
      },
      warning = function(war) return(war),
      error   = function(err) return(err)
      )
      
      if (!is(runResCo, "SummarizedExperiment")) {
        showModal(modalDialog(title = paste("Something went wrong: ", runResCo$message)))
      }
      validate({need(is(runResCo, "SummarizedExperiment"), message = paste0("Something went wrong: ", runResCo$message))})
      
      session$userData$FlomicsMultiAssay[[dataset]] <-  runResCo
      shiny::callModule(module  = module_runEnrichment, id = "custom_CoExpEnrichAnal", 
                        dataset = dataset, dom.select = dom.select,
                        list.source = "CoExpEnrichAnal",
                        rea.values = rea.values, local.rea.values = local.rea.values)
      
      rea.values[[dataset]]$coExpAnnot <- TRUE
    }
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep = ""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
}

############## functions ###############

# ----- check run enrichement analysis execution ------
#' @title .checkRunORAExecution
#' @noRd
#' @keywords internal
.checkRunORAExecution <- function(object.SE, db.type = NULL, param.list = NULL){
  

  if (!is.null(param.list$diffList)) {

    if (is.null(object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]])) return(TRUE)
    if (isFALSE(dplyr::setequal(names(object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$enrichResult), param.list$diffList)) ) return(TRUE)
    
    # check param
    if (isFALSE(dplyr::setequal(param.list$Domain, object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$list_args$Domain))) return(TRUE)
    if (param.list$pvalueCutoff != object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$list_args$pvalueCutoff) return(TRUE)
    
    switch (db.type,
            "GO" = {
              if(param.list$OrgDb        != object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$list_args$OrgDb)    return(TRUE)
              if(param.list$keyType      != object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$list_args$keyType)  return(TRUE)
            },
            "KEGG" = {
              if(param.list$keyType      != object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$list_args$keyType)  return(TRUE)
              if(param.list$organism     != object.SE@metadata[["DiffExpEnrichAnal"]][[db.type]]$list_args$organism) return(TRUE)
            }
    )
    
  }
  if (!is.null(param.list$CoexpList)) {
    
    if (is.null(object.SE@metadata[["CoExpEnrichAnal"]][[db.type]])) return(TRUE)
    if (isFALSE(dplyr::setequal(names(object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$enrichResult), param.list$CoexpList)) ) return(TRUE)
    
    # check param
    if (isFALSE(dplyr::setequal(param.list$Domain, object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$list_args$Domain))) return(TRUE)
    if (param.list$pvalueCutoff != object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$list_args$pvalueCutoff) return(TRUE)
    
    switch (db.type,
            "GO" = {
              if (param.list$OrgDb        != object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$list_args$OrgDb)   return(TRUE)
              if (param.list$keyType      != object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$list_args$keyType) return(TRUE)
            },
            "KEGG" = {
              if (param.list$keyType      != object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$list_args$keyType)  return(TRUE)
              if (param.list$organism     != object.SE@metadata[["CoExpEnrichAnal"]][[db.type]]$list_args$organism) return(TRUE)
            }
    )
    
  }
  
  return(FALSE)
}




