##########################################
# Part6 : ANNOTAION
##########################################

### main 

AnnotationEnrichmentClusterProfUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    fluidRow(
      box(title = span(tagList(icon("chart-pie"), " ", a("ClusterProfiler/", href="https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html"), 
                               a("Pathview",         href="https://bioconductor.org/packages/devel/bioc/vignettes/pathview/inst/doc/pathview.pdf"), "   " , tags$small("(Scroll down for instructions)"))),
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
    fluidRow(
      column(3, uiOutput(ns("AnnotParam"))),
      column(9, 
             uiOutput(ns("summary_diff")),
             uiOutput(ns("AnnotResultsCPR_diff")),
             uiOutput(ns("summary_coex")),
             uiOutput(ns("AnnotResultsCPR_coex")))
    )
  )
}

AnnotationEnrichmentClusterProf <- function(input, output, session, dataset, rea.values){
  
  ns <- session$ns
  
  local.rea.values <- reactiveValues(dataset.SE = NULL)
  local.rea.values$KEGG   <- FALSE 
  local.rea.values$custom <- FALSE 
  local.rea.values$GO     <- FALSE 
  
  
  # ---- settings: ----
  output$AnnotParam <- renderUI({
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis and validate your choices...")
    )
    
    ## gene lists
    ListNames.diff  <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName
    
    box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
        
        # select DEG list
        pickerInput(
          inputId = ns("GeneList.diff"), label = "Select DEG lists:", 
          choices = ListNames.diff, multiple = TRUE, selected = ListNames.diff,
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        # select clusters
        uiOutput(ns("GeneList.coseqCPRUI")),
        
        # select enrichDomain
        pickerInput(
          inputId = ns("dom.select"), label = "Select Ontology:",
          choices = c("custom", "GO", "KEGG"),
          selected = "GO",
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        # custum
        conditionalPanel(
          condition = paste0("input[\'", ns('dom.select'), "\'] == \'custom\'"),
          fileInput(inputId = ns("annotationFileCPR"), label = "Annotation file:",
                    accept  = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
          
          uiOutput(ns("selectColumnsFromCustomAnnot"))
        ),
        
        # KEGG 
        conditionalPanel(
          condition = paste0("input[\'", ns('dom.select'), "\'] == \'KEGG\'"),
          
          pickerInput(
            inputId = ns("keytype.kegg"), label = "Select id type:",
            choices = c("", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
            selected = "",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")),
          
          # type organism 
          textInput(inputId = ns("KEGG_org"), 
                    label = "Organism: (eg. ath)", value = "", 
                    width = NULL, placeholder = NULL)
        ),
        conditionalPanel(
          condition = paste0("input[\'", ns('dom.select'), "\'] == \'GO\'"),
          
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
          )
        ),
        
        ## alpha threshold
        numericInput(inputId = ns("pValue"), label="Adjusted p-value threshold:", value = 0.01 , min = 0, max = 1, step = 0.01),
        
        # run enrichment
        actionButton(inputId = ns("runEnrich_CPR"), label = "Run")
    )
    
  })
  
  # ---- if custom annoattion ----
  observeEvent(input$annotationFileCPR$datapath, {
    
    output$selectColumnsFromCustomAnnot <- renderUI({
      
      annotation <- fread(file = input$annotationFileCPR$datapath, sep="\t", header = TRUE)
      
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
  output$GeneList.coseqCPRUI <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnal == FALSE) return()
    
    ListNames.coseq <- names(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$CoExpAnal[["clusters"]])
    
    pickerInput(
      inputId = ns("GeneList.coseq"), label = "Select Clusters :",
      choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
    )
  })
  
  
  ##================================= RUN ====================================##
  
  # ---- run Annotation  ----
  observeEvent(input$runEnrich_CPR, {
    
    library(clusterProfiler)
    library(pathview)
    library(enrichplot)
    library(org.At.tair.db)
    library(AnnotationDbi)
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Annot", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    # # reset everything
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata[["DiffExpEnrichAnal"]][[input$dom.select]] <- NULL
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata[["CoExpEnrichAnal"]][[input$dom.select]]   <- NULL
    
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    local.rea.values[[input$dom.select]] <- FALSE
    
    local.rea.values$dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    # ---- Checks: ----
    
    # check list of genes
    if(length(c(input$GeneList.diff, input$GeneList.coseq )) == 0){
      showModal(modalDialog( title = "Error message", "Please select at least 1 variable list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff, input$GeneList.coseq)) != 0, message = "Please select at least 1 variable list") 
    })
    
    list_args <- list()
    Domain <- NULL
    annotation2 <- NULL
    col_domain_arg <- NULL
    
    switch (input$dom.select,
            #==== DB : GO ====
            "GO" = {
              
              # check param
              ## databases
              if(input$db.select == ""){
                showModal(modalDialog(title = "Error message", "Select databases..."))
              }
              validate({ 
                need(input$db.select != "", "Select databases...") 
              })
              library(input$db.select, character.only = TRUE)
              
              ## ID type
              accepted_keytypes <- AnnotationDbi::keytypes(get(input$db.select))
              
              if(!input$keytype.go %in% accepted_keytypes){
                showModal(modalDialog(title = paste("Keytype must be one of: ", paste(accepted_keytypes, collapse = ", "), sep = " ")))
              }
              validate({ 
                need(input$keytype.go %in% accepted_keytypes, message = paste("Keytype must be one of: ", paste(accepted_keytypes, collapse = ", "), sep = " ")) 
              })
              
              # set list args
              list_args[["OrgDb"]]   <- input$db.select
              list_args[["keyType"]] <- input$keytype.go
              
              Domain <- input$ont.select
              if (input$ont.select == "ALL") Domain <- c("MF", "BP", "CC")
              
            },
            #==== DB : KEGG ====
            "KEGG" = {
              
              # check param
              ## ID key
              if(!input$keytype.kegg %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")){
                showModal(modalDialog( title = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot"))
              }
              validate({ 
                need(input$keytype.kegg %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"), message = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot") 
              })
              
              ## code organism KEGG_org
              if(input$KEGG_org == ""){
                showModal(modalDialog( title = "Error message", "set organism kegg code."))
              }
              validate({ 
                need(input$KEGG_org != "", message = "Set organism kegg code.") 
              })              
              
              # set list args
              list_args[["organism"]] <- input$KEGG_org
              list_args[["keyType"]]  <- input$keytype.kegg
              
            },
            #==== custom annotation ====
            "custom" = {
              
              # check param
              ## if custom annotation file 
              if(is.null(input$annotationFileCPR$datapath)){
                showModal(modalDialog(title = "Error message", "load custom annotation file."))
              }
              # Check some custom elements and annotation file
              if(input$col_geneName == "" || input$col_termID == ""){
                showModal(modalDialog( title = "Error message", "Please choose columns names for the gene names/ID and the ontology terms ID!"))
              }
              validate({ 
                need(input$col_geneName != "" && input$col_termID != "", message = "Please choose columns names for the gene names/ID and the ontology terms ID!")
              })
              # Check if geneName correspond to variable list
              
              #======
              annotation <- fread(file = input$annotationFileCPR$datapath, sep="\t", header = TRUE)
              list.char.rm <- c(".", " ", "")
              annotation2 <- list()
              
              # if column with domain exist
              if(input$col_domain != "") {
                
                # remove column with domain == c(".", " ", "")
                domain.rm <- intersect(annotation[[input$col_domain]], list.char.rm)
                if (length(domain.rm) != 0){
                  annotation <- dplyr::filter(annotation, ! get(input$col_domain) %in% list.char.rm)
                }
                
                annotation2[["domain"]] <- annotation[[input$col_domain]]
                Domain <- unique(annotation2$domain)
                col_domain_arg <- "domain"
              }
              
              annotation2[["gene"]] <- annotation[[input$col_geneName]]
              annotation2[["term"]] <- annotation[[input$col_termID]]
              if(input$col_termName != "") annotation2[["name"]] <- annotation[[input$col_termName]]
              
              annotation2 <- data.frame(annotation2)
              
              # print(head(annotation2))
              # list_args[["annotation"]] <- annotation2
            }
    )
    
    list_args[["pvalueCutoff"]] <- input$pValue
    
    print(paste("# 11- Enrichment Analysis...", dataset))
    
    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#
    
    ## run annotation
    dom.select <- input$dom.select
    
    # ---- Annotation on diff results: ----  
    if(length(input$GeneList.diff) != 0){
      
      # run annotation
      local.rea.values$dataset.SE <- runAnnotationEnrichment_CPR(local.rea.values$dataset.SE, list_args = list_args, from = "DiffExpAnal",
                                                                 annot = annotation2,
                                                                 dom.select = dom.select, Domain = Domain,
                                                                 col_domain = col_domain_arg) 
      
      rea.values[[dataset]]$diffAnnot <- TRUE
    }
    
    # ---- Annotation on COEXP results: ----
    if(length(input$GeneList.coseq) != 0){
      
      local.rea.values$dataset.SE <- runAnnotationEnrichment_CPR(local.rea.values$dataset.SE, list_args = list_args, from = "CoExpAnal",
                                                                 annot = annotation2,
                                                                 dom.select = dom.select, Domain = Domain,
                                                                 col_domain = col_domain_arg) 
      
      rea.values[[dataset]]$coExpAnnot <- TRUE
    }
    
    local.rea.values[[input$dom.select]] <- TRUE
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]] <- local.rea.values$dataset.SE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#
  })
  
  ####================================= DISPLAY ====================================####
  
  ### DIFF summary results 
  output$summary_diff <- renderUI({
    
    if(rea.values[[dataset]]$diffAnnot == FALSE || local.rea.values[[input$dom.select]] == FALSE) return()
    
    results <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata[["DiffExpEnrichAnal"]][[input$dom.select]]
    
    if(is.null(results[["summary"]])){
      
      fluidRow(
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = "Summary from differential expression analysis",
            "There is no results for enrichment analysis! Check geneID"))
    }
    else{
      fluidRow(
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", title = "Summary from differential expression analysis",
            DT::renderDataTable({ DT::datatable(results[["summary"]], rownames = FALSE, options = list(pageLength = 6, dom = 'tip')) })))
    }
  })
  
  ### coseq summary results 
  output$summary_coex <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnnot == FALSE || local.rea.values[[input$dom.select]] == FALSE) return()
    
    results <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata[["CoExpEnrichAnal"]][[input$dom.select]]
    
    if(is.null(results[["summary"]])){
      
      fluidRow(
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", 
            title = paste0(input$dom.select," : summary from co-expression analysis"),
            "There is no results for enrichment analysis! Check geneID"))
      
    }
    else{   
      fluidRow(
        hr(),
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", 
            title = paste0(input$dom.select," : Summary from co-expression analysis"),
            DT::renderDataTable({ DT::datatable(results[["summary"]], rownames = FALSE, options = list(pageLength = 6, dom = 'tip')) })))
    }
    
  })
  
  ### DIFF display results
  
  #---- Figure Code for diff results: ----
  output$AnnotResultsCPR_diff <- renderUI({
    
    if(rea.values[[dataset]]$diffAnnot == FALSE || is.null(local.rea.values$dataset.SE@metadata[["DiffExpEnrichAnal"]][[input$dom.select]][["summary"]])) return()
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    results <- dataset.SE@metadata[["DiffExpEnrichAnal"]][[input$dom.select]]
    
    ## log2FC of the genes (put colors on graphs)
    log2FC_lists <- lapply(input$GeneList.diff, function(listname){
      vect <- dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[listname]][["logFC"]]
      names(vect) <- rownames(dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[listname]])
      return(vect)
    })
    names(log2FC_lists) <- input$GeneList.diff
    
    # foreach genes list selected (contrast)
    lapply(names(results[["enrichResult"]]), function(listname) {
      
      # annot.nbr <- lapply(results[["enrichResult"]][[listname]], function(enrichResult.ont){
      #   if(!is.null(enrichResult.ont)) table(enrichResult.ont@result$p.adjust < 0.1)["TRUE"] }) %>% unlist() # %>% sum(na.rm = TRUE)
      # 
      # # pseudo_title <- paste(annot.nbr, names(results[["enrichResult"]][[listname]])) %>%
      # #   paste(collapse = "; ")
      # 
      # if(sum(annot.nbr, na.rm = TRUE) == 0){}
      
      data <- results[["enrichResult"]][[listname]]
      
      if(length(data) != 0){
        if(sum(unlist(results$summary[which(results$summary$Contrast == listname),-1]), na.rm = TRUE) == 0){
          
          fluidRow(
            box(width = 12, title = paste0(listname, " : 0 enriched terms found"), status = "danger")
          )
        }
        else{
          
          log2FC_vect <- log2FC_lists[[listname]]
          
          tabPanel.list <- list(
            #tabsetPanel(
            # ---- Tab Panel : dotPlot : ----
            tabPanel("DotPlot",
                     
                     verticalLayout(
                       renderPlot({
                         
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over")]])
                         Categories <- dataPlot@result$Description[1:NbtoPlot]
                         if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                         
                         enrichplot::dotplot(dataPlot, showCategory = Categories)
                         
                       })
                     )
            ),
            # ---- Tab Panel : heatplot : ----
            tabPanel("Heatplot",
                     verticalLayout(
                       renderPlot({
                         
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over")]])
                         Categories <- dataPlot@result$Description[1:NbtoPlot]
                         if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                         if(length(Categories) == 0){ # TODO improve this, it doesn't work
                           renderText(expr = {
                             "There is no result to display \n
                                 - The mapping was unsuccessful \n
                                 - The number of enriched terms found is 0, please refer to the overview panel to check this information \n
                                 - You searched for an expression that is not present in the first enriched terms, you can try to increase the number of terms to display to see if there is a change \n
                                 - You tried to display 0 results"
                           })
                         }else{
                           suppressMessages(print(# delete warnings for scale fill replacement
                             enrichplot::heatplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect) +
                               labs(fill="log2FC") +
                               scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                               theme(axis.text.y = element_text(size = 10))
                           ))
                         }
                       }),
                     )
                     
            ) ,
            # ---- Tab Panel : cnetplot : ----
            tabPanel("Cnetplot",
                     verticalLayout(
                       renderPlot({
                         
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         dataTab <- dataPlot@result[dataPlot@result$p.adjust <  results$list_args$pvalueCutoff, ]
                         NbtoPlot <- min(nrow(dataTab),input[[paste0(listname, "-top.over")]])
                         Categories <- dataTab$Description[1:NbtoPlot]
                         if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                         
                         node_label_arg <- "none"
                         if(input[[paste0(listname, "-genesLabels_cnet")]] && input[[paste0(listname, "-termsLabels_cnet")]]){
                           node_label_arg <- "all"
                         }else if(input[[paste0(listname, "-genesLabels_cnet")]] && !input[[paste0(listname, "-termsLabels_cnet")]]){
                           node_label_arg <- "gene"
                         }else if(input[[paste0(listname, "-termsLabels_cnet")]] && !input[[paste0(listname, "-genesLabels_cnet")]]){
                           node_label_arg <- "category"
                         }
                         
                         suppressMessages(print( # delete warnings for scale fill replacement
                           enrichplot::cnetplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, node_label = node_label_arg) +
                             guides(colour=guide_colourbar(title = "log2FC")) +
                             scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
                         ))
                         
                       }),
                       fluidRow(
                         column(3,
                                checkboxInput(inputId = ns(paste0(listname, "-genesLabels_cnet")), label = "Genes Labels", value = FALSE),
                                checkboxInput(inputId = ns(paste0(listname, "-termsLabels_cnet")), label = "Terms Labels", value = TRUE)
                         )
                       )
                     )
            ),
            # ---- Tab Panel : Results table : ----
            tabPanel("Result Table",
                     # hr(),
                     verticalLayout(
                       
                       tags$br(),
                       DT::renderDataTable({
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         
                         DT::datatable(dataPlot@result[dataPlot@result$p.adjust <  results$list_args$pvalueCutoff,], # sometimes it prints non significant results...
                                       rownames = FALSE,
                                       options = list(pageLength = 5, lengthMenu = c(5, 10, 15, 20), scrollX = T))
                         
                       })
                     )
            )
            
          )# TabsetPanel
          
          # ---- Tab Panel : only for KEGG, pathview : ----
          if(input$dom.select == "KEGG")
            tabPanel.list <- c(tabPanel.list,
                               list(
                                 tabPanel("Pathview results",
                                          fluidRow(
                                            column(2,
                                                   selectInput(
                                                     inputId = ns(paste0(listname, "-MAP.sel")), label = "Select map:",
                                                     choices = sort(data[[1]]@result$ID[data[[1]]@result$p.adjust< results$list_args$pvalueCutoff]), multiple = FALSE, selectize = FALSE,
                                                     size = 5
                                                   ),
                                            ),
                                          ),
                                          fluidRow(
                                            column(12,
                                                   renderUI({
                                                     
                                                     # From the browseKEGG function:
                                                     link_to_map <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
                                                                           input[[paste0(listname, "-MAP.sel")]],
                                                                           "/",
                                                                           data[[1]][input[[paste0(listname, "-MAP.sel")]], "geneID"])
                                                     
                                                     a(href=link_to_map, "Link to interactive map online")
                                                     
                                                   })),
                                            column(12,
                                                   
                                                   renderPlot({
                                                     see_pathview(gene.data = log2FC_vect,
                                                                  pathway.id = input[[paste0(listname, "-MAP.sel")]],
                                                                  species = input$KEGG_org,
                                                                  gene.idtype = input$keytype,
                                                                  map.symbol = FALSE,
                                                                  same.layer = FALSE,
                                                                  low = list(gene = "blue"),
                                                                  mid = list(gene = "gray"),
                                                                  high = list(gene = "red"),
                                                                  na.col = "transparent"
                                                                  # cex = 1 # too much
                                                     )
                                                   }, res = 300, width = 1000, height = 1000),
                                            )
                                          )
                                 )
                               )
            )
          
          # display results
          fluidRow(
            
            box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success", title = listname, #paste0(listname, " (", pseudo_title, ")"),
                
                fluidRow(
                  column(4, 
                         radioButtons(inputId = ns(paste0(listname, "-domain")), label="Domain",
                                      choices = names(data), selected = names(data)[1])),
                  column(4,
                         numericInput(inputId = ns(paste0(listname, "-top.over")), label="Top terms:" , value=15 ,
                                      min = 1, max=10000, step = 5)), # max code en dur : pas bien
                  column(4,
                         textInput(inputId = ns(paste0(listname, "-grep")), label="Search Expression"))
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
  
  #---- Figure Code for diff results: ----
  output$AnnotResultsCPR_coex <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnnot == FALSE || is.null(local.rea.values$dataset.SE@metadata[["DiffExpEnrichAnal"]][[input$dom.select]][["summary"]])) return()
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    results <- dataset.SE@metadata[["CoExpEnrichAnal"]][[input$dom.select]]
    
    
    ## log2FC of the genes (put colors on graphs)
    log2FC_lists <- lapply(input$GeneList.diff, function(listname){
      vect <- dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[listname]][["logFC"]]
      names(vect) <- rownames(dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[listname]])
      return(vect)
    })
    names(log2FC_lists) <- input$GeneList.diff
    
    # foreach genes list selected (contrast)
    lapply(names(results[["enrichResult"]]), function(listname) {
      
      # annot.nbr <- lapply(results[["enrichResult"]][[listname]], function(enrichResult.ont){
      #   if(!is.null(enrichResult.ont)) table(enrichResult.ont@result$p.adjust < 0.1)["TRUE"] }) %>% unlist() # %>% sum(na.rm = TRUE)
      # 
      # pseudo_title <- paste(annot.nbr, names(results[["enrichResult"]][[listname]])) %>%
      #   paste(collapse = "; ")
      # 
      # if(sum(annot.nbr, na.rm = TRUE) == 0){}
      
      
      
      data <- results[["enrichResult"]][[listname]]
      
      if(length(data) != 0){
        if(sum(unlist(results$summary[which(results$summary$Contrast == listname),-1]), na.rm = TRUE) == 0){
          
          fluidRow(
            box(width = 12, title = paste0(listname, " : 0 enriched terms found"), status = "danger")
          )
        }
        else{
          
          log2FC_vect <- log2FC_lists[[listname]]
          
          tabPanel.list <- list(
            #tabsetPanel(
            
            # ---- Tab Panel : dotPlot : ----
            tabPanel("DotPlot",
                     
                     verticalLayout(
                       renderPlot({
                         
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over")]])
                         Categories <- dataPlot@result$Description[1:NbtoPlot]
                         if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                         
                         enrichplot::dotplot(dataPlot, showCategory = Categories)
                         
                       }),
                     ),
                     
            ),
            # ---- Tab Panel : heatplot : ----
            tabPanel("Heatplot",
                     verticalLayout(
                       renderPlot({
                         
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over")]])
                         Categories <- dataPlot@result$Description[1:NbtoPlot]
                         if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                         if(length(Categories) == 0){ # TODO improve this, it doesn't work
                           renderText(expr = {
                             "There is no result to display \n
                                 - The mapping was unsuccessful \n
                                 - The number of enriched terms found is 0, please refer to the overview panel to check this information \n
                                 - You searched for an expression that is not present in the first enriched terms, you can try to increase the number of terms to display to see if there is a change \n
                                 - You tried to display 0 results"
                           })
                         }else{
                           suppressMessages(print(# delete warnings for scale fill replacement
                             enrichplot::heatplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect) +
                               labs(fill="log2FC") +
                               scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                               theme(axis.text.y = element_text(size = 10))
                           ))
                         }
                       }),
                     ),
                     
            ) ,
            # ---- Tab Panel : cnetplot : ----
            tabPanel("Cnetplot",
                     verticalLayout(
                       renderPlot({
                         
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         dataTab <- dataPlot@result[dataPlot@result$p.adjust <  results$list_args$pvalueCutoff, ]
                         NbtoPlot <- min(nrow(dataTab),input[[paste0(listname, "-top.over")]])
                         Categories <- dataTab$Description[1:NbtoPlot]
                         if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                         
                         node_label_arg <- "none"
                         if(input[[paste0(listname, "-genesLabels_cnet")]] && input[[paste0(listname, "-termsLabels_cnet")]]){
                           node_label_arg <- "all"
                         }else if(input[[paste0(listname, "-genesLabels_cnet")]] && !input[[paste0(listname, "-termsLabels_cnet")]]){
                           node_label_arg <- "gene"
                         }else if(input[[paste0(listname, "-termsLabels_cnet")]] && !input[[paste0(listname, "-genesLabels_cnet")]]){
                           node_label_arg <- "category"
                         }
                         
                         suppressMessages(print( # delete warnings for scale fill replacement
                           enrichplot::cnetplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, node_label = node_label_arg) +
                             guides(colour=guide_colourbar(title = "log2FC")) +
                             scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
                         ))
                         
                       }),
                       fluidRow(
                         column(3,
                                checkboxInput(inputId = ns(paste0(listname, "-genesLabels_cnet")), label = "Genes Labels", value = FALSE),
                                checkboxInput(inputId = ns(paste0(listname, "-termsLabels_cnet")), label = "Terms Labels", value = TRUE)
                         ),
                       )
                     )
            ),# ---- Tab Panel : Results table : ----
            tabPanel("Result Table",
                     # hr(),
                     verticalLayout(
                       
                       tags$br(),
                       DT::renderDataTable({
                         dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                         
                         DT::datatable(dataPlot@result[dataPlot@result$p.adjust <  results$list_args$pvalueCutoff,], # sometimes it prints non significant results...
                                       rownames = FALSE,
                                       options = list( pageLength = 5,
                                                       lengthMenu = c(5, 10, 15, 20), scrollX = T))
                       }),
                     ),
            )
          )# TabsetPanel
          
          # ---- Tab Panel : only for KEGG, pathview : ----
          if(input$dom.select == "KEGG")
            tabPanel.list <- c(tabPanel.list,
                               list(
                                 tabPanel("Pathview results",
                                          fluidRow(
                                            column(2,
                                                   selectInput(
                                                     inputId = ns(paste0(listname, "-MAP.sel")), label = "Select map:",
                                                     choices = sort(data[[1]]@result$ID[data[[1]]@result$p.adjust< results$list_args$pvalueCutoff]), multiple = FALSE, selectize = FALSE,
                                                     size = 5
                                                   ),
                                            ),
                                          ),
                                          fluidRow(
                                            column(12,
                                                   renderUI({
                                                     
                                                     # From the browseKEGG function:
                                                     link_to_map <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
                                                                           input[[paste0(listname, "-MAP.sel")]],
                                                                           "/",
                                                                           data[[1]][input[[paste0(listname, "-MAP.sel")]], "geneID"])
                                                     
                                                     a(href=link_to_map, "Link to interactive map online")
                                                     
                                                   })),
                                            column(12,
                                                   
                                                   renderPlot({
                                                     see_pathview(gene.data = log2FC_vect,
                                                                  pathway.id = input[[paste0(listname, "-MAP.sel")]],
                                                                  species = input$KEGG_org,
                                                                  gene.idtype = input$keytype,
                                                                  map.symbol = FALSE,
                                                                  same.layer = FALSE,
                                                                  low = list(gene = "blue"),
                                                                  mid = list(gene = "gray"),
                                                                  high = list(gene = "red"),
                                                                  na.col = "transparent"
                                                                  # cex = 1 # too much
                                                     )
                                                   }, res = 300, width = 1000, height = 1000),
                                            )
                                          )
                                 )
                               )
            )
          
          # display results
          fluidRow(
            
            box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success", title = listname, #paste0(listname, " (", pseudo_title, ")"),
                
                fluidRow(
                  column(4, 
                         radioButtons(inputId = ns(paste0(listname, "-domain")), label="Domain",
                                      choices = names(data), selected = names(data)[1])),
                  column(4,
                         numericInput(inputId = ns(paste0(listname, "-top.over")), label="Top terms:" , value=15 ,
                                      min = 1, max=10000, step = 5)), # max code en dur : pas bien
                  column(4,
                         textInput(inputId = ns(paste0(listname, "-grep")), label="Search Expression"))
                ),
                fluidRow(
                  column(width = 12,
                         do.call(what = tabsetPanel, args = tabPanel.list))
                )
            ) # box
          ) # fluidrow
        }
      }
    })# lapply
    # )
    
  })
  
}

######## ANNOTATION CLUSTERPROFILER #########

# Code from: https://stackoverflow.com/questions/60141841/how-to-get-pathview-plot-displayed-directly-rather-than-saving-as-a-file-in-r
# It deletes every file created by pathview
see_pathview <- function(...){
  
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = TRUE)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  img <- png::readPNG(filename)
  # grid::grid.raster(img, width = dim(img)[1]*2, height = dim(img)[2]*2) # does not work
  grid::grid.raster(img)
  nam <- str_split(filename, "[.]")
  invisible(file.remove(filename))
  invisible(file.remove(paste0(nam[[1]][1], ".xml")))
  invisible(file.remove(paste0(nam[[1]][1], ".png")))
  return()
}