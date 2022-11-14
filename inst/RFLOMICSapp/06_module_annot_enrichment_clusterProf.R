AnnotationEnrichmentClusterProfUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    fluidRow(
      box(title = span(tagList(icon("chart-pie"), "   Enrichment analysis ", tags$small("(Scroll down for instructions)"))),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
          p("..."))
    ),
    ## parameters + input data
    fluidRow(
      column(3, uiOutput(ns("AnnotParamCPRUI")),
             uiOutput(ns("AnnotParamCPRUI2")),
             uiOutput(ns("AnnotParamCPRUI3"))
      ),
      column(9, uiOutput(ns("AnnotDiffResultsCPR"))),
      column(9, uiOutput(ns("AnnotCoExpResultsCPR")))
    )
    
  )
}

AnnotationEnrichmentClusterProf <- function(input, output, session, dataset, rea.values){
  
  local.rea.values <- reactiveValues(dataset.SE = NULL)
  local.rea.values$settings_ok <- NULL
  
  # ---- First set of settings: ----
  output$AnnotParamCPRUI <- renderUI({
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis and validate your choices...")
    )
    
    ## gene lists
    ListNames.diff  <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName
    
    box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
        
        # select DEG list
        pickerInput(
          inputId = session$ns("GeneList.diff"), label = "Select DEG lists:", 
          choices = ListNames.diff, multiple = TRUE, selected = ListNames.diff,
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        # select clusters
        uiOutput(session$ns("GeneList.coseqCPRUI")),
        
        # select enrichDomain
        pickerInput(
          inputId = session$ns("dom.select"), label = "Select Domain:",
          choices = c("custom", "GO", "GO:BP", "GO:CC", "GO:MF", "KEGG"),
          selected = "custom",
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        ## alpha threshold
        numericInput(inputId = session$ns("pValue"), label="Adjusted p-value threshold:", value = 0.01 , min = 0, max = 0.1, step = 0.01),
        
        actionButton(inputId = session$ns("settings_ok"), "set"),
        
    )
    
  })
  
  
  # ---- Second set of settings: ----
  local.rea.values$settings_ok <- NULL
  observeEvent(input$settings_ok, {local.rea.values$settings_ok <- NULL; local.rea.values$settings_ok <- TRUE})
  
  output$AnnotParamCPRUI2 <- renderUI({
    if(is.null(local.rea.values$settings_ok)) return()
    
    if(input$dom.select == "custom"){
      box(title = NULL, width = 14, status = "warning",  headerBorder = FALSE,
          
          # Select annotation file
          hr(),
          fileInput(inputId = session$ns("annotationFileCPR"), label = "Annotation file:", 
                    accept  = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
          
          actionButton(inputId = session$ns("settings_custom_ok"), label = "OK")
      )
    }else if(input$dom.select == "KEGG"){
      box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
          # select keytype
          pickerInput(
            inputId = session$ns("keytype"), label = "Select KeyType:",
            choices = c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
            selected = "kegg",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")),
          
          # type organism
          textInput(inputId = session$ns("KEGG_org"), 
                    label = "Organism", value = "ath", 
                    width = NULL, placeholder = NULL),
          
          actionButton(inputId = session$ns("runEnrich_CPR"), label = "Run")
      )
    }else if(length(grep("GO", input$dom.select))!=0){
      box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
          # select database
          pickerInput(
            inputId = session$ns("db.select"), label = "Select Data Base:",
            choices = c("", dir(.libPaths())[grep("org.*db", dir(.libPaths()))]),
            selected = "",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          
          # select keytype
          pickerInput(
            inputId = session$ns("keytype"), label = "Select KeyType:",
            choices = c("SYMBOL", "ENTREZID", "TAIR", "ENSEMBL", "PMID"),
            selected = "",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          actionButton(inputId = session$ns("runEnrich_CPR"), label = "Run")
      )
    }
  })
  
  # ---- Third set of settings: custom only ----
  local.rea.values$settings_custom_ok <- NULL
  observeEvent(input$settings_ok, {local.rea.values$settings_custom_ok <- NULL}) # erase everytime there's a change in the first panel of settings
  # observeEvent(input$settings_custom_ok, {local.rea.values$settings_custom_ok <- NULL; local.rea.values$settings_custom_ok <- TRUE})
  
  # Check if the annotation file is ok
  observeEvent(input$settings_custom_ok, {
    local.rea.values$settings_custom_ok <- NULL

    if(is.null(input$annotationFileCPR$datapath)){
      showModal(modalDialog(title = "Error message", "need annotation file."))
    }

    annotation <- fread(file = input$annotationFileCPR$datapath, sep="\t", header = TRUE)

    # check if file is not empty
    if(dim(annotation)[1] == 0){showModal(modalDialog( title = "Error message", "Empty file"))}
    validate({need(dim(annotation)[1] != 0, message="Empty file")  })

    # check if file is at least two columns
    if(dim(annotation)[2] < 2){ showModal(modalDialog(title = "Need at least two columns (genes names and corresponding ontology term)"))}
    validate({need(dim(annotation)[2] >= 2, message="Need at least two columns (genes names and corresponding ontology term)")})

    local.rea.values$settings_custom_ok <- TRUE
  })
  
  output$AnnotParamCPRUI3 <- renderUI({
    if(is.null(local.rea.values$settings_custom_ok)) return()
    
    annotation <- fread(file = input$annotationFileCPR$datapath, sep="\t", header = TRUE)
    
    box(title = span(tagList(icon("sliders"), "  ", "Chose columns names")), width = 14, status = "warning",
        
        # Select the right columns for the analysis
        pickerInput(inputId = session$ns("col_geneName"), label = "Genes Name:",
                    choices = c("", colnames(annotation)),
                    selected = ""),
        pickerInput(inputId = session$ns("col_termID"), label = "Terms IDs:",
                    choices = c("", colnames(annotation)),
                    selected = ""),
        pickerInput(inputId = session$ns("col_termName"), label = "Term Names:",
                    choices = c("", colnames(annotation)),
                    selected = ""),
        pickerInput(inputId = session$ns("col_domain"), label = "Domain/Ontology:",
                    choices = c("", colnames(annotation)),
                    selected = ""),
        
        actionButton(inputId = session$ns("runEnrich_CPR"), label = "Run")
    )
  })
  
  # ---- Coseq UI if there are results ----
  output$GeneList.coseqCPRUI <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnal == FALSE) return()
    
    ListNames.coseq <- names(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$CoExpAnal[["clusters"]])
    
    pickerInput(
      inputId = session$ns("GeneList.coseq"), label = "Select Clusters :",
      choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
    )
  })
  
  # ---- Annotation function ----
  observeEvent(input$runEnrich_CPR, {
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Annot", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    local.rea.values$dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    
    # reset everything
    local.rea.values$list_arg <- list()
    local.rea.values$func_to_use <- NULL
    local.rea.values$results_enrich <- NULL
    local.rea.values$list_arg$keyType <- NULL
    local.rea.values$Domains <- NULL
    local.rea.values$domain <- NULL
    annotation <- NULL
    
    # ---- Checks: ----
    
    isGoAnnotation <- length(grep("GO", input$dom.select))!=0
    print(isGoAnnotation) # TODO Delete
    
    print(input$dom.select) # TODO Delete
    
    # load if needed
    if(isGoAnnotation) library(input$db.select, character.only = TRUE)
    if(input$dom.select == "custom")   annotation <- fread(file = input$annotationFileCPR$datapath, sep="\t", header = TRUE)
    
    # check list of genes
    if(length(c(input$GeneList.diff, input$GeneList.coseq )) == 0){
      showModal(modalDialog( title = "Error message", "Please select at least 1 gene list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff, input$GeneList.coseq)) != 0, message = "Please select at least 1 gene list") 
    })
    
    print("I validated everything until now") # TODO delete
    
    # check keytypes (GO and KEGG enrichment)
    if(input$dom.select == "KEGG"){
      if(!input$keytype %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")){
        showModal(modalDialog( title = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot"))
      }
      validate({ 
        need(input$keytype %in% c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"), message = "KEGG keytype must be one of kegg, ncbi-geneid, ncbi-proteinid or uniprot") 
      })
    }else if(isGoAnnotation){
      accepted_keytypes <- AnnotationDbi::keytypes(get(input$db.select))
      print(accepted_keytypes) # TODO delete
      if(!input$keytype %in% accepted_keytypes){
        showModal(modalDialog(title = paste("Keytype must be one of: ", paste(accepted_keytypes, collapse = ", "), sep = " ")))
      }
      validate({ 
        need(input$keytype %in% accepted_keytypes, message = paste("Keytype must be one of: ", paste(accepted_keytypes, collapse = ", "), sep = " ")) 
      })
    }
 
    print("I'm here, checked the keytypes!") # TODO delete
    
    # Chose function parameters for annotation enrichment (three cases)
    if(input$dom.select == "custom"){
      local.rea.values$func_to_use <- "enricher"
      
      message("Performing enrichment analysis with custom annotation file")
      
      annotation <- fread(file = input$annotationFileCPR$datapath, sep="\t", header = TRUE)
      
      if(input$col_geneName == "" || input$col_termID == ""){
        showModal(modalDialog( title = "Error message", "Please chose columns names for the gene names/ID and the ontology terms ID!"))
      }
      
      local.rea.values$list_arg$universe <- names(local.rea.values$dataset.SE)
      local.rea.values$list_arg$pvalueCutoff <- input$pValue
      local.rea.values$list_arg$qvalueCutoff <- 1 # no threshold on qvalue
      local.rea.values$list_arg$minGSSize <- 10 # default in clusterprofiler
      local.rea.values$list_arg$maxGSSize <- 500 # default in clusterprofiler
      
      # Handling domains/ontologies in custom files
      if(input$col_domain == ""){
        annotation <- annotation %>% mutate("Domain" = "unique_domain")
        local.rea.values$Domains <- unique(annotation %>% dplyr::select(Domain))
      }else{
        local.rea.values$Domains <-  unlist(unique(annotation %>% dplyr::select(matches(input$col_domain))))
        if(length(grep("", local.rea.values$Domains))!=0)  local.rea.values$Domains <- local.rea.values$Domains[-which(local.rea.values$Domains=="")]
      }  
      
    }else{
      
      local.rea.values$domain <- input$dom.select
      
      if(length(grep('GO:', input$dom.select))!=0){
        local.rea.values$domain <- str_split(input$dom.select, ":")[[1]][1]
        local.rea.values$Domains <- str_split(input$dom.select, ":")[[1]][2]  # unique value
      }
      
      if(isGoAnnotation){
        # library(input$db.select, character.only = TRUE)
        local.rea.values$list_arg$OrgDb <- input$db.select
        if(length(grep(":", input$dom.select)) == 0){
          local.rea.values$Domains <- c("MF", "CC", "BP")
        }
      }
      
      if(length(grep('KEGG', input$dom.select))!=0){
        local.rea.values$list_arg$organism <- input$KEGG_org
        local.rea.values$Domains <- "KEGG" # unique value
      }
      
      local.rea.values$func_to_use <- paste0("enrich", local.rea.values$domain)
      
      local.rea.values$list_arg$universe <- names(local.rea.values$dataset.SE)
      local.rea.values$list_arg$keyType <- input$keytype
      local.rea.values$list_arg$pvalueCutoff <- input$pValue
      local.rea.values$list_arg$qvalueCutoff <- 1 # no threshold on qvalue
      local.rea.values$list_arg$minGSSize <- 10 # default in clusterprofiler
      local.rea.values$list_arg$maxGSSize <- 500 # default in clusterprofiler
      
    }
    
    print(paste("# 11- Enrichment Analysis...", dataset))
    
    #---- progress bar ----#
    progress$inc(1/5, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#
    
    ## run annotation
    
    # ---- Annotation on diff results: ----  
    if(length(input$GeneList.diff) != 0){
      
      ## log2FC of the genes (put colors on graphs) 
      local.rea.values$log2FC_lists <- lapply(input$GeneList.diff, function(listname){
        vect <- local.rea.values$dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[listname]][["logFC"]]
        names(vect) <- rownames(local.rea.values$dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[listname]])
        return(vect)
      })
      names(local.rea.values$log2FC_lists) <- input$GeneList.diff
      
      local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]] <- runAnnotationEnrichment_CPR(local.rea.values$dataset.SE,
                                                                                                         func_to_use = local.rea.values$func_to_use,
                                                                                                         ListNames = input$GeneList.diff,
                                                                                                         list_args = local.rea.values$list_arg,
                                                                                                         from = "DiffExpAnal",
                                                                                                         Domains = local.rea.values$Domains,
                                                                                                         dom.select = input$dom.select,
                                                                                                         col_termID = ifelse(input$dom.select == "custom", input$col_termID, ""),
                                                                                                         col_geneName = ifelse(input$dom.select == "custom", input$col_geneName, ""),
                                                                                                         col_termName = ifelse(input$dom.select == "custom", input$col_termName, ""),
                                                                                                         col_domain = ifelse(input$dom.select == "custom", input$col_domain, ""),
                                                                                                         annotationPath = input$annotationFileCPR$datapath)
      
      local.rea.values$resAnnot <- NULL
      
      rea.values[[dataset]]$diffAnnot <- TRUE
    }
    
    # ---- Annotation on COEXP results: ----  
    if(length(input$GeneList.coseq) != 0){
      
      local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]] <- runAnnotationEnrichment_CPR(local.rea.values$dataset.SE,
                                                                                                       func_to_use = local.rea.values$func_to_use,
                                                                                                       ListNames = input$GeneList.coseq,
                                                                                                       list_args = local.rea.values$list_arg,
                                                                                                       from = "CoExpAnal",
                                                                                                       Domains = local.rea.values$Domains,
                                                                                                       dom.select = input$dom.select,
                                                                                                       col_termID = ifelse(input$dom.select == "custom", input$col_termID, ""),
                                                                                                       col_geneName = ifelse(input$dom.select == "custom", input$col_geneName, ""),
                                                                                                       col_termName = ifelse(input$dom.select == "custom", input$col_termName, ""),
                                                                                                       col_domain = ifelse(input$dom.select == "custom", input$col_domain, ""),
                                                                                                       annotationPath = input$annotationFileCPR$datapath)
      
      rea.values[[dataset]]$coExpAnnot <- TRUE
    }
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]] <- local.rea.values$dataset.SE
    
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#
  })
  
  # ---- Figures Code : ----
  ## print results
  output$AnnotDiffResultsCPR <- renderUI({
    
    if(rea.values[[dataset]]$diffAnnot == FALSE) return()
    
    # check if result is empty
    if(is.null(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]])){
      
      showModal(modalDialog( title = "Error message", "There is no correspondence between the feature lists and the annotation file"))
    }
    validate({ 
      need(!is.null(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]]), message="There is no correspondence between the feature lists and the annotation file") 
    })
    
    # foreach gene list selected (contrast)
    lapply(names(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]]), function(listname) {
      
      # if result is empty
      if(is.null(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]][[listname]])){
        
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
            "There is no correspondence between the feature list and the annotation file!"
        )
        
      }else{
        
        # display result for each list
        data <- local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]][[listname]]
        log2FC_vect <- local.rea.values$log2FC_lists[[listname]]
        
        fluidRow(
          
          box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
              
              tabsetPanel(
                # ---- Tab Panel : dotPlot : ----
                tabPanel("Result Table",
                         # hr(),
                         verticalLayout(
                           
                           tags$br(),
                           DT::renderDataTable({
                             dataPlot <- data[[input[[paste0(listname, "-domain_DT")]]]]
                             
                             DT::datatable(dataPlot@result,
                                           rownames = FALSE,
                                           options = list(rownames = FALSE, pageLength = 5, lengthMenu = c(5, 10, 15, 20), scrollX = T))
                             
                           }),
                           fluidRow(
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain_DT")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                           )),
                         
                         
                ),
                # ---- Tab Panel : dotPlot : ----
                tabPanel("DotPlot",
                         
                         verticalLayout(
                           renderPlot({ 
                             
                             dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                             NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over")]])  
                             Categories <- dataPlot@result$Description[1:NbtoPlot]
                             if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                             
                             dotplot(dataPlot, showCategory = Categories)
                             
                           }),
                           fluidRow(
                             column(3,
                                    numericInput(inputId = session$ns(paste0(listname, "-top.over")), label="Top terms:" , value=15 , 
                                                 min = 0, max=10000, step = 5)), # max code en dur : pas bien
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                             column(4,
                                    textInput(inputId = session$ns(paste0(listname, "-grep")), label="Search Expression")),
                           )),
                         
                ), 
                # ---- Tab Panel : heatplot : ----
                tabPanel("Heatplot",
                         verticalLayout(
                           renderPlot({ 
                             
                             dataPlot <- data[[input[[paste0(listname, "-domain_heat")]]]]
                             NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over_heat")]])  
                             Categories <- dataPlot@result$Description[1:NbtoPlot]
                             if(input[[paste0(listname, "-grep_heat")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep_heat")]]), toupper(Categories))]
                             
                             heatplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect) + 
                               labs(fill="log2FC") + 
                               scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 
                             
                           }),
                           fluidRow(
                             column(3,
                                    numericInput(inputId = session$ns(paste0(listname, "-top.over_heat")), label="Top terms:" , value=15 , 
                                                 min = 0, max=10000, step = 5)), # max code en dur : pas bien
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain_heat")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                             column(4,
                                    textInput(inputId = session$ns(paste0(listname, "-grep_heat")), label="Search Expression")),
                           )),   
                         
                ) ,
                # ---- Tab Panel : cnetplot : ----
                tabPanel("Cnetplot",
                         verticalLayout(
                           renderPlot({ 
                             
                             dataPlot <- data[[input[[paste0(listname, "-domain_cnet")]]]]
                             dataTab <- dataPlot@result[dataPlot@result$p.adjust < input$pValue, ]
                             NbtoPlot <- min(nrow(dataTab),input[[paste0(listname, "-top.over_cnet")]])  
                             Categories <- dataTab$Description[1:NbtoPlot]
                             if(input[[paste0(listname, "-grep_cnet")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep_cnet")]]), toupper(Categories))]
                             
                             node_label_arg <- "none"
                             if(input[[paste0(listname, "-genesLabels_cnet")]] && input[[paste0(listname, "-termsLabels_cnet")]]){
                               node_label_arg <- "all"
                             }else if(input[[paste0(listname, "-genesLabels_cnet")]] && !input[[paste0(listname, "-termsLabels_cnet")]]){
                               node_label_arg <- "gene"
                             }else if(input[[paste0(listname, "-termsLabels_cnet")]] && !input[[paste0(listname, "-genesLabels_cnet")]]){
                               node_label_arg <- "category"
                             }
                             
                             cnetplot(dataPlot, showCategory = Categories, foldChange = log2FC_vect, node_label = node_label_arg) + 
                               guides(colour=guide_colourbar(title = "log2FC")) +
                               scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 
                             
                           }),
                           fluidRow(
                             column(3,
                                    numericInput(inputId = session$ns(paste0(listname, "-top.over_cnet")), label="Top terms:" , value=15 , 
                                                 min = 0, max=10000, step = 5)), # max code en dur : pas bien
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain_cnet")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                             column(1,
                                    checkboxInput(inputId = session$ns(paste0(listname, "-genesLabels_cnet")), label = "Genes Labels", value = TRUE),
                                    checkboxInput(inputId = session$ns(paste0(listname, "-termsLabels_cnet")), label = "Terms Labels", value = TRUE)
                             ),
                             column(4,
                                    textInput(inputId = session$ns(paste0(listname, "-grep_cnet")), label="Search Expression")),
                           )
                         )),   
                
                # ),  
              )# TabsetPanel
          ) # box
        ) # fluidrow
      } # if/else
    })# lapply
    # ) 
  }) 
  
  
  # ---- Figure Code for Coexpr results: ----
  output$AnnotCoExpResultsCPR <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnnot == FALSE) return()
    
    #dataset.SE <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]
    
    
    # check if result is empty
    if(is.null(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]])){
      
      showModal(modalDialog( title = "Error message", "No enrichment found"))
    }
    validate({ 
      need(!is.null(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]]), message="No enrichment found") 
    })
    
    # foreach gene list selected (contrast)
    lapply(names(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]]), function(listname) {
      
      # if result is empty
      if(is.null(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]][[listname]])){
        
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
            "No enrichment found"
        )
        
      }else{
        # display result for each list
        data <- local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]][[listname]]
        log2FC_vect <- local.rea.values$log2FC_lists[[listname]]
        
        fluidRow(
          
          box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
              
              tabsetPanel(
                # ---- Tab Panel : dotPlot : ----
                tabPanel("Result Table",
                         # hr(),
                         verticalLayout(
                           
                           tags$br(),
                           DT::renderDataTable({
                             dataPlot <- data[[input[[paste0(listname, "-domain_DT")]]]]
                             
                             DT::datatable(dataPlot@result,
                                           rownames = FALSE,
                                           options = list(rownames = FALSE, pageLength = 5, lengthMenu = c(5, 10, 15, 20), scrollX = T))
                             
                           }),
                           fluidRow(
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain_DT")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                           )),
                         
                         
                ),
                # ---- Tab Panel : dotPlot : ----
                tabPanel("DotPlot",
                         
                         verticalLayout(
                           renderPlot({ 
                             
                             dataPlot <- data[[input[[paste0(listname, "-domain")]]]]
                             NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over")]])  
                             Categories <- dataPlot@result$Description[1:NbtoPlot]
                             if(input[[paste0(listname, "-grep")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep")]]), toupper(Categories))]
                             
                             dotplot(dataPlot, showCategory = Categories)
                             
                           }),
                           fluidRow(
                             column(3,
                                    numericInput(inputId = session$ns(paste0(listname, "-top.over")), label="Top terms:" , value=15 , 
                                                 min = 0, max=10000, step = 5)), # max code en dur : pas bien
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                             column(4,
                                    textInput(inputId = session$ns(paste0(listname, "-grep")), label="Search Expression")),
                           )),
                         
                ), 
                # ---- Tab Panel : heatplot : ----
                tabPanel("Heatplot",
                         verticalLayout(
                           renderPlot({ 
                             
                             dataPlot <- data[[input[[paste0(listname, "-domain_heat")]]]]
                             NbtoPlot <- min(nrow(dataPlot@result),input[[paste0(listname, "-top.over_heat")]])  
                             Categories <- dataPlot@result$Description[1:NbtoPlot]
                             if(input[[paste0(listname, "-grep_heat")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep_heat")]]), toupper(Categories))]
                             
                             heatplot(dataPlot, showCategory = Categories) 
                             
                           }),
                           fluidRow(
                             column(3,
                                    numericInput(inputId = session$ns(paste0(listname, "-top.over_heat")), label="Top terms:" , value=15 , 
                                                 min = 0, max=10000, step = 5)), # max code en dur : pas bien
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain_heat")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                             column(4,
                                    textInput(inputId = session$ns(paste0(listname, "-grep_heat")), label="Search Expression")),
                           )),   
                         
                ) ,
                # ---- Tab Panel : cnetplot : ----
                tabPanel("Cnetplot",
                         verticalLayout(
                           renderPlot({ 
                             
                             dataPlot <- data[[input[[paste0(listname, "-domain_cnet")]]]]
                             dataTab <- dataPlot@result[dataPlot@result$p.adjust < input$pValue, ]
                             NbtoPlot <- min(nrow(dataTab),input[[paste0(listname, "-top.over_cnet")]])  
                             Categories <- dataTab$Description[1:NbtoPlot]
                             if(input[[paste0(listname, "-grep_cnet")]]!="") Categories <- Categories[grep(toupper(input[[paste0(listname, "-grep_cnet")]]), toupper(Categories))]
                             
                             node_label_arg <- "none"
                             if(input[[paste0(listname, "-genesLabels_cnet")]] && input[[paste0(listname, "-termsLabels_cnet")]]){
                               node_label_arg <- "all"
                             }else if(input[[paste0(listname, "-genesLabels_cnet")]] && !input[[paste0(listname, "-termsLabels_cnet")]]){
                               node_label_arg <- "gene"
                             }else if(input[[paste0(listname, "-termsLabels_cnet")]] && !input[[paste0(listname, "-genesLabels_cnet")]]){
                               node_label_arg <- "category"
                             }
                             
                             cnetplot(dataPlot, showCategory = Categories, node_label = node_label_arg) 
                             
                           }),
                           fluidRow(
                             column(3,
                                    numericInput(inputId = session$ns(paste0(listname, "-top.over_cnet")), label="Top terms:" , value=15 , 
                                                 min = 0, max=10000, step = 5)), # max code en dur : pas bien
                             column(4,
                                    radioButtons(inputId = session$ns(paste0(listname, "-domain_cnet")), label="Domain",
                                                 choices = names(data), selected = names(data)[1], inline = FALSE, width = 1.5)),
                             column(1,
                                    checkboxInput(inputId = session$ns(paste0(listname, "-genesLabels_cnet")), label = "Genes Labels", value = FALSE),
                                    checkboxInput(inputId = session$ns(paste0(listname, "-termsLabels_cnet")), label = "Terms Labels", value = TRUE)
                             ),
                             column(4,
                                    textInput(inputId = session$ns(paste0(listname, "-grep_cnet")), label="Search Expression")),
                           )
                         )),   
                
                # ),  
              )# TabsetPanel
          ) # box
        ) # fluidrow
      } # if/else
    })
  })
  
}

