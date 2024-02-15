##########################################
# module 06 - ANNOTATION & ENRICHMENT
##########################################
#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput 
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom colourpicker colourInput
#' @importFrom magrittr "%>%"

# ---- Module enrichment : main -----

#' @keywords internal
.modEnrichmentUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    fluidRow(
      box(title = span(tagList(icon("chart-pie"), " ", 
                               a("ClusterProfiler/", 
                                 href = "https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html"), 
                               a("Pathview",         
                                 href = "https://bioconductor.org/packages/devel/bioc/vignettes/pathview/inst/doc/pathview.pdf"), 
                               "   " , 
                               tags$small("(Scroll down for instructions)"))),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
          div(  
            p("Analyses in this module are conducted using the clusterprofiler R-package. 
              If you have more questions or interest in this package,
              please check the associated paper or the online vignette at
              https://yulab-smu.top/biomedical-knowledge-mining-book/index.html."),
            p(""),
            h4(tags$span("Parameters set up:", style = "color:orange")),
            p("Choose the lists of omics features you want to run the enrichment for.
              Default option selects all the available lists (contrasts or co-expression clusters)."),
            p("Then choose the ontology you want to refer to for the analysis.
              Multiple choices are not allowed. 
            If you select custom, you'll have to enter an annotation file with at least two columns : 
              the names of the entity (same as the rownames of your dataset) and 
            an id for an ontology term (eg: GO:0030198). 
              It can also contains a column for a more explicit name for the term 
            (eg: extracellular matrix organization) 
              and a column for specifying the domain (eg: MF, BP or CC). 
              An enrichment analysis will be ran on each specified domain."),
            p("You will have to specify the names of the columns after validating the annotation file."),
            p("If you choose GO, the three GO:BP, GO:MF and GO:CC will be analyzed. 
              You can chose to only analyze one of them by selecting the wanted ontology domain).
              It requires to indicate an R-package for the database, in the form of org.*db.
              (eg: org.At.tair.db for Arabidopsis thaliana), 
              and to specify which type of identifier is used in the data 
              (eg: TAIR for arabidopsis)."),
            p("For KEGG analysis, only four identifiers are possible. 
              Please check if the rownames correspond to one of them
              (eg: TAIR is also kegg identifiers)."),
            p("KEGG analysis uses an API to search for pathways online,
              it requires to have access to an internet connection."),
            p("Set the adjusted pvalue threshold. 
              Only results below this threshold will be displayed."),
            
            h4(tags$span("Outputs:", style = "color:orange")),
            p("For each list, either contrast results or co-expression cluster,
              multiple table and plots are displayed."),
            p("- Overview: shows all results for all ontology domain 
              (usefull when multiple domains)
              for all contrasts or all clusters. The blue line
              indicates the current results"),
            p("- Results table: for the current list, 
              all terms that passed the adjusted pvalue threshold, 
              ordered by increasing adjusted pvalue. You can change the domain 
              at the bottom of the table.
              You can also order the table with the columns or search 
              for a particular expression"),
            p("- Dotplot: for each domain, shows the 15 (default) 
              first terms by adj. pvalue. 
              You can also change the domain or search for a particular expression. 
              When searching for the expression, 
              you might want to increase the number of terms to consider 
              (if there is enough that passed the threshold)"),
            p("- Heatplot: for each domain, shows the 15 (default) first terms 
              and the features that are both part of the list and the pathway 
              in the form of a heatmap. 
            For the contrasts lists, colors indicate the log2FC of the features 
            as found in the differential analysis.
            You can change the domain and search for a particular expression 
              (adjusting the number of terms to consider if you want to check 
            further on the list of terms)"),
            p("- cnetplot: for each domain, shows the 15 (default) 
              first terms and the features that are both part of the list and 
              the pathway in the form of a network. 
            As for the heatplot, only the contrasts list have colors, according 
            to the log2FC of each features.
              Default only shows the terms labels, you can turn on the features
            names as well (it can be unreadable). 
            You can also search for a particular expression."),
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
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, 
           "Please run a differential analysis and validate your choices.")
    )
    
    omicsType <- getOmicsTypes(session$userData$FlomicsMultiAssay[[dataset]])
    
    tabPanel.list <- list(tabPanel(title = "Custom annotation", 
                                   .modEnrichmentDBUI(id = ns("custom"))))
    
    if (omicsType %in% c("RNAseq", "proteomics")) {
      
      tabPanel.list <- c(list(
        tabPanel(title = "Gene Onlology database", 
                 .modEnrichmentDBUI(id = ns("GO"))),
        tabPanel(title = "KEGG database", 
                 .modEnrichmentDBUI(id = ns("KEGG")))),
        tabPanel.list )
    }
    do.call(what = tabsetPanel, args = tabPanel.list)
  })
  
  # SERVEUR
  callModule(module = .modEnrichmentDB, id = "GO", 
             dataset = dataset, database = "GO", rea.values)
  callModule(module = .modEnrichmentDB, id = "KEGG", 
             dataset = dataset, database = "KEGG", rea.values)
  callModule(module = .modEnrichmentDB, id = "custom", 
             dataset = dataset, database = "custom", rea.values)
  
}

# ---- Module enrichment : per database ----

#' @title .modEnrichmentDBUI
#' @keywords internal
#' @noRd
.modEnrichmentDBUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
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
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table fread
#' @keywords internal
#' @noRd
.modEnrichmentDB <- function(input, output, session, 
                             dataset, database, rea.values){
  
  ns <- session$ns
  local.rea.values <- reactiveValues(dataPathAnnot = NULL)
  
  # UI
  
  ## setting
  output$ParamDB_UI <- renderUI({
    
    switch(database,
           "GO"     = { uiOutput(ns("AnnotParamGO_UI")) },
           "KEGG"   = { uiOutput(ns("AnnotParamKEGG_UI")) },
           "custom" = { uiOutput(ns("AnnotParamCustom_UI")) }
    )
    
  })
  
  ## display results of enrichment analysis for selected database
  output$tabsetRes_UI <- renderUI({
    
    ### for diff analysis results
    tabPanel.list <- list(
      tabPanel(title = "Enrichment from differential expression analysis", 
               br(),
               .modRunEnrichmentUI(id = ns("DiffExpEnrichAnal")))
    )
    
    ### for coexp analysis results if exist
    if (rea.values[[dataset]]$coExpAnal) { 
      tabPanel.list <- c(tabPanel.list,
                         list(tabPanel(title = "Enrichment from co-expression analysis",
                                       br(),
                                       .modRunEnrichmentUI(id = ns("CoExpEnrichAnal")))
                         )
      )
    }
    
    do.call(what = tabsetPanel, args = tabPanel.list )
  })
  
  ## GO settings
  output$AnnotParamGO_UI <- .annotParamGO(session, rea.values, dataset)
  
  ## KEGG settings
  output$AnnotParamKEGG_UI <- .annotParamKEGG(session, rea.values, dataset)  
  
  ## custom settings
  output$AnnotParamCustom_UI <- .annotParamCustom(session, rea.values, dataset) 
  
  ### select columns
  observeEvent(input$annotationFileCPR, {
    
    local.rea.values$dataPathAnnot <- NULL
    local.rea.values$dataPathAnnot <- input$annotationFileCPR$datapath
    output$selectColumnsCustom <- .annotFileColumns(session, local.rea.values,
                                                    dataset)
    
  })
  
  ### use example
  observeEvent(input$useExampleFile, {
    
    local.rea.values$dataPathAnnot <- NULL
    local.rea.values$dataPathAnnot <- paste0(system.file(package = "RFLOMICS"), 
                                             "/ExamplesFiles/ecoseed/AT_GOterm_EnsemblPlants.txt")
    output$selectColumnsCustom <- .annotExFileColumns(session, local.rea.values,
                                                      dataset)
    
  })
  
  # SERVEUR
  callModule(module     = .modRunEnrichment, 
             id         = "DiffExpEnrichAnal", 
             dataset    = dataset,
             database   = database, 
             listSource = "DiffExpEnrichAnal",
             rea.values = rea.values, 
             input2     = input, 
             local.rea.values = local.rea.values)
  
  callModule(module     = .modRunEnrichment, 
             id         = "CoExpEnrichAnal", 
             dataset    = dataset,
             database   = database, 
             listSource = "CoExpEnrichAnal",
             rea.values = rea.values, 
             input2     = input, 
             local.rea.values = local.rea.values)
  
}

# ---- module enrichment : run per analysis type ----

#' @title .modRunEnrichmentUI
#' @keywords internal
#' @noRd
.modRunEnrichmentUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("settings")),
    uiOutput(ns("enrichSummary")),
    uiOutput(ns("AnnotResults"))
  )
}

#' @title .modRunEnrichment
#' @keywords internal
#' @importFrom DT renderDataTable datatable
#' @importFrom grid unit
#' @importFrom RCurl url.exists
#' @importFrom AnnotationDbi keytypes
#' @importMethodsFrom ComplexHeatmap draw
#' @noRd
.modRunEnrichment <- function(input, output, session, dataset, database, 
                              listSource, rea.values, input2, local.rea.values){
  ns <- session$ns
  
  switch(listSource,
         "DiffExpEnrichAnal" = {
           title <- "Summary"
           fromAnnot <- "diffAnnot"
           from <- "DiffExpAnal"
         },
         "CoExpEnrichAnal" = {
           title <- "Summary"
           fromAnnot <- "coExpAnnot"
           from <- "CoExpAnal"
         }
  )
  
  # UI
  ## setting
  output$settings <- .annotSettings(session, rea.values, dataset, listSource)
  
  ## display results
  output$enrichSummary <- .outEnrichSummary(session, rea.values, input,
                                            fromAnnot, title, from, listSource,
                                            dataset, database)
  output$AnnotResults <- .outAnnotResults(session, input, input2,
                                          rea.values, fromAnnot,
                                          dataset, from, listSource, database)
  
  # SERVER
  
  ## run enrichment
  # ---- run Annotation  ----
  observeEvent(input$run, {
    
    #---- progress bar ----#
    progress <- Progress$new()
    progress$set(message = "Run annotation enrichment", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "Checking everything")
    #-----------------------#
    
    # check list of genes
    condition <- length(input$listToAnnot) != 0
    messCond <- "Please select at least one contrast name."
    if (!condition) {
      showModal(modalDialog(title = "No list selected", messCond))
    }
    validate({need(condition, message = messCond)})
    
    ### GO checks ###
    if (database == "GO") {
      
      condition <- input2$db.select != ""
      messCond <- "Select a database for the analysis. If none is suggested in 
      the list, please check and install the appropriate org.db package
      from Bioconductor."
      # check parameters
      if (!condition) {
        showModal(modalDialog(title = "No database selected", messCond))
      }
      validate({need(condition, messCond)})
      
      textToEval <- paste0("keytypes(", input2$db.select, "::", 
                           input2$db.select, ")")
      accepted_keytypes <- eval(parse(text = textToEval))
      
      # Check keytypes for GO
      condition <- input2$keytype.go %in% accepted_keytypes
      messCond <- paste("Keytype must be one of: ",
                        paste(accepted_keytypes, collapse = ", "),
                        sep = " ")
      if (!condition) {
        showModal(modalDialog(title = "Error in argument", messCond))
      }
      validate({need(condition, message = messCond)})
      
      # set list args
      database <- "GO"
      list_args <- list()
      list_args[["OrgDb"]]        <- input2$db.select
      list_args[["keyType"]]      <- input2$keytype.go
      list_args[["pvalueCutoff"]] <- input2$pValue_GO
      
      domain <- input2$ont.select
      if (input2$ont.select == "ALL") domain <- c("MF", "BP", "CC")
      
      paramList <- list(
        database     = database,
        list_args    = list_args,
        domain       = domain
      )
    }
    
    ### KEGG checks ###
    if (database == "KEGG") {
      
      # check param
      ## ID key
      condition <- input2$keyTypeKEGG %in% c("kegg", "ncbi-geneid",
                                             "ncib-proteinid", "uniprot")
      messCond  <- "Please chose an appriopriate
      KEGG keytype. It must be one of kegg, ncbi-geneid, 
          ncbi-proteinid or uniprot"
      
      if (!condition) {
        showModal(modalDialog(title = "Error in keytype", messCond))
      }
      validate({need(condition, message = messCond) })
      
      ## code organism KEGG
      condition <- input2$organism != ""
      messCond  <- "It seems you forgot to enter the KEGG organism code 
      (three or more letters, eg: hsa for human or 
      ath for arabidopsis thaliana)"
      if (!condition) {
        showModal(modalDialog(title = "Missing organism code", messCond))
      }
      validate({need(condition, message = messCond)})
      
      # set list args
      database <- "KEGG"
      list_args <- list()
      list_args[["organism"]] <- input2$organism
      list_args[["keyType"]]  <- input2$keyTypeKEGG
      list_args[["pvalueCutoff"]] <- input2$pValueKEGG
      
      # prevent multiple execution
      paramList <- list(
        list_args    = list_args,
        database     = database,
        domain       = "no-domain"
      )
    }
    
    ### Custom checks ###
    if (database == "custom") {
      
      # check param
      ## if custom annotation file
      condition <- !is.null(local.rea.values$dataPathAnnot)
      messCond  <- "Please load a custom annotation file."
      if (!condition) {
        showModal(modalDialog(title = "Missing annotation file", messCond))
      }
      validate({need(condition,  message = messCond)})
      
      # Check some custom elements and annotation file
      condition <- input2$col_geneName != "" && input2$col_termID != ""
      messCond <- "Please choose the columns names for the omics features  
                               names/ID and the ontology terms ID."
      if (!condition) {
        showModal(modalDialog(title = "Missing mandatory settings", messCond))
      }
      validate({need(condition, message = messCond)})
      
      # set param
      domain <- NULL
      annotation2 <- NULL
      col_domain_arg <- NULL
      list_args <- list()
      list_args[["pvalueCutoff"]] <- input2$pValue_custom
      
      annotation <- fread(file = local.rea.values$dataPathAnnot,
                          sep = "\t", header = TRUE)
      listCharRm <- c(".", " ", "")
      annotation2 <- list()
      
      # if column with domain exist
      if (input2$col_domain != "") {
        
        # remove column with domain == c(".", " ", "")
        domainRm <- intersect(annotation[[input2$col_domain]], listCharRm)
        if (length(domainRm) != 0) {
          annotation <- filter(annotation,
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
      paramList <- list(
        list_args    = list_args,
        database     = "custom",
        domain       = domain,
        annot        = annotation2,
        col_domain   = col_domain_arg)
    }
    
    # prevent multiple execution
    #if (.checkRunORAExecution(session$userData$FlomicsMultiAssay[[dataset]],
    # database, paramList) == FALSE) return()
    
    # ---- Annotation on diff results: ----
    # run analysis
    message("# 11- ", database, 
            " Enrichment Analysis of ", listSource ,
            " lists...", dataset)
    
    #---- progress bar ----#
    progress$inc(5/10, detail = paste("Run analyses ", 50, "%", sep = ""))
    #----------------------#
    
    session$userData$FlomicsMultiAssay[[dataset]]@metadata[[listSource]][[database]] <- NULL
    rea.values[[dataset]][[fromAnnot]] <- FALSE
    
    # run annotation diff analysis
    runRes <- tryCatch({
      paramList <- c(paramList, 
                     list(object = session$userData$FlomicsMultiAssay[[dataset]], 
                          from = from, nameList = input$listToAnnot))
      
      do.call(runAnnotationEnrichment, paramList)
    },
    warning = function(war) return(war),
    error   = function(err) return(err)
    )
    
    condition <- is(runRes, "RflomicsSE")
    if (!condition) {
      showModal(modalDialog(title = paste("Something went wrong: ", 
                                          runRes$message)))
    }
    validate({need(condition,
                   message = paste0("Something went wrong: ", 
                                    runRes$message))})
    
    session$userData$FlomicsMultiAssay[[dataset]] <-  runRes
    rea.values[[dataset]][[fromAnnot]] <- TRUE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("All done", 100, "%", sep = ""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
  
}

# ---- check run enrichement analysis execution (à adapter) ----
#' @title .checkRunORAExecution
#' @noRd
#' @keywords internal
#' @importFrom dplyr setequal
.checkRunORAExecution <- function(object.SE, dbType = NULL, paramList = NULL){
  
  
  if (!is.null(paramList$diffList)) {
    
    if (is.null(object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]])) return(TRUE)
    if (isFALSE(setequal(names(object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$enrichResult), paramList$diffList)) ) return(TRUE)
    
    # check param
    if (isFALSE(setequal(paramList$domain, object.SE@metadata[["DiffExpEnrichAnal"]][[dbType]]$list_args$domain))) return(TRUE)
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
    if (isFALSE(setequal(names(object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$enrichResult), paramList$CoexpList)) ) return(TRUE)
    
    # check param
    if (isFALSE(setequal(paramList$domain, object.SE@metadata[["CoExpEnrichAnal"]][[dbType]]$list_args$domain))) return(TRUE)
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



# ----- explain content ----
# functions for popify/bsshiny content


# ---- UI param functions ----
#' @noRd
#' @keywords internal
.annotSettings <- function(session, rea.values, dataset, 
                           listSource){
  
  ns <- session$ns
  
  renderUI({
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, 
           "Please run a differential analysis and validate the results before 
           using this panel.")
    )
    
    ### lists
    ListNames <- switch(listSource,
                        "DiffExpEnrichAnal" = rea.values[[dataset]]$DiffValidContrast$contrastName,
                        "CoExpEnrichAnal"   = rea.values[[dataset]]$CoExpClusterNames)
    
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
      column(width = 6, 
             br(),
             actionButton(inputId = ns("run"), 
                          label = "Run over-representation analysis"))
    )
  })
}


#' @noRd
#' @keywords internal
.annotParamCustom <- function(session, rea.values, dataset){
  
  ns <- session$ns
  
  renderUI({
    
    if (rea.values[[dataset]]$diffValid == FALSE) return()
    
    list(
      
      fileInput(inputId = ns("annotationFileCPR"), 
                label = popify(
                  actionLink(inputId = "customAnnotFilepop",
                             label = "Annotation file"),
                  title = "Annotation File format",
                  content = paste0("<p>Required format:</p>",
                                   "a tsv or csv file with at least two columns, ",
                                   "one for omic features names and one for", 
                                   " their associated terms."),
                  trigger = "click"
                ),
                accept  = c("text/csv", 
                            "text/comma-separated-values,text/plain", 
                            ".csv")), 
      
      if (rea.values$exampleData) {
        popify(bsButton(inputId = ns("useExampleFile"),
                        label = "Use example annotation file",
                        style = "primary", size = "default",
                        type = "action"),
               title = "Use example file for ecoseed data",
               content = paste0("<p> You are conducting an analysis using the ",
                                "example dataset ecossed. You can ",
                                "run an annotation enrichment ",
                                "using the example annotation file, which is an ",
                                "extract of GO terms for Arabidopsis thaliana",
                                " genes from plant ensembl, using TAIR ids"),
               placement = "top", trigger = "hover")
      },
      
      uiOutput(ns("selectColumnsCustom")),
      
      ## alpha threshold
      numericInput(inputId = ns("pValue_custom"), 
                   label = "Adjusted p-value", 
                   value = 0.01 , min = 0, max = 1, step = 0.01)
    )
  })
}

#' @noRd
#' @keywords internal
.annotFileColumns <- function(session, local.rea.values, 
                              dataset){
  
  ns <- session$ns
  dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
  varLabel <- omicsDic(dataSE)$variableName
  varLabel <- paste0(toupper(substr(varLabel, 1,1)),
                     substr(varLabel, 2, nchar(varLabel)))
  
  renderUI({
    
    annotation <- fread(file = local.rea.values$dataPathAnnot, 
                        sep = "\t", header = TRUE)
    
    column(width = 12,
           br(),
           hr(style = "border-top: 1px solid #000000;"),
           popify(
             actionLink(inputId = ns(paste0("colSelAnnot")),
                        label = "Help", 
                        icon = icon("question-circle")), 
             title = "Custom annotation",
             content = paste0("<p>Chose the appropriate columns names. ",
                              "They are automatically derived from your file. ",
                              "Items with a star are mandatory. </p>",
                              "<p>- Term Name is used in case of a term name ",
                              "that is different from its id ",
                              "(eg GO:0009409 is the id, response to ",
                              "cold is the term name).</p>",
                              "- Domain is used to differentiate databases or ",
                              "ontologies in the same file (eg: BP, CC and ",
                              "MF would be indicated in the domain column)."),
             trigger = "click", placement = "top"
           ),
           
           # Select the right columns for the analysis
           pickerInput(inputId = ns("col_geneName"),
                       label = paste0(varLabel, " ID: *"),
                       choices = c("", colnames(annotation)),
                       selected = ""),
           pickerInput(inputId = ns("col_termID"), 
                       label = "Terms IDs: *",
                       choices = c("", colnames(annotation)),
                       selected = ""),
           pickerInput(inputId = ns("col_termName"), 
                       label = "Term Names:",
                       choices = c("", colnames(annotation)),
                       selected = ""),
           pickerInput(inputId = ns("col_domain"), 
                       label = "Domain:",
                       choices = c("", colnames(annotation)),
                       selected = ""),
           hr(style = "border-top: 1px solid #000000;"),
    )
  })
}  

#' @noRd
#' @keywords internal
.annotExFileColumns <- function(session, local.rea.values,
                                dataset){
  ns <- session$ns
  dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
  varLabel <- omicsDic(dataSE)$variableName
  varLabel <- paste0(toupper(substr(varLabel, 1,1)),
                     substr(varLabel, 2, nchar(varLabel)))
  
  renderUI({
    
    # annotation <- fread(file = local.rea.values$dataPathAnnot, sep = "\t",
    #  header = TRUE)
    
    column(width = 12,
           br(),
           hr(style = "border-top: 1px solid #000000;"),
           
           # Select the right columns for the analysis
           pickerInput(inputId = ns("col_geneName"), 
                       label = paste0(varLabel, " ID: *"),
                       choices = "Gene stable ID",
                       selected = "Gene stable ID"),
           pickerInput(inputId = ns("col_termID"), 
                       label = "Terms IDs: *",
                       choices = "GO term accession",
                       selected = "GO term accession"),
           pickerInput(inputId = ns("col_termName"), 
                       label = "Term Names:",
                       choices = "GO term name",
                       selected = "GO term name"),
           pickerInput(inputId = ns("col_domain"), 
                       label = "Domain:",
                       choices = "GO domain",
                       selected = "GO domain"),
           hr(style = "border-top: 1px solid #000000;"),
    )
  })
}

#' @noRd
#' @keywords internal
.annotParamGO <- function(session, rea.values, dataset){
  
  ns <- session$ns
  
  renderUI({
    
    if (rea.values[[dataset]]$diffValid == FALSE) return()
    
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
                   label = "Adjusted p-value",
                   value = 0.01 , min = 0, max = 1, step = 0.01)
    )
    
    #}
  })
}  

#' @noRd
#' @keywords internal
.annotParamKEGG <- function(session, rea.values, dataset){
  ns <- session$ns
  renderUI({
    
    if (rea.values[[dataset]]$diffValid == FALSE) return()
    
    list(
      pickerInput(
        inputId = ns("keyTypeKEGG"), 
        label = "Select id type:",
        choices = c("", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
        selected = "",
        options = list(`actions-box` = TRUE, size = 10, 
                       `selected-text-format` = "count > 3")),
      
      # type organism 
      textInput(inputId = ns("organism"), 
                label = popify(actionLink(inputId = ns("colSelAnnot"),
                                          label = "Organism code:"),
                               title = "KEGG Organism code",
                               content = paste0("Three or more letters ",
                                                "indicating the species. ",
                                                "For example: ",
                                                "<ul><li>hsa for <i>Homo sapiens</i>,</li> ",
                                                "<li>ath for <i>Arabidopsis thaliana</i>.</li></ul>")
                ),
                value = "", 
                width = NULL, placeholder = NULL),
      
      ## alpha threshold
      numericInput(inputId = ns("pValueKEGG"), 
                   label = "Adjusted p-value",
                   value = 0.01 , min = 0, max = 1, step = 0.01)
      
      
    )
  })
}

# ---- main output functions -----
#' @noRd
#' @keywords internal
.outAnnotResults <- function(session, input, input2,
                             rea.values, fromAnnot,
                             dataset, from, listSource, database){
  
  ns <- session$ns
  
  renderUI({
    
    dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    if (rea.values[[dataset]][[fromAnnot]] == FALSE ||
        is.null(sumORA(dataSE, from = listSource, database = database))) { 
      return()
    }
    
    #foreach genes list selected (contrast)
    lapply(names(getEnrichRes(dataSE, from = listSource, database = database)), function(listname) {
      if (length(getEnrichRes(dataSE, from = listSource, database = database, contrastName = listname)) != 0) {
        if (sum(unlist(sumORA(dataSE, from = listSource, database = database)[-1]), na.rm = TRUE) == 0) {
          
          fluidRow(
            box(width = 12, 
                title = paste0(listname, ": 0 enriched terms found"), 
                status = "danger")
          )
        }
        else{
          
          choices <- names(getEnrichRes(dataSE, from = listSource,
                                        contrastName = listname,
                                        database = database))
          # ---- TabPanel plots ----
          tabPanel.list <- list(
            tabPanel("DotPlot",
                     br(),
                     .outDotPlot(input, dataSE, listname, 
                                 from, "dotplot", database)
            ),
            tabPanel("Heatplot",
                     br(),
                     .outHeatPlot(input, dataSE, listname, 
                                  from, "heatplot", database)
            ),
            tabPanel("Cnetplot",
                     br(),
                     .outCnetPlot(session, input, dataSE, listname, 
                                  from, "cnetplot", database)
            ),
            # ---- Tab Panel : Results table : ----
            tabPanel("Result Table",
                     br(),
                     verticalLayout(
                       br(),
                       renderDataTable({
                         
                         dataPlot <- getEnrichRes(object = dataSE,
                                                  from = listSource,
                                                  contrastName = listname,
                                                  database = database,
                                                  domain = input[[paste0(listname, "-domain")]])
                         
                         pvalue <- getEnrichPvalue(dataSE,
                                                   from = listSource, 
                                                   database = database)
                         datPlot <- dataPlot@result[dataPlot@result$p.adjust <  pvalue,]
                         datPlot$pvalue <- round(datPlot$pvalue, 3)
                         datPlot$p.adjust <- round(datPlot$p.adjust, 3)
                         datPlot$qvalue <- round(datPlot$qvalue, 3)
                         datPlot$geneID <- unlist(lapply(datPlot$geneID, 
                                                         FUN = function(longString){
                           return(gsub("/", ", ", longString))
                         }))
                         
                         datatable(datPlot,
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
                                   contrastName = listname, 
                                   from = listSource,
                                   database = "KEGG")[["no-domain"]]@result
            
            pvalue <- getEnrichPvalue(dataSE,
                                      from = listSource,
                                      database = database)
            mapChoices <- sort(data$ID[data$p.adjust < pvalue])
            
            tabPanel.list <- c(tabPanel.list,
                               list(
                                 .outPathview(session, input,
                                              input2, data, dataSE,
                                              listSource, listname, mapChoices)
                               )
            )
          }
          
          ### 
 
          # display results
          fluidRow(
            box(width = 12, solidHeader = TRUE, collapsible = TRUE,
                collapsed = TRUE, status = "success", title = listname,
                
                fluidRow(
                  column(3,
                         # choice of domain
                         radioButtons(ns(paste0(listname, "-domain")),
                                      label = "Domain",
                                      choices = choices,
                                      selected = choices[1])),
                  column(3,
                         # number of term
                         numericInput(ns(paste0(listname, "-top.over")),
                                      label = "Top terms:", value = 15 ,
                                      min = 1, max = 10000, step = 5)), 
                  column(3,
                         # search term or gene
                         textInput(ns(paste0(listname, "-grep")),
                                   label = "Search Expression")),
                  column(1,
                  ),
                  .popoverHelp(label = "")
                ),
                fluidRow(
                  # display all tabsets
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

.popoverHelp <- function(label = ""){
  
  renderUI({
    column(1,
           br(),  
           tags$style(HTML(
             "#classPop .popover {
              text-align: left;
              border-color: black;
              background-color: #d9edf7;
              width: 600px;
              max-width: 700px;
              color: black;
              opacity: 1;}
              "
           )),
           div(id = "classPop",
               span(label, 
                    popify(h4(icon("question-circle"), "Help"), 
                           title = "Parameters for enrichment results plots", 
                           content = paste0("For each domain indicated ",
                                            "(typically BP, CC and MF for GO ,",
                                            "enrichment), three plots and the result ",
                                            "table are displayed. For more ",
                                            "information on the plots, please ",
                                            "check the package vignette. ",
                                            "<p> You can chose the number of terms ",
                                            "to display, ordered by increasing ",
                                            "adjusted pvalue.</p>",
                                            "<p> Search expression will allow you ",
                                            "to display only the top terms containing ",
                                            "the regular expression of interest. ",
                                            "You can use regular expression patterns:",
                                            "<ul><li> <b>||</b> for <i>or</i> statement; </li>",
                                            "<li> <b>&</b> for <i>and</i> statement;</li>",
                                            "</ul>"), trigger = "click", placement = "top"),
                    style = "color: #337ab7; border-color: #337ab7;")
           )
    )
  })
  
}

# ---- plots ----
#' @noRd
#' @keywords internal
.outHeatPlot <- function(input,
                         dataSE, 
                         listname, from, plotType, database){
  renderUI({
    outplot <- .doNotSpeak(plotClusterProfiler(dataSE,
                                               contrastName = listname,
                                               from = from,
                                               plotType = plotType,
                                               database = database,
                                               domain = input[[paste0(listname, "-domain")]],
                                               showCategory = input[[paste0(listname, "-top.over")]],
                                               searchExpr = input[[paste0(listname, "-grep")]]))
    
    
    if (is(outplot, "gg")) renderPlot({
      warnOpt <- getOption("warn")
      options(warn = -1) # avoid ggrepel warnings, or at least trying to
      suppressMessages(suppressWarnings(print(outplot)))
      options(warn = warnOpt)
    })
    else renderText({outplot$message})
  })
} 

#' @noRd
#' @keywords internal
.outDotPlot <- function(input, dataSE, listname, from, plotType, database){
  renderUI({
    outdot <- .doNotSpeak(plotClusterProfiler(dataSE,
                                              contrastName = listname,
                                              from = from,
                                              plotType = plotType,
                                              database = database,
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
  })
}

#' @noRd
#' @keywords internal
.outCnetPlot <- function(session, input, dataSE, listname, 
                         from, plotType, database){
  
  ns <- session$ns
  varLabel <- omicsDic(dataSE)$variableName
  varLabel <- paste0(toupper(substr(varLabel, 1,1)),
                     substr(varLabel, 2, nchar(varLabel)))
  
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
      
      outcnet <- .doNotSpeak(plotClusterProfiler(dataSE,
                                                 contrastName = listname,
                                                 from = from,
                                                 plotType = plotType,
                                                 database = database,
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
      
      column(2, checkboxInput(inputId = ns(paste0(listname, 
                                                  "-genesLabels_cnet")),
                              label = paste0(varLabel, " Labels"), 
                              value = FALSE)),
      column(2, checkboxInput(inputId = ns(paste0(listname, 
                                                  "-termsLabels_cnet")),
                              label = "Terms Labels", value = TRUE))
      
    )
  )
}

#' @noRd
#' @keywords internal
.outPathview <- function(session, input,
                         input2,
                         data,
                         dataSE,
                         listSource,
                         listname, mapChoices
){
  
  ns <- session$ns
  tabPanel("Pathview results",
           br(),
           fluidRow(
             column(3,
                    selectInput(
                      inputId = ns(paste0(listname, "-MAP.sel")), 
                      label = "Select map:",
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
                      if (url.exists(link_to_map)) {
                        
                        renderPlot({
                          plotKEGG(object = dataSE,
                                   contrastName = listname,
                                   from = listSource,
                                   pathway_id = input[[paste0(listname, "-MAP.sel")]],
                                   species = input2$organism,
                                   gene_idtype = input2$keyTypeKEGG
                          )}, res = 300, width = 1000, height = 1000)
                      } else {
                        renderText("Please check your connection. 
                                   It seems the URL does not exist, or you're not connected.")
                      }
                      
                    })
             )))
}

# ---- Summary ----
#' @noRd
#' @keywords internal
.outEnrichSummary <- function(session, rea.values, input,
                              fromAnnot, title, from, listSource,
                              dataset, database){
  
  ns <- session$ns
  
  renderUI({
    
    if (rea.values[[dataset]][[fromAnnot]] == FALSE ||
        is.null(sumORA(session$userData$FlomicsMultiAssay[[dataset]],
                       from     = listSource,
                       database = database)))
      return()
    
    results <- getEnrichRes(session$userData$FlomicsMultiAssay[[dataset]], 
                            listSource, database)
    
    if (is.null(sumORA(session$userData$FlomicsMultiAssay[[dataset]], 
                       listSource, database)))  {
      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,
            collapsed = TRUE, status = "warning", title = title,
            "There is no results for enrichment analysis! Check the ids."))
    } else {
      fluidRow(
        box(width = 12, solidHeader = TRUE, collapsible = TRUE,
            collapsed = TRUE, status = "primary", title = title,
            
            tabsetPanel(
              tabPanel(title = "table",
                       br(),
                       renderDataTable({
                         
                         datatable(
                           sumORA(session$userData$FlomicsMultiAssay[[dataset]],
                                  listSource, database),
                           rownames = FALSE,
                           options = list(pageLength = 6, dom = 'tip'))
                       })),
              tabPanel(title = "heatmap", 
                       br(),
                       .outCompResults(session, rea.values, input, 
                                       fromAnnot, 
                                       dataset, from,
                                       listSource, database))
            )
            
        )
        
      )
    }
  })
}

#' @description draw the heatmap for the enrichment summary
#' @noRd
#' @keywords internal
.outCompResults <- function(session, rea.values, input, 
                            fromAnnot, 
                            dataset, from,
                            listSource, database) {
  
  ns <- session$ns
  dataSE <- session$userData$FlomicsMultiAssay[[dataset]]
  renderUI({
    
    if (rea.values[[dataset]][[fromAnnot]] == FALSE || 
        is.null(sumORA(dataSE, listSource, database))) { 
      return()
    }
    
    # Possible domains of ontology:
    possDomain <-  unique(unlist(
      lapply(dataSE@metadata[[listSource]][[database]][["enrichResult"]], 
             FUN = function(x) names(x)))) 
    
    # display results
    verticalLayout(
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
                   outHeatmap <-  plotEnrichComp(dataSE,
                                                 from = from, 
                                                 database = database, 
                                                 domain = input[[paste0(database, "-compDomain")]],
                                                 matrixType = input[[paste0(database, "-compType")]])
                 },
                 warning = function(warn) warn,
                 error = function(err) err
                 )
                 
                 if (is(outHeatmap, "Heatmap")) {
                   renderPlot({
                     draw(outHeatmap, 
                          heatmap_legend_side = "top", 
                          padding = unit(5, "mm"),
                          gap = unit(2, "mm"))}, 
                     width = "auto", 
                     height = min(400 + nrow(outHeatmap@matrix)*50, 1000)
                   )
                 }else { renderText({outHeatmap$message}) }
               })
        ) # column
      )#,  # fluidrow
    )
  })
}
