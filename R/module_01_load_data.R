### ============================================================================
### [01_Load_Data] shiny modules
### ----------------------------------------------------------------------------

# ---- modLoadOmicsData UI ----
.modLoadOmicsDataUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      # load exp design
      # display tab
      box(width = 12, title = "Load experimental design", status = "warning", solidHeader = TRUE,
          # ExpDesign
          fluidRow(
            column(width = 12,
                   # 1- set project name
                   column(width = 3, textInput(inputId = ns("projectName"), label = "Project name")),
                   # 2- matrix count/abundance input
                   column(width = 7, fileInput(inputId = ns("Experimental.Design.file"), label = "Experimental design (tsv)",
                                               accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv"))),
                   column(width = 2, actionButton(inputId = ns("loadEx_metadata"), label = "Load Example Metadata", icon = shiny::icon("file-import"))),
            )
          ),
          fluidRow(
            column(width = 12,
                   uiOutput(ns("ExpDesignTable"))
            ),
            column(width = 12,
                   uiOutput(ns("selectSamplesUI"))),
          )
      )
    ),
    fluidRow( uiOutput(outputId = ns("LoadDataUI"))),
    fluidRow( column(width = 2, actionButton(inputId = ns("loadData"),"load Data")),
              column(width = 2, actionButton(inputId = ns("loadExData"),"load Example Data", 
                                             icon = shiny::icon("file-import")))),
    br(),
    
    fluidRow(
      uiOutput(ns("overViewUI"))
    )
    
  )
  
}

# ---- modLoadOmicsData SERVER ----
#' @importFrom stringr str_subset
.modLoadOmicsData <- function(input, output, session, rea.values){
  
  # ---- initialization ----
  local.rea.values <- reactiveValues(plots      = FALSE, 
                                     ExpDesign  = NULL, 
                                     #FactorList = NULL, 
                                     dataPath   = NULL, # chemin  
                                     omicsData  = NULL,
                                     omicsNames = NULL,
                                     omicsTypes = NULL)
  
  observe({
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    
    # rea.values$contrastMat  <- NULL
    session$userData$FlomicsMultiAssay <- NULL
    
    rea.values$validate.status <- 0
    
  })
  
  
  # ---- Load experimental design ----
  # load user own metadata file
  observeEvent(input$Experimental.Design.file, {
    rea.values$exampleData    <- FALSE
    local.rea.values$dataPath <- NULL
    local.rea.values$dataPath <- input$Experimental.Design.file$datapath
  })
  
  # load example metadata file
  observeEvent(input$loadEx_metadata, {
    rea.values$exampleData    <- TRUE
    local.rea.values$dataPath <- NULL
    local.rea.values$dataPath <- paste0(system.file(package = "RFLOMICS"), 
                                        "/ExamplesFiles/ecoseed/condition.txt")
  })
  
  # as soon as the "design file" has been loaded
  # => check file format
  # => display content
  # => set ref for each factor
  # => set type of factors (bio, batch, meta)
  observeEvent(local.rea.values$dataPath, {
    
    # reset
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    local.rea.values$plots <- FALSE
    
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    
    session$userData$FlomicsMultiAssay      <- NULL
    
    # read and check design file
    design.tt <- tryCatch(expr =  readExpDesign(file = local.rea.values$dataPath),
                          error = function(e) e, warning = function(w) w)
    
    if (!is.null(design.tt$message)) {
      
      showModal(modalDialog( title = "Error message", design.tt$message))
    }
    validate({ need(expr = is.null(design.tt$message), message = design.tt$message) })
    
    ExpDesign.tbl <- design.tt
    
    local.rea.values$ExpDesign <- ExpDesign.tbl
    
    ####### Display tab of exp design  ########
    
    # display table
    output$ExpDesignTable <- renderUI({
      box(width = 12, background = "light-blue", solidHeader = TRUE, collapsible = TRUE, 
          collapsed = TRUE, title = "Overview of experimental design table", 
          DT::renderDataTable( DT::datatable(data = ExpDesign.tbl,
                                             options = list( pageLength = 5, autoWidth = TRUE, dom = 'tp' )))
      )
    })
    
    #### order and select modality for each factor
    output$selectSamplesUI <- renderUI({
      
      # condition list
      box(width = 12, background = "green", # Valid colors are: red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
          tags$b("The list and order of conditions:"),
          fluidRow(
            lapply(names(ExpDesign.tbl), function(i) {
              
              column(width = round(12/length(names(ExpDesign.tbl))),
                     selectizeInput(inputId = session$ns(paste0("selectFactors.", i)),
                                    label   = tags$span(style = "color: black;",i) ,
                                    choices = levels(as.factor(ExpDesign.tbl[[i]])),
                                    selected= levels(as.factor(ExpDesign.tbl[[i]])), multiple = TRUE,
                                    options = NULL))
            }),
            uiOutput(session$ns("GetdFactorRef")))
      )
    })
    
    ####### set ref and type of each factor ########
    output$GetdFactorRef <- renderUI({
      
      # filter samples
      # filtering per conditions
      for(factor in names(ExpDesign.tbl)) {
        
        # if select 1 modality or 0 for a Factor we exclude this factor
        if(length(input[[paste0("selectFactors.", factor)]]) > 0){
          ExpDesign.tbl <- dplyr::filter(ExpDesign.tbl, get(factor) %in% input[[paste0("selectFactors.", factor)]])
          ExpDesign.tbl[[factor]] <- factor(ExpDesign.tbl[[factor]], levels = input[[paste0("selectFactors.", factor)]])
        }
        
        if(length(input[[paste0("selectFactors.", factor)]]) <= 1){
          ExpDesign.tbl <- dplyr::select(ExpDesign.tbl, -tidyselect::all_of(factor))
        }
      }
      
      local.rea.values$ExpDesign <- ExpDesign.tbl
      
      column(width = 12,
             
             # Construct the form to set the reference factor level
             tags$b("Set the type and the reference fo each design factor:"),
             fluidRow(
               lapply(names(ExpDesign.tbl), function(i){
                 
                 column(width= round(12/length(names(ExpDesign.tbl))),
                        selectInput(inputId = session$ns(paste0("dF.RefLevel.", i)),
                                    label   = tags$span(style="color: black;",i) ,
                                    choices = unique(ExpDesign.tbl[[i]])),
                        
                        radioButtons(session$ns(paste0("dF.Type.", i)), label=NULL , choices = c("Bio","batch", "Meta"),
                                     selected = "Bio", inline = FALSE, width = 2, choiceNames = NULL, choiceValues = NULL))
               })
             )
      )
    })
  })
  
  # ---- Load omics Data ----
  # display interface for load data
  output$LoadDataUI <- renderUI({
    
    box(width = 12, title = "Load omic data", status = "warning",  height = NULL, solidHeader = TRUE,
        #h5("[warning] dataset with omics (type == none) was igored..."),
        fluidRow(
          column(2,
                 # omic type
                 selectInput(inputId = session$ns('omicType1'), label='Omic type',
                             choices = c("None"="none", "RNAseq"="RNAseq", "Proteomics"="proteomics", "Metabolomics"="metabolomics"), selected = "none")
          ),
          column(6,
                 # matrix count/abundance input
                 fileInput(inputId = session$ns("data1"), "Dataset matrix (tsv)",
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          ),
          column(3,
                 # dataset Name
                 textInput(inputId = session$ns("DataName1"), label="Dataset name",
                           value="set1")
          )
        ),
        uiOutput(   outputId = session$ns("toAddData2")),
        actionButton(inputId = session$ns("addData"),   "Add data")
    )
  })
  
  # ---- observe Event add Data ----
  # as soon as the "add data" button has been clicked
  # => a new select/file Input was display
  # =>
  addDataNum <- 1
  dataName.vec  <- c()
  observeEvent(input$addData, {
    
    # add input select for new data
    addDataNum <<- addDataNum + 1
    output[[paste("toAddData", addDataNum, sep="")]] <- renderUI({
      list(
        fluidRow(
          column(2,
                 # omic type
                 selectInput(inputId = session$ns(paste0('omicType', addDataNum)), 
                             label = 'Omic type', choices = c("None" = "none", "RNAseq" = "RNAseq", "Proteomics" = "proteomics", "Metabolomics" = "metabolomics"), selected = "none")),
          column(6,
                 # matrix count/abundance input
                 fileInput(inputId = session$ns(paste0("data", addDataNum)), 
                           "Dataset matrix (tsv)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
          column(3,
                 # dataset Name
                 textInput(inputId = session$ns(paste0("DataName", addDataNum)), 
                           label = "Dataset name", value=paste0("set", as.character(addDataNum))))
        ),
        
        uiOutput(session$ns(paste("toAddData", addDataNum + 1, sep = "")))
      )
    })
  })
  
  # ---- observeEvent loadData | loadExData ----
  observeEvent(input$loadExData, {
    
    message("# 1- Load example data ...")
    
    local.rea.values$omicsData  <- NULL
    local.rea.values$omicsNames <- NULL
    local.rea.values$omicsTypes <- NULL
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel  <- NULL
    rea.values$datasetList    <- NULL
    rea.values$datasetDiff    <- NULL
    rea.values$datasetProcess <- NULL
    rea.values$Contrasts.Sel  <- NULL
    
    session$userData$FlomicsMultiAssay <- NULL
    
    local.rea.values$plots <- FALSE
    
    inputs <- list()
    rea.values$validate.status <- 0
    
    # RNASeq
    data.mat.tt <- tryCatch( readOmicsData(file = paste0(system.file(package = "RFLOMICS"),
                                                           "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")), 
                             error = function(e) e, warning = function(w) w)
    
    local.rea.values$omicsData  <- list("RNAseq.set1" = data.mat.tt)
    local.rea.values$omicsNames <- c("RNAseq.set1")
    local.rea.values$omicsTypes <- c("RNAseq")
    
    # Metabolomic
    data.mat.tt <- tryCatch( readOmicsData(file = paste0(system.file(package = "RFLOMICS"), 
                                                           "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
                             error = function(e) e, warning = function(w) w)
    
    local.rea.values$omicsData  <- c(local.rea.values$omicsData , list("metabolomics.set2" = data.mat.tt))
    local.rea.values$omicsNames <- c(local.rea.values$omicsNames, "metabolomics.set2")
    local.rea.values$omicsTypes <- c(local.rea.values$omicsTypes, "metabolomics")
    
    # proteomics
    data.mat.tt <- tryCatch( readOmicsData(file = paste0(system.file(package = "RFLOMICS"), 
                                                           "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")), 
                             error = function(e) e, warning = function(w) w)
    
    local.rea.values$omicsData  <- c(local.rea.values$omicsData , list("proteomics.set3" = data.mat.tt))
    local.rea.values$omicsNames <- c(local.rea.values$omicsNames, "proteomics.set3")
    local.rea.values$omicsTypes <- c(local.rea.values$omicsTypes, "proteomics")
    
    
    rea.values$exampleData <- TRUE
  })
  
  # ---- Load Data button observe ---- 
  observeEvent(input$loadData, {
    
    message("# 1- Load ", input$projectName ," data ...")
    
    local.rea.values$omicsData  <- NULL
    local.rea.values$omicsNames <- NULL
    local.rea.values$omicsTypes <- NULL
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel  <- NULL
    rea.values$datasetList    <- NULL
    rea.values$datasetDiff    <- NULL
    rea.values$datasetProcess <- NULL
    rea.values$Contrasts.Sel  <- NULL
    
    session$userData$FlomicsMultiAssay <- NULL
    
    inputs <- list()
    
    omicsData  <- list()
    omicsNames <- vector()
    omicsTypes <- vector()
    
    #### list of omic data laoded from interface
    rea.values$validate.status <- 0
    dataName.vec <- c()
    for (k in 1:addDataNum){
      
      if(input[[paste0("omicType", k)]] != "none"){
        
        ### omics type ###
        omicType <- input[[paste0("omicType", k)]]
        
        ### dataset name ### 
        # => check presence of dataname
        dataName.tmp <- gsub("[[:space:]]", "", input[[paste0("DataName", k)]])
        
        if(dataName.tmp == ""){
          showModal(modalDialog( title = "Error message", "Dataset names is required: dataset ", k ))
          rea.values$validate.status <- 1
        }
        dataName <- paste0(omicType, ".", dataName.tmp)
        dataName.vec <- c(dataName.vec, dataName)
        
        # => check duplicat dataset name
        if(any(duplicated(dataName.vec)) == TRUE){
          showModal(modalDialog( title = "Error message", 
                                 "Dataset names must be unique: dataset ", 
                                 (1:addDataNum)[duplicated(dataName.vec)] ))
          rea.values$validate.status <- 1
        }
        
        #### omics dataset
        # => check omics data
        if(is.null(input[[paste0("data", k)]])){
          showModal(modalDialog( title = "Error message",
                                 "omics dataset is required: dataset ", k ))
          rea.values$validate.status <- 1
        }
        validate({ need(expr = !is.null(input[[paste0("data", k)]]), message="error") })
        
        # => read data matrix
        dataFile <- input[[paste0("data", k)]]
        data.mat.tt <- tryCatch( readOmicsData(file = dataFile$datapath),
                                 error=function(e) e, warning=function(w) w)
        
        if(!is.null(data.mat.tt$message)){
          
          showModal(modalDialog( title = "Error message", data.mat.tt$message)) 
          rea.values$validate.status <- 1
        }
        validate({ need(expr = is.null(data.mat.tt$message), message=data.mat.tt$message) })
        
        data.mat <- data.mat.tt
        
        omicsData[[dataName]]  <- data.mat
        omicsNames <- c(omicsNames, dataName)
        omicsTypes <- c(omicsTypes, omicType)
        
        validate({ need(rea.values$validate.status == 0, message = "error") })
      }
    }
    
    local.rea.values$omicsData  <- omicsData
    local.rea.values$omicsNames <- omicsNames
    local.rea.values$omicsTypes <- omicsTypes

    rea.values$exampleData <- FALSE
    
  })
  
  # ---- observe Event load Data ----
  # as soon as the "load" buttom has been clicked
  # => create ExpDesign object
  # => create flomicsMultiAssay object
  # => upsetR
  observeEvent(local.rea.values$omicsData, {
    
    ### load Design
    
    # reset objects and UI
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    rea.values$Contrasts.Sel <- NULL
    
    session$userData$FlomicsMultiAssay <- NULL
    
    local.rea.values$plots <- FALSE
    # #updateTabItems(session, "tabs", selected = "importData")
    
    ### check omicsData # no reason to check for null?
    if (is.null(local.rea.values$omicsData)) {
      showModal(modalDialog(title = "Error message", "Please load data"))
    }
    validate({ need(!is.null(local.rea.values$omicsData), message = "Please load data") })
    
    if (length(local.rea.values$omicsData) == 0) {
      showModal(modalDialog(title = "Error message", "Please load data"))
    }
    validate({ need(length(local.rea.values$omicsData) > 0, message = "Please load data") })
    
    ### check project name
    if (input$projectName == "") {
      showModal(modalDialog(title = "Error message", "project name is required"))
    }
    validate({ need(input$projectName != "", message="project name is required") })
    
    if (is.null(local.rea.values$dataPath)) {
      showModal(modalDialog(title = "Error message", "Experimental Design is required"))
    }
    validate({ need(! is.null(local.rea.values$dataPath), message = "Experimental Design is required") })
    
    ExpDesign.tbl <- local.rea.values$ExpDesign
    
    ### Get the Type and ref of the factors that the users enter in the form
    dF.Type.dFac <- vector()
    dF.List.Name <- vector()
    dF.List.ref  <- vector()
    
    rea.values$validate.status <- 0
    
    for(dFac in names(ExpDesign.tbl)){
      
      # list of type of factors (bio or batch)
      dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
      # list of level reference of factors
      dF.List.ref[dFac]  <- input[[paste0("dF.RefLevel.",dFac)]]
    }
    
    # check number of factor bio
    if(! length(stringr::str_subset(dF.Type.dFac, "Bio")) %in% 1:3){
      showModal(modalDialog(title = "Error message", "1 to 3 bio factor(s)")) }
    
    # check number of factor batch
    if(! length(stringr::str_subset(dF.Type.dFac, "batch")) %in% 1:2){
      showModal(modalDialog(title = "Error message", "at least 1 batch factor (max = 2)"))}
    
    validate({ need((length(stringr::str_subset(dF.Type.dFac, "Bio"  )) %in% 1:3) &
                      (length(stringr::str_subset(dF.Type.dFac, "batch")) %in% 1:2), message="") })
    
    #### load omics data
    #### list of omic data laoded from interface
    rea.values$validate.status <- 0
    
    FlomicsMultiAssay.try <- tryCatch( 
      createRflomicsMAE(projectName = input$projectName, 
                        omicsData   = local.rea.values$omicsData,
                        omicsNames  = local.rea.values$omicsNames,
                        omicsTypes  = local.rea.values$omicsTypes,
                        ExpDesign   = ExpDesign.tbl,
                        factorRef   = data.frame(factorName = names(dF.List.ref),
                                                 factorRef   = dF.List.ref,
                                                 factorType  = dF.Type.dFac)),
      error = function(e) e, warning = function(w) w)
    
    if(!is.null(FlomicsMultiAssay.try$message)) {
      showModal(modalDialog( title = "Error message", FlomicsMultiAssay.try$message))
      rea.values$validate.status <- 1
    }
    validate({ need(is.null(FlomicsMultiAssay.try$message), message="error") })
    
    session$userData$FlomicsMultiAssay <- FlomicsMultiAssay.try
    
    # => completeness
    #### check experimental design : experimental design must be a complete and balanced.
    print(paste0("#    => Design Completeness Check..."))
    
    # plot : ok
    local.rea.values$plots <- TRUE

    rea.values$datasetList <- session$userData$FlomicsMultiAssay@metadata$omicList
    
    rea.values$loadData <- TRUE
    
  }, ignoreInit = TRUE)
  
  ##########################################
  # (3) check data
  ##########################################
  
  # upset of all data
  output$overViewUI <- renderUI({
    
    if (local.rea.values$plots == FALSE) return()
    
    box(width = 12, status = "warning", title = "Data Overview", solidHeader = TRUE,
        
        tabsetPanel(
          tabPanel(title = "Samples",
                   h5(" Overview of the input omic data. Each color represents a distinct dataset, with their respective samples on the x-axis and the number of features on the y-axis. It illustrates the overlap between samples across datasets."),
                   renderPlot( isolate({  plotDataOverview(session$userData$FlomicsMultiAssay) }))
          ),
          tabPanel(title = "Conditions",
                   h5(" Number of Datasets per Condition. Each axis represents a distinct biological factor, and each cell value signifies the count of datasets associated with that specific condition."),
                   renderPlot( plotConditionsOverview(session$userData$FlomicsMultiAssay) )
          )
        )
    )
  })
  
  return(input)
}


