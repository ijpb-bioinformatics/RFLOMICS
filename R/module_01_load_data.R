### ============================================================================
### [01_Load_Data] shiny modules
### ----------------------------------------------------------------------------

#' @importFrom stringr str_subset
#' @importFrom dplyr filter select
#' @importFrom DT renderDataTable datatable
#' @importFrom tidyselect all_of
#' @importFrom stringr str_subset
#' @importFrom shinyBS bsTooltip popify


# ---- modCreateRflomicsObject UI ----
.modLoadOmicDataUI <- function(id) {
    ns <- NS(id)
    
    tagList(uiOutput(outputId = ns("overViewUI")))
}

# ---- modCreateRflomicsObject SERVER ----
.modLoadOmicData <-
    function(input,
             output,
             session,
             rea.values,
             local.rea.values) {
        message("[RFLOMICS] # 01- Load data: ", local.rea.values$projectName)
        
        # create Rflomics object
        FlomicsMultiAssay.try <- tryCatch(
            createRflomicsMAE(
                projectName = local.rea.values$projectName,
                omicsData   = local.rea.values$omicsData,
                omicsNames  = local.rea.values$omicsNames,
                omicsTypes  = local.rea.values$omicsTypes,
                ExpDesign   = local.rea.values$ExpDesign,
                factorRef   = data.frame(
                    factorName = names(local.rea.values$dF.List.ref),
                    factorRef  = local.rea.values$dF.List.ref,
                    factorType = local.rea.values$dF.Type.dFac
                )
            ),
            error = function(e)
                e,
            warning = function(w)
                w
        )
        
        if (!is.null(FlomicsMultiAssay.try$message)) {
            showModal(modalDialog(title = "Error message", 
                                  FlomicsMultiAssay.try$message))
            rea.values$validate.status <- 1
        }
        validate({
            need(is.null(FlomicsMultiAssay.try$message), message = "error")
        })
        
        session$userData$FlomicsMultiAssay <- FlomicsMultiAssay.try
        
        # update rea values
        rea.values$datasetList <-
          metadata(session$userData$FlomicsMultiAssay)$omicList
        rea.values$loadData    <- TRUE
        local.rea.values$plots <- TRUE
        
        # data overview
        output$overViewUI <- renderUI({
            if (local.rea.values$plots == FALSE)
                return()
            
            box(
                width = 12,
                status = "warning",
                title = "Data Overview",
                solidHeader = TRUE,
                
                tabsetPanel(
                    tabPanel(
                        title = "Samples",
                        tags$br(),
                        tags$i(
                            "Overview of the input omic data. Each color 
                            represents distinct dataset, with their respective 
                            samples on the x-axis and the number of features 
                            on the y-axis. It illustrates the samples overlap 
                            across dataset."
                        ),
                        renderPlot(isolate({
                            plotDataOverview(session$userData$FlomicsMultiAssay)
                        }))
                    ),
                    tabPanel(
                        title = "Conditions",
                        tags$br(),
                        tags$i(
                            "Number of Datasets per Condition. Each axis 
                            represents a distinct biological factor, and each 
                            cell value signifies the count of datasets 
                            associated with that specific condition."
                        ),
                        renderPlot(
                            plotConditionsOverview(session$userData$FlomicsMultiAssay)
                        )
                    )
                )
            )
        })
        
    }

# ---- modLoadOmicsData UI ----
.modLoadDataUI <- function(id) {
    ns <- NS(id)
    
    tagList(
        fluidRow(
            # load exp design
            # display tab
            box(
                width = 12,
                title = "Load Experimental Design",
                status = "warning",
                solidHeader = TRUE,
                # ExpDesign
                fluidRow(column(
                    width = 12,
                    # 1- set project name
                    column(width = 4, textInput(
                        inputId = ns("projectName"), label = "Project name"
                    )),
                    # 2- matrix count/abundance input
                    column(
                        width = 8,
                        fileInput(
                            inputId = ns("Experimental.Design.file"),
                            label   = .addBSpopify(
                                label   = "Experimental design (tsv) ",
                                title   = .generateExample("design", 
                                                           title = TRUE),
                                content = .generateExample("design")
                            ),
                            accept  = c(
                                "text/csv",
                                "text/comma-separated-values, text/plain",
                                ".csv"
                            )
                        )
                    )
                )),
                fluidRow(uiOutput(ns(
                    "displayDesignUI"
                )), )
            )
        ),
        fluidRow(uiOutput(outputId = ns("LoadDataUI"))),
        fluidRow(column(
            width = 12,
            actionButton(inputId = ns("loadData"), "Load Data", class = "butt"),
            popify(
            actionButton(
              inputId = ns("loadEcoseedData"),
              label = "Load Ecoseed Data",
              icon = icon("file-import")
            ),
              title = "Use example file for ecoseed data",
              content = paste0("<p> This will load the example data: ",
                               "the experimental design file and three differents",
                               "type of omics dataset: ",
                               "RNAseq, proteomics, and metabolomics. ",
                               "Experimental factors settings will be defined:",
                               "the factor name, type, and the reference modality."
                              
              ),
              placement = "top",
              trigger = "hover"
            )
        )),
        br(),
        
        fluidRow(.modLoadOmicDataUI(id = ns("MAE")))
    )
}

# ---- modLoadOmicsData SERVER ----
.modLoadData <- function(input, output, session, rea.values) {
    omicTypes <- c(
        "None" = "none",
        "Transcriptomics (count)" = "RNAseq",
        "Proteomics" = "proteomics",
        "Metabolomics" = "metabolomics"
    )
    
    # ---- initialization ----
    local.rea.values <- reactiveValues(
        plots        = FALSE,
        ExpDesign    = NULL,
        ExpDesignOrg = NULL,
        addDataNum   = 1
    )
    
    # ---- Load experimental design ----
    # as soon as the "design file" has been loaded
    observeEvent(input$Experimental.Design.file, {
        # reset
        local.rea.values$ExpDesignOrg <- NULL
        local.rea.values$ExpDesign    <- NULL
        local.rea.values$plots        <- FALSE
        rea.values$exampleData        <- FALSE
        
        # read and check design file
        design.tt <-
            tryCatch(
                expr = readExpDesign(file = input$Experimental.Design.file$datapath),
                error = function(e)
                    e,
                warning = function(w)
                    w
            )
        
        if (!is.null(design.tt$message)) {
            showModal(modalDialog(title = "Error message", design.tt$message))
        }
        validate({
            need(expr = is.null(design.tt$message),
                 message = design.tt$message)
        })
        
        local.rea.values$ExpDesignOrg <- design.tt
    })
    
    # ---- Add new omic data ----
    # => a new select/file Input was display
    observeEvent(input$addData, {
        # add input select for new data
        addDataNum <- local.rea.values$addDataNum
        if (input[[paste0('omicType', addDataNum)]] == "none")
            return()
        
        addDataNum <- addDataNum + 1
        
        output[[paste("toAddData", addDataNum, sep = "")]] <-
            renderUI({
                list(fluidRow(
                    column(
                        4,
                        # omic type
                        selectInput(
                            inputId = session$ns(paste0('omicType', addDataNum)),
                            label = 'Omic type',
                            choices = omicTypes,
                            selected = "none"
                        )
                    ),
                    column(
                        6,
                        # matrix count/abundance input
                        fileInput(
                            inputId = session$ns(paste0("data", addDataNum)),
                            label  = "Dataset matrix (tsv)",
                            accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"
                            )
                        )
                    ),
                    column(
                        2,
                        # dataset Name
                        textInput(
                            inputId = session$ns(paste0("DataName", addDataNum)),
                            label = "Dataset name",
                            value = paste0("set", as.character(addDataNum))
                        )
                    )
                ),
                
                uiOutput(session$ns(
                    paste("toAddData", addDataNum + 1, sep = "")
                )))
            })
        
        local.rea.values$addDataNum <- addDataNum
    })
    
    # ---- Load Data button observe ----
    observeEvent(input$loadData, {
        rea.values$loadData        <- FALSE
        rea.values$model           <- FALSE
        rea.values$analysis        <- FALSE
        rea.values$Contrasts.Sel   <- NULL
        rea.values$datasetList     <- NULL
        rea.values$datasetDiff     <- NULL
        rea.values$datasetProcess  <- NULL
        session$userData$FlomicsMultiAssay <- NULL
        
        if (isTRUE(rea.values$exampleData)) {
            local.rea.values$ExpDesignOrg <- NULL
            local.rea.values$ExpDesign    <- NULL
            updateTextInput(session, inputId = "projectName", value = "")
        }
        rea.values$exampleData      <- FALSE
        local.rea.values$plots      <- FALSE
        rea.values$validate.status  <- 0
        
        # check project name
        if (input$projectName == "") {
            showModal(modalDialog(title = "Error message", 
                                  "project name is required"))
        }
        validate({
            need(input$projectName != "", 
                 message = "project name is required")
        })
        
        # check design input
        if (is.null(local.rea.values$ExpDesign)) {
            showModal(modalDialog(title = "Error message", 
                                  "Experimental Design is required"))
        }
        validate({
            need(!is.null(local.rea.values$ExpDesign), 
                 message = "Experimental Design is required")
        })
        local.rea.values$projectName  <- input$projectName
        
        # check design
        checkDesign <- .checkDesignInput(input, local.rea.values)
        local.rea.values$dF.List.ref  <- checkDesign$dF.List.ref
        local.rea.values$dF.Type.dFac <- checkDesign$dF.Type.dFac
        
        # check omic data
        checkData <-
            .checkOmicInput(input, local.rea.values, rea.values)
        local.rea.values$omicsData    <- checkData$omicsData
        local.rea.values$omicsNames   <- checkData$omicsNames
        local.rea.values$omicsTypes   <- checkData$omicsTypes
        
        # create Rflomics object and plot data over view
        callModule(module = .modLoadOmicData,
                   id = "MAE",
                   rea.values,
                   local.rea.values)
        
    }, ignoreInit = TRUE)
    
    # ---- Load example data ----
    # load user own metadata file
    observeEvent(input$loadEcoseedData, {
        rea.values$loadData        <- FALSE
        rea.values$model           <- FALSE
        rea.values$analysis        <- FALSE
        rea.values$Contrasts.Sel   <- NULL
        rea.values$datasetList     <- NULL
        rea.values$datasetDiff     <- NULL
        rea.values$datasetProcess  <- NULL
        session$userData$FlomicsMultiAssay <- NULL
        
        # reset
        local.rea.values$plots  <- FALSE
        rea.values$exampleData  <- TRUE
        
        exampleData <- .generateEcoseedExampleData()
        
        local.rea.values$dF.List.ref  <- exampleData$dF.List.ref
        local.rea.values$dF.Type.dFac <- exampleData$dF.Type.dFac
        local.rea.values$projectName  <- exampleData$projectName
        local.rea.values$ExpDesignOrg <- exampleData$ExpDesign
        local.rea.values$ExpDesign    <- exampleData$ExpDesign
        local.rea.values$omicsNames   <- exampleData$omicsNames
        local.rea.values$omicsTypes   <- exampleData$omicsTypes
        local.rea.values$omicsData    <- exampleData$omicsData
        
        updateTextInput(session,
                        inputId = "projectName",
                        value = exampleData$projectName)
        
        callModule(module = .modLoadOmicData,
                   id = "MAE",
                   rea.values,
                   local.rea.values)
        
    }, ignoreInit = TRUE)
    
    output$displayDesignUI <- renderUI({
        if (is.null(local.rea.values$ExpDesignOrg))
            return()
        
        fluidRow(column(width = 12,
                        uiOutput(session$ns(
                            "ExpDesignTable"
                        ))),
                 column(width = 12,
                        uiOutput(session$ns(
                            "dipslayFactors"
                        ))),
        )
    })
    
    # Display tab of exp design
    output$ExpDesignTable <- renderUI({
        box(
            width = 12,
            background = "light-blue",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = TRUE,
            title = span(
                "Experimental Design Overview ",
                tags$small("(Scroll down)")
            ),
            renderDataTable(
                datatable(
                    data = local.rea.values$ExpDesignOrg,
                    options = list(
                        pageLength = 5,
                        autoWidth = TRUE,
                        dom = 'tp'
                    )
                )
            )
        )
    })
    # order and select modality for each factor
    output$dipslayFactors <- renderUI({
        ExpDesign <- local.rea.values$ExpDesignOrg
        box(
            width = 12,
            background = "green",
            .addBSpopify(
                tags$b("Select and order the modalities of each factor"),
                "We can remove modalities and consequently the associated samples. Removing all modalities of a factor results in ignoring that factor.</p>"
            ),
            fluidRow(lapply(names(ExpDesign), function(i) {
                column(
                    width = round(12 / length(names(
                        ExpDesign
                    ))),
                    selectizeInput(
                        inputId = session$ns(paste0("selectFactors.", i)),
                        label   = tags$span(style = "color: black;", i) ,
                        choices = levels(as.factor(ExpDesign[[i]])),
                        selected = levels(as.factor(ExpDesign[[i]])),
                        multiple = TRUE,
                        options = NULL
                    )
                )
            }),
            uiOutput(session$ns(
                "GetdFactorRef"
            )))
        )
        
    })
    # set ref and type of each factor
    output$GetdFactorRef <- renderUI({
        if (is.null(local.rea.values$ExpDesignOrg))
            return()
        
        ExpDesign <- local.rea.values$ExpDesignOrg
        # filter samples
        # filtering per conditions
        for (factor in names(ExpDesign)) {
            # if select 1 modality or 0 for a Factor we exclude this factor
            if (length(input[[paste0("selectFactors.", factor)]]) > 0) {
                ExpDesign    <-
                    filter(ExpDesign, get(factor) %in% input[[paste0("selectFactors.", factor)]])
                ExpDesign[[factor]] <-
                    factor(ExpDesign[[factor]], levels = input[[paste0("selectFactors.", factor)]])
            }
            
            if (length(input[[paste0("selectFactors.", factor)]]) <= 1) {
                ExpDesign <- select(ExpDesign, -all_of(factor))
            }
        }
        local.rea.values$ExpDesign <- ExpDesign
        
        radioButtons.choices <- rep("Bio", ncol(ExpDesign))
        names(radioButtons.choices) <- names(ExpDesign)
        
        if (isTRUE(rea.values$exampleData))
            radioButtons.choices <- local.rea.values$dF.Type.dFac
        
        # dispaly updated ui for selecting factors
        column(
            width = 12,
            # Construct the form to set the reference factor level
            .addBSpopify(
                tags$b("Set the reference and the type of each factor"),
                "It is mandatory to have at least one biological factor and one batch factor. Rflomics accepts 2 batch factors and supports up to 3 biological factors. If we have more than 3 biological factors, the remainder must be defined as metadata factors."
            ),
            fluidRow(lapply(names(ExpDesign), function(i) {
                column(
                    width = round(12 / length(names(
                        ExpDesign
                    ))),
                    selectInput(
                        inputId = session$ns(paste0("dF.RefLevel.", i)),
                        label   = tags$span(style = "color: black;", i) ,
                        choices = levels(ExpDesign[[i]])
                    ),
                    
                    radioButtons(
                        session$ns(paste0("dF.Type.", i)),
                        label = NULL ,
                        inline = FALSE,
                        width = 2,
                        selected = radioButtons.choices[i],
                        choiceNames = c("bio", "batch", "meta"),
                        choiceValues = c("Bio", "batch", "Meta")
                    )
                )
            }))
        )
    })
    
    # ---- load data ui ----
    # display interface for load data
    output$LoadDataUI <- renderUI({
        box(
            width = 12,
            title = "Load Omic Data",
            status = "warning",
            height = NULL,
            solidHeader = TRUE,
            fluidRow(
                column(
                    4,
                    # omic type
                    selectInput(
                        inputId = session$ns('omicType1'),
                        label = .addBSpopify(label = 'Omic type', 
                                             content = "Rflomics supports up to 3 types of omics and up to 10 datasets per omics type."),
                        choices = omicTypes,
                        selected = "none"
                    )
                ),
                column(
                    6,
                    # matrix count/abundance input
                    fileInput(
                        inputId = session$ns("data1"),
                        label   = .addBSpopify(
                            label   = 'Dataset matrix (tsv) ',
                            title   = .generateExample("matrix", title = TRUE),
                            content = .generateExample("matrix")
                        ),
                        accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv"
                        )
                    )
                ),
                column(
                    2,
                    # dataset Name
                    textInput(
                        inputId = session$ns("DataName1"),
                        label = "Dataset name",
                        value = "set1"
                    )
                )
            ),
            uiOutput(outputId = session$ns("toAddData2")),
            actionButton(inputId = session$ns("addData"),   "Add data", class = "butt")
        )
    })
    
    return(input)
}

# ---- internal function to module ----
.checkDesignInput <- function(input, local.rea.values) {
    # Get the Type and ref of the factors that the users enter in the form
    dF.Type.dFac <- vector()
    dF.List.Name <- vector()
    dF.List.ref  <- vector()
    for (dFac in names(local.rea.values$ExpDesign)) {
        # list of type of factors (bio or batch)
        dF.Type.dFac[dFac] <- input[[paste0("dF.Type.", dFac)]]
        # list of level reference of factors
        dF.List.ref[dFac]  <- input[[paste0("dF.RefLevel.", dFac)]]
    }
    
    # check number of factor bio
    if (!length(str_subset(dF.Type.dFac, "Bio")) %in% seq_len(3)) {
        showModal(modalDialog(title = "Error message", "1 to 3 bio factor(s)"))
    }
    
    # check number of factor batch
    if (!length(str_subset(dF.Type.dFac, "batch")) %in% c(1, 2)) {
        showModal(modalDialog(title = "Error message", 
                              "at least 1 batch factor (max = 2)"))
    }
    
    validate({
        need((length(str_subset(
            dF.Type.dFac, "Bio"
        )) %in% seq_len(3)) &
            (length(
                str_subset(dF.Type.dFac, "batch")
            ) %in% c(1, 2)), message = "")
    })
    
    return(
        list(
            dF.Type.dFac = dF.Type.dFac,
            dF.List.Name = dF.List.Name,
            dF.List.ref = dF.List.ref
        )
    )
}

.checkOmicInput <- function(input, local.rea.values, rea.values) {
    omicsData  <- list()
    omicsNames <- vector()
    omicsTypes <- vector()
    
    # get list of omic data laoded from interface
    dataName.vec <- c()
    for (k in seq_len(local.rea.values$addDataNum)) {
        if (input[[paste0("omicType", k)]] != "none") {
            ### omics type ###
            omicType <- input[[paste0("omicType", k)]]
            
            ### dataset name ###
            # => check presence of dataname
            dataName.tmp <-
                gsub("[[:space:]]", "", input[[paste0("DataName", k)]])
            
            if (dataName.tmp == "") {
                showModal(
                    modalDialog(title = "Error message", 
                                "Dataset names is required: dataset ", k)
                )
                rea.values$validate.status <- 1
            }
            dataName <- paste0(omicType, ".", dataName.tmp)
            dataName.vec <- c(dataName.vec, dataName)
            
            # => check duplicat dataset name
            if (any(duplicated(dataName.vec)) == TRUE) {
                showModal(
                    modalDialog(
                        title = "Error message",
                        "Dataset names must be unique: dataset ",
                        (seq_len(addDataNum))[duplicated(dataName.vec)]
                    )
                )
                rea.values$validate.status <- 1
            }
            
            #### omics dataset
            # => check omics data
            if (is.null(input[[paste0("data", k)]])) {
                showModal(
                    modalDialog(title = "Error message",
                                "omics dataset is required: dataset ", k)
                )
                rea.values$validate.status <- 1
            }
            validate({
                need(expr = !is.null(input[[paste0("data", k)]]),
                     message = "error")
            })
            
            # => read data matrix
            dataFile <- input[[paste0("data", k)]]
            data.mat.tt <-
                tryCatch(
                    readOmicsData(file = dataFile$datapath),
                    error = function(e)
                        e,
                    warning = function(w)
                        w
                )
            
            if (!is.null(data.mat.tt$message)) {
                showModal(modalDialog(title = "Error message", 
                                      data.mat.tt$message))
                rea.values$validate.status <- 1
            }
            validate({
                need(
                    expr = is.null(data.mat.tt$message),
                    message = data.mat.tt$message
                )
            })
            
            data.mat <- data.mat.tt
            
            omicsData[[dataName]]  <- data.mat
            omicsNames <- c(omicsNames, dataName)
            omicsTypes <- c(omicsTypes, omicType)
            
            validate({
                need(rea.values$validate.status == 0, message = "error")
            })
        }
    }
    
    # check omicsData # no reason to check for null?
    if (is.null(omicsData)) {
        showModal(modalDialog(title = "Error message", "Please load data"))
    }
    validate({
        need(!is.null(omicsData), message = "Please load data")
    })
    
    if (length(omicsData) == 0) {
        showModal(modalDialog(title = "Error message", "Please load data"))
    }
    validate({
        need(length(omicsData) > 0, message = "Please load data")
    })
    
    return(list(
        omicsData = omicsData,
        omicsNames = omicsNames,
        omicsTypes = omicsTypes
    ))
}

.addBSTooltip <- function(id, message = "") {
    bsTooltip(
        id = id,
        message,
        placement = "right",
        options = list(container = "body")
    )
}

.addBSpopify <- function(label = "",
                         content = "",
                         title = "",
                         color = "black",
                         placement = "right",
                         trigger = "click") {
    id <-
        paste0("id" , paste0(sample(letters, 4, replace = TRUE), collapse = ""))
    span(label,
         span(
             popify(
                 actionLink(id, icon("question-circle")),
                 title = title,
                 content = content,
                 trigger = trigger,
                 placement = placement
             ),
             style = paste0("color:", color)
         ))
    #tags$a(icon("question-circle"))
}

#' @keywords internal
#' @importFrom tidyr unite
#' @noRd
.generateExample <-
    function(what = c("design", "matrix", "annotation"),
             title = FALSE) {
        switch (what,
                "design" = {
                    table <- data.frame(
                        Sample   = c("indiv1", "indiv2", "indiv3"),
                        Genotype = c("Mutant1", "Mutant2", "Mutant1"),
                        Repeat   = c("rep1", "rep1", "rep2")
                    )
                    
                    res <-
                        c(
                            paste0("<b>", paste(names(table), collapse = "\t"), 
                                   "</b>"),
                            unite(table, "collapse", colnames(table), sep =
                                      " \t")$collapse
                        ) |>
                        paste(collapse = "<br>")
                    res   <-
                        paste0("Example:", "<pre>", res, "</pre>")
                    title.res <- "File contining experimental information and conditions for each sample, in tab-separated values (tsv) format."
                },
                "matrix" = {
                    table <- data.frame(
                        Genes   = c("gene1", "gene2", "gene3"),
                        Indiv1 = c(435, 400, 500),
                        Indiv2   = c(30, 0, 23)
                    )
                    
                    res <-
                        c(
                            paste0("<b>", paste(names(table), collapse = "\t"), "</b>"),
                            unite(table, "collapse", colnames(table), sep =
                                      "\t")$collapse
                        ) |>
                        paste(collapse = "<br>")
                    res <- paste0("Example:", "<pre>", res, "</pre>")
                    title.res <- "File containing experimental measurements (read count for transcripts, and abundance for proteins and metabolites)."
                })
        
        if (title == TRUE)
            return(title.res)
        return(res)
    }
