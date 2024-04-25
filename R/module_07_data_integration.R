##########################################
# module 07 : integration analyses
##########################################
#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinyBS popify
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch

# ---- Integration module ----
#' @title .modIntegrationAnalysisUI
#' @keywords internal
#' @noRd
.modIntegrationAnalysisUI <- function(id, method) {
    #name space for id
    ns <- NS(id)
    commontext1 <- "Select the datasets to integrate then select the features
        for each of them. In the feature selection panel, you can choose to run
         the analysis on all features: this type of option is fine if your data
         is not very high dimensional (<500). Otherwise, it is advised to 
        choose to run it only with the results of differential analysis. You 
        can refer to the overview plot below to see what you selected. Once 
        all is set, click on the Run preparation button and go the the Data 
        integration panel. "
    
    textExp <- switch(
        method,
        "mixOmics" = paste0(commontext1, 
                            "<p><b> MixOmics only runs on complete cases 
                            dataset, samples that are 
                            not present in all tables will be removed 
                            </b></p>"
        ),
        "MOFA" = paste0(commontext1, 
                        "<p><b> It is recommanded to run MOFA with at 
                        least 15 samples and a balanced number of features 
                        between tables. If you have less than 15, 
                        be cautious with the results interpretation. 
                        </b></p>"
        )
    )
    
    tagList(
        
        fluidRow(uiOutput(ns("overview"))),
        ##
        fluidRow(
            tabsetPanel(
                tabPanel(
                    "Dataset and variable selection",
                    br(),
                    box(
                        title = "Prepare data for integration",
                        width = 12,
                        status = "warning",
                        tags$style(
                            ".explain-p {
                            color: Gray;
                            text-justify: inter-word;
                            font-style: italic;
                            }"
                        ),
                        div(class = "explain-p", HTML(textExp)),
                        hr(),
                        fluidRow(uiOutput(ns("selectDataUI"))),
                        fluidRow(uiOutput(ns("selectVariablesUI")))
                    ),
                    uiOutput(ns("prepareDataUI"))
                ),
                tabPanel(
                    "Data Integration",
                    br(),
                    column(width = 3,
                           uiOutput(ns("ParamUI"))),
                    column(width = 9,
                           uiOutput(ns("ResultViewUI")))
                )
            )
        )
    )
}

#' @title .modIntegrationAnalysis
#' @importFrom DT renderDataTable datatable
#' @importFrom purrr reduce
#' @importFrom dplyr mutate across full_join
#' @keywords internal
#' @noRd
.modIntegrationAnalysis <- function(input, output, session,
                                    rea.values, method) {
    output$overview <- renderUI({
        switch(method,
               "mixOmics" = {
                   .mixOmicsText()
               },
               "MOFA" = {
                   .mofaText()
               })
    })
    
    local.rea.values <- reactiveValues(runintegration = FALSE,
                                       preparedObject = NULL)
    
    # if any preprocessing or validation of differential analysis is done
    observeEvent(rea.values$datasetProcess, {
        local.rea.values$runintegration    <- FALSE
        local.rea.values$preparedObject    <- NULL
        session$userData$FlomicsMultiAssay <- resetRflomicsMAE(
            object  = session$userData$FlomicsMultiAssay,
            multiAnalyses = c("IntegrationAnalysis")
        )
    })
    
    # select datasets to integrate
    output$selectDataUI <- .integrationSelectDataUI(session, rea.values,
                                                    input)
    
    # select method of variable reduction
    output$selectVariablesUI <-
        .integrationPrepareParamUI(session, input, rea.values)
    
    # over view of selected data after variable reduction
    output$prepareDataUI <- .integrationPrepareDataUI(session, input)
    
    # before MOFA integration
    observeEvent(input$run_prep, {
        # check: number of selected dataset min = 2
        condition <- length(input$selectData) > 1
        messCond <- "Please select at least 2 table."
        if (!condition) {
            showModal(modalDialog(title = "ERROR: ", messCond))
        }
        validate({
            need(condition, messCond)
        })
        
        local.rea.values$runintegration <- FALSE
        MAE.red <-
            session$userData$FlomicsMultiAssay[, , input$selectData]
        
        # check: sample covering
        samples.mat <- lapply(names(MAE.red), function(set) {
            tab <- data.frame("samples" = MAE.red[[set]]$samples,
                              "bin"     = rep(1, length(MAE.red[[set]]$samples)))
            names(tab) <- c("samples", set)
            return(tab)
        }) |>
            reduce(full_join, by = "samples") |>
            mutate(sum = rowSums(across(seq(2, (
                length(MAE.red) + 1
            ))), na.rm = TRUE))
        
        #---- progress bar ----#
        progress <- Progress$new()
        progress$set(message = "Prepare data for integration...", value = 0)
        on.exit(progress$close())
        #-----------------------#
        
        # check: if no common samples or no 100% overlapping
        if (length(unique(samples.mat$sum)) != 1) {
            condition <- any(unique(samples.mat$sum) != length(MAE.red))
            messCond <-
                "There is no common samples between tables. The integration
      step cannot be performed on completely unrelated tables in RFLOMICS"
            if (!condition) {
                showModal(modalDialog(title = "ERROR: ", messCond))
            }
            validate({
                need(condition, messCond)
            })
            
            if (method == "mixOmics")
                showModal(
                    modalDialog(
                        title = "WARNING: ",
                        "mixOmics is using the same samples across all tables.
             Samples not present in all tables will be removed."
                    )
                )
        } else if (unique(samples.mat$sum) != length(MAE.red)) {
            messCond <-
                "There is no common samples between tables. The integration
      step cannot be performed on completely unrelated tables in RFLOMICS"
            
            showModal(modalDialog(title = "ERROR: ", messCond))
            validate({
                need(unique(samples.mat$sum) == length(MAE.red), messCond)
            })
        }
        
        # check: MOFA : nb sample should be higher than 15
        if (method == "MOFA" &&
            length(unique(samples.mat$samples)) < 15)
            showModal(
                modalDialog(
                    title = "WARNING: ",
                    "MOFA recommends the use of more than 15 samples.
                    Be careful with the results interpretation."
                )
            )
        
        #---- progress bar ----#
        progress$inc(1 / 10, detail = "Checking everything")
        #-----------------------#
        
        # create list with variations to keep per table
        variableLists <- lapply(input$selectData, function(set) {
            switch(
                input[[paste0("selectmethode", set)]],
                "diff"  = getDEList(
                    object = session$userData$FlomicsMultiAssay[[set]],
                    contrasts = input[[paste0("selectContrast", set)]],
                    operation = input[[paste0("unionORintersect", set)]]
                ),
                "none"  = names(MAE.red[[set]])
            )
        })
        names(variableLists) <- input$selectData
        
        # check: table with nb of variables less then 5
        lowNbVarTab <- names(variableLists)[lengths(variableLists) < 5]
        
        condition <- length(lowNbVarTab) == 0
        messCond <-  paste0("number of variables is lower than 5 in
                        this(these) table(s): ",
                        lowNbVarTab)
        if (!condition) {
            showModal(modalDialog(title = "ERROR: ", messCond))
        }
        validate({
            need(condition, messCond)
        })
        
        #---- progress bar ----#
        progress$inc(5 / 10, detail = paste("Run analyses ", 50, "%", sep = ""))
        #-----------------------#
        
        # MAE to mixOmics object (list)
        message("# 8- prepare data for integration with ", method)
        local.rea.values$preparedObject <-  prepareForIntegration(
            object           = MAE.red,
            omicsNames       = input$selectData,
            rnaSeq_transfo   = input$RNAseqTransfo,
            variableLists    = variableLists,
            method           = method,
            cmd              = TRUE,
            silent           = TRUE
        )
        message("#   => Ready for integration")
        
        #---- progress bar ----#
        progress$inc(1, detail = paste("All done", 100, "%", sep = ""))
        #-----------------------#
    })
    
    ## list of parameters
    output$ParamUI <-
        .integrationMethodsParam(session, local.rea.values, method)
    
    ## observe the button run mixOmics
    observeEvent(input$run_integration, {
        #---- progress bar ----#
        progress <- Progress$new()
        progress$set(message = "Run integration", value = 0)
        on.exit(progress$close())
        progress$inc(1 / 10, detail = "in progress...")
        #----------------------#
        
        local.rea.values$runintegration <- FALSE
        
        # check: settings
        if (method == "mixOmics") {
            # Check: selection of response variables (at least one)
            condition <- !is.null(input$MO_selectedResponse)
            messCond <- "To run mixOmics, please select at least
      one response variable."
            if (!condition) {
                showModal(modalDialog(title = "Error message", messCond))
            }
            validate({
                need(condition, messCond)
            })
        }
        
        #---- progress bar ----#
        progress$inc(1 / 10, detail = paste("Checks ", 10, "%", sep = ""))
        #----------------------#
        
        # Prepare for integration run
        list_args <- list(
            object         = session$userData$FlomicsMultiAssay,
            preparedObject = local.rea.values$preparedObject,
            method         = method,
            scale_views    = input$scale_views,
            cmd = TRUE,
            silent = TRUE
        )
        
        list_args <- switch(method,
                            "mixOmics" = {
                                c(
                                    list_args,
                                    list(
                                        link_datasets = 0.5,
                                        link_response = 1,
                                        ncomp = input$MO_ncomp,
                                        sparsity = input$MO_sparsity,
                                        cases_to_try = input$MO_cases_to_try,
                                        selectedResponse = input$MO_selectedResponse
                                    )
                                )
                            },
                            "MOFA" = {
                                c(
                                    list_args,
                                    list(
                                        # maxiter = input$MOFA_maxiter,
                                        maxiter = 1000,
                                        num_factors = input$MOFA_numfactor
                                    )
                                )
                            })
        
        # Run the analysis
        message("# 8- integration Analysis with ", method)
        tryRomics <- NULL
        tryRomics <-  tryCatch(
            expr    = do.call(
                getFromNamespace("runOmicsIntegration", ns = "RFLOMICS"),
                list_args), 
            error   = function(err) err
        )
        
        #errors
        condition <- is(tryRomics, "RflomicsMAE") || !is(tryRomics, "simpleError")
        messCond <- "Something went wrong during the integration, 
        please try again with different parameters."
        if (!condition) {
            showModal(modalDialog(title = "ERROR: ", 
                                  paste0(messCond, tryRomics$message)))
        }
        validate(need(condition, paste0(messCond, tryRomics$message)))
        
        mess <- getMOFA(tryRomics, onlyResults = FALSE)$MOFA_messages
        condition <- length(mess) == 0
        if (!condition) {
            messCond <- paste(unlist(mess), collapse = "<br>")
            messCond <- gsub("prepare_mofa(.*)", "", messCond)
            showModal(modalDialog(title = "Warnings and messages from MOFA run: ", 
                                  HTML(messCond)))
        }
        
        session$userData$FlomicsMultiAssay <- tryRomics
        rm(tryRomics)
        
        listSelection <- list()
        listSelection <- lapply(input$selectData, FUN = function(set) {
            c("method"    = input[[paste0("selectmethode", set)]],
              "operation" =
                  switch(input[[paste0("selectmethode", set)]],
                         "diff" = {input[[paste0("unionORintersect", set)]]},
                         "none")
            )
        })
        names(listSelection) <- input$selectData
        metadata(session$userData$FlomicsMultiAssay)$IntegrationAnalysis[[method]]$settings$selectionMethod <-
            listSelection
        
        #---- progress bar ----#
        progress$inc(1, detail = paste("Finished ", 100, "%", sep = ""))
        #----------------------#
        
        local.rea.values$runintegration <- TRUE
    })
    
    output$ResultViewUI <- renderUI({
        switch(method,
               "mixOmics" = {
                   .modMixOmicsResultViewUI(id = session$ns("mixOmics_res"))
               },
               "MOFA" = {
                   .modMOFAResultViewUI(id = session$ns("MOFA_res"))
               })
    })
    
    switch(method,
           "mixOmics" = {
               callModule(
                   module  = .modMixOmicsResultView,
                   id = "mixOmics_res",
                   rea.values = rea.values,
                   local.rea.values = local.rea.values
               )
           },
           "MOFA" = {
               callModule(
                   module  = .modMOFAResultView,
                   id = "MOFA_res",
                   rea.values = rea.values,
                   local.rea.values = local.rea.values
               )
           })
}

# ---------- mixOmics results ----------


#' @title .modMixOmicsResultViewUI
#' @keywords internal
#' @noRd
.modMixOmicsResultViewUI <- function(id) {
    #name space for id
    ns <- NS(id)
    
    tagList(uiOutput(ns("resultsUI")))
}

#' @title .modMixOmicsResultView
#' @importFrom DT renderDataTable datatable
#' @importFrom mixOmics plotIndiv plotVar plotLoadings network cimDiablo
#' @keywords internal
#' @noRd
.modMixOmicsResultView <-
    function(input,
             output,
             session,
             rea.values,
             local.rea.values) {
        ## Output results
        output$resultsUI <- renderUI({
            if (isFALSE(local.rea.values$runintegration))
                return()
            
            settings <-
                getMixOmicsSettings(session$userData$FlomicsMultiAssay)
            
            lapply(settings$selectedResponse, function(Response) {
                Data_res <- getMixOmics(session$userData$FlomicsMultiAssay,
                                        response = Response)
                
                box(width = 14,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    status = "success",
                    title = Response,
                    
                    tabsetPanel(
                        tabPanel(
                            "Overview", 
                            .outMOOverview(Data_res, settings)),
                        # tabPanel(
                        #     "Explained Variance",
                        #     .outMOexplainedVar(session, Response)
                        # ),
                        tabPanel(
                            "Individuals",
                            .outMOIndividuals(session, input, settings,
                                              Response, Data_res)
                        ),
                        tabPanel(
                            "Features",
                            .outMOFeatures(session, input, settings,
                                           Response, Data_res)
                        ),
                        tabPanel(
                            "Loadings",
                            .outMOLoadings(session, input, settings,
                                           Response, Data_res)
                        ),
                        tabPanel(
                            "CimPlot",
                            .outMOCimPlot(session, input, settings,
                                          Response, Data_res)
                        ),
                    ) # tabsetpanel
                ) #box
            }) # lapply
        }) #renderui
        
    }

# ---------- MOFA results ----------
#' @title .modMOFAResultViewUI
#' @keywords internal
#' @noRd
.modMOFAResultViewUI <- function(id) {
    #name space for id
    ns <- NS(id)
    
    tagList(uiOutput(ns("resultsUI")))
}

#' @title .modMOFAResultView
#' @importFrom MOFA2 plot_data_overview plot_factor_cor plot_variance_explained
#' plot_weights views_names plot_factor get_dimensions plot_data_heatmap
#' @importFrom ggplot2 ggtitle
#' @importFrom DT renderDataTable datatable
#' @importFrom dplyr select
#' @importFrom grid grid.draw
#' @importFrom ggpubr ggarrange
#' @keywords internal
#' @noRd
.modMOFAResultView <- function(input,
                               output,
                               session,
                               rea.values,
                               local.rea.values) {
    output$resultsUI <- renderUI({
        if (isFALSE(local.rea.values$runintegration))
            return()
        
        resMOFA <- getMOFA(session$userData$FlomicsMultiAssay)
        
        box(
            width = 14,
            solidHeader = TRUE,
            status = "success",
            title = "MOFA results",
            
            tabsetPanel(
                # tabPanel("Overview",
                #          renderPlot(plot_data_overview(resMOFA) +
                #                         ggtitle("Data Overview"))
                # ),
                tabPanel("Factors Correlation",
                         .outMOFAFactorsCor(resMOFA)),
                tabPanel("Explained Variance",
                         .outMOFAexplainedVar(resMOFA)),
                tabPanel(
                    "Weights Plot",
                    .outMOFAWeightPlot(session, resMOFA, input)
                ),
                tabPanel(
                    "Weights table",
                    .outMOFAWeightTable(session, resMOFA, input)
                ),
                tabPanel(
                    "Factor Plots",
                    .outMOFAFactorsPlot(session, resMOFA, input)
                ),
                tabPanel("Relations",
                         .outMOFARelations(resMOFA)),
                tabPanel("Heatmap",
                         .outMOFAHeatmap(session, resMOFA, input)),
            )# tabsetpanel
        )#box
    })
    
}


# ---- UI parameter settings output ----

#' @noRd
#' @keywords internal
.integrationSelectDataUI <- function(session, rea.values,
                                     input){
    
    ns <- session$ns
    contentExp <- paste0("Only the pre-processed datasets will appear here.",
                         " If any is missing, go to the corresponding",
                         " tabset under Omics Analysis and run the",
                         " pre-processing step.")
    renderUI({
        tagList(
            column(
                width = 8,
                # select data to integrate
                checkboxGroupInput(
                    inputId  = ns("selectData"),
                    label    = h4("Select datasets for integration",
                                  popify(actionLink(
                                      inputId = ns("datSel"),
                                      label = "",
                                      icon = icon("question-circle")
                                  ), title = "Help",
                                  content = contentExp,
                                  trigger = "click", placement = "right"
                                  )),
                    choices  = rea.values$datasetProcess,
                    selected = rea.values$datasetProcess, 
                    inline = TRUE
                ), hr()),
            column(
                width = 4,
                actionButton(ns("run_prep"), 
                             label = "Run preparation", 
                             class = "btn-lg", 
                             styleclass = "primary")
            ),
        ) 
    })
}

#' @noRd
#' @keywords internal
.integrationPrepareParamUI <- function(session, input, rea.values) {
    ns <- session$ns
    
    renderUI({
        # for each dataset selected
        lapply(input$selectData, function(set) {
            
            if (set %in% rea.values$datasetDiff) {
                ValidContrasts <-
                    getValidContrasts(session$userData$FlomicsMultiAssay[[set]])
                ListNames.diff <- ValidContrasts$tag
                names(ListNames.diff) <- ValidContrasts$contrastName
                
            }else{
                ValidContrasts <- NULL
                ListNames.diff <- NULL
            }
            
            SelectTypeChoices        <- c('none', 'diff')
            names(SelectTypeChoices) <- c('none', 
                                          'from differential analysis')
            if (is.null(ValidContrasts)) {
                SelectTypeChoices <- SelectTypeChoices[c(1)]
            }
            
            box(
                title = set,
                width = 3,
                status = "primary", solidHeader = TRUE,
                # choose the method to select variable : none, diff, other?
                radioButtons(
                    inputId = ns(paste0("selectmethode", set)),
                    label   = .addBSpopify(
                        label = "Select type of variable selection",
                        title = "", 
                        content = paste0("Type of selection depends on",
                                         " the analyses performed before.",
                                         " Do not forget to validate your",
                                         " differential analysis results!"),
                        trigger = "click", placement = "right"),
                    choices = SelectTypeChoices,
                    selected = 'none',
                    inline = FALSE
                ),
                
                # if methode == diff -> dispay list of contrasts
                conditionalPanel(
                    condition = paste0("input[\'",
                                       ns(paste0(
                                           "selectmethode", set
                                       )),
                                       "\'] == \'diff\'"),
                    
                    pickerInput(
                        inputId  = ns(paste0("selectContrast", set)),
                        label    = "Select contrasts",
                        choices  = ListNames.diff,
                        options  = list(
                            `actions-box` = TRUE,
                            size = 10,
                            `selected-text-format` = "count > 3"
                        ),
                        multiple = TRUE,
                        selected = ListNames.diff
                    ),
                    
                    radioButtons(
                        inputId = ns(paste0("unionORintersect", set)),
                        label = NULL,
                        choiceNames = c("union", "intersection"),
                        choiceValues = c("union", "intersection"),
                        selected = "union",
                        inline = TRUE
                    )
                ),
                if (getOmicsTypes(session$userData$FlomicsMultiAssay[[set]]) == "RNAseq") {
                    selectInput(
                        inputId  = ns("RNAseqTransfo"),
                        label    = .addBSpopify(
                            label = "RNAseq transformation",
                            title = "", 
                            content = paste0("This parameter is fixed",
                                             " and cannot be changed.",
                                             " RNAseq data will automatically",
                                             " be transformed using",
                                             " limma::voom."),
                            trigger = "click", placement = "right"),
                        choices  = c("limma (voom)"),
                        selected = "limma (voom)"
                    )
                }
            )
        })
    })
}
#
#' @noRd
#' @keywords internal
.integrationPrepareDataUI <- function(session, input) {
    
    renderUI({
        if (length(input$selectData) == 0) return()
        if (is.null(input$selectData)) return()
        if (is.null(input[[paste0("selectmethode", input$selectData[1])]])) return()
        
        MAE2Integrate <- subRflomicsMAE(session$userData$FlomicsMultiAssay, input$selectData)
        
        for(set in input$selectData){
            
            if(input[[paste0("selectmethode", set)]] == "diff"){
                
                variable.to.keep <- getDEList(
                    object = session$userData$FlomicsMultiAssay[[set]],
                    contrasts = input[[paste0("selectContrast", set)]],
                    operation = input[[paste0("unionORintersect", set)]]
                )
                
                MAE2Integrate[[set]] <- 
                    session$userData$FlomicsMultiAssay[[set]][variable.to.keep]
                
            }
            else{
                MAE2Integrate[[set]] <- 
                    session$userData$FlomicsMultiAssay[[set]]
            }
        }
        
        textExp <- "This graph represents the dataset you will use in 
        the integration (tables and samples). 
        Gray areas represent missing samples. "
        
        box(
            title = "Overview",
            width = 12,
            status = "warning",
            solidHeader = FALSE,
            tags$style(
                ".explain-p {
                            color: Gray;
                            text-justify: inter-word;
                            font-style: italic;
                            }"
            ),
            div(class = "explain-p", 
                HTML(textExp)),
            hr(),
            renderPlot(
                plotDataOverview(MAE2Integrate,
                                 omicNames = input$selectData) +
                    theme(axis.text.y = element_text(size = 12),
                          axis.text.x = element_text(size = 10, 
                                                     angle = 45, 
                                                     hjust = 1)) 
            )
        )
    })
}

# ---- Settings methods ----
# 
#' @noRd
#' @keywords internal
.integrationMethodsParam <-
    function(session, local.rea.values, method) {
        ns <- session$ns
        
        renderUI({
            if (is.null(local.rea.values$preparedObject)) {
                return(renderText({
                    "Select the data in the previous panel to access the settings for integration."}))
            }
            
            bioFacts <- getBioFactors(session$userData$FlomicsMultiAssay)
            box(
                title = span(tagList(icon("sliders"), "  ", "Settings")),
                width = 14,
                status = "warning",
                hr(),
                span("Most of these settings are the defaults settings in the
                     method. If you don't know what to set, 
                     just run the analysis."),
                hr(),
                switch(
                    method,
                    "mixOmics" = {
                        list(
                            checkboxInput(
                                inputId = ns("scale_views"),
                                label = .addBSpopify(
                                    label = "Scale Datasets",
                                    title = "", 
                                    content = 
                                        paste0("If TRUE, each table will be transformed",
                                               " such that the global variance of the table is 1"),
                                    trigger = "click", placement = "right"),
                                value = TRUE,
                                width = NULL
                            ),
                            checkboxInput(
                                inputId = ns("MO_sparsity"),
                                label = 
                                    .addBSpopify(
                                        label = "Sparse analysis",
                                        title = "", 
                                        content = 
                                            paste0("If TRUE, sparse version of the mixOmics functions",
                                                   " will be used (block.splsda). The results of the",
                                                   " analyses will have a feature selection step",
                                                   " for each of the component and each of the tables"),
                                        trigger = "click", placement = "right"),
                                value = FALSE,
                                width = NULL
                            ),
                            numericInput(
                                inputId = ns("MO_ncomp"),
                                label = 
                                    .addBSpopify(
                                        label = "Components",
                                        title = "", 
                                        content = 
                                            paste0("Number of components to search for.",
                                                   " Equivalent to the components in the PCA."),
                                        trigger = "click", placement = "right"),
                                value = 5,
                                min = 1,
                                max = 20
                            ),
                            numericInput(
                                inputId = ns("MO_cases_to_try"),
                                label = 
                                    .addBSpopify(
                                        label = "Tuning cases",
                                        title = "", 
                                        content = paste0("The number of tuning cases will determine ",
                                                         "the number of features selection to try ",
                                                         "when performing a sparse analysis."), 
                                        # determine le nombre de cas à essayer pour le tuning.
                                        # Chaque cas ajoute un nombre de variables au précédent
                                        # Nombre de variable ajouté : nvar table/ncasestotry
                                        trigger = "click", placement = "right"),
                                value = 5,
                                min = 1,
                                max = 100
                            ),
                            checkboxGroupInput(
                                inputId  = ns("MO_selectedResponse"),
                                label    = 
                                    .addBSpopify(
                                        label = "Select response variables",
                                        title = "", 
                                        content = paste0("On which feature to perform the analysis?", 
                                                         " You can select as many response variables as you want, the analysis is performed on each of them separately.",
                                                         " The same parameters are applied for all of them."), 
                                        trigger = "click", placement = "right"),
                                choices  = bioFacts,
                                selected = bioFacts
                            )
                        )
                    },
                    "MOFA" = {
                        list(
                            selectInput(
                                ns("scale_views"),
                                label = .addBSpopify(
                                    label = "Scale views",
                                    title = "", 
                                    content = 
                                        paste0("If TRUE, each table will be transformed",
                                               " such that the global variance of the table is 1"),
                                    trigger = "click", placement = "right"),
                                choices  = c(FALSE, TRUE),
                                selected = TRUE
                            ),
                            numericInput(
                                inputId = ns("MOFA_numfactor"),
                                label = .addBSpopify(
                                    label = "Factors:",
                                    title = "", 
                                    content = 
                                        paste0("Number of components to search for.",
                                               " Roughly equivalent to the components in the PCA."),
                                    trigger = "click", placement = "right"),
                                value = 10,
                                min = 5,
                                max = 15
                            )
                            # column(
                            #     12,
                            #     numericInput(
                            #         inputId = ns("MOFA_maxiter"),
                            #         label = "Max iteration:",
                            #         value = 1000,
                            #         min = 1000,
                            #         max = 1000
                            #     )
                            # )
                        )
                    }),
                actionButton(
                    ns("run_integration"), "Run Analysis"
                )
            )
        })
    }



# ---- Plots MixOmics ----
#' @noRd
#' @keywords internal
#' @importFrom DT datatable
.outMOOverview <- function(Data_res, settings) {
    tagList(
        tags$style(
            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
        ),
        br(),
        div(class = "explain-p", 
            HTML(paste0("This overview present the number of observations and the number of features per dataset.",
                        " If a sparse analysis was performed, for each component and dataset, the best number of features is also printed."))),
        hr(),
        column(6 , renderDataTable({
            if (is.list(Data_res$X)) {
                # block.plsda or block.splsda
                df <- t(vapply(Data_res$X, dim, c(1, 1)))
                colnames(df) <- c("Ind", "Features")
                
                if (settings$sparsity) {
                    df <- cbind(df, do.call("rbind", Data_res$keepX))
                    colnames(df)[!colnames(df) %in% c("Ind", "Features")] <-
                        paste("Comp", seq_len(length(Data_res$keepX[[1]])))
                }
                
            } else {
                # plsda or splsda
                df <-  data.frame("Ind" = nrow(Data_res$X),
                                  "Features" = ncol(Data_res$X))
                if (settings$sparsity) {
                    df <- cbind(df, do.call("cbind", as.list(Data_res$keepX)))
                    colnames(df)[!colnames(df) %in% c("Ind", "Features")] <-
                        paste("Comp", seq_len(length(Data_res$keepX)))
                }
            }
            
            datatable(t(df))
            
        })))
}

#' @noRd
#' @keywords internal
.outMOexplainedVar <- function(session, Response) {
    column(12, renderPlot({
        plotMOVarExp(session$userData$FlomicsMultiAssay,
                     selectedResponse = Response)
    }))
}

#' @noRd
#' @keywords internal
#' @importFrom mixOmics plotIndiv
.outMOIndividuals <- function(session,
                              input,
                              settings,
                              Response,
                              Data_res) {
    ns <- session$ns
    
    textExplained <- paste0("These graphs represent the samples coordinates projected onto the components found by mixOmics.",
                            " If <b>Ellipses</b> is turned on, ellipses are added to the graphs according to the response variable.",
                            " When toggled, ellipses force all windows to be on the same scale. <br>", 
                            "Samples are colored according to ",
                            Response, 
                            ". It is best when the samples are grouping into response modalities.",
                            " If it is not the case, the analysis was maybe not successful."
    )
    
    renderUI({
        tagList(
            tags$style(
                ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
            ),
            br(),
            div(class = "explain-p", HTML(textExplained)),
            hr(),
            column(
                1,
                checkboxInput(
                    inputId = ns(paste0(Response, "ellipse_choice")),
                    label = "Ellipses",
                    value = FALSE,
                    width = NULL
                ),
                numericInput(
                    inputId = ns(paste0(Response, "ind_comp_choice_1")),
                    label = "Comp x:",
                    min = 1,
                    max =  settings$ncomp,
                    value = 1,
                    step = 1
                ),
                numericInput(
                    inputId = ns(paste0(Response, "ind_comp_choice_2")),
                    label = "Comp y:",
                    min = 1,
                    max =  settings$ncomp,
                    value = 2,
                    step = 1
                )
            ),
            column(11 ,
                   renderPlot(
                       plotIndiv(
                           Data_res,
                           comp = c(input[[paste0(Response, "ind_comp_choice_1")]],
                                    input[[paste0(Response, "ind_comp_choice_2")]]),
                           ellipse = input[[paste0(Response, "ellipse_choice")]],
                           legend = TRUE
                       ), height = 500)
            )
        )
    })
    
}

#' @noRd
#' @keywords internal
#' @importFrom mixOmics plotVar
.outMOFeatures <- function(session,
                           input,
                           settings,
                           Response,
                           Data_res) {
    ns <- session$ns
    moFeatText <- paste0("For each block (datatable), these graph represent",
                         ifelse(input$MO_sparsity, " the selected features", " all variables"),
                         " in a correlation circle.",
                         " The interpretation is the same as the ones from a PCA.",
                         " Important features (high loadings) are on the edge of the circle.")
    
    tagList(
        tags$style(
            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
        ),
        br(),
        div(class = "explain-p", HTML(moFeatText)),
        hr(),
        column(
            1,
            checkboxInput(
                inputId = ns("overlap"),
                label = "Overlap",
                value = FALSE,
                width = NULL
            ),
            numericInput(
                inputId = ns(paste0(Response, "var_comp_choice_1")),
                label = "Comp x:",
                min = 1,
                max =  settings$ncomp,
                value = 1,
                step = 1
            ),
            numericInput(
                inputId = ns(paste0(Response, "var_comp_choice_2")),
                label = "Comp y:",
                min = 1,
                max =  settings$ncomp,
                value = 2,
                step = 1
            )
        ),
        column(11 , renderPlot(
            plotVar(
                Data_res,
                comp = c(input[[paste0(Response, "var_comp_choice_1")]],
                         input[[paste0(Response, "var_comp_choice_2")]]),
                overlap = input$overlap,
                legend = TRUE
            )
        )))
}

#' @noRd
#' @keywords internal
#' @importFrom mixOmics plotLoadings
.outMOLoadings <-
    function(session,
             input,
             settings,
             Response,
             Data_res) {
        
        ns <- session$ns
        
        explanMess <-  paste0("For each block (datatable), these graphs represent the ",
                              # input[[paste0(Response, "Load_ndisplay")]],  
                              " first features in terms of loadings score (coefficients) for the selected component, by ascending order",
                              # input[[paste0(Response, "Load_comp_choice")]],
                              ". The higher (absolute) the score, the more important the variable is in the component construction. ",
                              " This means the highest-scores features also help the most to discriminate between categories of the response variable.")
        renderUI({
            tagList(
                tags$style(
                    ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
                ),
                br(),
                div(class = "explain-p", HTML(explanMess)),
                hr(),
                column(
                    1,
                    numericInput(
                        inputId = ns(paste0(Response, "Load_comp_choice")),
                        label = "Component:",
                        min = 1,
                        max =  settings$ncomp,
                        value = 1,
                        step = 1
                    ),
                    numericInput(
                        inputId = ns(paste0(Response, "Load_ndisplay")),
                        label = "Number of features to display:",
                        min = 1,
                        max =  ifelse(is.list(Data_res$X),
                                      max(vapply(
                                          Data_res$X, ncol, c(1)
                                      )),
                                      ncol(Data_res$X)),
                        value = 25,
                        step = 1
                    ),
                ),
                column(11 , renderPlot(
                    plotLoadings(Data_res,
                                 comp = input[[paste0(Response, "Load_comp_choice")]],
                                 ndisplay = input[[paste0(Response, "Load_ndisplay")]])
                )))
            
        })
    }

#' @noRd
#' @keywords internal
.outMOCimPlot <- function(session,
                          input,
                          settings,
                          Response,
                          Data_res) {
    ns <- session$ns
    moText <- paste0("This is a heatmap generated by the mixOmics::cimDiablo function.",
                     " All datatables are mixed together in this plot and the samples and features are clustered using euclidean distance and complete agregation method.",
                     " Features are centered but not scaled. ")
    
    if (is(Data_res, "block.splsda")) {
        tagList(
            tags$style(
                ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
            ),
            br(),
            div(class = "explain-p", HTML(moText)),
            hr(),
            column(
                1,
                numericInput(
                    inputId = ns(paste0(Response, "cimComp")),
                    label = "Comp",
                    min = 1,
                    max = settings$ncomp,
                    value = 1,
                    step = 1
                )
            ),
            column(12,
                   renderPlot(
                       cimDiablo(
                           Data_res,
                           legend.position = "bottomleft",
                           size.legend = 0.8,
                           comp = input[[paste0(Response, "cimComp")]]
                       )
                   )))
    } else{
        renderText({
            "This plot is only available for sparse multi-block discriminant analysis results."
        })
    }
}

# ---- Plots MOFA -----
#' @noRd
#' @keywords internal
.outMOFAFactorsCor <- function(resMOFA){
    mofaText <- paste0("This graph represents the correlation between calculated factors.",
                       " Factors are supposed to be orthogonal (correlation close to 0 between them).",
                       " You must check that no correlation outside of the diagonal is greater than 0.5 (absolute value).",
                       " <b>If there are factors that are strongly correlated, the results are not reliable.</b>")
    tagList(
        tags$style(
            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
        ),
        br(),
        div(class = "explain-p", HTML(mofaText)),
        hr(),      
        renderPlot(plot_factor_cor(resMOFA))
    )
}


#' @noRd
#' @keywords internal
.outMOFAexplainedVar <- function(resMOFA) {
    mofaText <- paste0("These two graphs represent the cumulative explained variance per dataset and its decomposition by factor and dataset.",
                       " This is the amount of variance within the dataset that is explained by the model (factors).",
                       " The second graph is interprated as follows: 'Factor k explains xx% of the total variance of dataset X and yy% of the total variance of dataset Y.'"
    )
    
    tagList(
        tags$style(
            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
        ),
        br(),
        div(class = "explain-p", HTML(mofaText)),
        hr(),
        column(6, renderPlot({
            g1 <- plot_variance_explained(resMOFA, plot_total = TRUE)[[2]]
            g1 + ggtitle("Total explained variance per omic data")
        })),
        column(6,
               renderPlot(
                   plot_variance_explained(resMOFA, x = "view", y = "factor") +
                       ggtitle("Explained variance by factors and omic data")
               )))
}

#' @noRd
#' @keywords internal
.outMOFAWeightPlot <- function(session, resMOFA, input) {
    ns <- session$ns
    mofaText <- paste0("Weights are an indicator of the contribution of each feature to the factors.",
                       " Features with the highest weight (absolute value) are the most related to the factor. <br>",
                       " In this graphs, the first features (by absolute weight value) are rank and presented by factor and dataset.",
                       " Weight scaling scales the factors so that the weights are all between -1 and 1." ,
                       " Weights should not be compared between datasets.")
    
    renderUI({
        tagList(
            tags$style(
                ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
            ),
            br(),
            div(class = "explain-p", HTML(mofaText)),
            hr(),
            
            column(
                3,
                sliderInput(
                    inputId = ns("WeightsPlot_Factors"),
                    label = 'Factors:',
                    min = 1,
                    max = resMOFA@dimensions$K,
                    value = c(1, 2),
                    step = 1
                )
            ),
            column(
                2,
                numericInput(
                    inputId = ns("nfeat_WeightsPlot"),
                    label = "Features:",
                    min = 1,
                    max = 500,
                    value = 10,
                    # default in MOFA function.
                )
            ),
            column(
                1,
                checkboxInput(
                    inputId = ns("scale_WeightsPlot"),
                    label = "Scale Weights",
                    value = FALSE,
                    width = NULL
                )
            ),
            
            column(12,
                   renderPlot({
                       ggplot_list <- list()
                       ggplot_list <- lapply(
                           seq(
                               min(input$WeightsPlot_Factors),
                               max(input$WeightsPlot_Factors)
                           ),
                           FUN = function(i) {
                               res_inter <- list()
                               res_inter <- lapply(
                                   views_names(resMOFA),
                                   FUN = function(vname) {
                                       res_inter[[length(res_inter) + 1]] <-
                                           plot_weights(
                                               resMOFA,
                                               view = vname,
                                               factors = i,
                                               nfeatures = input$nfeat_WeightsPlot,
                                               scale = input$scale_WeightsPlot
                                           ) +
                                           ggtitle(paste0(vname, " - Factor ", i))
                                       
                                       return(res_inter)
                                   }
                               )
                               
                               return(unlist(res_inter, recursive = FALSE))
                           }
                       )
                       ggplot_list <-
                           unlist(ggplot_list, recursive = FALSE)
                       nrowC <- max(input$WeightsPlot_Factors) -
                           min(input$WeightsPlot_Factors) + 1
                       
                       ggarrange(
                           plotlist = ggplot_list,
                           ncol = length(views_names(resMOFA)),
                           nrow = nrowC
                       )
                   }, execOnResize = TRUE)))
    })
}

#' @noRd
#' @keywords internal
.outMOFAWeightTable <- function(session, resMOFA, input) {
    ns <- session$ns
    mofaText <- paste0("Weights are an indicator of the contribution of each feature to the factors.",
                       " Features with the highest weight (absolute value) are strongly related to the factor. <br>",
                       " This table presents all the weights for the selected factors.",
                       " You can search for a specific feature."
    )
    
    verticalLayout(
        tags$style(
            ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
        ),
        br(),
        div(class = "explain-p", HTML(mofaText)),
        hr(),
        
        
        column(
            12,
            sliderInput(
                inputId = ns("Factors_select_MOFA"),
                label = 'Factors:',
                min = 1,
                max = resMOFA@dimensions$K,
                value = c(1, 2),
                step = 1
            )
        ),
        column(
            12,
            DT::renderDataTable({
                resTable <- get_weights(
                    object = resMOFA,
                    views = "all",
                    factors = seq(
                        min(input$Factors_select_MOFA),
                        max(input$Factors_select_MOFA)
                    ),
                    abs = FALSE,
                    scale = FALSE,
                    as.data.frame = TRUE
                )
                
                datatable(
                    resTable,
                    extensions = 'Buttons',
                    options = list(
                        dom = 'lfrtipB',
                        rownames = FALSE,
                        pageLength = 10,
                        buttons = c('csv', 'excel'),
                        lengthMenu = list(c(10, 25, 50, -1),
                                          c(10, 25, 50, "All"))
                    )
                )
            }, server = FALSE)
        ))
}

#' @noRd
#' @keywords internal
.outMOFAFactorsPlot <- function(session, resMOFA, input) {
    ns <- session$ns
    
    mofaText <- paste0("These graphs represent the projection of the samples/observation onto factors.",
                       " They allow you to search for a biological interpretation of the factors. <br>",
                       " You can adjust the <b>color</b> of the samples and their <b>shape</b> to help caracterize a factor.",
                       " You can also group them according to a specific response variable, and add violins or box plots to the graph.",
                       " Scaling the factors places all coordinates between -1 and 1."
    )
    
    moreChoices <- colnames(resMOFA@samples_metadata |>
                                select(!c(sample, group)))
    
    renderUI({
        tagList(
            tags$style(
                ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
            ),
            br(),
            div(class = "explain-p", HTML(mofaText)),
            hr(),
            fluidRow(
                column(2, sliderInput(
                    inputId = ns("factors_choices_MOFA"),
                    label = 'Factors:',
                    min = 1,
                    max = resMOFA@dimensions$K,
                    value = c(1, 2),
                    step = 1)),
                column(2, radioButtons(
                    inputId = ns("color_by_MOFA"),
                    label = "Color:",
                    choices = c('none', moreChoices),
                    selected = "none")),
                column(2, radioButtons(
                    inputId = ns("shape_by_MOFA"),
                    label = "Shape:",
                    choices = c('none', moreChoices),
                    selected = "none")),
                column(2, radioButtons(
                    inputId = ns("group_by_MOFA"),
                    label = "Group by:",
                    choices = c('none', moreChoices),
                    selected = "none")),
                column(
                    2,
                    verticalLayout(
                        checkboxInput(
                            inputId = ns("add_violin_MOFA"),
                            label = "Violin",
                            value = FALSE,
                            width = NULL
                        ),
                        checkboxInput(
                            inputId = ns("add_boxplot_MOFA"),
                            label = "Boxplot",
                            value = FALSE,
                            width = NULL
                        ),
                        checkboxInput(
                            inputId = ns("scale_scatter_MOFA"),
                            label = "Scale factors",
                            value = FALSE,
                            width = NULL
                        )
                    ))),
            column(12, renderPlot({
                color_by_par <- input$color_by_MOFA
                group_by_par <- input$group_by_MOFA
                shape_by_par <- input$shape_by_MOFA
                if (input$color_by_MOFA == "none")
                    color_by_par <- "group"
                if (input$group_by_MOFA == "none")
                    group_by_par <- "group"
                if (input$shape_by_MOFA == "none")
                    shape_by_par <- NULL
                dodge_par <- FALSE
                if (any(input$add_violin_MOFA, input$add_boxplot_MOFA)) {
                    dodge_par <- TRUE
                }
                
                plot_factor(
                    resMOFA,
                    factors = seq(
                        min(input$factors_choices_MOFA),
                        max(input$factors_choices_MOFA)
                    ),
                    color_by = color_by_par,
                    group_by = group_by_par,
                    shape_by = shape_by_par,
                    legend = TRUE,
                    add_violin = input$add_violin_MOFA,
                    violin_alpha = 0.25,
                    add_boxplot = input$add_boxplot_MOFA,
                    boxplot_alpha = 0.25,
                    dodge = dodge_par,
                    scale = input$scale_scatter_MOFA,
                    dot_size = 3
                )
            })))
    })
}

#' @noRd
#' @keywords internal
.outMOFAHeatmap <- function(session, resMOFA, input) {
    ns <- session$ns
    mofaText <- paste0("This heatmap represents the first features, in terms of absolute weight value, of a factor, for a specific dataset.",
                       " Samples and features are clustered using euclidean distance and complete aggregation method. <br>",
                       " You can add <b>annotations</b>: choosing a specific response variable will add a color bar at the top of the heatmap, colored according to the modality of each sample.",
                       " The <b>denoise</b> options will print the weights and factors product (the model results without the unexplained part) for this factor and this dataset. ")
    
    renderUI({
        verticalLayout(
            tags$style(
                ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
            ),
            br(),
            div(class = "explain-p", HTML(mofaText)),
            hr(),
            fluidRow(
                # buttons - choices for heatmap
                column(
                    1,
                    numericInput(
                        inputId = ns("factor_choice_heatmap_MOFA"),
                        label = "Factor:",
                        min = 1,
                        max =  resMOFA@dimensions$K,
                        value = 1,
                        step = 1
                    )
                ),
                column(1,),
                column(
                    2,
                    radioButtons(
                        inputId = ns("view_heatmap_MOFA"),
                        label = "Data:",
                        choices = views_names(resMOFA),
                        selected = views_names(resMOFA)[1]
                    )
                ),
                column(1,),
                column(
                    2,
                    numericInput(
                        inputId = ns("nfeat_choice_heatmap_MOFA"),
                        label = "Features:",
                        min = 1,
                        max = 500,
                        value = 50,
                        # default in MOFA function.
                    )
                ),
                column(
                    1,
                    checkboxInput(
                        inputId = ns("denoise_heatmap_MOFA"),
                        label = "Denoise",
                        value = FALSE,
                        width = NULL
                    )
                ),
                column(1,),
                column(
                    2,
                    radioButtons(
                        inputId = ns("annot_samples_MOFA"),
                        label = "Annotation:",
                        choices = c('none',
                                    colnames(
                                        resMOFA@samples_metadata |>
                                            select(!c(sample, group))
                                    )),
                        selected = "none"
                    )
                )
                
            ),
            fluidRow(column(12, # heatmap is too large with just renderPlot.
                            renderPlot({
                                annot_samples_values <- input$annot_samples_MOFA
                                if (input$annot_samples_MOFA == "none") {
                                    annot_samples_values <- NULL
                                }
                                
                                maxSlider <- get_dimensions(resMOFA)
                                maxSlider <-
                                    maxSlider$D[input$view_heatmap_MOFA][[1]]
                                
                                observeEvent(input$view_heatmap_MOFA, {
                                    updateSliderInput(session,
                                                      inputId = "nfeat_choice_heatmap_MOFA",
                                                      max = maxSlider)
                                })
                                
                                res_heatmap <- plot_data_heatmap(
                                    resMOFA,
                                    factor = input$factor_choice_heatmap_MOFA,
                                    view = input$view_heatmap_MOFA,
                                    features = input$nfeat_choice_heatmap_MOFA,
                                    denoise = input$denoise_heatmap_MOFA,
                                    annotation_samples = annot_samples_values
                                )
                                
                                grid.draw(res_heatmap)
                            })))
        )
    })
}


#' @noRd
#' @keywords internal
#' @importFrom stats aov
.outMOFARelations <- function(resMOFA) {
    
    mofaText <- paste0("Each of the cells is the pvalue result from an ANOVA between the factor and the response variable.",
                       " <b>No correction is applied.</b>",
                       " This table is to help you detect a biological relation between a factor and a variable that you may hae observed at the 'Factor Plots' panel. ")
    
    renderUI({
        verticalLayout(
            tags$style(
                ".explain-p {
                    color: Gray;
                    text-justify: inter-word;
                    font-style: italic;
                  }"
            ),
            br(),
            div(class = "explain-p", HTML(mofaText)),
            hr(),
            renderTable({
                factors   <- get_factors(resMOFA)
                ExpDesign <- resMOFA@samples_metadata
                ExpDesign$group  <- NULL
                ExpDesign$sample <- NULL
                
                res_aov <- lapply(
                    seq_len(ncol(ExpDesign)),
                    FUN = function(i) {
                        aov1 <- aov(factors$group1  ~ ExpDesign[, i])
                        unlist(lapply(
                            summary(aov1),
                            FUN = function(list_res) {
                                list_res[["Pr(>F)"]][[1]]
                            }
                        ))
                    }
                )
                names(res_aov) <- colnames(ExpDesign)
                
                res_res <- do.call("rbind", res_aov)
                colnames(res_res) <- gsub("Response ", "", colnames(res_res))
                
                res_res
                
            }, striped = TRUE, rownames = TRUE)
            
        )})
}

# ----- Beginning text ----

#' @noRd
#' @keywords internal
.mixOmicsText <- function() {
    box(
        title = span(tagList(
            icon('chart-line'),
            "   ",
            a("mixOmics", href = "http://mixomics.org/"),
            tags$small("")
        )),
        solidHeader = TRUE,
        status = "warning",
        width = 12,
        collapsible = TRUE,
        collapsed = TRUE,
        div(
            h4(tags$span("(Sparse) Multi-Block Discriminant Analysis", 
                         style = "color:orange")),
            p("This module uses the mixOmics package, which is dedicated to the 
        integration of multiple datasets from the same experimental design 
        (same samples)."),
        
        p("In RFLOMICS, as many data sets as you like, as long as the samples 
          overlap between tables. "),
        
        p("Analysis is performed using the block.(s)plsda function, based on 
          the sparsity parameter, which is a linear discriminant analysis 
          adapted to multiple data: axes, 
          similar to those in a PCA, are calculated according to a given
          response variable, leading to the best separation between 
          response modalities."),
        
        p("The use of multiple response variables in the same analysis 
          is not currently permitted. "),
        
        p("Sparse analysis  allows you to select interesting 
          features that discriminate between response modalities, 
          but requires tuning, i.e. running several models to 
          find the best one. This may take some time. 
      "),
        )
    )
}

#' @noRd
#' @keywords internal
.mofaText <- function() {
    box(
        title =
            span(tagList(
                icon('chart-line'),
                "   ",
                a("MOFA+", href = "https://biofam.github.io/MOFA2/"),
                tags$small("")
            )),
        solidHeader = TRUE,
        status = "warning",
        width = 12,
        collapsible = TRUE,
        collapsed = TRUE,
        div(
            h4(tags$span("Unsupervised integration analysis:", 
                         style = "color:orange")),
            p("This module uses the MOFA2 package, provinding a general 
              framework for the unsupervised integration analysis 
              of multi-omics dataset."), 
            
            p("Given a set of tables, it calculates the best axes as linear 
          combination of
         features (similar to those of a PCA). "), 
         
         p("t is based on the assumption that each dataset X can 
           be written as a linear model of the form: X = WF+e,
           where W is the weight matrix (different for each dataset), 
           F the factor coefficients (identical for each dataset) and 
           the unexplained part of the model (residuals). "), 
         
         p("Factors are constructed with each dataset, enabling the contribution 
           of each dataset to the overall model to be explored.")
        ),
    )
}
