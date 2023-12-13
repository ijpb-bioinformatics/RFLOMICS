##########################################
# module 08 : mixOmics
##########################################

IntegrationAnalysis_moduleUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      uiOutput(ns("overview"))),
    ##
    fluidRow(
      tabsetPanel(
        # ---- Tab panel Overview ----
        tabPanel("dataset and variable selection",
                 br(),
                 box(title = "Prepare data for integration", width = 12, status = "warning",
                     fluidRow(
                       uiOutput(ns("selectDataUI"))
                     ),
                     fluidRow(
                       uiOutput(ns("selectVariablesUI"))
                     )
                 ),
                 uiOutput(ns("prepareDataUI"))
        ),
        tabPanel("data integration",
                 br(),
                 column(width = 3,
                        uiOutput(ns("ParamUI"))),
                 column(width = 9, 
                        uiOutput(ns("ResultViewUI"))))
      ),
    )
  )
}

#' @importFrom mixOmics plotVar plotIndiv network circosPlot cimDiablo
#' @importFrom DT renderDataTable datatable
IntegrationAnalysis_module <- function(input, output, session, rea.values, method){
  
  output$overview <- renderUI({
    
    switch (method,
            "mixOmics" = {
              box(title = span(tagList(icon('chart-line'), "   ",a("mixOmics", href="http://mixomics.org/"), tags$small("(Scroll down for instructions)")  )),
                  solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
                  div(
                    h4(tags$span("Blocks settings:", style = "color:orange")),
                    p("You can choose which omics data you want to analyze together. It is required to have selected at least two to run an analysis."),
                    p("As of now, rflomics only allows you to run an analysis on filtered tables, taking into account differential analysis performed previously.
            You can choose which contrasts you want to take the DE genes from. You can select multiple ones (defaults select all the contrasts) 
              and choose to perform the analysis on the union or intersection of the DE lists."),
                    p("RNASeq data, given in the form of counts, will be processed using limma::voom transformation. 
              If you have indicated a batch effect when loading your data, it will be corrected in all datatables using limma::removebatcheffect before running mixOmics, for each table."),
                    p("Link between tables and response is set to 0.5 automatically."),
                    
                    h4(tags$span("Analysis settings:", style = "color:orange")),
                    p("- Scale Datasets: in each table, scale every feature to unit variance"),
                    p("- Components: number of components to be computed, default is 5"),
                    p("- Sparse Analysis: if checked, function block.splsda is run, a variable selection is performed for each component and each response variable"),
                    p("- Tuning cases: if \"sparse analysis\" is checked, to select the relevant features for each component, tuning has to be performed. Tuning cases determines the number of feature selection
              to try and decide on. A little warning here: tuning cases are applied on each table and component, and all combinaisions are tested (for example: two datatables and five tuning cases will 
              make 5*5 cases for each component to test), it can be quite long."),
                    
                    h4(tags$span("Response variables:", style = "color:orange")),
                    p("You can select as many response variables as you want, the analysis is performed on each of them separately. The same parameters are applied for all of them."),
                    
                    h4(tags$span("Outputs:", style = "color:orange")),
                    p("- Data Overview: shows how many samples and omic features are left per table after applying all filters. If you chose a sparse analysis, for each component, the number of
              selected features will be displayd."),
                    p("- Explained variance: two graphs are displayed in this section. 
              The first graph is showing the total explained variance per omic table. 
              The second one is a detailed version, showing the percentage of explained variance per feature per omic data."), 
                    p("- Individuals plot: coordinates of individuals on each component for each table. When toggled, ellipses force all windows to be on the same scale."),
                    p("- Features plot: similar to a pca correlation plot, coordinates of the (selected) features on the correlation circles for all datasets. "),
                    p("- Loadings: coefficients for the (selected) features for each dataset, by ascending order. Default is showing the first 25 entities for each block, on the first component."),
                    p("- Networks: similarity networks computed on selected components factors. It only shows the between table correlation, not the intra-table correlaction!"),
                    p("- CircosPlot: only available for sparse analyses. 
                Similarity networks computed on selected components factors. It only shows the between table correlation, not the intra-table correlaction!"),
                    p("- CimPlot: only available for sparse analyses. 
                For each component, shows the heatmap of selected features for all tables."),
                    
                  ))
            },
            "MOFA" = {
              box(title = span(tagList(icon('chart-line'), "   ", a("MOFA+", href = "https://biofam.github.io/MOFA2/"),
                                       tags$small("(Scroll down for instructions)")  )),
                  solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
                  div(      
                    h4(tags$span("Parameters set up:", style = "color:orange")),
                    p("You can choose which omics data you want to analyze together. It is required to have selected at least two to run an analysis."),
                    p("As of now, rflomics only allows you to run an analysis on filtered tables, taking into account differential analysis performed previously.
            You can choose which contrast you want to take the DE genes from. You can select multiple ones (defaults select all the contrasts) 
              and choose to perform the analysis on the union or intersection of the DE lists."),
                    p("RNASeq data, given in the form of counts, will be processed using limma::voom transformation. 
              If you have indicated a batch effect when loading your data, it will be corrected in all datatables using limma::removebatcheffect before running MOFA."),
                    
                    h4(tags$span("Outputs:", style = "color:orange")),
                    p("- Data Overview: shows how many samples and omic features are left per table after applying all filters. These are the data on which MOFA is performed."),
                    p("- Factors correlation: shows the correlation between factors. It is essential that the factors are decorrelated between each other. If there is a correlation >0.5 in absolute value, 
              the results of the analysis cannot be trusted."),
                    p("- Explained variance: two graphs are displayed in this section. 
            The first graph is showing the total explained variance per omic table. 
              The second one is a detailed version, showing the percentage of explained variance per feature per omic data."), 
                    # The percentage of explained variance cannot be compared between views."),
                    p("- Weights plot: displays the features weights, for selected factors, for each omic. 
              You can scale the weights by checking the \"scale weight\" box. This will ensure weights are comprised between -1 and 1."),
                    p("- Weights table: table of weights per selected factors."),
                    p("- Heatmap: shows the data for the selected number of features, per factor and per omic. Original data is displayed by default. 
              You can use the denoise option to show the reconstruction using factors and weights computed by MOFA2. 
              The annotation parameter is used to place a color annotation on the sample, based on a metadata feature." ),
                  ))
            }
    )
  })
  
  local.rea.values <- reactiveValues(runintegration = FALSE, preparedObject = NULL)
  
  observeEvent(rea.values$datasetProcess, {
    local.rea.values$runintegration    <- FALSE
    local.rea.values$preparedObject    <- NULL
    session$userData$FlomicsMultiAssay <- resetFlomicsMultiAssay(object  = session$userData$FlomicsMultiAssay, 
                                                                 results = c("IntegrationAnalysis"))
  })
  
  # select datasets to integrate
  output$selectDataUI <- renderUI({
    
    tagList(
      column(width = 8,
             # select data to integrate
             checkboxGroupInput(inputId  = session$ns("selectData"),
                                label    = "Select datasets to integrate together",
                                choices  = rea.values$datasetProcess,
                                selected = rea.values$datasetProcess)
      ),
      column(width = 4,
             actionButton(session$ns("run_prep"), label = "Let's go ;)")
      )
    )
  })
  
  # select methode of variable reduction
  output$selectVariablesUI <- renderUI({
    
    # for each dataset select 
    lapply(input$selectData, function(set){
      
      ValidContrasts        <- getValidContrasts(session$userData$FlomicsMultiAssay[[set]])
      ListNames.diff        <- ValidContrasts$tag
      names(ListNames.diff) <- ValidContrasts$contrastName
      SelectTypeChoices        <- c('none', 'diff')
      names(SelectTypeChoices) <- c('none', 'from diff analysis')
      
      if(is.null(ValidContrasts)) SelectTypeChoices <- SelectTypeChoices[c(1)]
      
      box(title = set, width = 4, background = "green",
          # choose the method to select variable : none, diff, other?
          radioButtons(inputId = session$ns(paste0("selectmethode", set)), 
                       label   = 'Select type of variable selection', 
                       choices = SelectTypeChoices, 
                       selected = 'none', inline = FALSE),
          
          # if methode == diff -> dispay list of contrasts
          conditionalPanel(
            
            condition = paste0("input[\'", session$ns(paste0("selectmethode", set)), "\'] == \'diff\'"),
            
            pickerInput(
              inputId  = session$ns(paste0("selectContrast", set)),
              label    = "Select contrasts",
              choices  = ListNames.diff,
              options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
              multiple = TRUE,
              selected = ListNames.diff),
            
            radioButtons(inputId = session$ns(paste0("unionORintersect", set)), label = NULL, 
                         choiceNames = c("union", "intersection"), 
                         choiceValues = c(TRUE, FALSE),
                         selected = TRUE, inline = TRUE)
          ),
          if (getOmicsTypes(session$userData$FlomicsMultiAssay[[set]]) == "RNAseq"){
            selectInput(inputId  = session$ns("RNAseqTransfo"),
                        label    = "RNAseq transformation",
                        choices  = c("limma (voom)"),
                        selected = "limma (voom)")
          }
      )
    })
  })
  
  # over view of selected data after variable reduction
  output$prepareDataUI <- renderUI({
    
    if(length(input$selectData) == 0) return()
    
    
    list.SE <- lapply(input$selectData, function(set){
      switch (input[[paste0("selectmethode", set)]],
              "diff" = {
                variable.to.keep <- getDE(object = session$userData$FlomicsMultiAssay[[set]], 
                                          contrast = input[[paste0("selectContrast", set)]], 
                                          union = input[[paste0("unionORintersect", set)]])$DEF
                session$userData$FlomicsMultiAssay[[set]][variable.to.keep]
              },
              {
                session$userData$FlomicsMultiAssay[[set]]
              }
      )
    })
    names(list.SE) <- input$selectData
    
    MAE2Integrate <- MultiAssayExperiment(experiments = list.SE,
                                          colData     = colData(session$userData$FlomicsMultiAssay),
                                          sampleMap   = sampleMap(session$userData$FlomicsMultiAssay),
                                          metadata    = metadata(session$userData$FlomicsMultiAssay))
    
    box(title = "", width = 12, status = "warning",
        renderPlot(Datasets_overview_plot(MAE2Integrate, dataset.list = input$selectData))
    )
  })
  
  # befor go to MOFA integration
  observeEvent(input$run_prep, {
    
    # check number of selected dataset min = 2
    if (length(input$selectData) < 1) {
      showModal(modalDialog(title = "ERROR : ", "Please select at least 1 table."))
    }
    validate({ 
      need(length(input$selectData) >= 1, "Please select at least 1 table.") 
    })
    
    local.rea.values$runintegration <- FALSE
    session$userData$FlomicsMultiAssay@metadata[[method]] <- NULL
    MAE.red <- session$userData$FlomicsMultiAssay[,, input$selectData]
    
    # check sample covering
    samples.mat <- lapply(names(MAE.red), function(set){
      tab <- data.frame("samples" = MAE.red[[set]]$samples,
                        "bin"     = rep(1, length(MAE.red[[set]]$samples)))
      names(tab) <- c("samples", set)
      return(tab)
    }) |> 
      purrr::reduce(dplyr::full_join, by = "samples") |> 
      dplyr::mutate(sum = rowSums(dplyr::across(2:(length(MAE.red)+1)), na.rm = TRUE))
    
    # if no common samples or no 100% overlapping
    if(length(unique(samples.mat$sum)) != 1){
      
      if(!any(unique(samples.mat$sum) != length(MAE.red))){
        showModal(modalDialog(title = "ERROR : ", "No commun samples between tables."))
      }
      validate({
        need(any(unique(samples.mat$sum) != length(MAE.red)), "") 
      })
      
      if(method == "mixOmics")
        showModal(modalDialog(title = "WARNING : ", "mixOmics recommends using the same samples across all tables. These samples will be removed from the other tables."))
    }
    else if(unique(samples.mat$sum) != length(MAE.red)){
      
      showModal(modalDialog(title = "ERROR : ", "No commun samples between tables."))
      
      validate({
        need(unique(samples.mat$sum) == length(MAE.red), "") 
      })
    }
    
    # MOFA : nb sample should be higher then 15
    if(method == "MOFA" && length(unique(samples.mat$samples)) < 15) 
      showModal(modalDialog(title = "WARNING : ", "MOFA recommand the use of more of 15 samples."))
    
    # creat list with variations to keep per table
    variableLists <- lapply(input$selectData, function(set){
      switch (input[[paste0("selectmethode", set)]],
              "diff"  = getDE(object = session$userData$FlomicsMultiAssay[[set]], 
                              contrast = input[[paste0("selectContrast", set)]], 
                              union = input[[paste0("unionORintersect", set)]])$DEF,
              "none"  = names(MAE.red[[set]])
      )
    })
    names(variableLists) <- input$selectData
    
    # table with nb of variables less then 5
    lowNbVarTab <- names(variableLists)[lengths(variableLists) < 5]
    
    if(length(lowNbVarTab) != 0){
      showModal(modalDialog(title = "EROOR : ", paste0("number of variations is lower then 5 in this(these) table(s): ", lowNbVarTab)))
    }
    validate({ 
      need(length(lowNbVarTab) == 0, "") 
    })
    
    # MAE to mixOmics object (list)
    print(paste0("# 8- prepare data for integration with ", method))
    local.rea.values$preparedObject <-  prepareForIntegration(
      object           = MAE.red,
      omicsToIntegrate = input$selectData,
      rnaSeq_transfo   = input$RNAseqTransfo,
      variableLists    = variableLists,
      method           = method,
      cmd              = TRUE, 
      silent           = TRUE)
    
  })
  
  ## list of parameters  
  output$ParamUI <- renderUI({
    
    # # Reinitialize if needed when validating another datatable in differential analysis. 
    # observeEvent(rea.values$datasetDiff, {
    #   local.rea.values$runMixOmics   <- FALSE
    #   if (!is.null(getMixOmics(session$userData$FlomicsMultiAssay))) {
    #     session$userData$FlomicsMultiAssay <- setMixOmics(session$userData$FlomicsMultiAssay, NULL)
    #   } 
    # })
    
    #listOfContrast <- getSelectedContrasts(session$userData$FlomicsMultiAssay)$contrastName
    
    if(is.null(local.rea.values$preparedObject)) return()
    
    box(title = span(tagList(icon("sliders"), "  ", "Settings")), width = 12, status = "warning",
        switch (method,
                "mixOmics" = {
                  
                  # set param in interface
                  list(
                    
                    column(12, checkboxInput(inputId = session$ns("scale_views"), label = "Scale Datasets",
                                             value = TRUE, width = NULL)),
                    column(12, checkboxInput(inputId = session$ns("MO_sparsity"), label = "Sparse analysis", 
                                             value = FALSE, width = NULL)),
                    column(12, numericInput(inputId = session$ns("MO_ncomp"), label = "Components", 
                                            value = 5, min = 1, max = 20)),
                    column(12, numericInput(inputId = session$ns("MO_cases_to_try"), label = "Tuning cases",
                                            value = 5, min = 1, max = 100)),
                    column(12,
                           checkboxGroupInput(
                             inputId  = session$ns("MO_selectedResponse"),
                             label    = "Select response variables",
                             choices  = bioFactors(session$userData$FlomicsMultiAssay),
                             selected = bioFactors(session$userData$FlomicsMultiAssay)))
                    
                  )
                },
                "MOFA" = {
                  
                  list(
                    column(12,
                           selectInput(session$ns("scale_views"),
                                       label    = "Scale views",
                                       choices  = c(FALSE, TRUE),
                                       selected = TRUE)),
                    column(12,
                           numericInput(inputId = session$ns("MOFA_numfactor"), label = "Num factors:", 
                                        value = 10, min = 5, max = 15)),
                    
                    column(12,
                           numericInput(inputId = session$ns("MOFA_maxiter"), label = "Max iteration:", 
                                        value = 1000, min = 1000, max = 1000))
                  )
                }
        ),
        column(12, actionButton(session$ns("run_integration"), "Run Analysis"))
    )
  })
  
  ## observe the button run mixOmics
  observeEvent(input$run_integration, {
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run integration", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    local.rea.values$runintegration <- FALSE
    session$userData$FlomicsMultiAssay@metadata[[method]] <- NULL
    
    # check setting
    if(method == "mixOmics"){
      
      # Check selection of response variables (at least one)
      if (is.null(input$MO_selectedResponse)) {
        showModal(modalDialog(title = "Error message", "To run mixOmics, please select at least one response variable"))
      }
      validate({ 
        need(!is.null(input$MO_selectedResponse), "To run mixOmics, please select at least one response variable") 
      })
    }
    # if(method == "MOFA"){
    #   
    #   
    # }
    
    
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Checks ", 10, "%", sep = ""))
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
    
    list_args <- switch (method,
                         "mixOmics" = {c(list_args, 
                                         list(
                                           link_datasets = 0.5,
                                           link_response = 1,
                                           ncomp = input$MO_ncomp,   
                                           sparsity = input$MO_sparsity,
                                           cases_to_try = input$MO_cases_to_try,
                                           selectedResponse = input$MO_selectedResponse))},
                         "MOFA" = {c(list_args,
                                     list(
                                       maxiter = input$MOFA_maxiter,
                                       num_factors = input$MOFA_numfactor
                                     ))}
    )
    
    # Run the analysis
    print(paste0("# 8- integration Analysis with ", method))
    session$userData$FlomicsMultiAssay <- do.call(getFromNamespace("runIntegration", ns = "RFLOMICS"), list_args)
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Finished ", 100,"%", sep = ""))
    #----------------------#
    
    local.rea.values$runintegration <- TRUE
  })
  
  output$ResultViewUI <- renderUI({
    
    
    switch (method,
            "mixOmics" = {
              mixOmics_Result_View_UI(id = session$ns("mixOmics_res"))
            },
            "MOFA" = {
              MOFA_Result_View_UI(id = session$ns("MOFA_res"))
            }
    )
  })
  
  # input$sparsity.  input$ncomp.  input$overlap. input$selectData.
  
  
  switch (method,
          "mixOmics" = {
            shiny::callModule(module  = mixOmics_Result_View, id = "mixOmics_res", rea.values = rea.values, local.rea.values = local.rea.values)
          },
          "MOFA" = {
            shiny::callModule(module  = MOFA_Result_View, id = "MOFA_res", rea.values = rea.values, local.rea.values = local.rea.values)
          }
  )
}

# ---------- mixOmics results ----------
mixOmics_Result_View_UI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    fluidRow(
      uiOutput(ns("resultsUI")))
  )
}

mixOmics_Result_View <- function(input, output, session, rea.values, local.rea.values){
  
  ## Output results
  output$resultsUI <- renderUI({
    
    if (isFALSE(local.rea.values$runintegration)) return()
    
    setting <- getMixOmicsSetting(session$userData$FlomicsMultiAssay)
    
    #lapply(names(getMixOmics(session$userData$FlomicsMultiAssay)), function(listname) { 
    lapply(setting$selectedResponse, function(Response) { 
      
      Data_res <- getMixOmics(session$userData$FlomicsMultiAssay, response = Response)
      
      fluidRow(
        
        box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success", title = Response,
            
            tabsetPanel(
              # ---- Tab panel Overview ----
              tabPanel("Overview",
                       column(6 , DT::renderDataTable({
                         
                         df <- t(sapply(Data_res$X, dim))
                         colnames(df) <- c("Ind", "Features")
                         
                         if (setting$sparsity) {
                           df <- cbind(df, do.call("rbind", Data_res$keepX))
                           colnames(df)[!colnames(df) %in% c("Ind", "Features")] <- paste("Comp", 1:length(Data_res$keepX[[1]]))
                         }
                         
                         t(df) %>% DT::datatable()
                         
                       })),
                       
              ),
              # ---- Tab panel Explained Variance ----
              tabPanel("Explained Variance",
                       column(12, renderPlot({
                         
                         plot_MO_varExp(session$userData$FlomicsMultiAssay, 
                                        selectedResponse = Response)
                       })),
                       
              ),
              # ---- Tab panel Individuals ----
              tabPanel("Individuals",
                       column(1,
                              checkboxInput(inputId = session$ns(paste0(Response, "ellipse_choice")), label = "Ellipses", value = FALSE, width = NULL),
                              numericInput(inputId = session$ns(paste0(Response, "ind_comp_choice_1")),
                                           label = "Comp x:",
                                           min = 1,
                                           max =  setting$ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(Response, "ind_comp_choice_2")),
                                           label = "Comp y:",
                                           min = 1,
                                           max =  setting$ncomp,
                                           value = 2, step = 1)
                       ),
                       column(11 , renderPlot(suppressWarnings(
                         mixOmics::plotIndiv(Data_res, 
                                             comp = c(input[[paste0(Response, "ind_comp_choice_1")]], input[[paste0(Response, "ind_comp_choice_2")]]),
                                             ellipse = input[[paste0(Response, "ellipse_choice")]],
                                             legend = TRUE))))
              ),
              # ---- Tab panel Features ----
              tabPanel("Features",
                       column(1,
                              checkboxInput(inputId = session$ns("overlap"), label = "Overlap", value = FALSE, width = NULL),
                              numericInput(inputId = session$ns(paste0(Response, "var_comp_choice_1")),
                                           label = "Comp x:",
                                           min = 1,
                                           max =  setting$ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(Response, "var_comp_choice_2")),
                                           label = "Comp y:",
                                           min = 1,
                                           max =  setting$ncomp,
                                           value = 2, step = 1)
                       ),
                       column(11 , renderPlot(mixOmics::plotVar(Data_res, 
                                                                comp = c(input[[paste0(Response, "var_comp_choice_1")]], input[[paste0(Response, "var_comp_choice_2")]]),
                                                                overlap = input$overlap,
                                                                legend = TRUE)))
              ),     
              # ---- Tab panel Loadings ----
              tabPanel("Loadings",
                       column(1,
                              numericInput(inputId = session$ns(paste0(Response, "Load_comp_choice")),
                                           label = "Component:",
                                           min = 1,
                                           max =  setting$ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(Response, "Load_ndisplay")),
                                           label = "Number of features to display:",
                                           min = 1,
                                           max =  max(sapply(Data_res$X, ncol)),
                                           value = 25, step = 1),
                       ),
                       column(11 , renderPlot(mixOmics::plotLoadings(Data_res, 
                                                                     comp = input[[paste0(Response, "Load_comp_choice")]],
                                                                     ndisplay = input[[paste0(Response, "Load_ndisplay")]])))
              ),   
              tabPanel("Networks",
                       column(1, numericInput(inputId = session$ns(paste0(Response, "Network_cutoff")),
                                              label = "Cutoff:",
                                              min = 0,
                                              max =  1,
                                              value = 0.9, step = 0.05)),
                       column(11 , 
                              renderUI({
                                outN <- .doNotPlot(mixOmics::network(mat = Data_res, 
                                                                     blocks = 1:length(setting$selectData),
                                                                     cutoff = input[[paste0(Response, "Network_cutoff")]], 
                                                                     shape.node = rep("rectangle", length(setting$selectData)))) 
                                
                                if (is(outN, "simpleError")){
                                  renderText({outN$message})
                                } else {
                                  renderPlot(
                                    .doNotSpeak(
                                      mixOmics::network(mat = Data_res, 
                                                        blocks = 1:length(setting$selectData),
                                                        cutoff = input[[paste0(Response, "Network_cutoff")]], 
                                                        shape.node = rep("rectangle", length(setting$selectData)))
                                    ))
                                }
                                
                              })
                              
                       )
              ), 
              # ---- Tab  Panel CircosPlot & cimPlot ----
              tabPanel("CircosPlot",
                       if (is(Data_res, "block.splsda") || is(Data_res, "block.plsda")) {
                         tagList(
                           fluidRow(
                             column(1, numericInput(inputId = session$ns(paste0(Response, "Circos_cutoff")),
                                                    label = "Cutoff:",
                                                    min = 0,
                                                    max =  1,
                                                    value = 0.9, step = 0.05))),
                           fluidRow(
                             column(11, 
                                    
                                    renderUI({
                                      out <- .doNotPlot(mixOmics::circosPlot(Data_res, cutoff = input[[paste0(Response, "Circos_cutoff")]]))
                                      
                                      if (is(out, "simpleWarning")) {
                                        renderText({out$message})
                                      } else {
                                        renderPlot(
                                          .doNotSpeak(
                                            mixOmics::circosPlot(Data_res, cutoff = input[[paste0(Response, "Circos_cutoff")]])
                                          )
                                        )
                                      }
                                      
                                    })
                                    
                             ))) # column # fluidrow
                       } else {
                         renderText({print("This plot is only available for multi-block discriminant analysis results.")})
                       }
              ),
              tabPanel("CimPlot",
                       if (is(Data_res, "block.splsda")) {
                         fluidRow(
                           column(1,
                                  numericInput(inputId = session$ns(paste0(Response, "cimComp")),
                                               label = "Comp",
                                               min = 1,
                                               max = setting$ncomp,
                                               value = 1, step = 1)),
                           column(12, 
                                  
                                  renderPlot(mixOmics::cimDiablo(Data_res, 
                                                                 legend.position = "bottomleft",
                                                                 size.legend = 0.8,
                                                                 comp = input[[paste0(Response, "cimComp")]])))
                         )
                       }else{
                         renderText({print("This plot is only available for sparse multi-block discriminant analysis results.")})
                       }
                       
              ),
            ) # tabsetpanel
        ) #box
      ) # fluidrow
    }) # lapply
  }) #renderui
}

# ---------- MOFA results ----------
MOFA_Result_View_UI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    fluidRow(
      uiOutput(ns("resultsUI")))
  )
}

MOFA_Result_View <- function(input, output, session, rea.values, local.rea.values){
  
  output$resultsUI <- renderUI({
    
    if (isFALSE(local.rea.values$runintegration)) return()
    
    resMOFA <- getMOFA(session$userData$FlomicsMultiAssay)
    
    colorPal_choices <- c("Greens", "Purples", "Oranges", "Reds", "Greys", "Blues")
    
    box(width = 14, solidHeader = TRUE, status = "success",
        title = "MOFA results",
        
        tabsetPanel(
          
          ###
          # ---- Tab panel Factors Overview ----
          tabPanel("Overview", 
                   column(12,
                          renderPlot(plot_data_overview(resMOFA) + ggtitle("Data Overview"))),
          ),
          ### 
          # ---- Tab panel Factors Correlation ----
          tabPanel("Factors Correlation", 
                   renderPlot(plot_factor_cor(resMOFA))
          ),
          ### 
          # ---- Tab panel Explained Variance ----
          tabPanel("Explained Variance", 
                   fluidRow(
                     column(6, renderPlot({
                       g1 <- plot_variance_explained(resMOFA, plot_total = TRUE)[[2]]
                       g1 + ggtitle("Total explained variance per omic data")
                     })),
                     column(6, renderPlot(plot_variance_explained(resMOFA, x = "view", y = "factor") +
                                            ggtitle("Explained variance by factors and omic data")))
                   )),
          ### 
          # ---- Tab panel Weights Plot ----
          tabPanel("Weights Plot", 
                   fluidRow(
                     column(3, sliderInput(inputId = session$ns("WeightsPlot_Factors_select"),
                                           label = 'Factors:',
                                           min = 1, 
                                           max = resMOFA@dimensions$K, 
                                           value = 1:2, step = 1)
                     ),
                     column(2, numericInput(inputId = session$ns("nfeat_choice_WeightsPlot"),
                                            label = "Features:",
                                            min = 1,
                                            max = 500,
                                            value = 10, # default in MOFA function.
                     )),
                     column(1,
                            checkboxInput(inputId = session$ns("scale_choice_WeightsPlot"), label = "Scale Weights", value = FALSE, width = NULL)
                     ),
                   ),
                   fluidRow(
                     column(12,
                            renderPlot({
                              
                              ggplot_list <- list()
                              ggplot_list <- lapply(min(input$WeightsPlot_Factors_select):max(input$WeightsPlot_Factors_select), FUN = function(i) {
                                
                                res_inter <- list()
                                res_inter <- lapply(views_names(resMOFA), FUN = function(vname) {
                                  
                                  res_inter[[length(res_inter) + 1]] <- plot_weights(resMOFA,
                                                                                     view = vname,
                                                                                     factors = i,
                                                                                     nfeatures = input$nfeat_choice_WeightsPlot,
                                                                                     scale = input$scale_choice_WeightsPlot) 
                                  res_inter[[length(res_inter)]] <- res_inter[[length(res_inter)]] + 
                                    ggtitle(paste0(vname, " - Factor ", i))
                                  
                                  return(res_inter)
                                })
                                
                                return(unlist(res_inter, recursive = FALSE))
                              })
                              ggplot_list <- unlist(ggplot_list, recursive = FALSE)
                              
                              ggpubr::ggarrange(plotlist = ggplot_list,
                                                ncol = length(MOFA2::views_names(resMOFA)),
                                                nrow = (max(input$WeightsPlot_Factors_select) - min(input$WeightsPlot_Factors_select) + 1))
                              
                            }, execOnResize = TRUE)
                            
                     ))
          ),
          
          ### 
          # ---- Tab panel Weights table ----
          tabPanel("Weights table",
                   
                   fluidRow(
                     column(11, sliderInput(inputId = session$ns("Factors_select_MOFA"),
                                            label = 'Factors:',
                                            min = 1, 
                                            max = resMOFA@dimensions$K,  
                                            value = 1:2, step = 1)) 
                   ),
                   fluidRow(
                     column(11, 
                            DT::renderDataTable({
                              resTable <- get_weights(resMOFA, views = "all", factors = min(input$Factors_select_MOFA):max(input$Factors_select_MOFA),
                                                      abs = FALSE, scale = FALSE, as.data.frame = TRUE)
                              
                              resTable %>% DT::datatable(extensions = 'Buttons',
                                                         options = list(dom = 'lfrtipB',
                                                                        rownames = FALSE,
                                                                        pageLength = 10,
                                                                        buttons = c('csv', 'excel'),
                                                                        lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All"))))
                            }, server = FALSE))
                   )),
          # ---- Tab panel Factor Plots ----
          tabPanel("Factor Plots",
                   
                   fluidRow(
                     column(2, sliderInput(inputId = session$ns("factors_choices_MOFA"),
                                           label = 'Factors:',
                                           min = 1, 
                                           max = resMOFA@dimensions$K, # pas forcement l'input, MOFA peut decider d'en enlever. 
                                           value = 1:2, step = 1)),
                     column(2, radioButtons(inputId = session$ns("color_by_MOFA"),
                                            label = "Color:",
                                            choices = c('none', 
                                                        colnames(resMOFA@samples_metadata %>% dplyr::select(!c(sample, group)))),
                                            selected = "none")
                     ),                           
                     column(2, radioButtons(inputId = session$ns("shape_by_MOFA"),
                                            label = "Shape:",
                                            choices = c('none', 
                                                        colnames(resMOFA@samples_metadata %>% dplyr::select(!c(sample, group)))),
                                            selected = "none")
                     ),                                  
                     column(2, radioButtons(inputId = session$ns("group_by_MOFA"),
                                            label = "Group by:",
                                            choices = c('none', 
                                                        colnames(resMOFA@samples_metadata %>% dplyr::select(!c(sample, group)))),
                                            selected = "none")
                     ),                           
                     column(1,
                            fluidRow(checkboxInput(inputId = session$ns("add_violin_MOFA"), label = "Violin", value = FALSE, width = NULL)),
                            fluidRow(checkboxInput(inputId = session$ns("add_boxplot_MOFA"), label = "Boxplot", value = FALSE, width = NULL)),
                            fluidRow(checkboxInput(inputId = session$ns("scale_scatter_MOFA"), label = "Scale factors", value = FALSE, width = NULL))
                     )),            
                   fluidRow(
                     column(11, 
                            renderPlot({
                              
                              color_by_par <- input$color_by_MOFA
                              group_by_par <- input$group_by_MOFA
                              shape_by_par <- input$shape_by_MOFA
                              if (input$color_by_MOFA == "none") color_by_par <- "group"
                              if (input$group_by_MOFA == "none") group_by_par <- "group"
                              if (input$shape_by_MOFA == "none") shape_by_par <- NULL
                              dodge_par <- FALSE
                              if (any(input$add_violin_MOFA, input$add_boxplot_MOFA)) dodge_par = TRUE
                              
                              
                              plot_factor(resMOFA,
                                          factors = min(input$factors_choices_MOFA):max(input$factors_choices_MOFA), 
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
                                          dot_size = 3) 
                            })
                     ))
          ),
          # ---- Tab panel ANOVA posteriori ----
          tabPanel("Relations",
                   
                   renderTable({
                     
                     factors   <- MOFA2::get_factors(resMOFA)
                     ExpDesign <- resMOFA@samples_metadata
                     ExpDesign$group  <- NULL 
                     ExpDesign$sample <- NULL

                     res_aov <- lapply(1:ncol(ExpDesign), FUN = function(i){
                       aov1 <- aov(factors$group1  ~ ExpDesign[,i])
                       sapply(summary(aov1), FUN = function(list_res) list_res[["Pr(>F)"]][[1]])
                     })
                     names(res_aov) <- colnames(ExpDesign)
                     
                     res_res <- do.call("rbind", res_aov)
                     colnames(res_res) <- gsub("Response ", "", colnames(res_res))
                     
                     res_res
                     
                   }, striped = TRUE, rownames = TRUE)
                   
          ),
          
          # ---- Tab panel Heatmap ----
          tabPanel("Heatmap",
                   
                   fluidRow(# buttons - choices for heatmap
                     
                     column(1, numericInput(inputId = session$ns("factor_choice_heatmap_MOFA"),
                                            label = "Factor:",
                                            min = 1,
                                            max =  resMOFA@dimensions$K,
                                            value = 1, step = 1)),
                     column(1,),
                     column(2, radioButtons(inputId = session$ns("view_choice_heatmap_MOFA"),
                                            label = "Data:",
                                            choices = views_names(resMOFA),
                                            selected = views_names(resMOFA)[1]
                     )),
                     column(1,),
                     column(2, numericInput(inputId = session$ns("nfeat_choice_heatmap_MOFA"),
                                            label = "Features:",
                                            min = 1,
                                            max = 500,
                                            value = 50, # default in MOFA function.
                     )),
                     column(1,
                            checkboxInput(inputId = session$ns("denoise_choice_heatmap_MOFA"), label = "Denoise", value = FALSE, width = NULL)
                     ),
                     column(1,),
                     column(2, radioButtons(inputId = session$ns("annot_samples_MOFA"),
                                            label = "Annotation:",
                                            choices = c('none', 
                                                        colnames(resMOFA@samples_metadata %>% 
                                                                   dplyr::select(!c(sample, group)))),
                                            selected = "none"))
                     
                   ),
                   fluidRow(
                     # if(input$annot_samples == "none") input$annot_samples <- NULL
                     
                     column(12, # heatmap is too large with just renderPlot. 
                            renderPlot({
                              annot_samples_values <- input$annot_samples_MOFA
                              if (input$annot_samples_MOFA == "none") annot_samples_values <- NULL
                              
                              observeEvent(input$view_choice_heatmap_MOFA, {
                                updateSliderInput(session, 
                                                  inputId = "nfeat_choice_heatmap_MOFA", 
                                                  max = get_dimensions(resMOFA)$D[input$view_choice_heatmap_MOFA][[1]])
                              })
                              
                              res_heatmap <- plot_data_heatmap(resMOFA, 
                                                               factor = input$factor_choice_heatmap_MOFA, 
                                                               view = input$view_choice_heatmap_MOFA, 
                                                               features = input$nfeat_choice_heatmap_MOFA,
                                                               denoise = input$denoise_choice_heatmap_MOFA,
                                                               annotation_samples = annot_samples_values)
                              grid::grid.draw(res_heatmap)
                            }))
                   ) # fluidrow heatmap
          ), # tabpanel heatmap
          tabPanel("Network",
                   
                   fluidRow(# buttons - choices for network
                     
                     column(1, numericInput(inputId = session$ns("factor_choice_network"),
                                            label = "Factor:",
                                            min = 1,
                                            max =  resMOFA@dimensions$K,
                                            value = 1, step = 1)),
                     column(2, numericInput(inputId = session$ns("abs_weight_network"),
                                            label = "Absolute Weight:",
                                            min = 0.05,
                                            max = 1,
                                            value = 0.8, 
                                            step = 0.05
                     )),                   
                     column(2, numericInput(inputId = session$ns("abs_min_cor_network"),
                                            label = "Minimum absolute \n correlation to display:",
                                            min = 0.05,
                                            max = 1,
                                            value = 0.75, 
                                            step = 0.05
                     )),
                     column(2, radioButtons(inputId = session$ns("network_layout"),
                                            label = "Layout:",
                                            choices = c("Spring", "Circle", "Circle + omics"),
                                            selected = "Circle + omics"
                     )),
                     column(1,  fluidRow(colourInput(inputId = session$ns("posCol"),
                                                     label = "Positive Edges Color:",
                                                     value = "red",
                                                     showColour = c("background"),
                                                     palette = c("square"),
                                                     allowedCols = NULL,
                                                     allowTransparent = FALSE,
                                                     returnName = FALSE,
                                                     closeOnClick = FALSE)),
                            fluidRow(colourInput(inputId = session$ns("negCol"),
                                                 label = "Negative Edges Color:",
                                                 value = "green",
                                                 showColour = c("background"),
                                                 palette = c("square"),
                                                 allowedCols = NULL,
                                                 allowTransparent = FALSE,
                                                 returnName = FALSE,
                                                 closeOnClick = FALSE)),
                     ),
                     column(2,  lapply(1:length(views_names(resMOFA)), function(i) {
                       selectInput(inputId = session$ns(paste0("colors_", views_names(resMOFA)[i])),
                                   label = views_names(resMOFA)[i],
                                   choices = colorPal_choices,
                                   multiple = FALSE, 
                                   selected = colorPal_choices[i]
                       )})
                     )),
                   
                   fluidRow(
                     column(12,
                            renderPlot({
                              
                              colors_list <- lapply(views_names(resMOFA), FUN = function(nam) input[[paste0("colors_", nam)]])
                              names(colors_list) <- views_names(resMOFA)
                              
                              MOFA_cor_network(resMOFA = resMOFA, 
                                               factor_choice = input$factor_choice_network,
                                               abs_weight_network = input$abs_weight_network, 
                                               network_layout = input$network_layout,
                                               omics_colors = colors_list,
                                               posCol = input$posCol,
                                               negCol = input$negCol
                              )
                              
                            }) # renderplot
                     )# column
                   ) # Fluidrow network
          ) # tabpanel network,         
        )# tabsetpanel 
    )#box
  })
  
}