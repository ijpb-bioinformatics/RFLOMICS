##########################################
# module 08 : mixOmics
##########################################

MixOmics_settingUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      box(title = span(tagList(icon('chart-line'), "   ",a("MixOmics", href="http://mixomics.org/"), tags$small("(Scroll down for instructions)")  )),
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
            
          ))),
    ##
    fluidRow(
      column(3, uiOutput(ns("MixOmics_ParamUI"))),
      column(9, uiOutput(ns("MOResultViewUI"))))
  )
}

MixOmics_setting <- function(input, output, session, rea.values){
  
  local.rea.values <- reactiveValues(runMixOmics = FALSE) # init local reactive values
  
  ## list of parameters  
  output$MixOmics_ParamUI <- renderUI({
    
    # Parameters to put
    
    listOfContrast <- session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName
    
    # set param in interface
    tagList(
      
      column(width = 12,
             ## Input parameters (For blocks/X)
             box(title = span(tagList(icon("sliders"), "  ", "Setting Blocks")), width = 12, status = "warning",
                 
                 # Select lists of dataset to integrate
                 fluidRow(
                   column(12,
                          
                          pickerInput(
                            inputId  = session$ns("MO_selectedData"),
                            label    = "Select dataset",
                            choices  = rea.values$datasetDiff,
                            options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            multiple = TRUE,
                            selected = rea.values$datasetDiff))),
                 
                 fluidRow(
                   column(12,
                          
                          pickerInput(
                            inputId  = session$ns("MO_selectedContrast"),
                            label    = "Select contrasts",
                            choices  = listOfContrast,
                            options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            multiple = TRUE,
                            selected = listOfContrast))),
                 
                 # select mode of feature filtering
                 fluidRow(
                   column(12,
                          radioButtons(inputId  = session$ns("MO_filtMode"), 
                                       label    = "Select type of filtering" ,
                                       choices  = c("union", "intersection"),
                                       selected = "union", inline = FALSE)),
                 ),
                 # set parameters
                 fluidRow(
                   column(12,
                          selectInput(inputId  = session$ns("MO_RNAseqTransfo"),
                                      label    = "RNAseq transfo",
                                      choices  = c("limma (voom)"),
                                      selected = "limma (voom)")))
                 
             ),
             box(title = span(tagList(icon("sliders"), "  ", "Settings for analysis")), width = 12, status = "warning",
                 
                 fluidRow(
                   column(12,
                          column(6, checkboxInput(inputId = session$ns("MO_scale_views"), label = "Scale Datasets", value = TRUE, width = NULL)),
                          column(6, checkboxInput(inputId = session$ns("MO_sparsity"), label = "Sparse analysis", value = FALSE, width = NULL)),
                          column(6, numericInput(inputId = session$ns("MO_ncomp"), label = "Components", value = 5, min = 1, max= 20)),
                          column(6, numericInput(inputId = session$ns("MO_cases_to_try"), label = "Tuning cases", value = 5, min = 1, max= 100)),
                          # column(6, numericInput(inputId = session$ns("link_datasets"), label = "Link between datasets", value = 1, min = 0, max= 1)),
                          # column(6, numericInput(inputId = session$ns("link_response"), label = "Link to response", value = 1, min = 0, max= 1))
                   ),
                   column(12,
                          checkboxGroupInput(
                            inputId  = session$ns("MO_selectedResponse"),
                            label    = "Select response variables",
                            choices  = c(colnames(colData(session$userData$FlomicsMultiAssay))),
                            selected = colnames(colData(session$userData$FlomicsMultiAssay))[1]))
                 ),
                 fluidRow(
                   column(12, actionButton(session$ns("runMixOmics"), "Run Analysis"))
                 )
             ),
      ) # column 3
    )# taglist
    
  })
  
  ## observe the button run mixOmics
  observeEvent(input$runMixOmics, {
    
    library(mixOmics)
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run MixOmics", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    print("# 8- MixOmics Analysis")
    
    local.rea.values$runMixOmics   <- FALSE
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Checks ", 10, "%", sep = ""))
    #----------------------#
    
    # Check selection of response variables (at least one)
    if (is.null(input$MO_selectedResponse)) {
      showModal(modalDialog(title = "Error message", "To run MixOmics, please select at least one response variable"))
    }
    validate({ 
      need(!is.null(input$MO_selectedResponse), "To run MixOmics, please select at least one response variable") 
    })
    
    # check number of tables (at least two)
    if (length(input$MO_selectedData) < 2) {
      showModal(modalDialog(title = "Error message", "To run a multi-omic analysis, please select at least two tables"))
    }
    validate({ 
      need(length(input$MO_selectedData) >= 2, "To run a multi-omic analysis, please select at least two tables") 
    })
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Preparing object ", 20, "%", sep = ""))
    #----------------------#
    
    # Prepare for MixOmics run  
    
    list_args_MO <- list(
      object = session$userData$FlomicsMultiAssay,
      omicsToIntegrate = paste0(input$MO_selectedData, ".filtred"),
      rnaSeq_transfo = input$MO_RNAseqTransfo,
      choice = "DE", 
      contrasts_names = input$MO_selectedContrast,
      type = input$MO_filtMode,
      group = NULL,
      method = "MixOmics",
      scale_views = input$MO_scale_views,
      ncomp = input$MO_ncomp,   
      link_datasets = 0.5,
      link_response = 1,
      sparsity = input$MO_sparsity,
      cases_to_try = input$MO_cases_to_try,
      selectedResponse = input$MO_selectedResponse,
      cmd = TRUE, 
      silent = TRUE
      
    )
    
    # ---- progress bar ----#
    progress$inc(1/10, detail = paste("Running MixOmics ", 30, "%", sep = ""))
    # ----------------------#
    
    # Run the analysis
    # print("#     =>Running MixOmics")
    
    session$userData$FlomicsMultiAssay <- do.call(getFromNamespace("integrationWrapper", ns = "RFLOMICS"), list_args_MO)
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Finished ", 100,"%", sep = ""))
    #----------------------#
    
    local.rea.values$runMixOmics <- TRUE
    
    
  })
  
  ## Output results
  output$MOResultViewUI <- renderUI({
    
    if (!local.rea.values$runMixOmics) return()
    
    lapply(names(session$userData$FlomicsMultiAssay@metadata[["mixOmics"]]), function(listname) { 
      
      # MAE@metadata$mixOmics$Genotype$MixOmics_results
      
      Data_res <- session$userData$FlomicsMultiAssay@metadata[["mixOmics"]][[listname]]$MixOmics_results
      
      fluidRow(
        
        box(width = 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success", title = listname,
            
            tabsetPanel(
              # ---- Tab panel Overview ----
              tabPanel("Overview",
                       column(6 , DT::renderDataTable({
                         
                         df <- t(sapply(Data_res$X, dim))
                         colnames(df) <- c("Ind", "Features")
                         
                         if (input$MO_sparsity) {
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
                                                  selectedResponse = listname)
                       })),
                       
              ),
              # ---- Tab panel Individuals ----
              tabPanel("Individuals",
                       column(1,
                              checkboxInput(inputId = session$ns(paste0(listname, "ellipse_choice")), label = "Ellipses", value = FALSE, width = NULL),
                              numericInput(inputId = session$ns(paste0(listname, "ind_comp_choice_1")),
                                           label = "Comp x:",
                                           min = 1,
                                           max =  input$MO_ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(listname, "ind_comp_choice_2")),
                                           label = "Comp y:",
                                           min = 1,
                                           max =  input$MO_ncomp,
                                           value = 2, step = 1)
                       ),
                       column(11 , renderPlot(suppressWarnings(
                         mixOmics::plotIndiv(Data_res, 
                                             comp = c(input[[paste0(listname, "ind_comp_choice_1")]], input[[paste0(listname, "ind_comp_choice_2")]]),
                                             ellipse = input[[paste0(listname, "ellipse_choice")]],
                                             legend = TRUE))))
              ),
              # ---- Tab panel Features ----
              tabPanel("Features",
                       column(1,
                              checkboxInput(inputId = session$ns("overlap"), label = "Overlap", value = FALSE, width = NULL),
                              numericInput(inputId = session$ns(paste0(listname, "var_comp_choice_1")),
                                           label = "Comp x:",
                                           min = 1,
                                           max =  input$MO_ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(listname, "var_comp_choice_2")),
                                           label = "Comp y:",
                                           min = 1,
                                           max =  input$MO_ncomp,
                                           value = 2, step = 1)
                       ),
                       column(11 , renderPlot(mixOmics::plotVar(Data_res, 
                                                                comp = c(input[[paste0(listname, "var_comp_choice_1")]], input[[paste0(listname, "var_comp_choice_2")]]),
                                                                overlap = input$overlap,
                                                                legend = TRUE)))
              ),     
              # ---- Tab panel Loadings ----
              tabPanel("Loadings",
                       column(1,
                              numericInput(inputId = session$ns(paste0(listname, "Load_comp_choice")),
                                           label = "Component:",
                                           min = 1,
                                           max =  input$MO_ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(listname, "Load_ndisplay")),
                                           label = "Number of features to display:",
                                           min = 1,
                                           max =  max(sapply(Data_res$X, ncol)),
                                           value = 25, step = 1),
                       ),
                       column(11 , renderPlot(mixOmics::plotLoadings(Data_res, 
                                                                     comp = input[[paste0(listname, "Load_comp_choice")]],
                                                                     ndisplay = input[[paste0(listname, "Load_ndisplay")]])))
              ),   
              # ---- Tab panel Tuning ----
              #   tabPanel("Tuning",
              #            # TODO 
              #            column(12 , renderPlot(mixOmics::plot(local.rea.values$tuning_res)))
              #   )), # plot.tune ne fait pas partie du package dans ma version ?!
              tabPanel("Networks",
                       # TODO Prevoir bouton pour comp selectionnee
                       # Comp fonctionne plus ?!
                       column(1, numericInput(inputId = session$ns(paste0(listname, "Network_cutoff")),
                                              label = "Cutoff:",
                                              min = 0,
                                              max =  1,
                                              value = 0.9, step = 0.05)),
                       column(11 , renderPlot(mixOmics::network(mat = Data_res, 
                                                                # comp = 1:2, 
                                                                blocks = 1:length(input$MO_selectedData),
                                                                cutoff = input[[paste0(listname, "Network_cutoff")]], 
                                                                shape.node = rep("rectangle", length(input$MO_selectedData)))))
              ), 
              # ---- Tab  Panel CircosPlot & cimPlot ----
              tabPanel("CircosPlot",
                       if (is(Data_res, "block.splsda")) {
                         fluidRow(
                           column(1, numericInput(inputId = session$ns(paste0(listname, "Circos_cutoff")),
                                                  label = "Cutoff:",
                                                  min = 0,
                                                  max =  1,
                                                  value = 0.9, step = 0.05)),
                           column(11, renderPlot(mixOmics::circosPlot(Data_res,
                                                                      cutoff = input[[paste0(listname, "Circos_cutoff")]])))
                         )
                       }else{
                         renderText({print("This plot is only available for sparse multi-block discriminant analysis results.")})
                       }
              ),
              tabPanel("CimPlot",
                       if (is(Data_res, "block.splsda")) {
                         print("=> Rendering cimPlot, be patient!")
                         fluidRow(
                           column(1,
                                  numericInput(inputId = session$ns(paste0(listname, "cimComp")),
                                               label = "Comp",
                                               min = 1,
                                               max = input$MO_ncomp,
                                               value = 1, step = 1)),
                           column(12, 
                                  
                                  renderPlot(mixOmics::cimDiablo(Data_res, 
                                                                 legend.position = "bottomleft",
                                                                 size.legend = 0.8,
                                                                 comp = input[[paste0(listname, "cimComp")]])))
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