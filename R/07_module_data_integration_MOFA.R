##########################################
# module 07 : MOFA
##########################################


MOFA_settingUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
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
          ))),
    ### parametres 
    fluidRow(
      column(3, uiOutput(ns("MOFA_ParamUI"))),
      column(9, uiOutput(ns("ResultViewUI"))))
  )
}

MOFA_setting <- function(input, output, session, rea.values){
  
  local.rea.values <- reactiveValues(runMOFA = FALSE) 
  
  # list of parameters  
  output$MOFA_ParamUI <- renderUI({
    
    listOfContrast <- getSelectedContrasts(session$userData$FlomicsMultiAssay)$contrastName
    
    # set param in interface
    tagList(
      
      ## Input parameters
      box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
          
          # Select lists of dataset to integrat
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("MOFA_selectedData"),
                     label    = "Select dataset",
                     choices  = rea.values$datasetDiff,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = rea.values$datasetDiff))),
          
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("MOFA_selectedContrasts"),
                     label    = "Select contrasts",
                     choices  = listOfContrast
                     ))),
          
          # select mode of feature filtering
          fluidRow(
            column(12,
                   
                   radioButtons(inputId = session$ns("MOFA_filtMode"), label = NULL ,
                                choices = c("union","intersection"),
                                selected = "union", inline = FALSE))),
          
          # set parameters
          fluidRow(
            column(12,
                   selectInput(session$ns("MOFA_RNAseqTransfo"),
                               label    = "RNAseq transfo",
                               choices  = c("limma (voom)"),
                               selected = "limma (voom)"))),
          
          fluidRow(
            column(12,
                   selectInput(session$ns("MOFA_scaleViews"),
                               label    = "Scale views",
                               choices  = c("FALSE", "TRUE"),
                               selected = "TRUE"))),
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("MOFA_numfactor"), label = "Num factors:", 
                                value = 10, min = 5, max = 15))),
          
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("MOFA_maxiter"), label = "Max iteration:", 
                                value = 1000, min = 1000, max = 1000))),
          
          fluidRow(
            column(4, actionButton(session$ns("runMOFA"),"Run")))
      ))
  })
  
  observeEvent(input$runMOFA, {
    
    # TODO put everything into one metadata list slot
    # TODO reinitialize everything when the person runs mofa.
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run MOFA2", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    print("# 7- MOFA Analysis")
    
    local.rea.values$runMOFA   <- FALSE
    
    # untrainedMOFA <- NULL 
    # resMOFA <- NULL
    
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Checks ", 10, "%", sep = ""))
    #----------------------#
    
    # check nbr dataset to integrate
    # if less then 2 -> error message
    if (length(input$MOFA_selectedData) < 2) {
      showModal(modalDialog( title = "Error message", "MOFA needs at least 2 datasets!"))
    }
    validate({
      need(length(input$MOFA_selectedData) >= 2, message = "MOFA needs at least 2 datasets!")
    })
    
    
    # check nbr of contrast 
    # if less than 1 -> error message
    if (length(input$MOFA_selectedContrasts) == 0) {
      showModal(modalDialog( title = "Error message", "Select at least one contast!"))
    }
    validate({
      need(length(input$MOFA_selectedContrasts) != 0, message = "Select at least one contast!")
    })
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Preparing object ", 20, "%", sep = ""))
    #----------------------#
    
    list_args_prepare_MOFA <- list(
      object = session$userData$FlomicsMultiAssay,
      omicsToIntegrate = input$MOFA_selectedData,
      rnaSeq_transfo = input$MOFA_RNAseqTransfo,
      choice = "DE", 
      contrasts_names = input$MOFA_selectedContrasts,
      type = input$MOFA_filtMode,
      group = NULL, 
      method = "MOFA",
      scale_views = as.logical(input$MOFA_scaleViews),
      maxiter = input$MOFA_maxiter,
      num_factors = input$MOFA_numfactor,
      cmd = TRUE, 
      silent = TRUE
    )
    
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Running MOFA ", 30, "%", sep = ""))
    #----------------------#
    
    session$userData$FlomicsMultiAssay <- do.call(getFromNamespace("integrationWrapper", ns = "RFLOMICS"), list_args_prepare_MOFA)

    #### TODO Try to catch MOFA2 warnings and put them on the interface. DOES NOT WORK. 
    # test <- run_MOFA_analysis(session$userData$FlomicsMultiAssay@metadata[["MOFA_untrained"]],
    #                           scale_views = FALSE,
    #                           MOFA_maxiter = 1000,
    #                           num_factors = 5)
    # local.rea.values <- list()
    # local.rea.values$warnings <- NULL
    # if(!is.null(names(warnings()))){
    #   local.rea.values$warnings <- names(warnings())
    # }else{local.rea.values$warnings <- "no warnings captured"}
    # output <- list()
    # output$warnings <- renderText({local.rea.values$warnings})
    #### End of catchning warnings.
    
    local.rea.values$runMOFA   <- TRUE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Finished ", 100,"%", sep = ""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
  
  output$ResultViewUI <- renderUI({
    
    if (local.rea.values$runMOFA == FALSE) return()
    
    resMOFA <- session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_results"]]
    
    # plot_height <- function() { # does not work ?
    #   local.rea.values$plotHeight <- length(input$WeightsPlot_Factors_select)*10
    #   return(local.rea.values$plotHeight)
    # }
    
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
                     
                     factors <- get_factors(resMOFA)
                     
                     res_aov <- lapply(1:ncol(session$userData$FlomicsMultiAssay@metadata$design@ExpDesign), FUN = function(i){
                       aov1 <- aov(factors$group1  ~ session$userData$FlomicsMultiAssay@metadata$design@ExpDesign[,i])
                       sapply(summary(aov1), FUN = function(list_res) list_res[["Pr(>F)"]][[1]])
                     })
                     names(res_aov) <- colnames(session$userData$FlomicsMultiAssay@metadata$design@ExpDesign)
                     
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




