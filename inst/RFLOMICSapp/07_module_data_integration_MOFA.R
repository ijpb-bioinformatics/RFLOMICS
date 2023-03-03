##########################################
# module 07 : MOFA
##########################################


MOFA_settingUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      box(title = span(tagList(icon('chart-line'), "   ",a("MOFA+", href="https://biofam.github.io/MOFA2/"), tags$small("(Scroll down for instructions)")  )),
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

# tags$a(href="www.rstudio.com", "Click here!")

MOFA_setting <- function(input, output, session, rea.values){
  
  local.rea.values <- reactiveValues(runMOFA = FALSE) # init local reactive values
  
  # list of parameters  
  output$MOFA_ParamUI <- renderUI({
    
    listOfContrast <- session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName
    # set param in interface
    tagList(
      
      ## Input parameters
      box(title = span(tagList(icon("sliders"), "  ", "Setting")), width = 14, status = "warning",
          
          # Select lists of dataset to integrat
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("MOFA_selectedData"),
                     label    = "Select dataset:",
                     choices  = rea.values$datasetDiff,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = rea.values$datasetDiff))),
          
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("MOFA_selectedContrasts"),
                     label    = "Select contrast:",
                     choices  = listOfContrast,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = listOfContrast))),
          
          # select mode of feature filtering
          fluidRow(
            column(12,
                   
                   radioButtons(inputId = session$ns("MOFA_filtMode"), label=NULL ,
                                choices = c("union","intersection"),
                                selected = "union", inline = FALSE))),
          
          # set parameters
          fluidRow(
            column(12,
                   selectInput(session$ns("MOFA_RNAseqTransfo"),
                               label    = "RNAseq transfo :",
                               choices  = c("limma (voom)"),
                               selected = "limma (voom)"))),
          
          fluidRow(
            column(12,
                   selectInput(session$ns("MOFA_scaleViews"),
                               label    = "Scale views :",
                               choices  = c("FALSE", "TRUE"),
                               selected = "TRUE"))),
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("MOFA_numfactor"), label="Num factors:", value=10, min = 5, max=15))),
          
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("MOFA_maxiter"), label="Max iteration:", value=1000, min = 1000, max=1000))),
          
          fluidRow(
            column(4, actionButton(session$ns("runMOFA"),"Run"))) ##### ACTION BUTTON
      ))
  })
  
  observeEvent(input$runMOFA, {
    
    library(MOFA2) 
    
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
    preparedMOFA <- NULL
    listResMOFA <- NULL
    
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_selected_contrasts"]] <- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_selected_filter"]] <- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_untrained"]] <- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_warnings"]] <- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_results"]] <- NULL
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Checks ", 10, "%", sep=""))
    #----------------------#
    
    # check nbr dataset to integrate
    # if less then 2 -> error message
    if(length(input$MOFA_selectedData) < 2){
      showModal(modalDialog( title = "Error message", "MOFA needs at least 2 datasets!"))
    }
    validate({
      need(length(input$MOFA_selectedData) >= 2, message="MOFA needs at least 2 datasets!")
    })
    
    
    # check nbr of contrast 
    # if less than 1 -> error message
    if(length(input$MOFA_selectedContrasts) == 0){
      showModal(modalDialog( title = "Error message", "Select at least one contast!"))
    }
    validate({
      need(length(input$MOFA_selectedContrasts) != 0, message="Select at least one contast!")
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
      method = "MOFA"
    )
    
    preparedMOFA <- do.call("prepareForIntegration", list_args_prepare_MOFA)
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Running MOFA ", 30, "%", sep = ""))
    #----------------------#
    
    list_run_MOFA <- list(
      object = preparedMOFA,
      scale_views = as.logical(input$MOFA_scaleViews),
      maxiter = input$MOFA_maxiter,
      num_factors = input$MOFA_numfactor
    )
    
    listResMOFA <- do.call("run_MOFA_analysis", list_run_MOFA) 
    
    # untrainedMOFA <- listResMOFA$MOFAObject.untrained
    # resMOFA <- listResMOFA$MOFAObject.trained
    
    
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
    
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_selected_contrasts"]] <- input$MOFA_selectedContrasts
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_selected_filter"]] <- input$MOFA_filtMode
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_untrained"]] <- listResMOFA$MOFAObject.untrained
    # session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_warnings"]] <- warnings
    session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_results"]] <- listResMOFA$MOFAObject.trained
    
    local.rea.values$runMOFA   <- TRUE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Finished ", 100,"%", sep = ""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
  
  output$ResultViewUI <- renderUI({
    
    if (local.rea.values$runMOFA == FALSE) return()
    
    resMOFA <- session$userData$FlomicsMultiAssay@metadata[["MOFA"]][["MOFA_results"]]
    
    plot_height <- function() { # does not work ?
      local.rea.values$plotHeight <- length(input$WeightsPlot_Factors_select)*10
      return(local.rea.values$plotHeight)
    }
    
    colorPal_choices <- c("Greens", "Purples", "Oranges", "Reds", "Greys", "Blues")
    
    box(width=14, solidHeader = TRUE, status = "success",
        title = "MOFA results",
        
        tabsetPanel(
          
          ###
          # ---- Tab panel Factors Overview ----
          tabPanel("Overview", 
                   column(12,
                          renderPlot(MOFA2::plot_data_overview(resMOFA) + ggtitle("Data Overview"))),
                   # column(6,
                   #        h4("Warnings:"),
                   #        renderText(expr ={
                   #          # if("no warning captured"%in%local.rea.values$warnings){
                   #          #   print(local.rea.values$warnings)
                   #          # }else{
                   #          # warnsms <- paste(local.rea.values$warnings, collapse = "\n")
                   #          # for(warns in local.rea.values$warnings) print(warns)
                   #          # HTML(paste(local.rea.values$warnings, collapse = '<br/>'))
                   #          HTML(paste(local.rea.values$warnings, collapse = ''))
                   #          # print(paste(local.rea.values$warnings, collapse = "\n"))
                   #          # for(i in 1:length(local.rea.values$warnings)){
                   #          #   cat(local.rea.values$warnings[i])
                   #          # }
                   #          # }
                   #        })
                   # )
          ),
          ### 
          # ---- Tab panel Factors Correlation ----
          tabPanel("Factors Correlation", 
                   renderPlot(MOFA2::plot_factor_cor(resMOFA))
          ),
          ### 
          # ---- Tab panel Explained Variance ----
          tabPanel("Explained Variance", 
                   fluidRow(
                     column(6, renderPlot({
                       g1 <- MOFA2::plot_variance_explained(resMOFA, plot_total = TRUE)[[2]] # compute both per factor and per table by default.
                       g1 + ggtitle("Total explained variance per omic data")
                     })),
                     column(6, renderPlot(MOFA2::plot_variance_explained(resMOFA, x = "view", y = "factor")+ ggtitle("Explained variance by factors and omic data")))
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
                              for(i in min(input$WeightsPlot_Factors_select):max(input$WeightsPlot_Factors_select)){
                                for(j in MOFA2::views_names(resMOFA)){
                                  ggplot_list[[length(ggplot_list)+1]] <- MOFA2::plot_weights(resMOFA,
                                                                                              view = j,
                                                                                              factor = i,
                                                                                              nfeatures = input$nfeat_choice_WeightsPlot,
                                                                                              scale = input$scale_choice_WeightsPlot) + ggtitle(paste0(j, " - Factor ", i))
                                }
                              }
                              
                              ggpubr::ggarrange(plotlist = ggplot_list,
                                                ncol = length(MOFA2::views_names(resMOFA)),
                                                nrow = (max(input$WeightsPlot_Factors_select)-min(input$WeightsPlot_Factors_select)+1))
                              
                            }, execOnResize = TRUE) # width = 60, height = plot_height(), 
                            
                     ))
          ),
          
          ### 
          # ---- Tab panel Weights table ----
          tabPanel("Weights table",
                   
                   fluidRow(
                     column(11, sliderInput(inputId = session$ns("Factors_select_MOFA"),
                                            label = 'Factors:',
                                            min = 1, 
                                            max = resMOFA@dimensions$K, # pas forcement l'input, MOFA peut decider d'en enlever. 
                                            value = 1:2, step = 1)) 
                   ),
                   fluidRow(
                     column(11, 
                            DT::renderDataTable({
                              resTable <- MOFA2::get_weights(resMOFA, views = "all", factors = min(input$Factors_select_MOFA):max(input$Factors_select_MOFA),
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
                              if(input$color_by_MOFA == "none") color_by_par <- "group"
                              if(input$group_by_MOFA == "none") group_by_par <- "group"
                              if(input$shape_by_MOFA == "none") shape_by_par <- NULL
                              dodge_par <- FALSE
                              if(any(input$add_violin_MOFA, input$add_boxplot_MOFA)) dodge_par = TRUE
                              
                              
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
                     
                     factors <- MOFA2::get_factors(resMOFA)
                     
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
                              if(input$annot_samples_MOFA == "none") annot_samples_values <- NULL
                              
                              observeEvent(input$view_choice_heatmap_MOFA, {
                                updateSliderInput(session, 
                                                  inputId = "nfeat_choice_heatmap_MOFA", 
                                                  max = get_dimensions(resMOFA)$D[input$view_choice_heatmap_MOFA][[1]])
                              })
                              
                              res_heatmap <- MOFA2::plot_data_heatmap(resMOFA, 
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
                     column(2,  lapply(1:length(MOFA2::views_names(resMOFA)), function(i) {
                       selectInput(inputId = session$ns(paste0("colors_", MOFA2::views_names(resMOFA)[i])),
                                   label = MOFA2::views_names(resMOFA)[i],
                                   choices = colorPal_choices,
                                   multiple = FALSE, 
                                   selected = colorPal_choices[i]
                       )})
                     )),
                   
                   fluidRow(
                     column(12,
                            renderPlot({
                              
                              # Correlation matrix is done on all ZW, not on the selected factor. 
                              data_reconst_list <- lapply(MOFA2::get_weights(resMOFA), FUN = function(mat){
                                MOFA2::get_factors(resMOFA)$group1 %*% t(mat)})
                              data_reconst <- do.call(cbind, data_reconst_list)
                              cor_mat <- stats::cor(data_reconst)
                              
                              features_metadata <- do.call(rbind, lapply(1:length(MOFA2::get_weights(resMOFA)), FUN = function(i){
                                mat_weights <- data.frame(MOFA2::get_weights(resMOFA, scale = TRUE)[[i]])
                                mat_weights$Table <- names(MOFA2::get_weights(resMOFA))[i]
                                return(mat_weights)
                              }))
                              
                              factor_selected <- paste0("Factor", input$factor_choice_network)
                              
                              feature_filtered <- features_metadata %>% 
                                tibble::rownames_to_column("EntityName") %>%
                                dplyr::mutate(F_selected = abs(get(factor_selected))) %>% 
                                dplyr::arrange(desc(abs(F_selected))) %>% 
                                dplyr::group_by(Table) %>% 
                                dplyr::filter(abs(F_selected)>input$abs_weight_network)
                              
                              if(nrow(feature_filtered)>0){
                                
                                omics_colors <- lapply(unique(feature_filtered$Table), FUN = function(omicTable){
                                  RColorBrewer::brewer.pal(name = input[[paste0("colors_", omicTable)]], n = 5)
                                })
                                names(omics_colors) <- unique(feature_filtered$Table)
                                
                                feature_filtered <- feature_filtered %>% dplyr::group_by(Table) %>%
                                  dplyr::mutate(Color = cut(F_selected, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
                                feature_filtered$Color2 <- sapply(1:nrow(feature_filtered), 
                                                                  FUN = function(i) omics_colors[[feature_filtered$Table[i]]][as.numeric(feature_filtered$Color[i])])
                                
                                # Layout
                                layout_arg <- tolower(input$network_layout)
                                if(tolower(layout_arg) == tolower("Circle + omics")){
                                  layout_arg <- "groups"
                                }
                                
                                # Network main graph
                                cor_display <- cor_mat[rownames(cor_mat) %in% feature_filtered$EntityName, colnames(cor_mat) %in% feature_filtered$EntityName]
                                
                                if(any(abs(cor_display[upper.tri(cor_display)])>=input$abs_min_cor_network)){
                                  qgraph_plot <- qgraph::qgraph(cor_display, minimum = input$abs_min_cor_network, 
                                                                cut = 0,
                                                                shape = "rectangle", labels = rownames(cor_display), vsize2 = 2, 
                                                                vsize = sapply(rownames(cor_display), nchar)*1.1,  layout = layout_arg,
                                                                esize = 2,
                                                                groups = gsub("[.]filtred", "", features_metadata$Table[match(rownames(cor_display), rownames(features_metadata))]),
                                                                posCol = input$posCol, negCol = input$negCol, 
                                                                details = FALSE,  legend = FALSE,
                                                                color = feature_filtered$Color2[match(rownames(cor_display), feature_filtered$EntityName)])
                                  qgraph_plot <- recordPlot()
                                  
                                  # Legend
                                  legend_matrix <- do.call("rbind", lapply(MOFA2::views_names(resMOFA), FUN = function(nam) RColorBrewer::brewer.pal(5, input[[paste0("colors_", nam)]])))
                                  colnames(legend_matrix) <- c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8, 1]")
                                  rownames(legend_matrix) <- gsub("[.]filtred|colors_", "", MOFA2::views_names(resMOFA))
                                  legend.reshape <- reshape2::melt(legend_matrix)
                                  
                                  gg.legend <-  ggplot2::ggplot(legend.reshape, ggplot2::aes(x = Var2, y = Var1)) + 
                                    ggplot2::geom_tile(fill = legend.reshape$value) + xlab("") + ylab("") + 
                                    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
                                                       axis.ticks.y = element_blank(),
                                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                                  
                                  gg.legend <- ggpubr::ggarrange(gg.legend, nrow = 3, ncol = 1) 
                                  
                                  # Actual plotting
                                  ggpubr::ggarrange(plotlist = list(qgraph_plot, gg.legend), nrow = 1, ncol = 2, widths = c(3, 1))
                                }else{
                                  renderText({print("There is nothing to plot. Please try to lower the absolute weight, the correlation threshold or change the factor.")})
                                }
                              }else{
                                renderText({print("There is nothing to plot. Please try to lower the absolute weight, the correlation threshold or change the factor.")})
                              }
                              
                              
                            }) # renderplot
                     )# column
                   ) # Fluidrow network
          ) # tabpanel network,         
        )# tabsetpanel 
    )#box
  })
}




