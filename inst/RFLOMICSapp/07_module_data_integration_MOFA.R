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
            p("You can chose which omics data you want to analyze together. It is required to have selected at least two to run an analysis."),
            p("As of now, rflomics only allows you to run an analysis on filtered tables, taking into account differential analysis performed previously.
            You can chose which contrast you want to take the DE genes from. You can select multiple ones (defaults select all the contrasts) 
              and chose to perform the analysis on the union or intersection of the DE lists."),
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
    ### parametres for Co-Exp
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
    
    # validate(
    #   need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis")
    # )
    
    # get good param :
    #listDataSet <- unlist(session$userData$FlomicsMultiAssay@metadata$omicList)
    #names(listDataSet) <- unlist(session$userData$FlomicsMultiAssay@metadata$omicList)
    
    
    listOfContrast <- session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName
    # set param in interface
    tagList(
      
      ## Input parameters
      box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
          
          # Select lists of dataset to integrat
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("selectedData"),
                     label    = "Setect dataset:",
                     choices  = rea.values$datasetDiff,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = rea.values$datasetDiff))),
          
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("selectedContrast"),
                     label    = "Setect contrast:",
                     choices  = listOfContrast,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = listOfContrast))),
          
          # select mode of feature filtering
          fluidRow(
            column(12,
                   
                   radioButtons(inputId = session$ns("filtMode"), label=NULL ,
                                choices = c("union","intersection"),
                                selected = "union", inline = FALSE))),
          
          # set parameters
          fluidRow(
            column(12,
                   selectInput(session$ns("RNAseqTransfo"),
                               label    = "RNAseq transfo :",
                               choices  = c("limma (voom)"),
                               selected = "limma (voom)"))),
          
          fluidRow(
            column(12,
                   selectInput(session$ns("scaleViews"),
                               label    = "Scale views :",
                               choices  = c("FALSE", "TRUE"),
                               selected = "TRUE"))),
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("numfactor"), label="num factors :", value=10, min = 5, max=15))),
          
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("maxiter"), label="Max iteration :", value=1000, min = 1000, max=1000))),
          
          fluidRow(
            column(4, actionButton(session$ns("runMOFA"),"Run"))) ##### ACTION BUTTON
      ))
  })
  
  
  
  observeEvent(input$runMOFA, {
    
    # TODO put everything into one metadata list slot
    # TODO reinitialize everything when the person runs mofa.
    
    local.rea.values$untrainedMOFA <- NULL # reactive ? # ADD 23/06
    local.rea.values$runMOFA   <- FALSE

    local.rea.values$resMOFA <- NULL
    local.rea.values$preparedMOFA <- NULL
    local.rea.values$listResMOFA <- NULL
    
    session$userData$FlomicsMultiAssay@metadata[["MOFA_selected_contrasts"]] <<- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA_selected_filter"]] <<- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA_untrained"]] <<- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA_warnings"]] <<- NULL
    session$userData$FlomicsMultiAssay@metadata[["MOFA_results"]] <<- NULL

    # check nbr dataset to integrate
    # if less then 2 -> error message
    if(length(input$selectedData) < 2){
      
      showModal(modalDialog( title = "Error message", "MOFA need at least 2 datasets!"))
    }
    validate({
      need(length(input$selectedData) >= 2, message="MOFA need at least 2 datasets!")
    })
    
    
    # check nbr of contrast 
    # if less then 1 -> error message
    if(length(input$selectedContrast) == 0){
      
      showModal(modalDialog( title = "Error message", "Select at least one contast!"))
    }
    validate({
      need(length(input$selectedContrast) != 0, message="Select at least one contast!")
    })
    
    # run MOFA
    print(input$selectedData)
    
    local.rea.values$preparedMOFA <- prepareMOFA(session$userData$FlomicsMultiAssay,
                                                  omicsToIntegrate = input$selectedData,
                                                  rnaSeq_transfo = input$RNAseqTransfo,
                                                  choice = "DE", 
                                                  contrasts_names = input$selectedContrast,
                                                  type = input$filtMode,
                                                  group = NULL)
    
    local.rea.values$listResMOFA <- run_MOFA_analysis(local.rea.values$preparedMOFA, 
                                                  scale_views = as.logical(input$scaleViews),
                                                  maxiter = input$maxiter,
                                                  num_factors = input$numfactor) 
    
    local.rea.values$untrainedMOFA <- local.rea.values$listResMOFA$MOFAObject.untrained
    local.rea.values$resMOFA <- local.rea.values$listResMOFA$MOFAObject.trained
    
    
    #### TODO Try to catch MOFA2 warnings and put them on the interface. DOES NOT WORK. 
    # test <- run_MOFA_analysis(session$userData$FlomicsMultiAssay@metadata[["MOFA_untrained"]],
    #                           scale_views = FALSE,
    #                           maxiter = 1000,
    #                           num_factors = 5)
    # local.rea.values <- list()
    local.rea.values$warnings <- NULL
    if(!is.null(names(warnings()))){
      local.rea.values$warnings <- names(warnings())
    }else{local.rea.values$warnings <- "no warnings captured"}
    # output <- list()
    # output$warnings <- renderText({local.rea.values$warnings})
    #### End of catchning warnings.
    
    session$userData$FlomicsMultiAssay@metadata[["MOFA_selected_contrasts"]] <- input$selectedContrast
    session$userData$FlomicsMultiAssay@metadata[["MOFA_selected_filter"]] <- input$filtMode
    session$userData$FlomicsMultiAssay@metadata[["MOFA_untrained"]] <- local.rea.values$untrainedMOFA
    session$userData$FlomicsMultiAssay@metadata[["MOFA_warnings"]] <- local.rea.values$warnings
    session$userData$FlomicsMultiAssay@metadata[["MOFA_results"]] <- local.rea.values$resMOFA
    
    local.rea.values$runMOFA   <- TRUE
    
  }, ignoreInit = TRUE)
  
  
  output$ResultViewUI <- renderUI({
    
    if (local.rea.values$runMOFA == FALSE) return()
    
    plot_height <- function() { # does not work ?
      local.rea.values$plotHeight <- length(input$WeightsPlot_Factors_select)*10
      return(local.rea.values$plotHeight)
    }
    
    box(width=14, solidHeader = TRUE, status = "warning",
        title = "MOFA results",
        
        tabsetPanel(
          
          ###
          tabPanel("Overview", 
                   column(6,
                          renderPlot(MOFA2::plot_data_overview(local.rea.values$resMOFA) + ggtitle("Data Overview"))),
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
          tabPanel("Factors Correlation", 
                   renderPlot(MOFA2::plot_factor_cor(local.rea.values$resMOFA))
          ),
          ### 
          tabPanel("Explained Variance", 
                   fluidRow(
                     column(6, renderPlot({
                       g1 <- MOFA2::plot_variance_explained(local.rea.values$resMOFA, plot_total = TRUE)[[2]] # compute both per factor and per table by default.
                       g1 + ggtitle("Total explained variance per omic data")
                     })),
                     column(6, renderPlot(MOFA2::plot_variance_explained(local.rea.values$resMOFA, x = "view", y = "factor")+ ggtitle("Explained variance by factors and omic data")))
                   )),
          ### 
          tabPanel("Weights Plot", 
                   fluidRow(
                     column(3, sliderInput(inputId = session$ns("WeightsPlot_Factors_select"),
                                           label = 'Factors:',
                                           min = 1, 
                                           max = local.rea.values$resMOFA@dimensions$K, 
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
                         for(j in views_names(local.rea.values$resMOFA)){
                           ggplot_list[[length(ggplot_list)+1]] <- plot_weights(local.rea.values$resMOFA,
                                                                                view = j,
                                                                                factor = i,
                                                                                nfeatures = input$nfeat_choice_WeightsPlot,
                                                                                scale = input$scale_choice_WeightsPlot) + ggtitle(paste0(j, " - Factor ", i))
                         }
                       }
                       
                       ggpubr::ggarrange(plotlist = ggplot_list,
                                         ncol = length(views_names(local.rea.values$resMOFA)),
                                         nrow = (max(input$WeightsPlot_Factors_select)-min(input$WeightsPlot_Factors_select)+1))
                       
                     }, execOnResize = TRUE) # width = 60, height = plot_height(), 
                     
                   ))
          ),
          
          ### 
          tabPanel("Weights table",
                   
                   fluidRow(
                     column(11, sliderInput(inputId = session$ns("Factors_select"),
                                            label = 'Factors:',
                                            min = 1, 
                                            max = local.rea.values$resMOFA@dimensions$K, # pas forcement l'input, MOFA peut decider d'en enlever. 
                                            value = 1:2, step = 1)) 
                   ),
                   fluidRow(
                     column(11, 
                            DT::renderDataTable({
                              resTable <- MOFA2::get_weights(local.rea.values$resMOFA, views = "all", factors = input$Factors_select, abs = FALSE, scale = FALSE, as.data.frame = TRUE)
                              
                              resTable %>% DT::datatable(extensions = 'Buttons',
                                                         options = list(dom = 'lfrtipB',
                                                                        rownames = FALSE,
                                                                        pageLength = 10,
                                                                        buttons = c('csv', 'excel'),
                                                                        lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All"))))
                            }))
                   )),
          tabPanel("Factor Plots",
                   
                   fluidRow(
                   column(2, sliderInput(inputId = session$ns("factors_choices"),
                                         label = 'Factors:',
                                         min = 1, 
                                         max = local.rea.values$resMOFA@dimensions$K, # pas forcement l'input, MOFA peut decider d'en enlever. 
                                         value = 1:2, step = 1)),
                   column(2, radioButtons(inputId = session$ns("color_by"),
                                          label = "Color:",
                                          choices = c('none', 
                                                      colnames(local.rea.values$resMOFA@samples_metadata %>% select(!c(sample, group)))),
                                          selected = "none")
                   ),                           
                   column(2, radioButtons(inputId = session$ns("shape_by"),
                                          label = "Shape:",
                                          choices = c('none', 
                                                      colnames(local.rea.values$resMOFA@samples_metadata %>% select(!c(sample, group)))),
                                          selected = "none")
                   ),                                  
                   column(2, radioButtons(inputId = session$ns("group_by"),
                                          label = "Group by:",
                                          choices = c('none', 
                                                      colnames(local.rea.values$resMOFA@samples_metadata %>% select(!c(sample, group)))),
                                          selected = "none")
                   ),                           
                   column(1,
                          fluidRow(checkboxInput(inputId = session$ns("add_violin"), label = "Violin", value = FALSE, width = NULL)),
                          fluidRow(checkboxInput(inputId = session$ns("add_boxplot"), label = "Boxplot", value = FALSE, width = NULL)),
                          fluidRow(checkboxInput(inputId = session$ns("scale_scatter"), label = "Scale factors", value = FALSE, width = NULL))
                   )),            
                   fluidRow(
                     column(11, 
                            renderPlot({
                              
                              color_by_par <- input$color_by
                              group_by_par <- input$group_by
                              shape_by_par <- input$shape_by
                              if(input$color_by == "none") color_by_par <- "group"
                              if(input$group_by == "none") group_by_par <- "group"
                              if(input$shape_by == "none") shape_by_par <- NULL
                              dodge_par <- FALSE
                              if(any(input$add_violin, input$add_boxplot)) dodge_par = TRUE
                              
                              
                              plot_factor(local.rea.values$resMOFA,
                                          factors = min(input$factors_choices):max(input$factors_choices), 
                                          color_by = color_by_par,
                                          group_by = group_by_par, 
                                          shape_by = shape_by_par,
                                          legend = TRUE,
                                          add_violin = input$add_violin,
                                          violin_alpha = 0.25, 
                                          add_boxplot = input$add_boxplot,
                                          boxplot_alpha = 0.25,
                                          dodge = dodge_par,
                                          scale = input$scale_scatter,
                                          dot_size = 3) 
                            })
                     ))
          ),
        tabPanel("Heatmap",
                   
                   fluidRow(# buttons - choices for heatmap
                     
                     column(1, numericInput(inputId = session$ns("factor_choice_heatmap"),
                                            label = "Factor:",
                                            min = 1,
                                            max =  local.rea.values$resMOFA@dimensions$K,
                                            value = 1, step = 1)),
                     column(1,),
                     column(2, radioButtons(inputId = session$ns("view_choice_heatmap"),
                                            label = "Data:",
                                            choices = views_names(local.rea.values$resMOFA),
                                            selected = views_names(local.rea.values$resMOFA)[1]
                     )),
                     column(1,),
                     column(2, numericInput(inputId = session$ns("nfeat_choice_heatmap"),
                                            label = "Features:",
                                            min = 1,
                                            max = 500,
                                            value = 50, # default in MOFA function.
                     )),
                     column(1,
                            checkboxInput(inputId = session$ns("denoise_choice_heatmap"), label = "Denoise", value = FALSE, width = NULL)
                     ),
                     column(1,),
                     column(2, radioButtons(inputId = session$ns("annot_samples"),
                                            label = "Annotation:",
                                            choices = c('none', 
                                                        colnames(local.rea.values$resMOFA@samples_metadata %>% select(!c(sample, group)))),
                                            selected = "none"))
                     
                   ),
                   fluidRow(
                     # if(input$annot_samples == "none") input$annot_samples <- NULL
                     
                     column(12, # heatmap is too large with just renderPlot. 
                            renderPlot({
                              annot_samples_values <- input$annot_samples
                              if(input$annot_samples == "none") annot_samples_values <- NULL
                              
                              observeEvent(input$view_choice_heatmap, {
                                updateSliderInput(inputId = "nfeat_choice_heatmap", max = get_dimensions(local.rea.values$resMOFA)$D[input$view_choice_heatmap][[1]])
                              })
                              
                              res_heatmap <- MOFA2::plot_data_heatmap(local.rea.values$resMOFA, 
                                                                      factor = input$factor_choice_heatmap, 
                                                                      view = input$view_choice_heatmap, 
                                                                      features = input$nfeat_choice_heatmap,
                                                                      denoise = input$denoise_choice_heatmap,
                                                                      annotation_samples = annot_samples_values)
                              grid::grid.draw(res_heatmap)
                            }))
                   ) # fluidrow heatmap
          ), # tabpanel heatmap
          tabPanel("Network",
                   p("Work in progress :)")
          ),         
        )# tabsetpanel 
    )#box
  })
}




