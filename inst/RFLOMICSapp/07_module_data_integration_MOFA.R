##########################################
# module 07 : MOFA
##########################################

MOFA_settingUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      box(title = span(tagList(icon('chart-line'), "   ",a("MOFA", href="https://biofam.github.io/MOFA2/"), tags$small("(Scroll down for instructions)")  )),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
          div(      
            h4(tags$span("Parameters set up:", style = "color:orange")),
            p("C'est la pour marquer des notes pour l'utilisateur")
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
    #listDataSet <- unlist(FlomicsMultiAssay@metadata$omicList)
    #names(listDataSet) <- unlist(FlomicsMultiAssay@metadata$omicList)
    
    
    listOfContrast <- FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName
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
    
    # local.rea.values <- reactiveValues(runMOFA = FALSE) #### ????
    
    # run MOFA
    print(input$selectedData)
    FlomicsMultiAssay@metadata[["MOFA_run"]] <- FALSE
    
    untrainedMOFA <- prepareMOFA(FlomicsMultiAssay,
                                 omicsToIntegrate = input$selectedData,
                                 rnaSeq_transfo = input$RNAseqTransfo,
                                 choice = "DE",
                                 contrast = input$selectedContrast, 
                                 group = NULL)
    FlomicsMultiAssay@metadata[["MOFA_results"]] <<- run_MOFA_analysis(untrainedMOFA, 
                                                                       scale_views = input$scaleViews,
                                                                       maxiter = input$maxiter,
                                                                       num_factors = input$numfactor) 
    local.rea.values$resMOFA <- FlomicsMultiAssay@metadata[["MOFA_results"]]
    
    FlomicsMultiAssay@metadata[["MOFA_run"]] <- TRUE
    local.rea.values$runMOFA   <- TRUE
    
  }, ignoreInit = TRUE)
  
  
  output$ResultViewUI <- renderUI({
    
    if (local.rea.values$runMOFA == FALSE) return()
    
    box(width=14, solidHeader = TRUE, status = "warning",
        title = "MOFA results",
        
        tabsetPanel(
          
          ###
          tabPanel("Overview", 
                   renderPlot(MOFA2::plot_data_overview(FlomicsMultiAssay@metadata[["MOFA_results"]]))
          ),
          ### 
          tabPanel("Factors Correlation", 
                   renderPlot(MOFA2::plot_factor_cor(FlomicsMultiAssay@metadata[["MOFA_results"]]))
          ),
          ### 
          tabPanel("Explained Variance", 
                   fluidRow(
                     column(6, renderPlot({
                       g1 <- MOFA2::plot_variance_explained(FlomicsMultiAssay@metadata[["MOFA_results"]], plot_total = TRUE)[[2]] # compute both per factor and per table by default.
                       g1 + ggtitle("Total explained variance per omic data")
                     })),
                     column(6, renderPlot(MOFA2::plot_variance_explained(FlomicsMultiAssay@metadata[["MOFA_results"]], x = "view", y = "factor")+ ggtitle("Explained variance by factors and omic data")))
                   )),
          ### 
          tabPanel("Weights Plot", 
                   fluidRow(
                     column(3, sliderInput(inputId = session$ns("WeightsPlot_Factors_select"),
                                           label = 'Factors:',
                                           min = 1, 
                                           max = FlomicsMultiAssay@metadata[["MOFA_results"]]@dimensions$K, 
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
                     renderPlot({
                       
                       # plot_weights(FlomicsMultiAssay@metadata[["MOFA_results"]], 
                       #              view = 1, 
                       #              factor = 1,
                       #              nfeatures = input$nfeat_choice_WeightsPlot,
                       #              scale = input$scale_choice_WeightsPlot)
                       
                       ggplot_list <- list()
                       for(i in min(input$WeightsPlot_Factors_select):max(input$WeightsPlot_Factors_select)){
                         for(j in views_names(FlomicsMultiAssay@metadata[["MOFA_results"]])){
                           ggplot_list[[length(ggplot_list)+1]] <- plot_weights(FlomicsMultiAssay@metadata[["MOFA_results"]],
                                                                                 view = j,
                                                                                 factor = i,
                                                                                 nfeatures = input$nfeat_choice_WeightsPlot,
                                                                                 scale = input$scale_choice_WeightsPlot) + ggtitle(paste0(j, " - Factor ", i))
                         }
                       }

                       ggpubr::ggarrange(plotlist = ggplot_list,
                                         ncol = length(views_names(FlomicsMultiAssay@metadata[["MOFA_results"]])),
                                         nrow = (max(input$WeightsPlot_Factors_select)-min(input$WeightsPlot_Factors_select)+1))
                       
                     })
                     
                   )
          ),
          
          ### 
          tabPanel("Weights table",
                   
                   fluidRow(
                     column(11, sliderInput(inputId = session$ns("Factors_select"),
                                            label = 'Factors:',
                                            min = 1, 
                                            max = FlomicsMultiAssay@metadata[["MOFA_results"]]@dimensions$K, # pas forcement l'input, MOFA peut decider d'en enlever. 
                                            value = c(1,2), step = 1)) 
                   ),
                   fluidRow(
                     column(11, 
                            DT::renderDataTable({
                              resTable <- MOFA2::get_weights(FlomicsMultiAssay@metadata[["MOFA_results"]], views = "all", factors = input$Factors_select, abs = FALSE, scale = FALSE, as.data.frame = TRUE)
                              
                              resTable %>% DT::datatable(extensions = 'Buttons',
                                                         options = list(dom = 'lfrtipB',
                                                                        rownames = FALSE,
                                                                        pageLength = 10,
                                                                        buttons = c('csv', 'excel'),
                                                                        lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All"))))
                            }))
                   )),
          tabPanel("Heatmap",
                   
                   fluidRow(# buttons - choices for heatmap
                     
                     column(1, numericInput(inputId = session$ns("factor_choice_heatmap"),
                                            label = "Factor:",
                                            min = 1,
                                            max =  FlomicsMultiAssay@metadata[["MOFA_results"]]@dimensions$K,
                                            value = 1, step = 1)),
                     column(3, radioButtons(inputId = session$ns("view_choice_heatmap"),
                                            label = "Data:",
                                            choices = views_names(FlomicsMultiAssay@metadata[["MOFA_results"]]),
                                            selected = views_names(FlomicsMultiAssay@metadata[["MOFA_results"]])[1]
                     )),
                     column(2, numericInput(inputId = session$ns("nfeat_choice_heatmap"),
                                            label = "Features:",
                                            min = 1,
                                            max = 500,
                                            value = 50, # default in MOFA function.
                     )),
                     column(1, radioButtons(inputId = session$ns("denoise_choice_heatmap"),
                                            label = "Denoise:",
                                            choices = c("TRUE", "FALSE"),
                                            selected = "FALSE")),
                     column(2, radioButtons(inputId = session$ns("annot_samples"),
                                            label = "Annotation:",
                                            choices = c('none', 
                                                        colnames(FlomicsMultiAssay@metadata$MOFA_results@samples_metadata %>% select(!c(sample, group)))),
                                            selected = "none"))
                     
                   ),
                   fluidRow(
                     # if(input$annot_samples == "none") input$annot_samples <- NULL
                     renderPlot({
                       annot_samples_values <- input$annot_samples
                       if(input$annot_samples == "none") annot_samples_values <- NULL
                       
                       observeEvent(input$view_choice_heatmap, {
                         updateSliderInput(inputId = "nfeat_choice_heatmap", max = get_dimensions(FlomicsMultiAssay@metadata$MOFA_results)$D[input$view_choice_heatmap][[1]])
                       })
                       
                       res_heatmap <- MOFA2::plot_data_heatmap(FlomicsMultiAssay@metadata[["MOFA_results"]], 
                                                               factor = input$factor_choice_heatmap, 
                                                               view = input$view_choice_heatmap, 
                                                               features = input$nfeat_choice_heatmap,
                                                               denoise = input$denoise_choice_heatmap,
                                                               annotation_samples = annot_samples_values)
                       grid::grid.draw(res_heatmap)
                     })
                   ) # fluidrow heatmap
          ) # tabpanel heatmap
        )# tabsetpanel 
    )#box
  })
}




