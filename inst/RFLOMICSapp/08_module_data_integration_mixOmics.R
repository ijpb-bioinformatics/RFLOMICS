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
            h4(tags$span("Parameters set up:", style = "color:orange")),
            p("This is where you have to put the parameters"),
            
            h4(tags$span("Outputs:", style = "color:orange")),
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
             box(title = span(tagList(icon("sliders-h"), "  ", "Setting Blocks")), width = 12, status = "warning",
                 
                 # Select lists of dataset to integrate
                 fluidRow(
                   column(12,
                          
                          pickerInput(
                            inputId  = session$ns("selectedData"),
                            label    = "Select dataset:",
                            choices  = rea.values$datasetDiff,
                            options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            multiple = TRUE,
                            selected = rea.values$datasetDiff))),
                 
                 fluidRow(
                   column(12,
                          
                          pickerInput(
                            inputId  = session$ns("selectedContrast"),
                            label    = "Select contrasts:",
                            choices  = listOfContrast,
                            options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            multiple = TRUE,
                            selected = listOfContrast))),
                 
                 # select mode of feature filtering
                 fluidRow(
                   column(12,
                          
                          radioButtons(inputId  = session$ns("filtMode"), 
                                       label    = "Select type of filtering:" ,
                                       choices  = c("union","intersection"),
                                       selected = "union", inline = FALSE))),
                 
                 # set parameters
                 fluidRow(
                   column(12,
                          selectInput(inputId  = session$ns("RNAseqTransfo"),
                                      label    = "RNAseq transfo :",
                                      choices  = c("limma (voom)"),
                                      selected = "limma (voom)")))
                 
             ),
             box(title = span(tagList(icon("sliders-h"), "  ", "Settings for analysis")), width = 12, status = "warning",
                 
                 fluidRow(
                   column(12,
                          column(5, checkboxInput(inputId = session$ns("scale_views"), label = "Scale Datasets", value = FALSE, width = NULL)),
                          column(5, checkboxInput(inputId = session$ns("sparsity"), label = "Sparse analysis", value = FALSE, width = NULL)),
                          numericInput(inputId = session$ns("ncomp"), label = "Number of component", value = 2, min = 1, max= 5),
                          numericInput(inputId = session$ns("cases_to_try"), label = "Tuning cases", value = 2, min = 1, max= 5),
                          numericInput(inputId = session$ns("link_datasets"), label = "Link between datasets", value = 1, min = 0, max= 1),
                          numericInput(inputId = session$ns("link_response"), label = "Link to response", value = 1, min = 0, max= 1)
                   ))
             ),
             box(title = span(tagList(icon("sliders-h"), "  ", "Response Variables")), width = 12, status = "warning",
                 
                 # Select lists of dataset to integrate
                 fluidRow(
                   # box(title = span(tagList(icon("sliders-h"), "  ", "Options")), width = 12, collapsible = TRUE, collapsed = TRUE,
                   #     div(      
                   #       h4(tags$span("Parameters set up:", style = "color:orange")),
                   #       p("Selecting one feature will set the analysis to discriminant analysis"),
                   #       p("Selecting more than one will set the analysis to non-discriminant: qualitative features will be turned into dummy variables and
                   #               the analysis is conducted on the result quantitative response matrix")
                   #     )
                   # ),
                   column(12,
                          checkboxGroupInput(
                            inputId  = session$ns("selectedResponse"),
                            label    = "Select response variables:",
                            choices  = c(colnames(colData(session$userData$FlomicsMultiAssay))),
                            # options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            selected = colnames(colData(session$userData$FlomicsMultiAssay))[1]))
                 ),
                 box(title = span(tagList(icon("sliders-h"), "  ", "Run Analysis")), width = 12, status = "warning",
                     fluidRow(
                       column(12, actionButton(session$ns("runMixOmics"),"Run Analysis")))
                 )
             ),
      ) # column 3
    )# taglist
    
  })
  
  ## observe the button run mixOmics
  observeEvent(input$runMixOmics, {
    
    library(mixOmics)
    
    local.rea.values$runMixOmics   <- FALSE
    local.rea.values$preparedMixOmics  <- NULL
    local.rea.values$resMixOmics <- NULL
    local.rea.values$dis_anal <- FALSE
    local.rea.values$Y <- NULL
    local.rea.values$Design_mat <- NULL
    local.rea.values$keepX <- NULL
    local.rea.values$tuning_res <- NULL
    # local.rea.values$multipleBlocks <- FALSE # TODO use ?
    # local.rea.values$tuning <- FALSE # TODO use ?
    
    
    ### TESTS TO DELETE
    # load("inst/ExamplesFiles/FlomicsMultiAssay.RData")
    # preparedMixOmics <- prepareForIntegration(FlomicsMultiAssay,
    #                                           omicsToIntegrate = c("proteomics.set1", "metabolomics.set2"),
    #                                           rnaSeq_transfo = "",
    #                                           choice = "DE",
    #                                           contrasts_names = c("(temperatureLow - temperatureElevated) in mean",
    #                                                               "(temperatureMedium - temperatureElevated) in mean"),
    #                                           type = "union",
    #                                           group = NULL,
    #                                           method = "MixOmics")
    
    ## choice of function according to user selection
    
    # selectedResponse = c("temperature", "imbibition") # TODO DELETE
    # Y = as.data.frame(preparedMixOmics$metadata) %>% dplyr::select(all_of(selectedResponse))  # TODO DELETE
    
    # local.rea.values <- list()
    # input <- list()
    # input$sparsity = TRUE
    # input$Tuning = TRUE
    # input$cases_to_try = 3
    # input$ncomp = 2
    # input$scale_views = TRUE
    # input$link_datasets=1
    # input$link_response = 1
    # input$selectedResponse = c("temperature", "imbibition")
    # # input$selectedResponse = c("temperature")
    # local.rea.values$preparedMixOmics <- preparedMixOmics
    # local.rea.values$Y = Y
    
    session$userData$FlomicsMultiAssay@metadata[["MixOmics_results"]] <<- NULL
    session$userData$FlomicsMultiAssay@metadata[["MixOmics_tuning_results"]] <<- NULL
    
    # Prepare for MixOmics run  
    local.rea.values$preparedMixOmics <- prepareForIntegration(session$userData$FlomicsMultiAssay,
                                                               omicsToIntegrate = input$selectedData,
                                                               rnaSeq_transfo = input$RNAseqTransfo,
                                                               choice = "DE", 
                                                               contrasts_names = input$selectedContrast,
                                                               type = input$filtMode,
                                                               group = NULL,
                                                               method = "MixOmics")
    
    
    # Run the analysis
    local.rea.values$MixOmics_res <- run_MixOmics_analysis(local.rea.values$preparedMixOmics,
                                                           scale_views = input$scale_views,
                                                           selectedResponse = input$selectedResponse,
                                                           ncomp = input$ncomp,
                                                           link_datasets = input$link_datasets,
                                                           link_resposne = input$link_response,
                                                           sparsity = input$sparsity,
                                                           cases_to_try = input$cases_to_try
    )
    
    # Store results
    session$userData$FlomicsMultiAssay@metadata[["MixOmics_tuning_results"]] <<- local.rea.values$MixOmics_res$tuning_res
    session$userData$FlomicsMultiAssay@metadata[["MixOmics_results"]] <<- local.rea.values$MixOmics_res$analysis_res
    
    local.rea.values$runMixOmics <- TRUE
    
  })
  
  ## Output results
  output$MOResultViewUI <- renderUI({
    
    if(local.rea.values$runMixOmics == FALSE) return()
    
    box(width=14, solidHeader = TRUE, status = "warning",
        title = "MixOmics results",
        
        tabsetPanel(
          # ---- Tab panel Overview ----
          tabPanel("Overview",
                   column(6 , DT::renderDataTable({
                     
                     df <- t(sapply(local.rea.values$MixOmics_res$analysis_res$X, dim))
                     colnames(df) <- c("Ind", "Features")
                     
                     if(input$sparsity){
                       df <- cbind(df, do.call("rbind", local.rea.values$MixOmics_res$analysis_res$keepX))
                       colnames(df)[!colnames(df) %in% c("Ind", "Features")] <- paste("Comp", 1:length(local.rea.values$MixOmics_res$analysis_res$keepX[[1]]))
                     }
                     
                     t(df) %>% DT::datatable()
                     
                   })),
                   
                   
                   
          ),
          # ---- Tab panel Explained Variance ----
          tabPanel("Explained Variance",
                   column(6, renderPlot({
                     dat_explained <- reshape2::melt(do.call("rbind", local.rea.values$MixOmics_res$analysis_res$prop_expl_var))
                     colnames(dat_explained) <- c("Dataset", "Component", "Percentage of explained variance")
                     dat_explained$`Percentage of explained variance` <- dat_explained$`Percentage of explained variance`*100
                     
                     print(dat_explained)
                     
                     dat_comb <- dat_explained %>% 
                       dplyr::group_by(Dataset) %>% 
                       dplyr::summarise("Cumulative Explained Variance" = sum(`Percentage of explained variance`))
                     
                     print(dat_comb)
                     
                     ggplot2::ggplot(dat_comb, aes(x = Dataset, y = `Cumulative Explained Variance`)) +
                       geom_col() + 
                       theme_classic() +
                       theme(
                         axis.text = element_text(size = 12),
                         axis.line = element_blank(),
                         axis.ticks =  element_blank(),
                         strip.text = element_text(size = 12),
                       ) + ylab("")   
                     
                     
                   })),
                   
                   column(6 , renderPlot({
                     # local.rea.values$MixOmics_res$analysis_res = TCGA.block.splsda
                     dat_explained <- reshape2::melt(do.call("rbind", local.rea.values$MixOmics_res$analysis_res$prop_expl_var))
                     colnames(dat_explained) <- c("Dataset", "Component", "Percentage of explained variance")
                     dat_explained$`Percentage of explained variance` <- dat_explained$`Percentage of explained variance`*100
                     
                     # Chunk of code to be cohesive with MOFA2::plot_explained_variance
                     ggplot2::ggplot(dat_explained, aes(x = Dataset, y = Component)) +
                       geom_tile(aes(fill = `Percentage of explained variance`)) + 
                       theme_classic() +
                       theme(
                         axis.text = element_text(size = 12),
                         axis.line = element_blank(),
                         axis.ticks =  element_blank(),
                         strip.text = element_text(size = 12),
                       ) + ylab("") +  
                       scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min(dat_explained$`Percentage of explained variance`),
                                                                                                      max(dat_explained$`Percentage of explained variance`)))
                     
                   })),
                   
          ),
          # ---- Tab panel Individuals ----
          tabPanel("Individuals",
                   # TODO Rep.space ? Ne pas laisser le choix, je pense ....
                   # TODO Le titre ?
                   column(3,
                          checkboxInput(inputId = session$ns("ellipse_choice"), label = "Ellipses", value = FALSE, width = NULL),
                          numericInput(inputId = session$ns("ind_comp_choice_1"),
                                       label = "Comp x:",
                                       min = 1,
                                       max =  input$ncomp,
                                       value = 1, step = 1),
                          numericInput(inputId = session$ns("ind_comp_choice_2"),
                                       label = "Comp y:",
                                       min = 1,
                                       max =  input$ncomp,
                                       value = 2, step = 1)
                   ),
                   column(9 , renderPlot(mixOmics::plotIndiv(local.rea.values$MixOmics_res$analysis_res, 
                                                             comp = c(input$ind_comp_choice_1, input$ind_comp_choice_2),
                                                             ellipse = input$ellipse_choice,
                                                             legend = TRUE)))
          ),
          # ---- Tab panel Features ----
          tabPanel("Features",
                   # TODO Permettre la selection du block affiche ? Bof ?   
                   column(3,
                          checkboxInput(inputId = session$ns("overlap"), label = "Overlap", value = FALSE, width = NULL),
                          numericInput(inputId = session$ns("var_comp_choice_1"),
                                       label = "Comp x:",
                                       min = 1,
                                       max =  input$ncomp,
                                       value = 1, step = 1),
                          numericInput(inputId = session$ns("var_comp_choice_2"),
                                       label = "Comp y:",
                                       min = 1,
                                       max =  input$ncomp,
                                       value = 2, step = 1)
                   ),
                   column(9 , renderPlot(mixOmics::plotVar(local.rea.values$MixOmics_res$analysis_res, 
                                                           comp = c(input$var_comp_choice_1, input$var_comp_choice_2),
                                                           overlap = input$overlap,
                                                           legend = TRUE)))
          ),     
          # ---- Tab panel Loadings ----
          tabPanel("Loadings",
                   # TODO Permettre la selection du block affiche ?
                   column(3,
                          numericInput(inputId = session$ns("Load_comp_choice"),
                                       label = "Component:",
                                       min = 1,
                                       max =  input$ncomp,
                                       value = 1, step = 1),
                          numericInput(inputId = session$ns("Load_ndisplay"),
                                       label = "Number of features to display:",
                                       min = 1,
                                       max =  max(sapply(local.rea.values$MixOmics_res$analysis_res$X, ncol)),
                                       value = 25, step = 1),
                   ),
                   column(9 , renderPlot(mixOmics::plotLoadings(local.rea.values$MixOmics_res$analysis_res, 
                                                                comp = input$Load_comp_choice,
                                                                ndisplay = input$Load_ndisplay)))
          ),   
          # ---- Tab panel Tuning ----
          #   tabPanel("Tuning",
          #            # TODO 
          #            column(12 , renderPlot(mixOmics::plot(local.rea.values$tuning_res)))
          #   )), # plot.tune ne fait pas partie du package dans ma version ?!
          tabPanel("Networks",
                   # TODO Prevoir bouton pour comp selectionnee
                   # Comp fonctionne plus ?!
                   column(3, numericInput(inputId = session$ns("Network_cutoff"),
                                          label = "Cutoff:",
                                          min = 0,
                                          max =  1,
                                          value = 0.75, step = 0.05)),
                   column(9 , renderPlot(mixOmics::network(mat = local.rea.values$MixOmics_res$analysis_res, 
                                                           # comp = 1:2, 
                                                           blocks = 1:length(input$selectedData),
                                                           cutoff = input$Network_cutoff, 
                                                           shape.node = c("rectangle", "rectangle"))))
          ), 
          # ---- Tab  Panel CircosPlot & cimPlot ----
          # conditionalPanel(condition = "is(local.real.values$MixOmics_res, 'block.splsda')",
          tabPanel("CircosPlot",
                   column(3, numericInput(inputId = session$ns("Circos_cutoff"),
                                          label = "Cutoff:",
                                          min = 0,
                                          max =  1,
                                          value = 0.75, step = 0.05)),
                   column(9, renderPlot(mixOmics::circosPlot(local.rea.values$MixOmics_res$analysis_res,
                                                             cutoff = input$Circos_cutoff))),
          ),
          tabPanel("CimPlot",
                   column(3,),
                   column(9, renderPlot(mixOmics::cimDiablo(local.rea.values$MixOmics_res$analysis_res))),
          ),
          # )
        ))
    
  })
  
}