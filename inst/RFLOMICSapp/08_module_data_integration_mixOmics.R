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
             box(title = span(tagList(icon("sliders"), "  ", "Setting Blocks")), width = 12, status = "warning",
                 
                 # Select lists of dataset to integrate
                 fluidRow(
                   column(12,
                          
                          pickerInput(
                            inputId  = session$ns("MOselectedData"),
                            label    = "Select dataset:",
                            choices  = rea.values$datasetDiff,
                            options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            multiple = TRUE,
                            selected = rea.values$datasetDiff))),
                 
                 fluidRow(
                   column(12,
                          
                          pickerInput(
                            inputId  = session$ns("MOselectedContrast"),
                            label    = "Select contrasts:",
                            choices  = listOfContrast,
                            options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                            multiple = TRUE,
                            selected = listOfContrast))),
                 
                 # select mode of feature filtering
                 fluidRow(
                   column(4,
                          radioButtons(inputId  = session$ns("MOfiltMode"), 
                                       label    = "Select type of filtering:" ,
                                       choices  = c("union", "intersection"),
                                       selected = "union", inline = FALSE)),
                   column(8,
                          verbatimTextOutput(session$ns("mergeValue_mixOmics")))
                 ),
                 # set parameters
                 fluidRow(
                   column(12,
                          selectInput(inputId  = session$ns("MORNAseqTransfo"),
                                      label    = "RNAseq transfo :",
                                      choices  = c("limma (voom)"),
                                      selected = "limma (voom)")))
                 
             ),
             box(title = span(tagList(icon("sliders"), "  ", "Settings for analysis")), width = 12, status = "warning",
                 
                 fluidRow(
                   column(12,
                          column(6, checkboxInput(inputId = session$ns("scale_views"), label = "Scale Datasets", value = FALSE, width = NULL)),
                          column(6, checkboxInput(inputId = session$ns("sparsity"), label = "Sparse analysis", value = FALSE, width = NULL)),
                          column(6, numericInput(inputId = session$ns("ncomp"), label = "Components", value = 2, min = 1, max= 5)),
                          column(6, numericInput(inputId = session$ns("cases_to_try"), label = "Tuning cases", value = 2, min = 1, max= 5)),
                          # column(6, numericInput(inputId = session$ns("link_datasets"), label = "Link between datasets", value = 1, min = 0, max= 1)),
                          # column(6, numericInput(inputId = session$ns("link_response"), label = "Link to response", value = 1, min = 0, max= 1))
                   ))
             ),
             box(title = span(tagList(icon("sliders"), "  ", "Response Variables")), width = 12, status = "warning",
                 
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
                 fluidRow(
                   column(12, actionButton(session$ns("runMixOmics"), "Run Analysis"))
                 )
             ),
      ) # column 3
    )# taglist
    
  })
  
  
  # TODO delete
  # input <- list()
  # input$MOselectedData <- c("metabolomics.set2.filtred", "proteomics.set1.filtred")
  # input$MOselectedContrast <- c("(imbibitionEI - imbibitionDS) in mean",            
  #                             "(imbibitionLI - imbibitionDS) in mean",            
  #                             "(imbibitionLI - imbibitionEI) in mean")
  
  # Entities list
  Ent_lists <- 
    lapply(input$MOselectedData, FUN = function(namData){
      lapply(session$userData$FlomicsMultiAssay@ExperimentList[[namData]]@metadata$DiffExpAnal$TopDEF[input$MOselectedContrast], 
             rownames)
    })
  
  # display number of selected entities
  output$mergeValue_mixOmics <- renderText({
    if(input$MOfiltMode == "union"){
      print(sum(sapply(Ent_lists, FUN = function(listDE){length(unique(unlist(listDE)))})))
    }else if(input$MOfiltMode == "intersection"){
      print(sum(sapply(Ent_lists, FUN = function(listDE){
        length(Reduce('intersect', listDE))
      })))
    }
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
    local.rea.values$preparedMixOmics  <- NULL
    local.rea.values$resMixOmics <- NULL
    session$userData$FlomicsMultiAssay@metadata[["mixOmics"]][["MixOmics_tuning_results"]] <- NULL
    session$userData$FlomicsMultiAssay@metadata[["mixOmics"]][["MixOmics_results"]] <- NULL
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Checks ", 10, "%", sep=""))
    #----------------------#
    
    # TODO missing checks in here !!
    
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Preparing object ", 20, "%", sep = ""))
    #----------------------#
    
    # Prepare for MixOmics run  
    print("#     =>Preparing data list")
    local.rea.values$preparedMixOmics <- prepareForIntegration(session$userData$FlomicsMultiAssay,
                                                               omicsToIntegrate = input$MOselectedData,
                                                               rnaSeq_transfo = input$MORNAseqTransfo,
                                                               choice = "DE", 
                                                               contrasts_names = input$MOselectedContrast,
                                                               type = input$MOfiltMode,
                                                               group = NULL,
                                                               method = "MixOmics")
    #---- progress bar ----#
    progress$inc(1/10, detail = paste("Running MixOmics ", 30, "%", sep = ""))
    #----------------------#
    
    # Run the analysis
    print("#     =>Running MixOmics")
    local.rea.values$MixOmics_res <- lapply(input$selectedResponse, 
                                            FUN = function(response_var){
                                              res_mixOmics <- run_MixOmics_analysis(local.rea.values$preparedMixOmics,
                                                                                    scale_views = input$scale_views,
                                                                                    selectedResponse = response_var,
                                                                                    ncomp = input$ncomp,
                                                                                    # link_datasets = input$link_datasets,
                                                                                    # link_resposne = input$link_response,      
                                                                                    link_datasets = 0.5,
                                                                                    link_resposne = 0.5,
                                                                                    sparsity = input$sparsity,
                                                                                    cases_to_try = input$cases_to_try)
                                              
                                              
                                              return(
                                                list(
                                                  "MixOmics_tuning_results" = res_mixOmics$tuning_res,
                                                  "MixOmics_results" = res_mixOmics$analysis_res
                                                )
                                              )
                                            })
    names(local.rea.values$MixOmics_res) <- input$selectedResponse
    
    # Store results
    # session$userData$FlomicsMultiAssay@metadata[["mixOmics"]][["MixOmics_tuning_results"]] <- local.rea.values$MixOmics_res$tuning_res
    # session$userData$FlomicsMultiAssay@metadata[["mixOmics"]][["MixOmics_results"]] <- local.rea.values$MixOmics_res$analysis_res
    session$userData$FlomicsMultiAssay@metadata[["mixOmics"]] <-  local.rea.values$MixOmics_res 

    local.rea.values$runMixOmics <- TRUE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Finished ", 100,"%", sep = ""))
    #----------------------#
  })
  
  ## Output results
  output$MOResultViewUI <- renderUI({
    
    if(!local.rea.values$runMixOmics) return()
    
    lapply(names(session$userData$FlomicsMultiAssay@metadata[["mixOmics"]]), function(listname) { # ADD
      
      Data_res <- session$userData$FlomicsMultiAssay@metadata[["mixOmics"]][[listname]]$MixOmics_results
      
      fluidRow(
        
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success", title = listname,
            
            tabsetPanel(
              # ---- Tab panel Overview ----
              tabPanel("Overview",
                       column(6 , DT::renderDataTable({
                         
                         df <- t(sapply(Data_res$X, dim))
                         colnames(df) <- c("Ind", "Features")
                         
                         if(input$sparsity){
                           df <- cbind(df, do.call("rbind", Data_res$keepX))
                           colnames(df)[!colnames(df) %in% c("Ind", "Features")] <- paste("Comp", 1:length(Data_res$keepX[[1]]))
                         }
                         
                         t(df) %>% DT::datatable()
                         
                       })),
                       
                       
                       
              ),
              # ---- Tab panel Explained Variance ----
              tabPanel("Explained Variance",
                       column(6, renderPlot({
                         dat_explained <- reshape2::melt(do.call("rbind", Data_res$prop_expl_var))
                         colnames(dat_explained) <- c("Dataset", "Component", "% of explained variance")
                         dat_explained$`% of explained variance` <- dat_explained$`% of explained variance`*100
                         
                         # print(dat_explained) # TODO delete
                         
                         dat_comb <- dat_explained %>% 
                           dplyr::group_by(Dataset) %>% 
                           dplyr::summarise("Cumulative Explained Variance" = sum(`% of explained variance`))
                         
                         # print(dat_comb) # TODO delete
                         
                         if(is(Data_res, "block.splsda") || is(Data_res, "block.plsda")){
                           dat_comb <- dat_comb %>% dplyr::filter(Dataset!="Y")
                         }
                         
                         ggplot2::ggplot(dat_comb, aes(x = Dataset, y = `Cumulative Explained Variance`)) +
                           geom_col(fill = "darkblue") + 
                           theme_classic() +
                           theme(
                             axis.text = element_text(size = 12),
                             axis.line = element_blank(),
                             axis.ticks =  element_blank(),
                             strip.text = element_text(size = 12),
                           ) + ylab("") + ggtitle("Cumulative explained variance")  
                         
                         
                       })),
                       
                       column(6 , renderPlot({
                         dat_explained <- reshape2::melt(do.call("rbind", Data_res$prop_expl_var))
                         colnames(dat_explained) <- c("Dataset", "Component", "% of explained variance")
                         dat_explained$`% of explained variance` <- dat_explained$`% of explained variance`*100
                         
                         # print(dat_explained) # TODO delete
                         
                         if(is(Data_res, "block.splsda") || is(Data_res, "block.plsda")){
                           dat_explained <- dat_explained %>% dplyr::filter(Dataset!="Y")
                         }
                         
                         # Chunk of code to be cohesive with MOFA2::plot_explained_variance
                         ggplot2::ggplot(dat_explained, aes(x = Dataset, y = Component)) +
                           geom_tile(aes(fill = `% of explained variance`)) + 
                           theme_classic() +
                           theme(
                             axis.text = element_text(size = 12),
                             axis.line = element_blank(),
                             axis.ticks =  element_blank(),
                             strip.text = element_text(size = 12),
                           ) + ylab("") +  
                           scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min(dat_explained$`% of explained variance`),
                                                                                                          max(dat_explained$`% of explained variance`))) + 
                           ggtitle("Percentage of explained variance \n per component per block")
                         
                       })),
                       
              ),
              # ---- Tab panel Individuals ----
              tabPanel("Individuals",
                       column(1,
                              checkboxInput(inputId = session$ns(paste0(listname, "ellipse_choice")), label = "Ellipses", value = FALSE, width = NULL),
                              numericInput(inputId = session$ns(paste0(listname, "ind_comp_choice_1")),
                                           label = "Comp x:",
                                           min = 1,
                                           max =  input$ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(listname, "ind_comp_choice_2")),
                                           label = "Comp y:",
                                           min = 1,
                                           max =  input$ncomp,
                                           value = 2, step = 1)
                       ),
                       column(11 , renderPlot(mixOmics::plotIndiv(Data_res, 
                                                                  comp = c(input[[paste0(listname, "ind_comp_choice_1")]], input[[paste0(listname, "ind_comp_choice_2")]]),
                                                                  ellipse = input[[paste0(listname, "ellipse_choice")]],
                                                                  legend = TRUE)))
              ),
              # ---- Tab panel Features ----
              tabPanel("Features",
                       column(1,
                              checkboxInput(inputId = session$ns("overlap"), label = "Overlap", value = FALSE, width = NULL),
                              numericInput(inputId = session$ns(paste0(listname, "var_comp_choice_1")),
                                           label = "Comp x:",
                                           min = 1,
                                           max =  input$ncomp,
                                           value = 1, step = 1),
                              numericInput(inputId = session$ns(paste0(listname, "var_comp_choice_2")),
                                           label = "Comp y:",
                                           min = 1,
                                           max =  input$ncomp,
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
                                           max =  input$ncomp,
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
                                                                blocks = 1:length(input$MOselectedData),
                                                                cutoff = input[[paste0(listname, "Network_cutoff")]], 
                                                                shape.node = rep("rectangle", length(input$MOselectedData)))))
              ), 
              # ---- Tab  Panel CircosPlot & cimPlot ----
              tabPanel("CircosPlot",
                       if(is(Data_res, "block.splsda")){
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
                       column(1,),
                       column(11, 
                              if(is(Data_res, "block.splsda")){
                                print("=> Rendering cimPlot, be patient!")
                                renderPlot(mixOmics::cimDiablo(Data_res))
                              }else{
                                renderText({print("This plot is only available for sparse multi-block discriminant analysis results.")})
                              }
                       ),
              ),
            ) # tabsetpanel
        ) #box
      ) # fluidrow
    }) # lapply
  }) #renderui
  
}