# ==============================================================================
# Common modules
# ==============================================================================

# ---- update radio button ----
UpdateRadioButtonsUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    radioButtons(inputId  = ns("Firstaxis"),
                 label    = "Choice of PCs :",
                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                 selected = 1, inline = TRUE),
    
    # select PCA axis 2 for plot
    radioButtons(inputId  = ns("Secondaxis"),
                 label    = "",
                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                 selected = 2, inline = TRUE)
  )
}


UpdateRadioButtons <- function(input, output, session){
  
  observeEvent(input$Firstaxis, {
    
    x <- input$Firstaxis
    # Can also set the label and select items
    choices <- c("PC1" = 1, "PC2" = 2, "PC3" = 3)
    updateRadioButtons(session, "Secondaxis",
                       choices = choices[-as.numeric(x)],
                       inline  = TRUE)
  })
  
}

RadioButtonsConditionUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    uiOutput(ns('condColor')),
  )
}

RadioButtonsCondition <- function(input, output, session, typeFact){
  
  # select factors for color PCA plot
  output$condColor <- renderUI({
    
    factors <- session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% typeFact]
    condition <- names(factors)
    
    if (! any(typeFact %in% "Meta")) 
      condition <- c("groups", condition)
    
    radioButtons(inputId = session$ns("condColorSelect"),
                 label = 'Levels:',
                 choices = condition,
                 selected = condition[1])
  })
}

# ---- summary of all analysed data ----
.modSingleOmicAnalysesSummaryUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    fluidPage(
      column(width = 12, 
             fluidRow( uiOutput(ns("overView"))),
             fluidRow( uiOutput(ns("DiffSummary"))),
             fluidRow( uiOutput(ns("CoExSummary")))
      )
    ))
}

#' @importFrom MOFA2 create_mofa_from_MultiAssayExperiment  plot_data_overview
#' @importFrom purrr reduce
#' @importFrom data.table data.table
.modSingleOmicAnalysesSummary <- function(input, output, session, rea.values){
  
  # over view of dataset dimensions after processing
  output$overView <- renderUI({
    
    if(is.null(rea.values$datasetProcess)) return()
    
    box(title = "Dataset overview after data processing", width = 12, status = "warning", 
        solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, 
        
        renderPlot({
          plotDataOverview(session$userData$FlomicsMultiAssay, omicNames = rea.values$datasetProcess, realSize = TRUE) 
        })
    )
  })
  
  # summary of diff analysis on all dataset
  output$DiffSummary <- renderUI({
    
    if(is.null(rea.values$datasetDiff)) return()
    
    box(title = "Summary of Differential expression analyses", width = 12, status = "warning", 
        solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, 
        # renderPlot({ 
        #   ggplot2::ggplot(data = summaryDiff.df) + 
        #     ggplot2::geom_col(ggplot2::aes(y=hypothesis, x=percent, label=Up_Down, fill=variable)) + 
        #     ggplot2::facet_grid(dataset~.) %>% dplyr::mutate(percent=Up_Down/All*100)
        # })
        tagList({
          tabPanel.list <- list(
            tabPanel(title = "diff results", 
                     renderPlot({
                       getDiffAnalysesSummary(session$userData$FlomicsMultiAssay, plot = TRUE)
                     })))
          
          p.list <- getAnnotAnalysesSummary(session$userData$FlomicsMultiAssay, 
                                            from="DiffExpEnrichAnal",
                                            matrixType = "presence")
          
          if(!is.null(rea.values$datasetDiffAnnot)){
            
            tabPanel.list <- 
              c(tabPanel.list,
                lapply(names(p.list), function(database){
                  
                  tabPanel(title = paste0("ORA results from ",database), 
                           fluidRow(
                             column(width = 12,
                                    radioButtons(inputId = session$ns(paste0(database, "-domain.diff")), 
                                                 label = "Domain",
                                                 choices = names(p.list[[database]]), selected = names(p.list[[database]])[1], 
                                                 inline = TRUE))
                           ),
                           fluidRow(
                             column(width = 12,
                                    renderPlot({ p.list[[database]][[input[[paste0(database, "-domain.diff")]]]] }) ))
                  )
                })
              )
          }
          
          do.call(what = tabsetPanel, args = tabPanel.list)
        })
    )
  })
  
  output$CoExSummary <- renderUI({
    
    if(is.null(rea.values$datasetCoEx)) return()
    
    box(title = "Summary of Co-expression analyses", width = 12, status = "warning", 
        solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, 
        
        tagList({
          tabPanel.list <- list(
            tabPanel(title = "CoExp results", 
                     renderPlot({
                       getCoExpAnalysesSummary(session$userData$FlomicsMultiAssay)
                     })
            ))
          
          p.list <- getAnnotAnalysesSummary(session$userData$FlomicsMultiAssay,
                                            from="CoExpEnrichAnal",
                                            matrixType = "presence")
          
          if(!is.null(rea.values$datasetCoExAnnot)){
            
            tabPanel.list <-
              c(tabPanel.list,
                lapply(names(p.list), function(database){
                  
                  tabPanel(title = paste0("ORA results from ",database),
                           fluidRow(
                             column(width = 12,
                                    radioButtons(inputId = session$ns(paste0(database, "-domain.coex")),
                                                 label = "Domain",
                                                 choices = names(p.list[[database]]), selected = names(p.list[[database]])[1],
                                                 inline = TRUE))
                           ),
                           fluidRow(
                             column(width = 12,
                                    renderPlot({ p.list[[database]][[input[[paste0(database, "-domain.coex")]]]] }) ))
                  )
                })
              )
          }
          
          do.call(what = tabsetPanel, args = tabPanel.list)
        })
    )
  })
  
}


