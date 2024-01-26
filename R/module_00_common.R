############# module 

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




###### summary of all analysed data ####
omics_data_analysis_summaryUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    column(width = 12, 
      fluidRow( uiOutput(ns("mofaPlot"))),
      fluidRow( uiOutput(ns("DiffSummary")))
    ))
}

#' @importFrom MOFA2 create_mofa_from_MultiAssayExperiment  plot_data_overview
#' @importFrom purrr reduce
#' @importFrom data.table data.table
omics_data_analysis_summary <- function(input, output, session, rea.values){
  
  # over view of dataset dimensions after processing
  output$mofaPlot <- renderUI({
    box(title = "Dataset overview after data processing", width = 12, status = "warning", 
        
        renderPlot({
          Datasets_overview_plot(session$userData$FlomicsMultiAssay, dataset.list = rea.values$datasetDiff, real.size = TRUE) 
        })
    )
  })
  
  # summary of diff analysis on all dataset
  output$DiffSummary <- renderUI({
    
    if(is.null(rea.values$datasetDiff)) return()
    
    summaryDiff.df <- lapply( rea.values$datasetDiff, function(dataset){
      as.data.frame(session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats) %>% 
        dplyr::mutate(dataset = dataset, hypothesis = rownames(.))
      
    }) %>% purrr::reduce(rbind) %>% reshape2::melt(id=c("dataset", "hypothesis", "All"), value.name = "Up_Down") %>%
      dplyr::mutate(percent=Up_Down/All*100)
    
    box(title = "Summary of Differential expression analysis", width = 12, status = "warning", 
        # renderPlot({ 
        #   ggplot2::ggplot(data = summaryDiff.df) + ggplot2::geom_col(ggplot2::aes(y=hypothesis, x=percent, label=Up_Down, fill=variable)) + 
        #     ggplot2::facet_grid(dataset~.) %>% dplyr::mutate(percent=Up_Down/All*100)
        # })
        
        renderPlot({
          ggplot2::ggplot(data = summaryDiff.df, ggplot2::aes(y=hypothesis, x=percent, fill=variable)) + 
            ggplot2::geom_col() +
            ggplot2::geom_text(ggplot2::aes(label=Up_Down), position = ggplot2::position_stack(vjust = 0.5)) +
            ggplot2::facet_grid(dataset~.) +
            ggplot2::scale_x_continuous(breaks = seq(0,100, 25), labels = paste0(seq(0,100, 25), "%")) + 
            ggplot2::labs(fill=NULL, x="")
        })
    )
  })
  
}


