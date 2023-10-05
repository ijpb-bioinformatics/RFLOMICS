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
    choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
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
                 label = 'Levels :',
                 choices = condition,
                 selected = condition[1])
  })
}




###### summary of all analysed data ####
omics_data_analysis_summaryUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box(title = "Summary of dataset analysis", width = 12, status = "warning", solidHeader = TRUE)),
    fluidRow(
      box(title = "Dataset processing", width = 12, status = "warning", 
          column(4, plotOutput(ns("procSummary"))),
          column(8, plotOutput(ns("mofaPlot"))))),
    fluidRow(
      box(title = "Differential analysiss", width = 12, status = "warning", plotOutput(ns("DiffSummary")))),
  )
}

#' @importFrom MOFA2 create_mofa_from_MultiAssayExperiment  plot_data_overview
#' @importFrom purrr reduce
#' @importFrom data.table data.table
omics_data_analysis_summary <- function(input, output, session, rea.values){

  output$procSummary <- renderPlot({ 
    
    summaryProcess.df <- lapply(rea.values$datasetDiff, function(dataset){
      
      c(paste(dataset, "filtred", sep = "."), 
        dim(session$userData$FlomicsMultiAssay[[dataset]]) - dim(session$userData$FlomicsMultiAssay[[paste(dataset, "filtred", sep = ".")]]))
    }) %>% purrr::reduce(rbind) %>% data.table::data.table()
    names(summaryProcess.df) <- c("dataset", "rm_entities", "rm_samples")
    
    ggplot2::ggplot(summaryProcess.df %>% reshape2::melt(id="dataset"), ggplot2::aes(x=variable, y=dataset)) + 
      ggplot2::geom_tile(aes(fill=dataset),color="white", lwd =1.5, linetype=1) + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
            panel.background = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position="none") + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::geom_text(ggplot2::aes(label=value)) 
    })
  
  
  output$mofaPlot <- renderPlot({ 

    MAE.filtred.Diff <- session$userData$FlomicsMultiAssay[,, paste(rea.values$datasetDiff, "filtred", sep = ".")]
    MAE.MOFA <- MOFA2::create_mofa_from_MultiAssayExperiment(MAE.filtred.Diff)
    MOFA2::plot_data_overview(MAE.MOFA)
  })
  
  output$DiffSummary <- renderPlot({ 
  
    
    summaryDiff.df <- lapply( paste(rea.values$datasetDiff, "filtred", sep = "."), function(dataset){
      lapply(names(session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats), function(contrast){  
        
        c(dataset, contrast, paste0(session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats[[contrast]][["gDEup"]], "/", 
                                    session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats[[contrast]][["gDEdown"]]) )
        
      }) %>% purrr::reduce(rbind)
      
    }) %>% purrr::reduce(rbind) %>% data.table::data.table()
    
    names(summaryDiff.df) <- c("dataset", "contrast", "value")
    
    
    ggplot2::ggplot(summaryDiff.df, ggplot2::aes(x=contrast, y=dataset)) + 
      ggplot2::geom_tile(ggplot2::aes(fill=dataset),color="white", lwd =1.5, linetype=1) + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
            panel.background = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position="none",
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + 
      ggplot2::ylab("") + ggplot2::xlab("") + ggplot2::geom_text(aes(label=value)) + ggplot2::ggtitle("DE entities (up/down)")
  })
  
}


