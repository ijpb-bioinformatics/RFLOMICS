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
                   numericInput(inputId = session$ns("numiter"), label="num factors :", value=10, min = 10, max=10))),
          
          fluidRow(
            column(12,
                   numericInput(inputId = session$ns("maxiter"), label="Max iteration :", value=1000, min = 1000, max=1000))),
          
          fluidRow(
            column(4, actionButton(session$ns("runMOFA"),"Run")))
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
    
    # run MOFA
    
    #!!!! lister les fonctions pour créer l'objet mofa
    
    # object.1 = dataPrepareForIntegration(object = object, omicsToIntegrate = c("RNAseq", "proteomics"), choice = "DE", contrast = c("H1", "H2"))
    # object.2 = MOFA_createObject(object.1, group = "temperature")
    # object.3 = MOFA_prepareObject(object.2, scale_views = TRUE, maxiter = 550, num_factors = 5) 
    
    print(input$selectedData)
    untrainedMOFA <- prepareMOFA(FlomicsMultiAssay,
                           omicsToIntegrate = input$selectedData,
                           rnaSeq_transfo = input$RNAseqTransfo,
                           choice = "DE",
                           contrast = input$selectedContrast, 
                           group = NULL)
    FlomicsMultiAssay@metadata[["MOFA_results"]] <<- run_MOFA_analysis(untrainedMOFA, 
                                                                     scale_views = input$scaleViews,
                                                                     maxiter = input$maxiter,
                                                                     num_factors = input$numiter) # numiter a changer !!!
    
    
  }, ignoreInit = TRUE)
  
  
  output$ResultViewUI <- renderUI({
    
    box(width=14, solidHeader = TRUE, status = "warning",
        title = "MOFA results",
        
        tabsetPanel(
          
          ###  ###
          tabPanel("view1", "? min de nbr de dataset à integrer?; "),
          #tabPanel("view1", MOFA2::plot_factor_cor(FlomicsMultiAssay@metadata$MOFA_results$MOFAobject)),
          ### 
          # tabPanel("view2", FlomicsMultiAssay@metadata$MOFA_results),
          
          ### 
          # tabPanel("view3", MOFA2::plot_factor_cor(FlomicsMultiAssay@metadata$MOFA_results)),
          
          ### 
          tabPanel("view2", "titi")
        )
    )
    
  })
}




