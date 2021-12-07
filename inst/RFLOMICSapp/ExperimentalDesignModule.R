

GLM_modelUI <- function(id){

  ns <- NS(id)

  tagList(

    fluidRow(
      column(width= 12, uiOutput(ns("SetModelFormula")))
    ),
    fluidRow(
      column(width= 12, uiOutput(ns("SetContrasts")))
    ),
    fluidRow(
      column(width= 12, verbatimTextOutput(ns("printContrast")))
    ),
    tags$br()
  )
}


GLM_model <- function(input, output, session, rea.values){
  
  
    # reactive value for reinitialisation of UIoutput
    local.rea.values <- reactiveValues()
    
    observe({
      local.rea.values$contrast <- FALSE
    })
  
    # Construct the form to select the model
    output$SetModelFormula <- renderUI({
      
      #if (rea.values$model == FALSE) return()
      
      validate(
        need(rea.values$model != FALSE, "Please load data")
      )
      
      box(status = "warning", width = 12, solidHeader = TRUE, title = "Select a model formulae",

          tags$i("Models are written from the simplest one to the most complete one. Only the two orders
          interaction terms between the biological factors will appear. Batch factor will never appear
           in interaction terms."),

          selectInput( inputId = session$ns("model.formulae"), label = "",
                       choices = rev(names(GetModelFormulae(Factors.Name=names(session$userData$Design@List.Factors),
                                                            Factors.Type=session$userData$Design@Factors.Type))),
                       selectize=FALSE,size=5),
          actionButton(session$ns("validModelFormula"),"Valid model choice")
      )
    })

    # as soon as the "valid model formulae" button has been clicked
    # => The model formulae is set and the interface to select the contrasts appear
    observeEvent(input$validModelFormula, {
      
      local.rea.values$contrast <- FALSE
      rea.values$analysis       <- FALSE
      rea.values$Contrasts.Sel   <- NULL 

      FlomicsMultiAssay <<- resetFlomicsMultiAssay(object=FlomicsMultiAssay, results=c("DiffExpAnal", "CoExpAnal", "EnrichAnal"))
      
      print("# 3- Choice of statistical model...")

      # => Set the model formulae
      session$userData$Design@Model.formula <- input$model.formulae
      print(paste0("#    => ", input$model.formulae))

      # => get list of expression contrast (hypothesis)
      session$userData$Design <- getExpressionContrast(object = session$userData$Design, model.formula = input$model.formulae)
      FlomicsMultiAssay@metadata$design <<- session$userData$Design
      
      local.rea.values$contrast <- TRUE
      
    })
    

    # => Get and Display all the contrasts
    #  => The contrasts have to be choosen
    output$SetContrasts <- renderUI({
      
      if (rea.values$model == FALSE) return()
      if (local.rea.values$contrast == FALSE) return()

      box(width=12, status = "warning", solidHeader = TRUE, title = "Select contrasts",
          
          tags$i("blabla..."),
          br(),
          br(),
          
          column(width = 12, 
                 lapply(names(session$userData$Design@Contrasts.List), function(contrastType) {
                   
                   vect        <- as.vector(session$userData$Design@Contrasts.List[[contrastType]]$contrast)
                   names(vect) <- as.vector(session$userData$Design@Contrasts.List[[contrastType]]$contrastName)
                   
                   # pickerInput(
                   #   inputId  = session$ns(paste0("ContrastType",contrastType)),
                   #   label    = tags$span(style="color: black;", paste0(contrastType, "contrasts" )),
                   #   choices  = vect,
                   #   options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                   #   multiple = TRUE,
                   #   selected = NULL)
                   
                   box(
                     checkboxGroupInput(inputId = session$ns(paste0("ContrastType",contrastType)),
                                        label   = paste0("Contrast type : ", contrastType), choices = vect)
                   )
                 })
          ),
          br(),
          actionButton(session$ns("validContrasts"),"Valid contrast(s) choice(s)")
      )
    })   

    # as soon as the "valid Contrasts" buttom has been clicked
    # => The selected contrasts are saved
    # => The load data item appear
    observeEvent(input$validContrasts, {

      print("# 4- Choice of contrasts...")
      
      rea.values$analysis <- FALSE
      
      # reset analysis

      FlomicsMultiAssay <<- resetFlomicsMultiAssay(object=FlomicsMultiAssay, results=c("DiffExpAnal", "CoExpAnal", "EnrichAnal"))

      #rea.values$validate.status <- 0

      #get list of selected contrast data frames with expression, name and type
      contrastList <- list()
      contrastList <- lapply(names(session$userData$Design@Contrasts.List), function(contrastType) {

          input[[paste0("ContrastType",contrastType)]]
      })
      #names(contrastList) <- names(Design@Contrasts.List)
      contrast.sel.vec <- contrastList %>% unlist()

      # check if user has selected the contrasts to test
      if(length(contrast.sel.vec) == 0){

        showModal(modalDialog( title = "Error message", "Please select the hypotheses to test."))
        #rea.values$validate.status <- 1
      }

      ## continue only if message is true
      validate({
        need(length(contrast.sel.vec) != 0, message="ok")
        })

      print(paste0("#    => ", contrast.sel.vec))

      # define all the coefficients of selected contrasts and return a contrast matrix with contrast sample name and associated coefficients
      session$userData$Design <- getContrastMatrix(object = session$userData$Design, contrastList = contrast.sel.vec)
      FlomicsMultiAssay@metadata$design <<- session$userData$Design
      
      rea.values$Contrasts.Sel <- session$userData$Design@Contrasts.Sel
      
      rea.values$analysis <- TRUE
      
      Design <<- session$userData$Design
      
      # # à supprimer à la fin du dev
      # output$printContrast <- renderPrint({
      # 
      #   #Design@Contrasts.Coeff
      #   A <- session$userData$Design@Contrasts.Coeff
      #   row.names(A) <- NULL
      #   print(A)
      # })
    }, ignoreInit = TRUE)

    return(input)
  }

