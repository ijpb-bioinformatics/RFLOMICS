

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


GLM_model <- function(input, output, session){

      # Construct the form to select the model
      output$SetModelFormula <- renderUI({
        box(status = "warning", width = 12,
            h4("Select a model formulae"),
            hr(),
            tags$i("Models are written from the simplest one to the most complete one. Only the two orders
            interaction terms between the biological factors will appear. Batch factor will never appear
             in interaction terms."),
            hr(),
            selectInput( inputId = session$ns("model.formulae"), label = "",
                         choices = rev(names(GetModelFormulae(Factors.Name=names(Design@List.Factors),
                                                              Factors.Type=Design@Factors.Type))),
                         selectize=FALSE,size=5),
            actionButton(session$ns("validModelFormula"),"Valid model choice")
        )
      })

      # as soon as the "valid model formulae" button has been clicked
      # => The model formulae is set and the interface to select the contrasts appear
      observeEvent(input$validModelFormula, {

        print("# 3- Choice of statistical model...")

        # => Set the model formulae
        #Design@Model.formula <<- input$model.formulae
        print(paste0("#    => ", input$model.formulae))

        # => get list of expression contrast (hypothesis)
        Design <<- getExpressionContrast(object = Design, model.formula = input$model.formulae)
        FlomicsMultiAssay@metadata$design <<- Design

        # => Set Model design matrix
        # => Get and Display all the contrasts


        #  => The contrasts have to be choosen
        output$SetContrasts <- renderUI({

          box(width=12, status = "warning", size=3,

            lapply(names(Design@Contrasts.List), function(contrastType) {

                  vect        <- as.vector(Design@Contrasts.List[[contrastType]]$contrast)
                  names(vect) <- as.vector(Design@Contrasts.List[[contrastType]]$contrastName)

                  box(
                    checkboxGroupInput(inputId = session$ns(paste0("ContrastType",contrastType)),
                                       label   = paste0("Contrast type : ", contrastType), choices = vect)
                  )
            }),
            column(width=4, actionButton(session$ns("validContrasts"),"Valid contrast(s) choice(s)"))
          )
        })
      })


      # as soon as the "valid Contrasts" buttom has been clicked
      # => The selected contrasts are saved
      # => The load data item appear
      observeEvent(input$validContrasts, {

        validate.status <<- 0

        #get list of selected contrast data frames with expression, name and type
        contrastList <- list()
        contrastList <- lapply(names(Design@Contrasts.List), function(contrastType) {

            input[[paste0("ContrastType",contrastType)]]
        })
        #names(contrastList) <- names(Design@Contrasts.List)
        contrast.sel.vec <- contrastList %>% unlist()

        # check if user has selected the contrasts to test
        if(length(contrast.sel.vec) == 0){

          showModal(modalDialog( title = "Error message", "Please select the hypotheses to test."))
          validate.status <<- 1
        }

        ## continue only if message is true
        validate({
          need(length(contrast.sel.vec) != 0, message="ok")
          })


        # define all the coefficients of selected contrasts and return a contrast matrix with contrast sample name and associated coefficients
        Design <<- getContrastMatrix(Design, contrastList = contrast.sel.vec)
        FlomicsMultiAssay@metadata$design <<- Design
        
        # à supprimer à la fin du dev
        output$printContrast <- renderPrint({

          #Design@Contrasts.Coeff
          A <- Design@Contrasts.Coeff
          row.names(A) <- NULL
          print(A)
        })
      })
      
      


    return(input)
  }

