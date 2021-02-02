ExperimentalDesignUI <- function(id){
  
  ns <- NS(id)
  
  tagList(  
    ### import Design file
    fluidRow(
      box(width = 9, status = "warning",
          
          column(width = 8,
                 # matrix count/abundance input
                 fileInput(inputId = ns("Experimental.Design.file"), label = "Import matrix of Experimental Design (txt)",
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                 actionButton(ns("loadExpDesign"),"load")
          )
      )
    ),
    tags$br(),
    
    ### table visualisation
    fluidRow(
      uiOutput(ns("ExpDesignTable"))
    ),
    tags$br(),
    
    ### level
    fluidRow(
      uiOutput(ns("GetdFactorRef"))
    ),
    tags$br(),
    
    ### completeness
    fluidRow(
      uiOutput(ns("Completeness"))
    )
  )
  
}

ExperimentalDesign <- function(input, output, session){
  
  ##########################################
  # Part2 : Define the Experimental design:
  #         -> load experimental plan
  #         -> the level of ref for each factor
  #         -> the formulae
  #         -> the model
  #         -> the contrasts
  ##########################################
  
  # as soon as the "load" button has been clicked
  #  => the loadExpDesign function is called and the experimental design item is printed
  observeEvent(input$loadExpDesign, {
    
        print("# 1- Load experimental design...")
        
        ### Experimental Design
        if(is.null(input$Experimental.Design.file)){
          showModal(modalDialog(title = "Error message", "Experimental Design is required"))
        }
        validate({ need(! is.null(input$Experimental.Design.file), message="Set a name") })
        
        ExpDesign.tbl <<- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1, sep = "\t")
        
        # display desgin table
        output$ExpDesignTable <- renderUI({
          
          box(width = 9, status = "warning",
              DT::renderDataTable( DT::datatable(ExpDesign.tbl) )
          )
        })
        
        ####### Set up Design model ########
        output$GetdFactorRef <- renderUI({
    
          box(status = "warning", width = 9, height = NULL,
    
              # Construct the form to enter the name of the factor
              h4("Enter a name for each design factor"),
              fluidRow(
                lapply(names(ExpDesign.tbl), function(i) {
                  box(width=3,  
                      textInput(session$ns(paste0("dF.Name.", i)), label=NULL , value = i, width = NULL, placeholder = NULL))})
              ),
              
              # Construct the form to set the reference factor level
              h4("Select the level of reference fo each design factor"),
              fluidRow(
                lapply(names(ExpDesign.tbl), function(i) {
                  box(width=3, 
                      selectInput(session$ns(paste0("dF.RefLevel.", i)), i, choices = levels(as.factor(ExpDesign.tbl[[i]])), selectize=FALSE, size=5))})
              ),
    
              # Construct the form to set the type of the factor (either biological or batch)
              h4("Select the type of the design factor"),
              fluidRow(
                lapply(names(ExpDesign.tbl), function(i) {
                  box(width=3, 
                      radioButtons(session$ns(paste0("dF.Type.", i)), label=NULL , choices = c("Bio","batch"), selected = "Bio", inline = FALSE, 
                                                                      width = 2, choiceNames = NULL, choiceValues = NULL))})
              ),
              actionButton(session$ns("ValidF"),"Valid factor set up")
          )
    
        })
  })
  
  # as soon as the "Valid factor set up" button has been clicked
  #  => The upvdateDesignFactors function is called
  #  => The interface to select the model formulae appear
  validate.status <<- 0
  observeEvent(input$ValidF, {
    print("# 2- Set design model...")

    # Get the Type, ref and the new name of the factors that the users enter in the form
    dF.Type.dFac<-vector()
    dF.List.Name<-vector()
    dF.List.ref <-vector()
    
    validate.status <<- 0
    for(dFac in names(ExpDesign.tbl)){
      
      # check if factor names are not empty
      if(input[[paste0("dF.Name.",dFac)]]==""){
        
        showModal(modalDialog( title = "Error message", "Empty name factor are not allowed" ))
        validate.status <<- 1
      }
      # list of names of factors
      dF.List.Name[dFac] <- input[[paste0("dF.Name.",dFac)]]
      # list of type of factors (bio or batch)
      dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
      # list of level reference of factors
      dF.List.ref[dFac]  <- input[[paste0("dF.RefLevel.",dFac)]]
      
    }
    
    validate({
      need(validate.status == 0, message="")
    })
        
    # check if factor names are unique
    if(length(names(ExpDesign.tbl)) != length(unique(dF.List.Name))){
      showModal(modalDialog( title = "Error message", "Factor names must be unique" ))
      validate.status <<- 1
    }

    validate({
      need(validate.status == 0, message="")
    })
    
    ## construct ExperimentalDesign object
    names(ExpDesign.tbl) <- dF.List.Name
    
    print(head(ExpDesign.tbl))
    print(dF.List.ref)
    print(dF.Type.dFac)
    Design <<- ExpDesign.constructor(ExpDesign = ExpDesign.tbl, refList = dF.List.ref, typeList = dF.Type.dFac)
    
    
    #### check experimental design : experimental design must be a complete and balanced.
    completeCheckRes <- CheckExpDesignCompleteness(Design)

    output$Completeness <- renderUI({

        box( status = "warning", width = 9,

             # print message
             renderText( completeCheckRes[["message"]][2] ),

             # plot of count per condition
             renderPlot( plotExperimentalDesign(completeCheckRes[["count"]] )))
    })



    ## error/warning message
    if(completeCheckRes[["message"]][1] == "false"){
      showModal(modalDialog(title = "Error message", completeCheckRes[["message"]][2]))
      validate.status <<- 1
    }

    if(completeCheckRes[["message"]][1] == "warning"){
      showModal(modalDialog( title = "Warning message", completeCheckRes[["message"]][2] ))
    }


    # continue only if message is true or warning
    validate({
      need(validate.status == 0 ,message="ok")
    })
    
  })

  return(input)
  

}


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

