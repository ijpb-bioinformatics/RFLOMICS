coverPageUI <- function(){
  
  fluidPage(
    
    fluidRow(style="display: flex; align-items: center;",
             column(12,
                    h1(tags$span("An interactive web application for multi-omics data analysis", style = "color:orange")),
                    tags$hr())),
    
    fluidRow(
      column(8,
             tags$div(style="text-align:justify",
                      p(tags$b("RFLOMICS")," is an R package and Shiny interface. 
                        It provides a guided, comprehensive, and reproducible analysis framework capable of handling multiple omics datasets. 
                        The interface offers great flexibility, allowing users to navigate between results exploration and visualization and 
                        between the different steps of the analyses."),
                      p("RFLOMICS code source is available at: https://forgemia.inra.fr/flomics/rflomics"), 
                        # tags$a(href="https://forgemia.inra.fr/flomics/rflomics", "gitLab")),
                      p("RFLOMICS user documentations are available at: https://flomics.pages.mia.inra.fr/rflomics/index.html"), 
                        # tags$a(href="https://flomics.pages.mia.inra.fr/rflomics/index.html", "Vignette")),
                      tags$hr(),
                      
                      h3(tags$span(tags$i("Features"), style = "color:orange")),
                      
                      p("RFLOMICS currently supports three types of omics (transcriptomics (counts), proteomics, and metabolomics),
                        along with multiple datasets per omics type (up to 10).
                        These datasets should be part of the same biological experiment."), 
                        # tags$a(href="https://flomics.pages.mia.inra.fr/rflomics/articles/RFLOMICS-input-data.html",
                      # ""),
                      tags$img(src = "/www/workflow.png", width = "600px", height = "400px" ),
                      
                      tags$hr(),
                      h3(tags$span(tags$i("Specifications"), style = "color:orange")),
                      tags$ul(
                        tags$li("Performs complete single and multi-omics project analysis,"),
                        tags$li("Supports multi-factorial experimental design,"),
                        tags$li("Guarantees the relevance of the used methods,"),
                        tags$li("Reduces the analysis time on the unavoidable steps,"),
                        tags$li("Ensures the reproducibility of omics analyses,"),
                        tags$li("Accessible via one and simple user-friendly interface")),
                      
                      tags$hr(),
                      h3(tags$span(tags$i("Datasets example"), style = "color:orange")),
                      tags$ul(
                        tags$li("Ecoseed dataset:
                    data have been provided by Loic Rajjou and Gwendal Cueff.
                    They are included in the inst/ExampleFiles/ecoseed directory of the package. Briefly,
                    A. thaliana's transcriptoms, proteoms and metaboloms have been obtained in the context of the study of seed
                    germination and vigor. In particular, authors were interested in the influence of
                    temperature (high, medium and low) and imbibition (Dry: DI, early imbibition: EI and late
                    imbibition: LI) on gene's expression."),
                        # tags$li(tags$a(href="ecoseed-report.html","An example of report")),
                      ),
                      tags$hr(),
                      h3(tags$span(tags$i("Contact and support"), style = "color:orange")),
                      tags$a(href="mailto:ijpb-bioinfo-team@inrae.fr","ijpb-bioinfo-team")
                      
             )),
      
      column(4, 
             tags$img(src = "/www/IJPB_logo.png",  width = "240px", height = "120px"),
             tags$hr(), 
             tags$img(src = "/www/INRAE_logo.png", width = "180px", height = "60px"),
             tags$hr(), 
             tags$img(src = "/www/catiSysmics_logo.png", width = "120px", height = "80px")
      )
      
    )
  )
}

GlossaryPageUI <- function(){
  
  GlossaryTable <- matrix(c("Dataset",     "?", "",
                            "Factor",      "?", "",
                            "Batch",       "or technical factor: technical impact on the experience (Replicate, scientist...), not related to the biological question behind the experiment", "",
                            "Bio",         "Biological factor: all the conditions that samples were subject to for the experiment (stress, temperature changes, genotype, time...)" , "",
                            "Meta",        "?",  "",
                            "Balanced design", "a balanced design is a an experimental design where all experimental group have the an equal number of subject observations",  "STATO_0000003",
                            "Complete Design","every combination of levels between the factors is present in the design","",
                             "modality",    "?",  "",
                            "reference",   "?",  "",
                            "Condition",   "?",  "",
                            "Sample",      "?",  "",
                            "Feature",     "?",  "",
                            "Contrast",    "?",  "",
                            "CPM","Count Per Million of mapped reads","",
                            "TMM","Trimmed Mean of M-values: is the weighted mean of log ratios between a sample and a reference, after excluding the most expressed genes and the genes with the largest log ratios","edgeR",
                            "Sample",      "?",  "",
                            "Average",     "?",  "",
                            "Interaction", "a model interaction effect term is a model term which accounts for variation explained by the combined effects of the factor levels of more than one (usually 2) independent variables","(STATO_0000469)"
                              ), byrow = TRUE, ncol = 3)
  
  colnames(GlossaryTable) <- c("Term", "Definition","Source")
  
  fluidRow(
    # load exp design
    # display tab
    box(width = 12, status = "warning", solidHeader = TRUE,
        title = (tags$i("You find on this page the definitions of the terms used in the analyses proposed by Rflomics.")), 
        
        renderDataTable(datatable(data = GlossaryTable))
    )
  )
}