coverPageUI <- function(){
  
  fluidPage(
    
    fluidRow(style="display: flex; align-items: center;",
             column(12,
                    div(
                      h2(tags$span("Introduction", style = "color:orange")),
                      hr(),
                      p(
                        tags$b("Rflomics is a (",tags$i("under devlopement"),") R package coupled with a shiny application"),
                        "dedicated to the management and analysis of multiple omics-datasets (RNAseq illumina, proteomic LC-MS/MS, … )
                    in the statistical framework of", tags$b("vertical integration")," of observations ",tags$i("(i.e. analysis of omics data
                    across experiments on the same indivuals)."),
                    tags$hr(), tags$img(src = "/www/Rflomics_features.png", width = "600px", height = "300px" )),
                    hr(),
                    p("For the moment, Rflomics can deal",
                      tags$b("with multi-factorial experiments"),"(up to 3 biological factors with a batch factor) and
                    helps to set up the statistical model and contrasts associated to the biological experiments. RNA,
                    proteomic and metabolomic datasets are then analysed with expert methods and parameters.", 
                    "It can allows the remote computing for time/cpu consuming tasks (",tags$a("clustermq package",href="#ref4"),")"),
                    hr(),
                    h2(tags$span("Aims", style = "color:orange")),
                    tags$ul(
                      tags$li("Guarantee the relevance of the used methods and parameters (RNAseq workflow: ", 
                              tags$a("DicoExpress",href="#ref1"),")"),
                      tags$li("Decrease the time spend to code robust analysis"),
                      tags$li("Allow experiment sharing and capitalizing"),
                      tags$li("Ensure the reproducibility of omics analysis")),
                    hr(),
                    h3(tags$span("Datasets example", style = "color:orange")),
                    tags$ul(
                      tags$li("Ecoseed dataset:
                    data have been provided by Loic Rajjou and Gwendal Cueff.
                    They are included in the inst/ExampleFiles/ecoseed directory of the package. Briefly,
                    A. thaliana's transcriptoms, proteoms and metaboloms have been obtained in the context of the study of seed
                    germination and vigor. In particular, authors were interested in the influence of
                    temperature (high, medium and low) and imbibition (Dry: DI, early imbibition: EI and late
                    imbibition: LI) on gene's expression."),
                    tags$li(tags$a(href="ecoseed-report.html","An example of report"),),
                    ),
                    hr(),
                    h3(tags$span("Vignettes", style = "color:orange")),
                    tags$ul(
                      tags$li(tags$a(href = "/docs/articles/RFLOMICS-input-data.html", "Input format")),
                      tags$li(tags$a(href = "/docs/articles/RFLOMICS.html", "Quick start"))),
                    hr(),
                    h3(tags$span("Contact and support", style = "color:orange")),
                    tags$a(href="mailto:ijpb-bioinfo-team@inrae.fr","ijpb-bioinfo-team")
                    )
             )
    )
    
  )
}

