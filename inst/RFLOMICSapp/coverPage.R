coverPage <- fluidPage(

  fluidRow(style="display: flex; align-items: center;",
           column(12,
                  div(
                    h2(tags$span("Introduction", style = "color:orange")),
                    hr(),
                    p(
                      tags$b("Rflomics is a (",tags$i("under devlopement"),") R package coupled with a shiny application"),
                      "dedicated to the management and analysis of multiple omics-datasets (RNAseq illumina, proteomic LC-MS/MS, … )
                    in the statistical framework of", tags$b("vertical integration")," of observations ",tags$i("(i.e. analysis of omics data
                    across experiments on the same indivuals)"),"see the figure below.",
                      tags$hr(), tags$img(src = "integration.png", width = "600px", height = "300px" )),
                    hr(),
                    p("For the moment, Rflomics can deal",
                      tags$b("with multi-factorial experiments"),"(up to 3 biological factors with a batch factor) and
                    helps to set up the statistical model and contrasts associated to the biological experiments. RNA,
                    proteomic and metabolomic datasets are then analysed with expert methods and parameters."),
                    hr(),
                    h2(tags$span("Aims", style = "color:orange")),
                    tags$ul(
                      tags$li("Guarantee the relevance of the used methods and parameters (RNAseq workflow: ", tags$a("DicoExpress",href="#ref1"),")"),
                      tags$li("Decrease the time spend to code robust analysis"),
                      tags$li("Allow experiment sharing and capitalizing"),
                      tags$li("Ensure the reproducibility of omics analysis")),
                    hr(),
                    h2(tags$span("Features", style = "color:orange")),
                    hr(),
                    p("",tags$img(src = "Rflomics_features.png", width = "600px", height = "400px"),
                      tags$hr(),
                      tags$ul(
                        tags$li( "It manages raw and processed (filtered, normalized) datasets (",tags$a("MultiAssayExperiment package",href="#ref2"),")"),
                        tags$li( "It can allows the remote computing for time/cpu consuming tasks (",tags$a("clustermq package",href="#ref4"),")")
                      ),""),
                    hr(),
                    h4("Incoming features"),
                    tags$ul(
                      tags$li("Start to include multi-omics integration methods (MCIA first ?)"),
                      tags$li("Imputation of missing values for proteomics and metabolomics data"),
                    ),
                  h2(tags$span("Input", style = "color:orange")),
                  tags$ul(
                    tags$li("RNAseq (illumina data sequencing technology): matrix of gene expression based on ",tags$b("raw read count"), "quantification. ",
                            tags$br(),"-> gene_id in row and individuals in column"),
                    tags$li("Proteomic (LC-MS/MS mass spectrometry): matrix of protein abundance based on ",tags$b("XIC quantification"), "(Extracted ion chromatograms) ",
                            tags$br(),"-> protein_id in row and individuals in column",
                            tags$br(),"-> preprocessed matrix are expected for the moment (NA imputation, filtering)"),
                    tags$li("Metabolomic (GC-MS mass spectrometry): matrix of metabolomic abundance based on ", tags$b("XIC quantification"), "(Extracted ion chromatograms) ",
                            tags$br(),"-> metabolite_id in row and individuals in column",
                            tags$br(),"-> preprocessed matrix are expected for the moment (NA imputation, filtering)"),
                    ),
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
                  p(tags$a(href="RFLOMICS-vignette-CommandLine.html","CommandLine-Vignette.html")),
                  hr(),
                  h3(tags$span("Contact and support", style = "color:orange")),
                  p("ijpb-bioinfo-team"),
                  hr(),
                  h3(tags$span("References", style = "color:orange")),
                  tags$ul(
                  tags$li("Paysant-Le Roux C. (2020), packageContrast_2020, unpublished"),
                  tags$li(tags$a(href="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=32426025&retmode=ref&cmd=prlinks",
                                   "Ilana L. et al. (2020), DiCoExpress")),
                  tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html","MultiAssayExperiment package")),
                  tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/edgeR.html","edgeR package")),
                  tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/limma.html","limma package")),
                  tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/coseq.html","coseq package")),
                  tags$li(tags$a(href="https://cran.r-project.org/web/packages/clustermq/index.html","clustermq package")),
                  ),
           )
  )
)

)


