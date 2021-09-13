coverPage <- fluidPage(

  fluidRow(style="display: flex; align-items: center;",
           column(12,
                  h1("Welcome"),
                  hr(),
                  div(
                    h3("Introduction"),
                  p(
                    tags$b("Rflomics is a R package coupled with a shiny application"),
                    "dedicated to the managment and analysis of multiple omics-datasets (RNAseq illumina, GC-MS, … )
                    in the statistical framework of", tags$b("vertical integration")," of observations."
                    ,tags$hr(), tags$img(src = "Integration.png", width = "500px", height = "300px")),
                  hr(),
                  p("Rflomics can deal",
                    tags$b("with multi-factorial experiments"),"(up to 3 biological factors with a technical factor) and
                    helps to set up models and contrasts associated to the biological experiment. Omics
                    data are analysed with the same statistical framework."),
                  hr(),
                  h3("Features"),
                  p("For each omic dataset:",tags$ul(
                    tags$li( "It performs data quality controls"),
                    tags$li( "It offers data normalization and filtering methods."),
                    tags$li( "It helps to decompose and explore biological and technical variability."),
                    tags$li( "It manages raw and processed (filtred, normalized) data (",tags$a("MultiAssayExperiment package",href="#ref2"),")"),
                    tags$li( "It does differential expression/abundance analysis (glm framework: ",tags$a("edgeR package",href="#ref2")," and", tags$a("limma package",href="#ref2"),")"),
                    tags$li( "It also performs co-expression/co-abundance analysis (gaussian mixture model: ",tags$a("coseq package",href="#ref2"),") of genes/proteins/metabolites"),
                    tags$li( "It does functional enrichment analysis of clusters (or list of differentially expressed) genes/proteins/metabolites"),
                    tags$li( "It can allows the remote computing for time/cpu consuming tasks (",tags$a("clustermq package",href="#ref4"),")")
                    ),""),
                  hr(),
                  h3("Aims:"),
                  tags$ul(
                    tags$li("Guarantee the relevance of the used methods and parameters (RNAseq workflow: ", tags$a("DicoExpress",href="#ref1"),")"),
                    tags$li("Decrease the time spend to code robust analysis"),
                    tags$li("Allow experiment sharing and capitalizing"),
                    tags$li("Ensure the reproducibility of omics analysis: FAIR code")),
                  hr(),
                  h3("Incoming Features"),
                  tags$ul(
                  tags$li("MOFA"),
                  tags$li("Imputation of missing values for proteomics and metabolomics data"),
                  ),
                  h3("References:"),
                  tags$ul(
                  h5("Experimental design setup: ",id="ref1"),
                  tags$ul(
                  tags$li("[Paysant-Le Roux C. (2020), packageContrast_2020, unpublished]")),
                  h5("RNAseq workflow and parameters: ",id="ref1"),
                  tags$ul(
                    tags$li(tags$a(href="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=32426025&retmode=ref&cmd=prlinks",
                  "[Ilana L. et al. (2020), DiCoExpress]"))),
                  h5("Data managment",id="ref2"),
                  tags$ul(
                    tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html","MultiAssayExperiment")),
                  ),
                  h5("Statistical analysis",id="ref3"),
                  tags$ul(
                    tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/edgeR.html","edgeR")),
                    tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/limma.html","limma")),
                    tags$li(tags$a(href="https://bioconductor.org/packages/release/bioc/html/coseq.html","coseq")),
                  ),
                  h5("Remote computing",id="ref4"),
                  tags$ul(
                    tags$li(tags$a(href="https://cran.r-project.org/web/packages/clustermq/index.html","clustermq")),
                  ),
                  ),
                  hr(),
                  h5("Input"),
                  tags$img(src = "Input1.png", width = "500px", height = "300px"),
                  hr(),
                  tags$img(src = "Input2.png", width = "500px", height = "300px"),
                  hr(),
                  h5("Datasets"),
                  tags$ul(
                    tags$li("ecoseed dataset: Data have been provided by Loic Rajjou and Gwendal Cueff.
                    They are included in the inst/ExampleFiles/ecoseed directory of the package. Briefly,
                    A. thaliana's transcriptomes have been obtained in the context of the study of seed
                    germination and vigor. In particular, the author were interested in the influence of
                    temperature (high, medium and low) and imbibition (Dry: DI, early imbibition: EI and late
                    imbibition: LI) on gene's expression."),
                    tags$li(tags$a(href="ecoseed-report.html","An example of report"),),
                  ),
                  hr(),
                  h3("Vignettes"),
                  p(tags$a(href="RFLOMICS-vignette-CommandLine.html","CommandLine-Vignette.html")),
                  hr(),
                  h3("Contact and support"),
                  p("ijpb-bioinfo-team"),
          )
  )
)
)



