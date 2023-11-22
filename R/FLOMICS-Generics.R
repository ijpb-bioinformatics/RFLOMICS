#' @import MultiAssayExperiment
#' @import SummarizedExperiment

# ---- Transformation and normalization: ----

methods::setGeneric(
  name = "RunNormalization",
  def  = function(object, NormMethod = NULL, ... ){standardGeneric("RunNormalization")}
)

methods::setGeneric(
  name = "TransformData",
  def  = function(object, 
                  transformMethod = NULL, 
                  modify_assay = FALSE, ...){standardGeneric("TransformData")}
)

methods::setGeneric(
  name = "FilterLowAbundance",
  def  = function(object, 
                  Filter_Strategy = "NbConditions", 
                  CPM_Cutoff = 5, ... ){standardGeneric("FilterLowAbundance")}
)

# ---- Diff Analysis ----

methods::setGeneric(
  name = "RunDiffAnalysis",
  def  = function(object,  
                  design = NULL,
                  Adj.pvalue.method="BH",
                  contrastList = NULL, 
                  DiffAnalysisMethod = NULL, 
                  Adj.pvalue.cutoff = 0.05, 
                  logFC.cutoff = 0, 
                  clustermq=FALSE, 
                  parallel = FALSE, 
                  nworkers = 1, ... ){standardGeneric("RunDiffAnalysis")}
)

methods::setGeneric(
  name = "FilterDiffAnalysis",
  def  = function(object, 
                  Adj.pvalue.cutoff = NULL, 
                  logFC.cutoff = NULL, ... ){standardGeneric("FilterDiffAnalysis")}
)

methods::setGeneric(
  name = "runDataProcessing",
  def  = function(object, ... ){standardGeneric("runDataProcessing")}
)


methods::setGeneric(
  name = "runSampleFiltering",
  def  = function(object, ... ){standardGeneric("runSampleFiltering")}
)

# ---- PCA ----

methods::setGeneric(
  name = "plotPCA",
  def  = function(object, PCA, PCs = c(1,2), condition = "groups", ... ){standardGeneric("plotPCA")}
)

methods::setGeneric(
  name = "RunPCA",
  def  = function(object, nbcp = 5, raw = FALSE, ... ){standardGeneric("RunPCA")}
)

# ---- Design ----

methods::setGeneric(
  name = "CheckExpDesign",
  def  = function(object, sampleList=NULL, ... ){standardGeneric("CheckExpDesign")}
)

methods::setGeneric(
  name = "CheckExpDesignCompleteness",
  def  = function(object, sampleList=NULL, ... ){standardGeneric("CheckExpDesignCompleteness")}
)

methods::setGeneric(
  name = "getContrastMatrix",
  def  = function(object, modelFormula, contrastList, ... ){standardGeneric("getContrastMatrix")}
)

methods::setGeneric(
  name = "Datasets_overview_plot",
  def  = function(object, ... ){standardGeneric("Datasets_overview_plot")}
)

# ---- Co Expression ----

methods::setGeneric(
  name = "runCoExpression",
  def  = function(object, 
                  K=2:20, 
                  replicates=5, 
                  nameList = NULL, 
                  merge = "union",
                  model = "Normal", 
                  GaussianModel = NULL,
                  transformation = NULL, 
                  normFactors = NULL, 
                  clustermq = FALSE,
                  meanFilterCutoff = NULL, 
                  scale = NULL,
                  silent = TRUE,
                  cmd = FALSE,  ... ){standardGeneric("runCoExpression")}
)

methods::setGeneric(
  name = "coseq.profile.plot",
  def  = function(object, ... ){standardGeneric("coseq.profile.plot")}
)

methods::setGeneric(
  name = "CoExpressionPlots",
  def  = function(object, ... ){standardGeneric("CoExpressionPlots")}
)


# ---- Diff Plots ----

methods::setGeneric(
  name = "DiffAnal.plot",
  def  = function(object, hypothesis, 
                  typeofplots = c("MA.plot", "volcano", "histogram"), ... ){standardGeneric("DiffAnal.plot")}
)


methods::setGeneric(
  name = "Data_Distribution_plot",
  def  = function(object, plot = "boxplot", raw = FALSE, ... ){standardGeneric("Data_Distribution_plot")}
)

methods::setGeneric(
  name = "Library_size_barplot.plot",
  def  = function(object, raw = FALSE, ... ){standardGeneric("Library_size_barplot.plot")}
)

methods::setGeneric(
  name = "heatmapPlot",
  def  = function(object, 
                  hypothesis, 
                  condition="none", 
                  title = "", 
                  annot_to_show = NULL, 
                  subset_list = NULL, 
                  draw_args = list(), 
                  heatmap_args = list(),
                  ...){standardGeneric("heatmapPlot")}
)

methods::setGeneric(
  name = "boxplot.DE.plot",
  def  = function(object, ... ){standardGeneric("boxplot.DE.plot")}
)

#---- ANNOTATION ----

methods::setGeneric(
  name = "runAnnotationEnrichment",
  def  = function(object,
                  nameList = NULL,
                  list_args = list(),
                  from = "DiffExpAnal",
                  dom.select = "custom",
                  Domain = "no-domain",
                  col_term = "term",
                  col_gene = "gene",
                  col_name = "name",
                  col_domain = NULL,
                  annot = NULL, ... ){standardGeneric("runAnnotationEnrichment")}
)

methods::setGeneric(
  name = "plotCPRKEGG",
  def  = function(object,
                  contrast,
                  pathway_id = NULL,
                  species = "ath",
                  gene_idtype = "kegg",
                  from = "DiffExpAnal",
                  pvalueCutoff = object@metadata$DiffExpEnrichAnal[["KEGG"]]$list_args$pvalueCutoff,
                  ...){standardGeneric("plotCPRKEGG")}
)

methods::setGeneric(
  name = "plotCPR",
  def  = function(object,
                  contrast,
                  ont,
                  Domain = NULL,
                  from = "DiffExpAnal",
                  type = "dotplot",
                  showCategory = 15,
                  searchExpr = "",
                  node_label = "all",
                  pvalueCutoff = object@metadata$DiffExpEnrichAnal[[ont]]$list_args$pvalueCutoff,
                  ...){standardGeneric("plotCPR")}
)

# ---- INTEGRATION ----

methods::setGeneric(
  name = "integrationWrapper",
  def  = function(object,
                  omicsToIntegrate = NULL,
                  rnaSeq_transfo = "limma (voom)",
                  choice = "DE",
                  contrasts_names = NULL,
                  type = "union",
                  group = NULL,
                  method = "MOFA",
                  scale_views = FALSE,
                  maxiter = 1000,
                  num_factors = 10,
                  selectedResponse = NULL,
                  ncomp = 2,
                  link_datasets = 1,
                  link_response = 1,
                  sparsity = FALSE,
                  cases_to_try = 5, 
                  silent = TRUE, 
                  cmd = FALSE, ... ){standardGeneric("integrationWrapper")}
)


methods::setGeneric(
  name = "prepareForIntegration",
  def  = function(object,
                  omicsToIntegrate = NULL,
                  rnaSeq_transfo = "limma (voom)",
                  choice = "DE",
                  contrasts_names = NULL,
                  type = "union",
                  group = NULL,
                  method = c("MOFA", "MixOmics"), ... ){standardGeneric("prepareForIntegration")}
)

# ---- RESET ----

methods::setGeneric(
  name = "resetFlomicsMultiAssay",
  def  = function(object, results, datasets = NULL, ... ){standardGeneric("resetFlomicsMultiAssay")}
)