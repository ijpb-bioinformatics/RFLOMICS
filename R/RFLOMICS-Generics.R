#' @import MultiAssayExperiment
#' @import SummarizedExperiment


################# set stat #################



methods::setGeneric(
  name = "CheckExpDesign",
  def  = function(object, sampleList=NULL, ... ){standardGeneric("CheckExpDesign")}
)

methods::setGeneric(
  name = "CheckExpDesignCompleteness",
  def  = function(object, sampleList=NULL, ... ){standardGeneric("CheckExpDesignCompleteness")}
)

# ---- GET / SET ----

methods::setGeneric(
  name = "getContrastMatrix",
  def  = function(object, modelFormula, contrastList, ... ){standardGeneric("getContrastMatrix")}
)


methods::setGeneric(
  name = "getExperimentNames",
  def  = function(object){standardGeneric("getExperimentNames")}
)

methods::setGeneric(
  name = "getExperimentTypes",
  def  = function(object, experimentNames = NULL, ...){standardGeneric("getExperimentTypes")}
)

# ---- plot ----

methods::setGeneric(
  name = "Datasets_overview_plot",
  def  = function(object, ... ){standardGeneric("Datasets_overview_plot")}
)


########### data Explore ###########

# ---- PCA ----

methods::setGeneric(
  name = "plotPCA",
  def  = function(object, PCA, PCs = c(1,2), condition = "groups", ... ){standardGeneric("plotPCA")}
)

methods::setGeneric(
  name = "RunPCA",
  def  = function(object, nbcp = 5, raw = FALSE, ... ){standardGeneric("RunPCA")}
)


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


methods::setGeneric(
  name = "runDataProcessing",
  def  = function(object, ... ){standardGeneric("runDataProcessing")}
)


methods::setGeneric(
  name = "runSampleFiltering",
  def  = function(object, ... ){standardGeneric("runSampleFiltering")}
)


################# Diff Analysis #################

methods::setGeneric(
  name = "runDiffAnalysis",
  def  = function(object,  
                  design = NULL,
                  Adj.pvalue.method="BH",
                  contrastList = NULL, 
                  DiffAnalysisMethod = NULL, 
                  Adj.pvalue.cutoff = 0.05, 
                  logFC.cutoff = 0, 
                  clustermq=FALSE, 
                  parallel = FALSE, 
                  nworkers = 1, ... ){standardGeneric("runDiffAnalysis")}
)

methods::setGeneric(
  name = "filterDiffAnalysis",
  def  = function(object, 
                  Adj.pvalue.cutoff = NULL, 
                  logFC.cutoff = NULL, ... ){standardGeneric("filterDiffAnalysis")}
)

# ----GET/SET---------

methods::setGeneric(
  name = "getDEList",
  def  = function(object, ... ){standardGeneric("getDEList")}
)

methods::setGeneric(
  name = "getDEMatrix",
  def  = function(object, ... ){standardGeneric("getDEMatrix")}
)

methods::setGeneric(
  name = "getDiffSetting",
  def  = function(object, ... ){standardGeneric("getDiffSetting")}
)


# ---- Plots ----

methods::setGeneric(
  name = "plotDiffAnalysis",
  def  = function(object, hypothesis, 
                  typeofplots = c("MA.plot", "volcano", "histogram"), ...){standardGeneric("plotDiffAnalysis")}
)

methods::setGeneric(
  name = "plotHeatmapDesign",
  def  = function(object, 
                  hypothesis, 
                  condition="none", 
                  title = "", 
                  annot_to_show = NULL, 
                  subset_list = NULL, 
                  draw_args = list(), 
                  heatmap_args = list(),
                  ...){standardGeneric("plotHeatmapDesign")}
)

methods::setGeneric(
  name = "plotBoxplotDE",
  def  = function(object, ... ){standardGeneric("plotBoxplotDE")}
)


################# Coexp Analysis #################

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

# -------- Coexp plot --------# 

methods::setGeneric(
  name = "plotCoExpressionProfile",
  def  = function(object, ... ){standardGeneric("plotCoExpressionProfile")}
)


methods::setGeneric(
  name = "plotCoExpression",
  def  = function(object, ... ){standardGeneric("plotCoExpression")}
)

methods::setGeneric(
  name = "plotCoseqContrasts",
  def  = function(object, ... ){standardGeneric("plotCoseqContrasts")}
)


# ---- Coexp GET/SET---------


methods::setGeneric(
  name = "getCoexpSetting",
  def  = function(object, ... ){standardGeneric("getCoexpSetting")}
)

methods::setGeneric(
  name = "getClusterEntities",
  def  = function(object, ... ){standardGeneric("getClusterEntities")}
)

methods::setGeneric(
  name = "getCoseqClusters",
  def  = function(object, ... ){standardGeneric("getCoseqClusters")}
)


methods::setGeneric(
  name = "Data_Distribution_plot",
  def  = function(object, plot = "boxplot", raw = FALSE, ... ){standardGeneric("Data_Distribution_plot")}
)

methods::setGeneric(
  name = "Library_size_barplot.plot",
  def  = function(object, raw = FALSE, ... ){standardGeneric("Library_size_barplot.plot")}
)


################# ANNOTATION #################

methods::setGeneric(
  name = "runAnnotationEnrichment",
  def  = function(object,
                  nameList = NULL,
                  list_args = list(
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 1,
                    # minGSSize = 10,
                    minGSSize = 3,
                    maxGSSize = 500,
                    universe = names(object)
                  ),
                  from = "DiffExpAnal",
                  database = "custom",
                  domain = "no-domain",
                  col_term = "term",
                  col_gene = "gene",
                  col_name = "name",
                  col_domain = NULL,
                  annot = NULL, ... ){standardGeneric("runAnnotationEnrichment")}
)

methods::setGeneric(
  name = "plotKEGG",
  def  = function(object,
                  contrastName,
                  pathway_id = NULL,
                  species = "ath",
                  gene_idtype = "kegg",
                  from = "DiffExpAnal",
                  pvalueCutoff = object@metadata$DiffExpEnrichAnal[["KEGG"]]$list_args$pvalueCutoff,
                  ...){standardGeneric("plotKEGG")}
)

methods::setGeneric(
  name = "plotClusterProfiler",
  def  = function(object,
                  contrastName,
                  database,
                  domain = NULL,
                  from = "DiffExpAnal",
                  plotType = "dotplot",
                  showCategory = 15,
                  searchExpr = "",
                  nodeLabel = "all",
                  pvalueCutoff = object@metadata$DiffExpEnrichAnal[[database]]$list_args$pvalueCutoff,
                  ...){standardGeneric("plotClusterProfiler")}
)


methods::setGeneric(
  name = "plotEnrichComp",
  def = function(object, 
                 from = "DiffExp", 
                 database = NULL, 
                 domain = "no-domain",
                 matrixType = "GeneRatio",
                 nClust = NULL,
                 ...
  ){standardGeneric("plotEnrichComp")}
)
methods::setGeneric(
  name = "getEnrichRes",
  def  = function(object,
                  contrastName = NULL,
                  from = "DiffExpEnrichAnal",
                  database = "GO",
                  domain = NULL, ...){standardGeneric("getEnrichRes")}
)


methods::setGeneric(
  name = "getEnrichPvalue",
  def  = function(object,
                  from = "DiffExpEnrichAnal",
                  database = "GO"){standardGeneric("getEnrichPvalue")}
)


methods::setGeneric(
  name = "sumORA",
  def  = function(object, 
                  from = "DiffExpEnrichAnal", 
                  database = NULL, 
                  contrastName = NULL){standardGeneric("sumORA")}
)

################# INTEGRATION #################

methods::setGeneric(
  name = "integrationWrapper",
  def  = function(object,
                  omicsNames = names(object),
                  rnaSeq_transfo = "limma (voom)",
                  selOpt = rep(list("DE"), length(omicsToIntegrate)),
                  type = rep(list("union"), length(selOpt)),
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
                  cmd = FALSE,
                  ... ){standardGeneric("integrationWrapper")}
)


methods::setGeneric(
  name = "prepareForIntegration",
  def  = function(object,
                  omicsNames = NULL,
                  rnaSeq_transfo = "limma (voom)",
                  variableLists = NULL,
                  group = NULL,
                  method = "MOFA",
                  cmd = FALSE, 
                  silent = TRUE, ... ){standardGeneric("prepareForIntegration")}
)


methods::setGeneric(
  name = "runOmicsIntegration",
  def  = function(object,
                  preparedObject = NULL,
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
                  cmd = FALSE,
                  ... ){standardGeneric("runOmicsIntegration")}
)


methods::setGeneric(
  name = "filterFeatures",
  def = function(object,
                 selOpt = rep(list("all"), length(object)),
                 type = rep(list("union"), length(selOpt))
  ){standardGeneric("filterFeatures")}
)


methods::setGeneric(
  name = "getMixOmics",
  def  = function(object, 
                  response = NULL,
                  onlyResults = TRUE){standardGeneric("getMixOmics")}
)

methods::setGeneric(
  name = "getMOFA",
  def  = function(object, onlyResults = TRUE){standardGeneric("getMOFA")}
)

methods::setGeneric(
  name = "getMOFASettings",
  def  = function(object){standardGeneric("getMOFASettings")}
)

methods::setGeneric(
  name = "setMOFA",
  def  = function(object, results = NULL){standardGeneric("setMOFA")}
)

methods::setGeneric(
  name = "setMixOmics",
  def  = function(object, results = NULL){standardGeneric("setMixOmics")}
)

methods::setGeneric(
  name = "sumMixOmics",
  def  = function(object, selectedResponse = NULL){
    standardGeneric("sumMixOmics")}
)

methods::setGeneric(
  name = ".getOneMORes",
  def  = function(object, selectedResponse){standardGeneric(".getOneMORes")}
)


methods::setGeneric(
  name = "getMixOmicsSettings",
  def  = function(object){standardGeneric("getMixOmicsSettings")}
)

# ---- RESET ----

methods::setGeneric(
  name = "resetFlomicsMultiAssay",
  def  = function(object, results, datasets = NULL, ... ){standardGeneric("resetFlomicsMultiAssay")}
)