#' @import MultiAssayExperiment
#' @import SummarizedExperiment


# ---- RflomicsMAE and RflomicsSE getteurs ----
methods::setGeneric(
  name = "subRflomicsMAE",
  def  = function(object, omicNames = NULL, ...){standardGeneric("subRflomicsMAE")}
)

methods::setGeneric(
  name = "getFactorNames",
  def  = function(object){standardGeneric("getFactorNames")}
)

methods::setGeneric(
  name = "getFactorTypes",
  def  = function(object){standardGeneric("getFactorTypes")}
)

methods::setGeneric(
  name = "getBioFactors",
  def  = function(object){standardGeneric("getBioFactors")}
)

methods::setGeneric(
  name = "getBatchFactors",
  def  = function(object){standardGeneric("getBatchFactors")}
)

methods::setGeneric(
  name = "getMetaFactors",
  def  = function(object){standardGeneric("getMetaFactors")}
)


methods::setGeneric(
  name = "bioFactors",
  def  = function(object){standardGeneric("bioFactors")}
)

methods::setGeneric(
  name = "batchFactors",
  def  = function(object){standardGeneric("batchFactors")}
)

methods::setGeneric(
  name = "metaFactors",
  def  = function(object){standardGeneric("metaFactors")}
)

methods::setGeneric(
  name = "getDesignMat",
  def  = function(object){standardGeneric("getDesignMat")}
)

methods::setGeneric(
  name = "getDatasetNames",
  def  = function(object){standardGeneric("getDatasetNames")}
)

methods::setGeneric(
  name = "getOmicsTypes",
  def  = function(object){standardGeneric("getOmicsTypes")}
)

methods::setGeneric(
  name = "getRflomicsSE",
  def  = function(object, datasetName = NULL, ...){standardGeneric("getRflomicsSE")}
)


methods::setGeneric(
  name = "getFactorModalities",
  def  = function(object, factorName, ...){standardGeneric("getFactorModalities")}
)
# ---- 02 ----

methods::setGeneric(
  name = "generateModelFormulae",
  def  = function(object, ... ){standardGeneric("generateModelFormulae")}
)
methods::setGeneric(
  name = "setModelFormula",
  def  = function(object, modelFormula = NULL, ... ){standardGeneric("setModelFormula")}
)
methods::setGeneric(
  name = "getModelFormula",
  def  = function(object, ... ){standardGeneric("getModelFormula")}
)
methods::setGeneric(
  name = "generateExpressionContrast",
  def  = function(object, ... ){standardGeneric("generateExpressionContrast")}
)
methods::setGeneric(
  name = "generateContrastMatrix",
  def  = function(object, ... ){standardGeneric("generateContrastMatrix")}
)
methods::setGeneric(
  name = "setContrastList",
  def  = function(object, contrastList=NULL, ... ){standardGeneric("setContrastList")}
)
methods::setGeneric(
  name = "setSelectedContrasts",
  def  = function(object, contrastList=NULL, ... ){standardGeneric("setSelectedContrasts")}
)
methods::setGeneric(
  name = "getSelectedContrasts",
  def  = function(object, contrastList=NULL, ... ){standardGeneric("getSelectedContrasts")}
)
methods::setGeneric(
  name = "setValidContrasts",
  def  = function(object, contrastList=NULL, ... ){standardGeneric("setValidContrasts")}
)
methods::setGeneric(
  name = "getValidContrasts",
  def  = function(object, contrastList=NULL, ... ){standardGeneric("getValidContrasts")}
)

methods::setGeneric(
  name = "plotConditionsOverview",
  def  = function(object, omicNames=NULL, ... ){standardGeneric("plotConditionsOverview")}
)

methods::setGeneric(
  name = "checkExpDesignCompleteness",
  def  = function(object, sampleList=NULL, ... ){standardGeneric("checkExpDesignCompleteness")}
)

methods::setGeneric(
  name = "plotExpDesignCompleteness",
  def  = function(object, sampleList=NULL, ... ){standardGeneric("plotExpDesignCompleteness")}
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
  name = "plotDataOverview",
  def  = function(object, omicNames=NULL, realSize=FALSE, ... ){standardGeneric("plotDataOverview")}
)


########### data Explore ###########

# ---- PCA ----

methods::setGeneric(
  name = "plotPCA",
  def  = function(object, raw=c("raw","norm"), axes = c(1,2), groupColor = "groups", ... ){standardGeneric("plotPCA")}
)

methods::setGeneric(
  name = "runOmicsPCA",
  def  = function(object, ncomp = 5, raw = FALSE , ... ){standardGeneric("runOmicsPCA")}
)


# ---- Transformation and normalization: ----

methods::setGeneric(
  name = "runNormalization",
  def  = function(object, normMethod = NULL, modifyAssay = FALSE, ... ){standardGeneric("runNormalization")}
)

methods::setGeneric(
  name = "runTransformData",
  def  = function(object, 
                  transformMethod = NULL, 
                  modifyAssay = FALSE, ...){standardGeneric("runTransformData")}
)

methods::setGeneric(
  name = "filterLowAbundance",
  def  = function(object, 
                  filterMethod= "CPM",
                  filterStrategy = "NbConditions", 
                  cpmCutoff = 5, ... ){standardGeneric("filterLowAbundance")}
)


methods::setGeneric(
  name = "runDataProcessing",
  def  = function(object, ... ){standardGeneric("runDataProcessing")}
)


methods::setGeneric(
  name = "runSampleFiltering",
  def  = function(object, samples=NULL, ... ){standardGeneric("runSampleFiltering")}
)

methods::setGeneric(
  name = "getTransSettings",
  def  = function(object, ... ){standardGeneric("getTransSettings")}
)

methods::setGeneric(
  name = "getTransSettings",
  def  = function(object, ... ){standardGeneric("getTransSettings")}
)

methods::setGeneric(
  name = "setTrans",
  def  = function(object, ... ){standardGeneric("setTrans")}
)

methods::setGeneric(
  name = "getNormSettings",
  def  = function(object, ... ){standardGeneric("getNormSettings")}
)

methods::setGeneric(
  name = "setNorm",
  def  = function(object, ... ){standardGeneric("setNorm")}
)

methods::setGeneric(
  name = "getFilterSettings",
  def  = function(object, ... ){standardGeneric("getFilterSettings")}
)


methods::setGeneric(
  name = "getFilteredFeatures",
  def  = function(object, ... ){standardGeneric("getFilteredFeatures")}
)


methods::setGeneric(
  name = "getCoeffNorm",
  def  = function(object, ... ){standardGeneric("getCoeffNorm")}
)

methods::setGeneric(
  name = "setCoeffNorm",
  def  = function(object, ... ){standardGeneric("setCoeffNorm")}
)


################# Diff Analysis #################

methods::setGeneric(
  name = "runDiffAnalysis",
  def  = function(object,  
                  p.adj.method="BH",
                  contrastList = NULL, 
                  method = NULL, 
                  p.adj.cutoff = 0.05, 
                  logFC.cutoff = 0, 
                  clustermq=FALSE, 
                  parallel = FALSE, 
                  nworkers = 1, ... ){standardGeneric("runDiffAnalysis")}
)

methods::setGeneric(
  name = "filterDiffAnalysis",
  def  = function(object, 
                  p.adj.cutoff = NULL, 
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
  name = "getDiffSettings",
  def  = function(object, ... ){standardGeneric("getDiffSettings")}
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
                  contrastName, 
                  splitFactor="none", 
                  title = "", 
                  annotNames = NULL, 
                  modalities = NULL, 
                  drawArgs = list(), 
                  heatmapArgs = list(),
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
                  contrastNames = NULL, 
                  merge="union",
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
  name = "plotDataDistribution",
  def  = function(object, plot = "boxplot", raw = FALSE, ... ){standardGeneric("plotDataDistribution")}
)

methods::setGeneric(
  name = "plotLibrarySize",
  def  = function(object, raw = FALSE, ... ){standardGeneric("plotLibrarySize")}
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