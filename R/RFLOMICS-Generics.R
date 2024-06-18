#' @import methods

setGeneric(
    name = "getCoExpAnalysesSummary",
    def  = function(object, ...) {
        standardGeneric("getCoExpAnalysesSummary")
    }
)
setGeneric(
    name = "getAnnotAnalysesSummary",
    def  = function(object,
                    from       = "DiffExpEnrichAnal",
                    matrixType = "presence",
                    ...) {
        standardGeneric("getAnnotAnalysesSummary")
    }
)
setGeneric(
    name = "getDiffAnalysesSummary",
    def  = function(object, plot = FALSE, ...) {
        standardGeneric("getDiffAnalysesSummary")
    }
)
setGeneric(
    name = "generateReport",
    def  = function(object,
                    fileName = NULL,
                    archiveName = NULL,
                    export = FALSE,
                    tmpDir = getwd(),
                    ...) {
        standardGeneric("generateReport")
    }
)

# ---- RflomicsMAE and RflomicsSE getteurs ----

setGeneric(
    name = "getProjectName",
    def  = function(object, omicNames = NULL, ...) {
        standardGeneric("getProjectName")
    }
)
setGeneric(
    name = "subRflomicsMAE",
    def  = function(object, omicNames = NULL, ...) {
        standardGeneric("subRflomicsMAE")
    }
)

setGeneric(
    name = "getFactorNames",
    def  = function(object) {
        standardGeneric("getFactorNames")
    }
)

setGeneric(
    name = "getFactorTypes",
    def  = function(object) {
        standardGeneric("getFactorTypes")
    }
)

setGeneric(
    name = "getBioFactors",
    def  = function(object) {
        standardGeneric("getBioFactors")
    }
)

setGeneric(
    name = "getBatchFactors",
    def  = function(object) {
        standardGeneric("getBatchFactors")
    }
)

setGeneric(
    name = "getMetaFactors",
    def  = function(object) {
        standardGeneric("getMetaFactors")
    }
)

setGeneric(
    name = "getDesignMat",
    def  = function(object) {
        standardGeneric("getDesignMat")
    }
)

setGeneric(
    name = "getDatasetNames",
    def  = function(object) {
        standardGeneric("getDatasetNames")
    }
)

setGeneric(
    name = "getOmicsTypes",
    def  = function(object) {
        standardGeneric("getOmicsTypes")
    }
)

setGeneric(
    name = "getRflomicsSE",
    def  = function(object, datasetName = NULL, ...) {
        standardGeneric("getRflomicsSE")
    }
)


setGeneric(
    name = "getFactorModalities",
    def  = function(object, factorName, ...) {
        standardGeneric("getFactorModalities")
    }
)
# ---- 02 ----

setGeneric(
    name = "generateModelFormulae",
    def  = function(object, ...) {
        standardGeneric("generateModelFormulae")
    }
)
setGeneric(
    name = "setModelFormula",
    def  = function(object, modelFormula = NULL, ...) {
        standardGeneric("setModelFormula")
    }
)
setGeneric(
    name = "getModelFormula",
    def  = function(object, ...) {
        standardGeneric("getModelFormula")
    }
)
setGeneric(
    name = "generateExpressionContrast",
    def  = function(object, ...) {
        standardGeneric("generateExpressionContrast")
    }
)
setGeneric(
    name = "generateContrastMatrix",
    def  = function(object, ...) {
        standardGeneric("generateContrastMatrix")
    }
)
setGeneric(
    name = "setContrastList",
    def  = function(object, contrastList = NULL, ...) {
        standardGeneric("setContrastList")
    }
)
setGeneric(
    name = "setSelectedContrasts",
    def  = function(object, contrastList = NULL, ...) {
        standardGeneric("setSelectedContrasts")
    }
)

setGeneric(
    name = "getSelectedContrasts",
    def  = function(object, contrastList = NULL, ...) {
        standardGeneric("getSelectedContrasts")
    }
)
setGeneric(
    name = "setValidContrasts",
    def  = function(object, contrastList = NULL, ...) {
        standardGeneric("setValidContrasts")
    }
)
setGeneric(
    name = "getValidContrasts",
    def  = function(object, contrastList = NULL, ...) {
        standardGeneric("getValidContrasts")
    }
)

setGeneric(
    name = "plotConditionsOverview",
    def  = function(object, omicNames = NULL, ...) {
        standardGeneric("plotConditionsOverview")
    }
)

setGeneric(
    name = "checkExpDesignCompleteness",
    def  = function(object, sampleList = NULL, ...) {
        standardGeneric("checkExpDesignCompleteness")
    }
)

setGeneric(
    name = "plotExpDesignCompleteness",
    def  = function(object, sampleList = NULL, ...) {
        standardGeneric("plotExpDesignCompleteness")
    }
)

# ---- GET / SET ----

setGeneric(
    name = "getContrastMatrix",
    def  = function(object, modelFormula, contrastList, ...) {
        standardGeneric("getContrastMatrix")
    }
)


setGeneric(
    name = "getExperimentNames",
    def  = function(object) {
        standardGeneric("getExperimentNames")
    }
)

setGeneric(
    name = "getExperimentTypes",
    def  = function(object, experimentNames = NULL, ...) {
        standardGeneric("getExperimentTypes")
    }
)

# ---- plot ----

setGeneric(
    name = "plotDataOverview",
    def  = function(object,
                    omicNames = NULL,
                    realSize = FALSE,
                    ...) {
        standardGeneric("plotDataOverview")
    }
)


########### data Explore ###########

# ---- PCA ----

setGeneric(
    name = "plotOmicsPCA",
    def  = function(object,
                    raw = c("raw", "norm"),
                    axes = c(1, 2),
                    groupColor = "groups",
                    ...) {
        standardGeneric("plotOmicsPCA")
    }
)

setGeneric(
    name = "runOmicsPCA",
    def  = function(object,
                    ncomp = 5,
                    raw = FALSE ,
                    ...) {
        standardGeneric("runOmicsPCA")
    }
)


# ---- Transformation and normalization: ----

setGeneric(
    name = "runNormalization",
    def  = function(object,
                    normMethod = NULL,
                    modifyAssay = FALSE,
                    ...) {
        standardGeneric("runNormalization")
    }
)

setGeneric(
    name = "runTransformData",
    def  = function(object,
                    transformMethod = NULL,
                    modifyAssay = FALSE,
                    ...) {
        standardGeneric("runTransformData")
    }
)

setGeneric(
    name = "filterLowAbundance",
    def  = function(object,
                    filterMethod = "CPM",
                    filterStrategy = "NbConditions",
                    cpmCutoff = 5,
                    ...) {
        standardGeneric("filterLowAbundance")
    }
)


setGeneric(
    name = "runDataProcessing",
    def  = function(object, ...) {
        standardGeneric("runDataProcessing")
    }
)


setGeneric(
    name = "runSampleFiltering",
    def  = function(object, samples = NULL, ...) {
        standardGeneric("runSampleFiltering")
    }
)

setGeneric(
    name = "getTransSettings",
    def  = function(object, ...) {
        standardGeneric("getTransSettings")
    }
)

setGeneric(
    name = "getTransSettings",
    def  = function(object, ...) {
        standardGeneric("getTransSettings")
    }
)

setGeneric(
    name = "setTrans",
    def  = function(object, ...) {
        standardGeneric("setTrans")
    }
)

setGeneric(
    name = "getNormSettings",
    def  = function(object, ...) {
        standardGeneric("getNormSettings")
    }
)

setGeneric(
    name = "setNorm",
    def  = function(object, ...) {
        standardGeneric("setNorm")
    }
)

setGeneric(
    name = "getFilterSettings",
    def  = function(object, ...) {
        standardGeneric("getFilterSettings")
    }
)


setGeneric(
    name = "getFilteredFeatures",
    def  = function(object, ...) {
        standardGeneric("getFilteredFeatures")
    }
)


setGeneric(
    name = "getCoeffNorm",
    def  = function(object, ...) {
        standardGeneric("getCoeffNorm")
    }
)

setGeneric(
    name = "setCoeffNorm",
    def  = function(object, ...) {
        standardGeneric("setCoeffNorm")
    }
)


################# Diff Analysis #################

setGeneric(
    name = "runDiffAnalysis",
    def  = function(object,
                    contrastList = NULL,
                    p.adj.method = "BH",
                    method = NULL,
                    p.adj.cutoff = 0.05,
                    logFC.cutoff = 0,
                    clustermq = FALSE,
                    ...) {
        standardGeneric("runDiffAnalysis")
    }
)

setGeneric(
    name = "filterDiffAnalysis",
    def  = function(object,
                    p.adj.cutoff = NULL,
                    logFC.cutoff = NULL,
                    ...) {
        standardGeneric("filterDiffAnalysis")
    }
)

# ----GET/SET---------

setGeneric(
    name = "getDEList",
    def  = function(object, contrasts = NULL, ...) {
        standardGeneric("getDEList")
    }
)

setGeneric(
    name = "getDEMatrix",
    def  = function(object, ...) {
        standardGeneric("getDEMatrix")
    }
)

setGeneric(
    name = "getDiffSettings",
    def  = function(object, ...) {
        standardGeneric("getDiffSettings")
    }
)


# ---- Plots ----

setGeneric(
    name = "plotDiffAnalysis",
    def  = function(object,
                    hypothesis,
                    typeofplots = c("MA.plot", "volcano", "histogram"),
                    ...) {
        standardGeneric("plotDiffAnalysis")
    }
)

setGeneric(
    name = "plotHeatmapDesign",
    def  = function(object,
                    contrastName,
                    splitFactor = "none",
                    title = "",
                    annotNames = NULL,
                    modalities = NULL,
                    drawArgs = list(),
                    heatmapArgs = list(),
                    ...) {
        standardGeneric("plotHeatmapDesign")
    }
)


setGeneric(
    name = "plotBoxplotDE",
    def  = function(object, ...) {
        standardGeneric("plotBoxplotDE")
    }
)


################# Coexp Analysis #################

setGeneric(
    name = "runCoExpression",
    def  = function(object,
                    K = 2:20,
                    replicates = 5,
                    contrastNames = NULL,
                    merge = "union",
                    model = "Normal",
                    GaussianModel = NULL,
                    transformation = NULL,
                    normFactors = NULL,
                    clustermq = FALSE,
                    meanFilterCutoff = NULL,
                    scale = NULL,
                    silent = TRUE,
                    cmd = FALSE,
                    ...) {
        standardGeneric("runCoExpression")
    }
)

# -------- Coexp plot --------#

setGeneric(
    name = "plotCoExpressionProfile",
    def  = function(object, ...) {
        standardGeneric("plotCoExpressionProfile")
    }
)


setGeneric(
    name = "plotCoExpression",
    def  = function(object, ...) {
        standardGeneric("plotCoExpression")
    }
)

setGeneric(
    name = "plotCoseqContrasts",
    def  = function(object, ...) {
        standardGeneric("plotCoseqContrasts")
    }
)


# ---- Coexp GET/SET---------


setGeneric(
    name = "getCoexpSettings",
    def  = function(object, ...) {
        standardGeneric("getCoexpSettings")
    }
)

setGeneric(
    name = "getClusterEntities",
    def  = function(object, ...) {
        standardGeneric("getClusterEntities")
    }
)

setGeneric(
    name = "getCoseqClusters",
    def  = function(object, ...) {
        standardGeneric("getCoseqClusters")
    }
)


setGeneric(
    name = "plotDataDistribution",
    def  = function(object,
                    plot = "boxplot",
                    raw = FALSE,
                    ...) {
        standardGeneric("plotDataDistribution")
    }
)

setGeneric(
    name = "plotLibrarySize",
    def  = function(object, raw = FALSE, ...) {
        standardGeneric("plotLibrarySize")
    }
)


################# ANNOTATION #################

setGeneric(
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
                    annot = NULL,
                    ...) {
        standardGeneric("runAnnotationEnrichment")
    }
)

setGeneric(
    name = "plotKEGG",
    def  = function(object,
                    contrastName,
                    pathway_id = NULL,
                    species = "ath",
                    gene_idtype = "kegg",
                    from = "DiffExpAnal",
                    ...) {
        standardGeneric("plotKEGG")
    }
)

setGeneric(
    name = "plotClusterProfiler",
    def  = function(object,
                    contrastName,
                    database,
                    domain = "no-domain",
                    from = "DiffExp",
                    plotType = "dotplot",
                    showCategory = 15,
                    searchExpr = "",
                    nodeLabel = "all",
                    p.adj.cutoff = NULL,
                    ...) {
        standardGeneric("plotClusterProfiler")
    }
)


setGeneric(
    name = "plotEnrichComp",
    def = function(object,
                   from = "DiffExp",
                   database = NULL,
                   domain = "no-domain",
                   matrixType = "FC",
                   nClust = NULL,
                   ...) {
        standardGeneric("plotEnrichComp")
    }
)
setGeneric(
    name = "getEnrichRes",
    def  = function(object,
                    contrastName = NULL,
                    from = "DiffExpEnrichAnal",
                    database = "GO",
                    domain = NULL,
                    ...) {
        standardGeneric("getEnrichRes")
    }
)


setGeneric(
    name = "getEnrichPvalue",
    def  = function(object,
                    from = "DiffExpEnrichAnal",
                    database = "GO") {
        standardGeneric("getEnrichPvalue")
    }
)

setGeneric(
    name = "getEnrichSettings",
    def  = function(object,
                    from = "DiffExpEnrichAnal",
                    database = "GO") {
        standardGeneric("getEnrichSettings")
    }
)

setGeneric(
    name = "sumORA",
    def  = function(object,
                    from = "DiffExpEnrichAnal",
                    database = NULL,
                    contrastName = NULL) {
        standardGeneric("sumORA")
    }
)

# setGeneric(
#     name = "setEnrichNull",
#     def = function(object,
#                    from = "DiffExp",
#                    database = NULL) {
#         standardGeneric("setEnrichNull")
#     }
# )

################# INTEGRATION #################

setGeneric(
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
                    cmd = FALSE,
                    ...) {
        standardGeneric("integrationWrapper")
    }
)


setGeneric(
    name = "prepareForIntegration",
    def  = function(object,
                    omicsNames = NULL,
                    rnaSeq_transfo = "limma (voom)",
                    variableLists = NULL,
                    group = NULL,
                    method = "MOFA",
                    transformData = TRUE,
                    cmd = FALSE,
                    ...) {
        standardGeneric("prepareForIntegration")
    }
)


setGeneric(
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
                    cmd = FALSE,
                    ...) {
        standardGeneric("runOmicsIntegration")
    }
)


setGeneric(
    name = "filterFeatures",
    def = function(object,
                   selOpt = rep(list("all"), length(object)),
                   type = rep(list("union"), length(selOpt))) {
        standardGeneric("filterFeatures")
    }
)


setGeneric(
    name = "getMixOmics",
    def  = function(object,
                    response = NULL,
                    onlyResults = TRUE) {
        standardGeneric("getMixOmics")
    }
)

setGeneric(
    name = "getMOFA",
    def  = function(object, onlyResults = TRUE) {
        standardGeneric("getMOFA")
    }
)

setGeneric(
    name = "getMOFASettings",
    def  = function(object) {
        standardGeneric("getMOFASettings")
    }
)

setGeneric(
    name = "setMOFA",
    def  = function(object, results = NULL) {
        standardGeneric("setMOFA")
    }
)

setGeneric(
    name = "setMixOmics",
    def  = function(object, results = NULL) {
        standardGeneric("setMixOmics")
    }
)

setGeneric(
    name = "sumMixOmics",
    def  = function(object, selectedResponse = NULL) {
        standardGeneric("sumMixOmics")
    }
)

setGeneric(
    name = ".getOneMORes",
    def  = function(object, selectedResponse) {
        standardGeneric(".getOneMORes")
    }
)


setGeneric(
    name = "getMixOmicsSettings",
    def  = function(object) {
        standardGeneric("getMixOmicsSettings")
    }
)

setGeneric(
    name = "plotMOVarExp",
    def  = function(object,
                    selectedResponse,
                    mode = NULL) {
        standardGeneric("plotMOVarExp")
    }
)

# ---- RESET ----

setGeneric(
    name = "resetRflomicsMAE",
    def  = function(object, 
                    singleAnalyses = NULL, 
                    multiAnalyses = NULL, 
                    datasetNames = NULL, ...) {
        standardGeneric("resetRflomicsMAE")
    }
)

setGeneric(
  name = "getAnalyzedDatasetNames",
  def  = function(object, analyses = NULL, ...) {
    standardGeneric("getAnalyzedDatasetNames")
  }
)


setGeneric(
  name = "getProcessedData",
  def  = function(object, SE.name, ...) {
    standardGeneric("getProcessedData")
  }
)
