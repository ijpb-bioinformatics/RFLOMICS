
methods::setGeneric(
  name = "mvQCdesign",
  def  = function(object, ... ){standardGeneric("mvQCdesign")}
)

methods::setGeneric(
  name = "mvQCdata",
  def  = function(object, ... ){standardGeneric("mvQCdata")}
)

# ---- Transformation and normalization: ----

methods::setGeneric(
  name = "RunNormalization",
  def  = function(object, NormMethod = NULL, ... ){standardGeneric("RunNormalization")}
)

methods::setGeneric(
  name = "TransformData",
  def  = function(object, 
                  transform_method = NULL, 
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
  def  = function(object, ... ){standardGeneric("RunDiffAnalysis")}
)

methods::setGeneric(
  name = "FilterDiffAnalysis",
  def  = function(object, ... ){standardGeneric("FilterDiffAnalysis")}
)

methods::setGeneric(
  name = "plotPCA",
  def  = function(object, ... ){standardGeneric("plotPCA")}
)





methods::setGeneric(
  name = "RunPCA",
  def  = function(object, ... ){standardGeneric("RunPCA")}
)

methods::setGeneric(
  name = "CheckExpDesignCompleteness",
  def  = function(object, ... ){standardGeneric("CheckExpDesignCompleteness")}
)

methods::setGeneric(
  name = "getExpressionContrast",
  def  = function(object, ... ){standardGeneric("getExpressionContrast")}
)


methods::setGeneric(
  name = "getContrastMatrix",
  def  = function(object, ... ){standardGeneric("getContrastMatrix")}
)

methods::setGeneric(
  name = "runAnnotationEnrichment",
  def  = function(object, ... ){standardGeneric("runAnnotationEnrichment")}
)

methods::setGeneric(
  name = "runCoExpression",
  def  = function(object, ... ){standardGeneric("runCoExpression")}
)

methods::setGeneric(
  name = "DiffAnal.plot",
  def  = function(object, ... ){standardGeneric("DiffAnal.plot")}
)


methods::setGeneric(
  name = "Data_Distribution_plot",
  def  = function(object, ... ){standardGeneric("Data_Distribution_plot")}
)

methods::setGeneric(
  name = "Library_size_barplot.plot",
  def  = function(object, ... ){standardGeneric("Library_size_barplot.plot")}
)

methods::setGeneric(
  name = "Enrichment.plot",
  def  = function(object, ... ){standardGeneric("Enrichment.plot")}
)

methods::setGeneric(
  name = "resetFlomicsMultiAssay",
  def  = function(object, ... ){standardGeneric("resetFlomicsMultiAssay")}
)


methods::setGeneric(
  name = "integrationWrapper",
  def  = function(object, ... ){standardGeneric("integrationWrapper")}
)


methods::setGeneric(
  name = "prepareForIntegration",
  def  = function(object, ... ){standardGeneric("prepareForIntegration")}
)

methods::setGeneric(
  name = "run_MOFA_analysis",
  def  = function(object, ... ){standardGeneric("run_MOFA_analysis")}

)

methods::setGeneric(
  name = "Datasets_overview_plot",
  def  = function(object, ... ){standardGeneric("Datasets_overview_plot")}
)

methods::setGeneric(
  name = "run_MixOmics_analysis",
  def  = function(object, ... ){standardGeneric("run_MixOmics_analysis")}
)

methods::setGeneric(
  name = "runAnnotationEnrichment_CPR",
  def  = function(object, ... ){standardGeneric("runAnnotationEnrichment_CPR")}
)

methods::setGeneric(
  name = "plot.CPRKEGG_Results",
  def  = function(object, ... ){standardGeneric("plot.CPRKEGG_Results")}
)

methods::setGeneric(
  name = "plot.CPR_Results",
  def  = function(object, ... ){standardGeneric("plot.CPR_Results")}
)

methods::setGeneric(
  name = "heatmap.plot",
  def  = function(object, ... ){standardGeneric("heatmap.plot")}
)

methods::setGeneric(
  name = "boxplot.DE.plot",
  def  = function(object, ... ){standardGeneric("boxplot.DE.plot")}
)

methods::setGeneric(
  name = "coseq.profile.plot",
  def  = function(object, ... ){standardGeneric("coseq.profile.plot")}
)