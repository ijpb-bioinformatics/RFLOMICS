
methods::setGeneric(
  name="mvQCdesign",
  def = function( object, ... ){standardGeneric("mvQCdesign")}
)

methods::setGeneric(
  name="mvQCdata",
  def = function( object, ... ){standardGeneric("mvQCdata")}
)

methods::setGeneric(
  name="RunNormalization",
  def = function( object, ... ){standardGeneric("RunNormalization")}
)


methods::setGeneric(
  name="RunDiffAnalysis",
  def = function( object, ... ){standardGeneric("RunDiffAnalysis")}
)

methods::setGeneric(
  name="FilterDiffAnalysis",
  def = function( object, ... ){standardGeneric("FilterDiffAnalysis")}
)

methods::setGeneric(
  name="plotPCA",
  def = function( object, ... ){standardGeneric("plotPCA")}
)

methods::setGeneric(
  name="TransformData",
  def = function( object, ... ){standardGeneric("TransformData")}
)

methods::setGeneric(
  name="FilterLowAbundance",
  def = function( object, ... ){standardGeneric("FilterLowAbundance")}
)

methods::setGeneric(
  name="RunPCA",
  def = function( object, ... ){standardGeneric("RunPCA")}
)

methods::setGeneric(
  name="CheckExpDesignCompleteness",
  def = function( object, ... ){standardGeneric("CheckExpDesignCompleteness")}
)

methods::setGeneric(
  name="getExpressionContrast",
  def = function( object, ... ){standardGeneric("getExpressionContrast")}
)


methods::setGeneric(
  name="getContrastMatrix",
  def = function( object, ... ){standardGeneric("getContrastMatrix")}
)

methods::setGeneric(
  name="runAnnotationEnrichment",
  def = function( object, ... ){standardGeneric("runAnnotationEnrichment")}
)

methods::setGeneric(
  name="runCoExpression",
  def = function( object, ... ){standardGeneric("runCoExpression")}
)

methods::setGeneric(
  name="DiffAnal.plot",
  def = function( object, ... ){standardGeneric("DiffAnal.plot")}
)


methods::setGeneric(
  name="Data_Distribution_plot",
  def = function( object, ... ){standardGeneric("Data_Distribution_plot")}
)

methods::setGeneric(
  name="Library_size_barplot.plot",
  def = function( object, ... ){standardGeneric("Library_size_barplot.plot")}
)

methods::setGeneric(
  name="Enrichment.plot",
  def = function( object, ... ){standardGeneric("Enrichment.plot")}
)

methods::setGeneric(
  name="resetFlomicsMultiAssay",
  def = function( object, ... ){standardGeneric("resetFlomicsMultiAssay")}
)

methods::setGeneric(
  name="prepareForIntegration",
  def = function( object, ... ){standardGeneric("prepareForIntegration")}
)

methods::setGeneric(
  name="run_MOFA_analysis",
  def = function( object, ... ){standardGeneric("run_MOFA_analysis")}

)

methods::setGeneric(
  name="run_MixOmics_analysis",
  def = function( object, ... ){standardGeneric("run_MixOmics_analysis")}
  
)
