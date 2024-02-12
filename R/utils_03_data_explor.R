
#' 
#' 
#' #' plotDistr
#' #'
#' #' @param abundances matrix or dataframe of feature/gene abundances/counts
#' #' @export
#' #' @importFrom ggplot2 geom_density xlab
#' #' @importFrom reshape2 melt
#' #' @noRd
#' plotDistr <- function(abundances, dataType, transform_method){
#'   
#'   
#'   value <- samples <- NULL
#'   
#'   switch(dataType,
#'          "RNAseq" = {
#'            pseudo_counts <- log2(abundances+1) %>% reshape2::melt()
#'            colnames(pseudo_counts) <- c("features", "samples", "value")
#'            x_lab <-"log2(feature abundances)"
#'          },
#'          "proteomics"={
#'            pseudo_counts <- abundances %>% reshape2::melt()
#'            colnames(pseudo_counts) <- c("features", "samples", "value")
#'            x_lab <- switch (transform_method,
#'                             "log2"= "log2(feature abundances)",
#'                             "none"="feature abundances")
#'          },
#'          "metabolomics"={
#'            pseudo_counts <- abundances %>% reshape2::melt()
#'            colnames(pseudo_counts) <- c("features", "samples", "value")
#'            x_lab <- switch (transform_method,
#'                             "log2"= "log2(feature abundances)",
#'                             "none"="feature abundances")
#'          }
#'   )
#'   
#'   p <- ggplot2::ggplot(pseudo_counts) + geom_density(aes(value, color=samples) ) + xlab(x_lab) +
#'     theme(legend.position='none')
#'   print(p)
#'   
#' }



######## INTERNAL - Transform the Data ###########

# .applyTransformation: apply the transformation method stored in object@metadata[["transform_method"]] and modify the assay.
#' @title apply_transformation
#'
#' @param object An object of class \link{RflomicsSE}
#' @keywords internal
#' @noRd
#'

.applyTransformation <- function(object) {
  if (is.null(getTransSettings(object)$method)) {
    stop("Expect transformation method.")
  }
  
  if (.isTransformed(object)) {
    warning("Data were already transformed before!")
  }
  
  transform_method <- getTransSettings(object)$method
  assayTransform <- assay(object, withDimnames = TRUE)
  validTransform <- TRUE
  
  
  switch(transform_method,
         "log1p" = {
           assay(object) <- log1p(assayTransform)
         },
         "log2" = {
           assay(object) <- log2(assayTransform + 1)
         },
         "log10" = {
           assay(object) <- log10(assayTransform + 1)
         },
         "squareroot" = {
           assay(object) <- sqrt(assayTransform)
         },
         "none" = {
           assay(object) <- assayTransform
         },
         {
           assay(object) <- assayTransform
           message("Could not recognize the transformation method. No transformation applied. Please check your parameters.")
           validTransform <- FALSE
         } # default is none
  )
  
  if (transform_method != "none" && validTransform) object@metadata[["DataProcessing"]][["Transformation"]][["transformed"]] <- TRUE
  
  return(object)
}

######## INTERNAL - Normalize the Data ###########

# .applyNorm: apply the normalization method stored in object@metadata[["Normalization"]] and modify the assay.
#' @title .applyNorm
#'
#' @param object An object of class \link{RflomicsSE}
#' @description apply the normalization to the assay. Usually, after the transformation,
#' unless in the case of counts RNASeq data (TMM), where log2 is the second step.
#' @keywords internal
#' @noRd
#'

.applyNorm <- function(object) {
  if (is.null(getNormSettings(object)$method)) {
    stop("Expects normalization method.")
  }
  
  if (.isNorm(object)) {
    warning("Data were already normalized before!")
  }
  
  norm_method <- getNormSettings(object)$method
  coefNorm <- getCoeffNorm(object)
  validNorm <- TRUE
  
  assayTransform <- assay(object)
  
  switch(norm_method,
         "median" = {
           assay(object) <- sweep(assayTransform, 2, coefNorm, "-")
         },
         "totalSum" = {
           assay(object) <- sweep(assayTransform, 2, coefNorm, "/")
         },
         "TMM" = {
           scales_factors <- coefNorm$norm.factors * coefNorm$lib.size
           assay(object) <- scale(assayTransform + 1, center = FALSE, scale = scales_factors)
         },
         "none" = {
           assay(object) <- assayTransform
         },
         { # default is none
           assay(object) <- assayTransform
           message("Could not recognize the normalization method. No normalization applied. Please check your parameters.")
           validNorm <- FALSE
         }
  )
  
  if (norm_method != "none" && validNorm) object@metadata[["DataProcessing"]][["Normalization"]][["normalized"]] <- TRUE
  
  return(object)
}

######## INTERNAL - Check, transform and normalize the data ###########

# checkTransNorm: check the data, transform them and normalize them.
#' @title .checkTransNorm
#'
#' @param object An object of class \link{RflomicsSE}
#' @description apply the normalization and the transformation stored into the metadata of the SE object.
#'  Applies TMM and log2 transformation for RNAseq data.
#' @keywords internal
#' @noRd
#'

.checkTransNorm <- function(object, raw = FALSE) {
  if (!is(object, "RflomicsSE")) stop("Object is not a RflomicsSE")
  
  # check things
  if (.checkNA(object)) stop("NA detected in the assay.")
  
  # No transformation (except for RNAseq, which expect counts...)
  if (raw) {
    if (.isTransformed(object)) warning("Your data are not raw (transformed)")
    if (.isNorm(object)) warning("Your data are not raw (normalized)")
    
    if (getOmicsTypes(object) == "RNAseq") {
      assay(object) <- log2(assay(object) + 1)
    }
  } else {
    # if RNAseq
    switch(getOmicsTypes(object),
           "RNAseq" = {
             # Really depends if TMM is the normalization or not.
             # Make it easier: force TMM and log2.
             if (.isTransformed(object)) stop("Expect untransformed RNAseq data at this point.")
             
             if (.isNorm(object)) {
               switch(getNormSettings(object)$method, 
                      "TMM" = {assay(object) <- log2(assay(object))}, # +1 in the apply_norm function
                      {message("RNAseq counts expects TMM normalization. Data were already normalized with another method.
                Skipping to the end without transforming or normalizing data.")}
               )
             } else {
               # Force "none" transformation.
               if (getTransSettings(object)$method != "none") {
                 message("RNAseq counts expects TMM normalization. Transformation is done after the normalization,
                  using 'none' as transform method. Data will be transformed using log2 after the normalization anyway")
                 
                 object <- setTrans(object, method = "none")
               }
               
               # Force TMM normalization
               if (getNormSettings(object)$method != "TMM") {
                 message("For RNAseq data (counts), only TMM applies for now. Forcing TMM normalization.")
                 object <- runNormalization(object, normMethod = "TMM")
               }
             }
             
             # Finally transforming the data.
             object <- .applyTransformation(object) # none
             object <- .applyNorm(object) # TMM
             assay(object) <- log2(assay(object)) # +1 in the .applyNorm function
             
           }, # end switch rnaseq
           { # default
             # in case any other omics type (does not expect counts)
             # transform and norm
             if (!.isTransformed(object)) object <- .applyTransformation(object)
             if (!.isNorm(object)) object <- .applyNorm(object)
           }
    )
  }
  
  # check things
  if (.checkNA(object)) stop("NA detected in the assay.")
  
  return(object)
}

######## INTERNAL - isNorm, isTransform, getNorm, getTransform ###########

#' @title isNorm, isTransform, getNorm, getTransform, setNorm, setTrans
#'
#' @param object An object of class \link{RflomicsSE}
#' @description get if an assay has been transformed or normalized.
#' @keywords internal
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @noRd
#'

.isTransformed <- function(object) {
  metadata(object)[["DataProcessing"]][["Transformation"]][["transformed"]]
}

.isNorm <- function(object) {
  metadata(object)[["DataProcessing"]][["Normalization"]][["normalized"]]
}


######## INTERNAL - CHECKS FUNCTIONS ###########

#  .checkNA: checks if there are NA/nan in the RflomicsSE assay
#' @title .checkNA
#'
#' @param object An object of class \link{RflomicsSE}
#' @return boolean. if TRUE, NA/nan are detected in the SE::assay.
#' @keywords internal
#' @noRd
#'
.checkNA <- function(object) {
  NA_detect <- ifelse(any(is.na(assay(object))), TRUE, FALSE)
  return(NA_detect)
}

######## INTERNAL

#' @title .tmmNormalization
#' Interface to the calcNormFactors function of the edgeR package  with the choosen TMM parameters as the normalization method
#' @param counts numeric matrix of read counts
#' @param groups vector or factor giving the experimental group/condition for each sample/library.
#' @return a data.frame with a row for each sample and columns group, lib.size and norm.factors containing the group labels, library sizes and normalization factors. Other columns can be optionally added to give more detailed sample information.
#' @keywords internal
#' @importFrom edgeR DGEList calcNormFactors
#' @noRd

.tmmNormalization <- function(counts, groups){
  dge <- edgeR::DGEList(counts=counts, group=groups)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  nf  <- dge$samples
  return(nf)
}