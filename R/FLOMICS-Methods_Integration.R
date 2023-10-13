######################## COMMON METHODS FOR OMICS INTEGRATION ########################


#' @title integrationWrapper
#' @description This function executes all the steps to ensure data integration from a \link{MultiAssayExperiment} object produced by FLOMICS.
#' @param object An object of class \link{MultiAssayExperiment}. It is expected the MAE object is produced by rflomics previous analyses, as it relies on their results.
#' @param omicsToIntegrate vector of characters strings, referring to the names of the filtered table in 'object@ExperimentList'.
#' @param rnaSeq_transfo character string, only supports 'limma (voom)' for now. Transformation of the rnaSeq data from counts to continuous data.
#' @param choice character. If choice is set to 'DE', filters the object to take only the DE omics using differential analysis results stored in object. If choice is different than DE, no filtering is applied.
#' @param contrasts_names contrasts names for the selection of DE entities.
#' @param type one of union or intersection.
#' @param group Not implemented yet in the interface. Useful for MOFA2 run.
#' @return a MultiAssayExperiment object.
#' @export
#' @rdname integrationWrapper
#' @exportMethod integrationWrapper
methods::setMethod(
  f = "integrationWrapper",
  signature = "MultiAssayExperiment",
  definition = function(object,
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
                        cmd = FALSE,
                        ...) {
    
    method <- switch(toupper(method),
                     "MIXOMICS" = "MixOmics",
                     "MOFA"  = "MOFA",
                     "MOFA2" = "MOFA",
                     "MOFA+" = "MOFA"
    )
    
    if (any(!omicsToIntegrate %in% names(object))) {
      stop("There are omics to integrate that are not names from the object")
    }
    
    # TODO handle NULL case better than that
    if (isTagName(object, contrasts_names)) 
      contrasts_names <- convertTagToContrast(object, contrasts_names)
    
    if (cmd) print("#     => Preparing for multi-omics analysis")
    
    preparedObject <- prepareForIntegration(
      object = object,
      omicsToIntegrate = omicsToIntegrate,
      rnaSeq_transfo = rnaSeq_transfo,
      choice = choice,
      contrasts_names = contrasts_names,
      type = type,
      group = group,
      method = method,
      cmd = cmd, 
      silent = silent
    )
    
    if (toupper(method) == "MOFA") {
      object@metadata[["MOFA"]] <- NULL
      
      if (cmd) print("#     => Running MOFA analysis")
      
      MOFA_run <- runMOFAAnalysis(
        object = preparedObject,
        scale_views = scale_views,
        maxiter = maxiter,
        num_factors = num_factors,
        silent = silent
      )
      
      
      object@metadata[["MOFA"]][["MOFA_results"]] <- MOFA_run$MOFAObject.trained
      object@metadata[["MOFA"]][["MOFA_untrained"]] <- MOFA_run$MOFAObject.untrained
      object@metadata[["MOFA"]][["MOFA_selected_filter"]] <- type
      object@metadata[["MOFA"]][["MOFA_selected_contrasts"]] <- contrasts_names
    } else if (toupper(method) == "MIXOMICS") {
      object@metadata[["mixOmics"]] <- NULL
      
      if (cmd) print("#     => Running mixOmics analysis")
      
      if (is.null(selectedResponse)) selectedResponse <- bioFactors(object)
      
      if (silent) {
        co <- capture.output({ 
          MixOmics_res <- lapply(selectedResponse,
                                 FUN = function(response_var) {
                                   res_mixOmics <- suppressWarnings(runMixOmicsAnalysis(
                                     object = preparedObject,
                                     selectedResponse = response_var,
                                     scale_views = scale_views,
                                     ncomp = ncomp,
                                     link_datasets = link_datasets,
                                     link_response = link_response,
                                     sparsity = sparsity,
                                     cases_to_try = cases_to_try
                                   ))
                                   
                                   return(
                                     list(
                                       "MixOmics_tuning_results" = res_mixOmics$tuning_res,
                                       "MixOmics_results"        = res_mixOmics$analysis_res
                                     )
                                   )
                                 }
          )
        })
        
      } else {
        MixOmics_res <- lapply(selectedResponse,
                               FUN = function(response_var) {
                                 res_mixOmics <- runMixOmicsAnalysis(
                                   object = preparedObject,
                                   selectedResponse = response_var,
                                   scale_views = scale_views,
                                   ncomp = ncomp,
                                   link_datasets = link_datasets,
                                   link_response = link_response,
                                   sparsity = sparsity,
                                   cases_to_try = cases_to_try
                                 )
                                 
                                 return(
                                   list(
                                     "MixOmics_tuning_results" = res_mixOmics$tuning_res,
                                     "MixOmics_results"        = res_mixOmics$analysis_res
                                   )
                                 )
                               }
        )
      }
      
      names(MixOmics_res) <- selectedResponse
      object@metadata[["mixOmics"]] <- MixOmics_res
    }
    
    return(object)
  }
)


#' @title prepareForIntegration
#' @description This function transforms a MultiAssayExperiment produced by rflomics into an untrained MOFA objects or a list to use for mixOmics. 
#' It checks for batch effect to correct them prior to the integration.
#' It also transforms RNASeq counts data into continuous data. This is the first step into the integration.
#' @param object An object of class \link{MultiAssayExperiment}. It is expected the MAE object is produced by rflomics previous analyses, as it relies on their results.
#' @param omicsToIntegrate vector of characters strings, referring to the names of the filtered table in 'object@ExperimentList'.
#' @param rnaSeq_transfo character string, only supports 'limma (voom)' for now. Transformation of the rnaSeq data from counts to continuous data.
#' @param choice character. If choice is set to 'DE', filters the object to take only the DE omics using differential analysis results stored in object. If choice is different than DE, no filtering is applied.
#' @param contrasts_names contrasts names for the selection of DE entities.
#' @param type one of union or intersection.
#' @param group Not implemented yet in the interface. Useful for MOFA2 run.
#' @return An untrained MOFA object or a list of dataset
#' @export
#' @exportMethod prepareForIntegration
#' @importFrom MOFA2 create_mofa
methods::setMethod(
  f = "prepareForIntegration",
  signature = "MultiAssayExperiment",
  definition = function(object,
                        omicsToIntegrate = NULL,
                        rnaSeq_transfo = "limma (voom)",
                        choice = "DE",
                        contrasts_names = NULL,
                        type = "union",
                        group = NULL,
                        method = c("MOFA", "MixOmics"),
                        cmd = FALSE, 
                        silent = TRUE) {
    if (is.null(omicsToIntegrate)) omicsToIntegrate <- names(object)
    
    # Checking for batch effects
    correct_batch <- FALSE
    ftypes <- getFactorTypes(object)
    
    if (any(ftypes == "batch")) {
      correct_batch <- TRUE
    }
    
    object <- object[, , omicsToIntegrate]
    
    # On each selected omics, according to its type, apply transformation if demanded.
    # Filter DE entities
    # TODO : add a possibility of choice for keeping every entity (small tables)
    # TODO change for lapply
    for (SEname in omicsToIntegrate) {
      SEobject <- object[[SEname]]
      omicsType <- getOmicsTypes(SEobject)
      
      list_args <- list(
        object = object,
        SEname = SEname,
        correctBatch = correct_batch,
        contrasts_names = contrasts_names,
        type = type,
        choice = choice,
        cmd = cmd
      )
      
      object <- switch(omicsType,
                       "RNAseq" = {
                         list_args$transformation <- rnaSeq_transfo
                         do.call("rnaseqRBETransform", list_args)
                       },
                       "proteomics" = do.call("RBETransform", list_args),
                       "metabolomics" = do.call("RBETransform", list_args)
      )
    }
    
    if (method == "MOFA") {    
      if (silent) {
        MOFAObject <- suppressMessages(
          suppressWarnings(create_mofa(object,
                                              group = group,
                                              extract_metadata = TRUE)))
      } else {
        MOFAObject <- create_mofa(object,
                                         group = group,
                                         extract_metadata = TRUE)
      }
      
      return(MOFAObject)
      
    } else if (method == "MixOmics") { 
      
      # Common samples names:
      nsamp <- nrow(colData(object))
      object <- intersectColumns(object)
      
      if (nsamp != nrow(colData(object))) {
        warning("Removing ", nsamp - nrow(colData(object)), " samples not present in every experiment.")
      }
      
      MixOmicsObject <- list(
        blocks   = lapply(object@ExperimentList, FUN = function(SE) t(assay(SE))),
        metadata = object@colData
      )
      
      MixOmicsObject$blocks <- lapply(MixOmicsObject$blocks, 
                                      FUN = function(mat) mat[match(rownames(mat), rownames(MixOmicsObject$metadata)), ])
      
      return(MixOmicsObject)
    }
  }
)


