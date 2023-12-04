######################## COMMON METHODS FOR OMICS INTEGRATION ########################


#' @title integrationWrapper
#' @description This function executes all the steps to ensure data integration from a \link{MultiAssayExperiment} object produced by FLOMICS.
#' @param object An object of class \link{MultiAssayExperiment}. It is expected the MAE object is produced by rflomics previous analyses, as it relies on their results.
#' @param omicsToIntegrate vector of characters strings, referring to the names of the filtered table in 'object@ExperimentList'.
#' @param rnaSeq_transfo character string, only supports 'limma (voom)' for now. Transformation of the rnaSeq data from counts to continuous data.
#' @param choice character. If choice is set to 'DE', filters the object to take only the DE omics using differential analysis results stored in object. If choice is different than DE, no filtering is applied.
#' @param variableLists list of variables to keep per dataset.
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
                        variableLists = NULL,
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
    
    
    # TODO should do the intersection, not stopping everything !!
    if (any(!omicsToIntegrate %in% names(object))) {
      stop("There are omics to integrate that are not names from the object")
    }
    
    
    # TODO handle NULL case better than that
    if (isTagName(object, contrasts_names)) 
      contrasts_names <- convertTagToContrast(object, contrasts_names)
    
    if (cmd) print("#     => Preparing for multi-omics analysis")
    
    preparedObject <- prepareForIntegration(object = object,
                                            omicsToIntegrate = omicsToIntegrate,
                                            rnaSeq_transfo = rnaSeq_transfo,
                                            choice = choice,
                                            variableLists = variableLists,
                                            type = type,
                                            group = group,
                                            method = method,
                                            cmd = cmd,
                                            silent = silent)
    
    if (cmd) print("#     => run data integration")
    object <- runIntegration(object = object,
                             preparedObject = preparedObject,
                             method = method,
                             scale_views = scale_views,
                             maxiter = maxiter,
                             num_factors = num_factors,
                             selectedResponse = selectedResponse,
                             ncomp = ncomp,
                             link_datasets = link_datasets,
                             link_response = link_response,
                             sparsity = sparsity,
                             cases_to_try = cases_to_try,
                             silent = TRUE, cmd = FALSE)
    
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
#' @param variableLists list of variables to keep per dataset.
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
                        variableLists = NULL,
                        group = NULL,
                        method = c("MOFA", "MixOmics"),
                        cmd = FALSE, 
                        silent = TRUE) {
    
    method <- switch(toupper(method),
                     "MIXOMICS" = "MixOmics",
                     "MOFA"  = "MOFA",
                     "MOFA2" = "MOFA",
                     "MOFA+" = "MOFA"
    )
    
    # if no omicsToIntegrate we keep all SE
    if (is.null(omicsToIntegrate)) omicsToIntegrate <- names(object)
    
    # # Checking if intersection is a possible type of selection
    # if (type == "intersection") {
    #   allrownames <- lapply(omicsToIntegrate, FUN = function(nam){
    #     tryCatch(opDEList(object, SE.name = nam, contrasts = contrasts_names, operation = type),
    #              error = function(e) e, 
    #              warning = function(w) w)
    #   })
    #   names(allrownames) <- omicsToIntegrate
    #   
    #   if (any(sapply(allrownames, class) != "character")) {
    #     probOmics <- names(allrownames)[sapply(allrownames, class) != "character"]
    #     stop("It seems there  is a problem with: ", 
    #          sapply(probOmics, FUN = function(namesError){
    #            paste("\n", namesError, " -- error message:\n", allrownames[[namesError]])}))
    #   }
    #   
    # } # does this chunk of code work?!
    
    object <- object[, , omicsToIntegrate]
    
    # Checking for batch effects
    correct_batch <- FALSE
    ftypes <- getFactorTypes(object)
    
    if (any(ftypes == "batch")) {
      correct_batch <- TRUE
    }
    
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
        variableNames = variableLists[[SEname]], # variableLists[[SEname]]
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
    
    # Check for duplicated features names across tables
    # MOFA add the entire view name if it's the case, MixOmics do not care. 
    # For visualization purpose and coherence, add .index at the end of duplicated variables.
    commonVarNames <- sum(duplicated(unlist(rownames(object))))
    if (commonVarNames > 0) {
      if (cmd) cat("#   => Duplicated features names across tables, changing names for integration\n")
      
      dupTab <- data.frame("dataTable" = rep(names(object), time = sapply(object@ExperimentList, nrow)),
                           "rownames" = unlist(rownames(object)), 
                           "dup" = duplicated(unlist(rownames(object))) + 
                             duplicated(unlist(rownames(object)), fromLast = TRUE))
      dupTab <- dupTab[dupTab$dup == 1,]
      omicstochange <- unique(dupTab$dataTable)
      
      res <- lapply(1:length(object), FUN = function(i){
        SE.object <- object[[names(object)[i]]]
        if (names(object)[i] %in% omicstochange) {
          rownames(SE.object) <- paste(rownames(SE.object), i, sep = ".")
        }
        return(SE.object)
      })
      names(res) <- names(object)
      
      object <- MultiAssayExperiment(experiments = res, 
                                     colData = colData(object), 
                                     sampleMap = sampleMap(object), 
                                     metadata = metadata(object))
    }
    
    # keep only columns corresponding to design factors (remove samples and groups)
    #groups <- object@colData$groups
    object@colData <- object@colData[c(bioFactors(object), batchFactors(object), metaFactors(object))]
    
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


#' @title runIntegration
#' @description This function executes all the steps to ensure data integration from a \link{MultiAssayExperiment} object produced by FLOMICS.
#' @param object An object of class \link{MultiAssayExperiment}. It is expected the MAE object is produced by rflomics previous analyses, as it relies on their results.
#' @param preparedObject An untrained MOFA object or a list of dataset.
#' @param type one of union or intersection.
#' @param group Not implemented yet in the interface. Useful for MOFA2 run.
#' @return a MultiAssayExperiment object.
#' @export
#' @rdname runIntegration
#' @exportMethod runIntegration
methods::setMethod(
  f = "runIntegration",
  signature = "MultiAssayExperiment",
  definition = function(object,
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
                        ...) {
    
    method <- switch(toupper(method),
                     "MIXOMICS" = "MixOmics",
                     "MOFA"  = "MOFA",
                     "MOFA2" = "MOFA",
                     "MOFA+" = "MOFA"
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
      
      object@metadata[["IntegrationAnalysis"]][["MOFA"]][["setting"]] <- list(scale_views = scale_views,
                                                                              maxiter     = maxiter,
                                                                              num_factors = num_factors,
                                                                              selectData  = names(preparedObject@data))
      
      object@metadata[["IntegrationAnalysis"]][["MOFA"]][["MOFA_results"]]   <- MOFA_run$MOFAObject.trained
      object@metadata[["IntegrationAnalysis"]][["MOFA"]][["MOFA_untrained"]] <- MOFA_run$MOFAObject.untrained
      
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
      object@metadata[["IntegrationAnalysis"]][["mixOmics"]] <- MixOmics_res
      object@metadata[["IntegrationAnalysis"]][["mixOmics"]][["setting"]] <- list(scale_views      = scale_views,
                                                                                  ncomp.           = ncomp,
                                                                                  sparsity         = sparsity,
                                                                                  cases_to_try     = cases_to_try,
                                                                                  selectedResponse = selectedResponse,
                                                                                  selectData  = names(preparedObject$blocks))
    }
    return(object)
  }
)
