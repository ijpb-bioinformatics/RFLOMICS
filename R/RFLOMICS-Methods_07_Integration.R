######################## METHODS FOR OMICS INTEGRATION ########################

# ---- Wrapper ----
#' @title Wrapper for integration of omics data using RFLOMICS
#' @description This function executes all the steps to ensure data integration
#' from a \link{RflomicsMAE} object produced by FLOMICS. It encapsulates the
#' three other functions: \code{\link{filterFeatures,RflomicsMAE-method}}, 
#' \code{\link{prepareForIntegration,RflomicsMAE-method}}
#' and \code{\link{runOmicsIntegration,RflomicsMAE-method}} otherwise necessary
#' to complete the integration.
#' @param object An object of class \link{RflomicsMAE}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param omicsNames vector of characters strings,
#' referring to the names of the filtered table in 'object@ExperimentList'.
#' @param rnaSeq_transfo character string, only supports 'limma (voom)'
#' for now.
#' Transformation of the rnaSeq data from counts to continuous data.
#' @param selOpt list of selection option for each dataset: one of none, DE 
#' or a vector of contrasts or cluster names.
#' @param type one of union or intersection. 
#' @param group Not implemented yet in the interface. Useful for MOFA2 run.
#' @param method one of MOFA or mixOmics
#' @param scale_views boolean. If TRUE, each dataset is scaled.
#' @param maxiter MOFA2 parameter. Number of maximum iteration to use.
#' @param num_factors MOFA2 parameter. Number of factor to compute. 
#' @param selectedResponse vector of character. Response variables for mixOmics
#' @param ncomp mixOmics parameter. Number of component to compute. 
#' @param link_datasets mixOmics parameter. Link between datasets in the design.
#' @param link_response mixOmics parameter. Link between dataset and response.
#' @param sparsity mixOmics parameter. If TRUE, uses block.splsda.
#' @param cases_to_try used for tuning when sparse analysis is TRUE. 
#' @param silent silence all functions. 
#' @param cmd used in the interface, print cmd lines.
#' @return a RflomicsMAE object. According to the method (MOFA or mixOmics),
#' the correct slot of metadata has been filled with the results and the
#' settings.
#' @noRd
#' @keywords internal
setMethod(
    f = "integrationWrapper",
    signature = "RflomicsMAE",
    definition = function(object,
                          omicsNames = names(object),
                          rnaSeq_transfo = "limma (voom)",
                          selOpt = rep(list("DE"), length(omicsNames)),
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
                          ...) {
        if (any(!omicsNames %in% names(object))) {
            stop("There are omics to integrate that are not names from the object")
        }
        
        if (cmd)
            message("#     => Preparing for multi-omics analysis")
        
        objectfilt <- filterFeatures(object = object,
                                     selOpt = selOpt,
                                     type = type)
        
        variableLists <- lapply(experiments(objectfilt), names)
        
        preparedObject <- prepareForIntegration(
            object = object,
            omicsNames = omicsNames,
            rnaSeq_transfo = rnaSeq_transfo,
            variableLists = variableLists,
            group = group,
            method = method,
            cmd = cmd,
            silent = silent
        )
        
        if (cmd)
            message("#     => run data integration")
        object <- runOmicsIntegration(
            object = object,
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
            silent = TRUE,
            cmd = FALSE
        )
        
        return(object)
    }
)

# ---- prepareForIntegration ----
#' @title Preparation step for integration
#' @description This function transforms a RflomicsMAE produced by rflomics
#' into an untrained MOFA object or a list to use for mixOmics.
#' It checks for batch effects to correct them before integration.
#' It also transforms RNASeq counts data into continuous data using 
#' \code{\link[limma]{voom}}.
#' This is the second step into the integration. It is usually preceded by
#' \link{filterFeatures} to extract the correct variables,
#'  and followed by \link{runOmicsIntegration}.
#' @param object An object of class \link{RflomicsMAE}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param omicsNames vector of characters strings,
#' referring to the names of the filtered table in 'object@ExperimentList'.
#' @param rnaSeq_transfo character string, only supports 'limma (voom)'
#' for now.
#' Transformation of the rnaSeq data from counts to continuous data.
#' @param variableLists list of variables to keep per dataset.
#' @param group Not implemented yet in the interface. Useful for MOFA2 run.
#' @param method one of MOFA or mixOmics. 
#' Method for which the object is prepared.
#' @param transformData boolean. 
#' Transform the data with the transform and normalization method?
#' Default is TRUE.
#' @param cmd used in the interface. Print cmd lines.
#' @param silent if TRUE, silence all functions.
#' @return An untrained MOFA object or a list of dataset
#' @exportMethod prepareForIntegration
#' @rdname prepareForIntegration
#' @importFrom MOFA2 create_mofa
#' @examples
#' MAEtest <- generateExample(annotation = FALSE, coexp = FALSE,
#'                            integration = FALSE)
#' mofaObj <- prepareForIntegration(MAEtest,
#'                 omicsNames = c("protetest", "metatest"),
#'                 variableLists = rownames(MAEtest),
#'                 method = "MOFA")
#'
#' mixOmicsList <- prepareForIntegration(MAEtest,
#'                   omicsNames = c("protetest", "metatest"),
#'                   variableLists = rownames(MAEtest),
#'                   method = "mixOmics")
#'
setMethod(
    f = "prepareForIntegration",
    signature = "RflomicsMAE",
    definition = function(object,
                          omicsNames = NULL,
                          rnaSeq_transfo = "limma (voom)",
                          variableLists = NULL,
                          group = NULL,
                          method = "MOFA",
                          transformData = TRUE,
                          cmd = FALSE,
                          silent = TRUE) {
        method <- switch(
            toupper(method),
            "MIXOMICS" = "MixOmics",
            "MOFA"  = "MOFA",
            "MOFA2" = "MOFA",
            "MOFA+" = "MOFA"
        )
        
        # if no omicsNames we keep all SE
        if (is.null(omicsNames))
            omicsNames <- names(object)
        
        object <- object[, , omicsNames]
        
        # Checking for batch effects
        correct_batch <- FALSE
        ftypes <- getFactorTypes(object)
        
        if (any(ftypes == "batch")) {
            correct_batch <- TRUE
        }
        
        # Transformation before anything else, except for RNAseq data.
        if (transformData) {
            for (SEname in omicsNames) {
                if (getOmicsTypes(object[[SEname]]) != "RNAseq") {
                    object[[SEname]] <- .checkTransNorm(object[[SEname]])
                }
            }
        }
        
        # On each selected omics, according to its type,
        # apply transformation if demanded.
        # Filter DE entities
        for (SEname in omicsNames) {
            SEobject <- object[[SEname]]
            omicsType <- getOmicsTypes(SEobject)
            
            list_args <- list(
                object = object,
                SEname = SEname,
                correctBatch = correct_batch,
                variableNames = variableLists[[SEname]],
                cmd = cmd
            )
            
            object <- switch(
                omicsType,
                "RNAseq" = {
                    list_args$transformation <- rnaSeq_transfo
                    do.call(".rnaseqRBETransform", list_args)
                },
                "proteomics" = do.call(".rbeTransform", list_args),
                "metabolomics" = do.call(".rbeTransform", list_args)
            )
        }
        
        # Check for duplicated features names across tables
        # MOFA add the entire view name if it's the case, MixOmics do not care.
        # For visualization purpose and coherence,
        # add .index at the end of duplicated variables.
        commonVarNames <- sum(duplicated(unlist(rownames(object))))
        if (commonVarNames > 0) {
            if (cmd) {
                message("#   => Duplicated features names across tables,
                changing names for integration")
            }
            
            dupTab <- data.frame(
                "dataTable" = rep(names(object),
                                  time = vapply(experiments(object), nrow, c(1))),
                "rownames" = unlist(rownames(object)),
                "dup" = duplicated(unlist(rownames(object))) +
                    duplicated(unlist(rownames(object)), fromLast = TRUE)
            )
            dupTab <- dupTab[which(dupTab$dup == 1), ]
            omicstochange <- unique(dupTab$dataTable)
            
            res <- lapply(
                seq_len(length(object)),
                FUN = function(i) {
                    SE.object <- object[[names(object)[i]]]
                    if (names(object)[i] %in% omicstochange) {
                        rownames(SE.object) <- paste(rownames(SE.object), i, sep = ".")
                    }
                    return(SE.object)
                }
            )
            names(res) <- names(object)
            
            object <- RflomicsMAE(
                experiments = res,
                colData   = colData(object),
                sampleMap = sampleMap(object),
                omicList    = metadata(object)$omicList, 
                projectName = getProjectName(object), 
                design      = metadata(object)$design
            )
            
            if (cmd) {
                message("#   => Done replacing features names")
            }
        }
        
        # keep only columns corresponding to design factors
        # (remove samples and groups)
        colData(object) <- colData(object)[c(getBioFactors(object),
                                             getBatchFactors(object),
                                             getMetaFactors(object))]
        
        if (method == "MOFA") {
            if (silent) {
                MOFAObject <- suppressMessages(suppressWarnings(
                    create_mofa(object,
                                groups = group,
                                extract_metadata = TRUE)
                ))
            } else {
                MOFAObject <- create_mofa(object,
                                          groups = group,
                                          extract_metadata = TRUE)
            }
            return(MOFAObject)
            
        } else if (method == "MixOmics") {
            # Common samples names:
            nsamp <- nrow(colData(object))
            object <- intersectColumns(object)
            
            if (nsamp != nrow(colData(object))) {
                warning("Removing ",
                        nsamp - nrow(colData(object)),
                        " samples not present in every experiment.")
            }
            
            MixOmicsObject <- list(blocks   = lapply(
                experiments(object),
                FUN = function(SE) {
                    t(assay(SE))
                }
            ), metadata = colData(object))
            
            MixOmicsObject$blocks <- lapply(
                MixOmicsObject$blocks,
                FUN = function(mat) {
                    mat[match(rownames(mat), rownames(MixOmicsObject$metadata)),]
                }
            )
            
            
            return(MixOmicsObject)
        }
    }
)




# ---- Select features to  keep for integration (MAE) ----
#
#' @title Feature selection in a Rflomics MAE
#' @description This function selects all the features to keep according to
#' user's choices on each omic data.
#' @param object An object of class \link{RflomicsMAE-class}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param selOpt list of vectors. Preferred named list with names corresponding
#' to the names of the experimentList in the object. For each Experiment, gives
#' the option for the filtering: either 'all', 'DE', 'none',
#' or a specific name of a
#' contrast or cluster (if coexpression results are available). Default is
#' taking all features for all experiment list. If the vector is named and an
#' Experiment is missing, no feature will be selected from it.
#' If the vector is not named, the selection will be applied in order of the
#' Experiments in the object.
#' @param type if selOpt is set to a specific set of contrasts or clusters,
#' indicates whether the selection is "union" or "intersection" of entities in
#' these sets.
#' @return a RflomicsMAE, filtered with only the corresponding features.
#' @rdname filterFeatures
#' @exportMethod filterFeatures
#' @examples
#' MAE <- generateExample(integration = FALSE, annotation = FALSE)
#'
#' selOpt = list("RNAtest" = c("cluster.1", "H3"), protetest = c("DE"))
#' MAE1 <- filterFeatures(MAE, selOpt)
#'
#' selOpt2 = list("RNAtest" = c("cluster.2", "H2", "H1"), protetest = c("DE"))
#' MAE2 <- filterFeatures(MAE, selOpt2,
#' type = c("RNAtest" = "intersection",
#' "protetest" = 'union'))
#'
setMethod(
    f = "filterFeatures",
    signature = "RflomicsMAE",
    definition = function(object,
                          selOpt = rep(list("all"), length(object)),
                          type = rep(list("union"), length(selOpt))) {
        # check if named vector
        if (is.null(names(selOpt))) {
            names(selOpt) <- names(object)[seq_along(selOpt)]
        }
        
        if (is.null(names(type))) {
            names(type) <- names(selOpt)[seq_along(selOpt)]
        }
        
        # Applying corresponding filtering
        res <- lapply(
            names(selOpt),
            FUN = function(nam) {
                SE.object <- object[[nam]]
                vectSel <- selOpt[[nam]]
                typeSel <- type[[nam]]
                
                if (is.null(typeSel))
                    typeSel <- "union"
                
                # intermediate list: list of features for each selectiontype
                # default: take all if typo somewhere.
                resInter <- lapply(
                    vectSel,
                    FUN = function(listSel) {
                        if (!listSel %in% c("all", "none", "DE")) {
                            originList <- .getOrigin(SE.object, listSel)
                        } else {
                            originList <- listSel
                        }
                        
                        switch(
                            originList,
                            "all" = {
                                names(SE.object)
                            },
                            "none" = {
                                NULL
                            },
                            "DE"  = {
                                getDEMatrix(object = SE.object)$DEF
                            },
                            "Contrast" = {
                                getDEList(object = SE.object, contrasts = listSel)
                            },
                            "Tag" = {
                                getDEList(object = SE.object, contrasts = listSel)
                                # TODO problem when only one selected
                            },
                            "CoexCluster" = {
                                getClusterEntities(SE.object, clusterName = listSel)
                            },
                            {
                                # Default: all
                                message("Cannot detect origin of ",
                                        listSel,
                                        " for ",
                                        nam ,
                                        " taking all features")
                                names(SE.object)
                            }
                        )
                        
                    }
                )
                
                # Union or intersection of all selected features
                filtKeep <- if (!is.null(typeSel)) {
                    switch(typeSel,
                           unique(unlist(resInter)),
                           # default
                           "intersection" = Reduce(intersect, resInter))
                }
                if (length(filtKeep) == 0) {
                    message("No feature to keep in ", nam,
                            ", it will be dropped.")
                }
                return(SE.object[filtKeep, ])
                
            }
        )
        names(res) <- names(selOpt)
        
        res <- res[lengths(res) > 0]
        
        return(
            RflomicsMAE(
                experiments = res,
                colData = colData(object),
                sampleMap = sampleMap(object),
                omicList    = metadata(object)$omicList, 
                projectName = getProjectName(object), 
                design      = metadata(object)$design
            )
        )
        
    }
)


# ---- Run Omics integration ----

#' @title runOmicsIntegration
#' @description Runs the integration according to the selected method (MOFA or
#' mixOmics) and the settings given by the user. Requires to have the correct
#' entry format in preparedObject before running.
#' @param object An object of class \link{RflomicsMAE}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param preparedObject An untrained MOFA object or a list of dataset.
#' Usually a result of \link{prepareForIntegration}.
#' @param method one of MOFA or mixOmics. 
#' Method for which the object is prepared.
#' @param scale_views boolean. If TRUE, scale each dataset to unit variance.
#' @param maxiter MOFA2 parameter. 
#' Number of max iteration (otherwise stop when converged.)
#' @param num_factors MOFA2 parameter. The number of factor to compute.
#' @param selectedResponse character vector, used for mixOmics. 
#' Response variables names for block.(s)plsda.
#' @param ncomp mixOmics parameter. Number of components to compute.
#' @param link_datasets mixOmics parameter. 
#' Link between datasets in the computation.
#' @param link_response mixOmics parameter. Link between dataset and response.
#' @param sparsity boolean. Used to determine which mixOmics function to apply (either
#' block.plsda if FALSE or block.splsda if TRUE).
#' @param cases_to_try integer. If sparsity is set to TRUE, then cases_to_try
#' is used to determine the number of sets of variables to test for tuning.
#' @param silent boolean. If TRUE, silence all functions.
#' @param cmd boolean. Used in the interface. If TRUE, print cmd in the console.
#' @param ... not in use at the moment
#' @return a RflomicsMAE object with the correct metadata slot filled with the
#' results and the settings.
#' @rdname runOmicsIntegration
#' @exportMethod runOmicsIntegration
#' @examples
#' # Generate MAE for test:
#' MAEtest <- generateExample(annotation = FALSE, coexp = FALSE,
#'                            integration = FALSE)
#'                            
#' # Prepare mofa object:
#' mofaObj <- prepareForIntegration(MAEtest,
#'                 omicsNames = c("protetest", "metatest"),
#'                 variableLists = rownames(MAEtest),
#'                 method = "MOFA")
#'
#' # Perform integration:
#' MAEtest <- runOmicsIntegration(MAEtest, mofaObj, method = "MOFA")
#' MOFA2::plot_variance_explained(getMOFA(MAEtest))
#'
setMethod(
    f = "runOmicsIntegration",
    signature = "RflomicsMAE",
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
        method <- switch(
            toupper(method),
            "MIXOMICS" = "MixOmics",
            "MOFA"  = "MOFA",
            "MOFA2" = "MOFA",
            "MOFA+" = "MOFA"
        )
        
        if (toupper(method) == "MOFA") {
            # object@metadata[["MOFA"]] <- NULL
            object <- setMOFA(object, NULL)
            
            if (cmd)
                message("#     => Running MOFA analysis")
            
            MOFA_run <- .runMOFAAnalysis(
                object = preparedObject,
                scale_views = scale_views,
                maxiter = maxiter,
                num_factors = num_factors,
                silent = silent
            )
            
            object <- setMOFA(
                object,
                list(
                    "MOFA_results" = MOFA_run$MOFAObject.trained,
                    "MOFA_untrained" = MOFA_run$MOFAObject.untrained,
                    "settings" = list(
                        scale_views = scale_views,
                        maxiter     = maxiter,
                        num_factors = num_factors,
                        selectData  = names(preparedObject@data)
                    )
                )
            )
            
        } else if (toupper(method) == "MIXOMICS") {
            # metadata(object)[["mixOmics"]] <- NULL
            object <- setMixOmics(object, NULL)
            
            if (cmd)
                message("#     => Running mixOmics analysis")
            
            if (is.null(selectedResponse))
                selectedResponse <- getBioFactors(object)
            
            if (silent) {
                co <- capture.output({
                    MixOmics_res <- lapply(
                        selectedResponse,
                        FUN = function(response_var) {
                            res_mixOmics <- suppressWarnings(
                                .runMixOmicsAnalysis(
                                    object = preparedObject,
                                    selectedResponse = response_var,
                                    scale_views = scale_views,
                                    ncomp = ncomp,
                                    link_datasets = link_datasets,
                                    link_response = link_response,
                                    sparsity = sparsity,
                                    cases_to_try = cases_to_try
                                )
                            )
                            
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
                MixOmics_res <- lapply(
                    selectedResponse,
                    FUN = function(response_var) {
                        res_mixOmics <- .runMixOmicsAnalysis(
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
            MixOmics_res$settings <- list(
                scale_views      = scale_views,
                ncomp            = ncomp,
                sparsity         = sparsity,
                cases_to_try     = cases_to_try,
                selectedResponse = selectedResponse,
                selectData  = names(preparedObject$blocks)
            )
            
            object <- setMixOmics(object, MixOmics_res)
            
        }
        return(object)
    }
)

# ----  Get a particular multi-omics result ----
#
#' @title Get a particular multi-omics result or settings.
#' @description
#' These methods are used to directly access the results of multi-omics 
#' analyses or their settings, usually stored in the metadata of the 
#' \link{RflomicsMAE} object. Setters are also available.
#'
#' @param object An object of class \link{RflomicsMAE}.
#' It is expected the MAE object is produced by rflomics previous analyses,
#' as it relies on their results.
#' @param response a character giving the response variable to access
#' specifically.
#' @param onlyResults default return only the MixOmics or MOFA2 results.
#' If you want to access all information of the integration,
#' set onlyResults to FALSE.
#' In MixOmics case, works only when response is specified.
#' @return
#' For getters:
#' in getMixOmics, if response is NULL,
#' then all the mixOmics results are returned.
#' Otherwise, it gives the particular mixOmics result.
#' For MOFA, returns the untrained object and the trained object as a list.
#'
#' For setters: always returns a \link{RflomicsMAE} object.
#'
#' @exportMethod getMixOmics
#' @rdname methods-for-integration
#' @examples
#'
#' MAEtest <- generateExample(annotation = FALSE, coexp = FALSE)
#'
#' # Access mixOmics results:
#' getMixOmics(MAEtest, response = "temperature")
#' getMixOmicsSettings(MAEtest)
#' mixOmics::plotIndiv(getMixOmics(MAEtest, response = "imbibition"))
#'
#' # Access MOFA2 results:
#' getMOFA(MAEtest)
#' getMOFASettings(MAEtest)
#' MOFA2::plot_variance_explained(getMOFA(MAEtest))
#'
setMethod(
    f = "getMixOmics",
    signature = "RflomicsMAE",
    definition = function(object,
                          response = NULL,
                          onlyResults = TRUE) {
        toreturn <- metadata(object)[["IntegrationAnalysis"]][["mixOmics"]]
        
        if (is.null(toreturn)) {
            return(toreturn)
        }
        
        if (!is.null(response)) {
            toreturn <- toreturn[[response]]
            if (onlyResults)
                toreturn <- toreturn$MixOmics_results
            return(toreturn)
        } else{
            return(toreturn)
        }
        
    }
)


#' @rdname methods-for-integration
#' @exportMethod getMOFA
setMethod(
    f = "getMOFA",
    signature = "RflomicsMAE",
    definition = function(object, onlyResults = TRUE) {
        toreturn <- metadata(object)[["IntegrationAnalysis"]][["MOFA"]]
        if (onlyResults && !is.null(toreturn)) {
            toreturn <- toreturn[["MOFA_results"]]
        }
        
        return(toreturn)
    }
)
# ---- Get integration setting ----

#' @exportMethod getMOFASettings
#' @rdname methods-for-integration
setMethod(
    f = "getMOFASettings",
    signature = "RflomicsMAE",
    definition = function(object) {
        return(getMOFA(object, onlyResults = FALSE)$settings)
    }
)

#' @exportMethod getMixOmicsSettings
#' @rdname methods-for-integration
setMethod(
    f = "getMixOmicsSettings",
    signature = "RflomicsMAE",
    definition = function(object) {
        return(metadata(object)[["IntegrationAnalysis"]][["mixOmics"]][["settings"]])
    }
)

# ---- Set Integration Results ----

#' @rdname methods-for-integration
#' @param results The MOFA or mixOmics results to set in the object. 
#' If null, set to NULL.
#' @exportMethod setMOFA
setMethod(
    f = "setMOFA",
    signature = "RflomicsMAE",
    definition = function(object, results = NULL) {
        metadata(object)[["IntegrationAnalysis"]][["MOFA"]] <- results
        return(object)
    }
)

#' @rdname methods-for-integration
#' @exportMethod setMixOmics
setMethod(
    f = "setMixOmics",
    signature = "RflomicsMAE",
    definition =  function(object, results = NULL) {
        metadata(object)[["IntegrationAnalysis"]][["mixOmics"]] <- results
        return(object)
    }
)



# ---- MixOmics summary ----

#' @title Get an overview of MixOmics integration results
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character.
#' Useful if MixOmics was run on several response variable.
#' If NULL, all variables are taken into account.
#' @return sumMixOmics: A data frame or a list of dataframe
#' (if selectedResponse is NULL) presenting the summary of mixOmics analyses.
#'
#' @rdname methods-for-integration
#' @exportMethod sumMixOmics

setMethod(
    f = "sumMixOmics",
    signature = "RflomicsMAE",
    definition =  function(object, selectedResponse = NULL) {
        if (is.null(metadata(object)$IntegrationAnalysis$mixOmics)) {
            stop("It seems this object has no mixOmics results.")
        }
        
        if (is.null(selectedResponse)) {
            posResponse <- names(metadata(object)$IntegrationAnalysis$mixOmics)
            posResponse <- posResponse[-which(posResponse == "settings")]
            res <- lapply(
                posResponse,
                FUN = function(selResponse) {
                    .getOneMORes(object, selectedResponse = selResponse)
                }
            )
            names(res) <- posResponse
            return(res)
        } else {
            .getOneMORes(object, selectedResponse = selectedResponse)
        }
    }
)

#' @title get one MixOmics result
#' @description Get an overview of MixOmics integration results for
#' a specific response variable.
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string.
#' @return A data frame.
#' @keywords internal
#' @noRd

setMethod(
    f = ".getOneMORes",
    signature = "RflomicsMAE",
    definition =   function(object, selectedResponse) {
        res <- getMixOmics(object, response = selectedResponse)
        # Data_res <- res$MixOmics_results
        Data_res <- res
        
        df <- t(sapply(Data_res$X, dim))
        colnames(df) <- c("Ind", "Features")
        
        if (!is.null(res$MixOmics_tuning_results)) {
            df <- cbind(df, do.call("rbind", Data_res$keepX))
            colnames(df)[!colnames(df) %in% c("Ind", "Features")] <-
                paste("Comp", seq_len(length(Data_res$keepX[[1]])))
        }
        
        return(df)
    }
)

# ----- MixOmics: plot variance explained ----

#' @title plotMOVarExp
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param selectedResponse a character string of the response variable to
#' consider
#' @param mode Can be NULL (default), "cumulative" or "comp".
#' Defines the type of graph to return
#' @return An object of class \link{RflomicsMAE}
#' @importFrom ggpubr ggarrange
#' @exportMethod plotMOVarExp
#'
setMethod(
    f = "plotMOVarExp",
    signature = "RflomicsMAE",
    definition =   function(object, selectedResponse, mode = NULL) {
        if (is.null(getMixOmics(object,
                                response = NULL,
                                onlyResults = TRUE))) {
            stop("It seems this object has no mixOmics results.")
        }
        if (is.null(getMixOmics(object,
                                response = selectedResponse,
                                onlyResults = TRUE))) {
            stop("It seems you didn't run MixOmics on this particular variable.")
        }
        
        Data_res <- getMixOmics(object,
                                response = selectedResponse,
                                onlyResults = TRUE)
        gg_return <- NULL
        
        if (is.null(mode)) {
            gg_return <- ggarrange(.plot_MO_1(Data_res),
                                   .plot_MO_2(Data_res),
                                   ncol = 2)
        }
        else if (tolower(mode) == "cumulative") {
            gg_return <- .plot_MO_1(Data_res)
        }
        else if (tolower(mode) == "comp") {
            gg_return <- .plot_MO_2(Data_res)
        }
        
        return(gg_return)
    }
)
