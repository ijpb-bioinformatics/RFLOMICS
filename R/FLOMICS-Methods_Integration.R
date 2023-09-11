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
                        cases_to_try = 5) {
    
    method <- switch(toupper(method),
      "MIXOMICS" = "MixOmics",
      "MOFA"  = "MOFA",
      "MOFA2" = "MOFA",
      "MOFA+" = "MOFA"
    )
    
    if (any(!omicsToIntegrate %in% names(object))) {
      stop("There are omics to integrate that are not names from the object")
    }
    
    if (RFLOMICS:::isTagName(object, contrasts_names)) 
      contrasts_names <- RFLOMICS:::convertTagToContrast(object, contrasts_names)

    preparedObject <- prepareForIntegration(
      object = object,
      omicsToIntegrate = omicsToIntegrate,
      rnaSeq_transfo = rnaSeq_transfo,
      choice = choice,
      contrasts_names = contrasts_names,
      type = type,
      group = group,
      method = method
    )

    if (toupper(method) == "MOFA") {
      object@metadata[["MOFA"]] <- NULL

      MOFA_run <- run_MOFA_analysis(
        object = preparedObject,
        scale_views = scale_views,
        maxiter = maxiter,
        num_factors = num_factors
      )

      object@metadata[["MOFA"]][["MOFA_results"]] <- MOFA_run$MOFAObject.trained
      object@metadata[["MOFA"]][["MOFA_untrained"]] <- MOFA_run$MOFAObject.untrained
      object@metadata[["MOFA"]][["MOFA_selected_filter"]] <- type
      object@metadata[["MOFA"]][["MOFA_selected_contrasts"]] <- contrasts_names
    } else if (toupper(method) == "MIXOMICS") {
      object@metadata[["mixOmics"]] <- NULL

      if (is.null(selectedResponse)) selectedResponse <- colnames(object@metadata$design@ExpDesign)
      
      MixOmics_res <- lapply(selectedResponse,
        FUN = function(response_var) {
          res_mixOmics <- run_MixOmics_analysis(
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
      names(MixOmics_res) <- selectedResponse

      object@metadata[["mixOmics"]] <- MixOmics_res
    }

    return(object)
  }
)

#' @title prepareForIntegration
#' @description This function transforms a MultiAssayExperiment produced by rflomics into an untrained MOFA objects or a list to use for mixOmics. It checks for batch effect to correct them prior to the integration.
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
                        method = c("MOFA", "MixOmics")) {
    if (is.null(omicsToIntegrate)) omicsToIntegrate <- names(object)

    # Checking for batch effects
    correct_batch <- FALSE
    ftypes <- getFactorTypes(object)

    if (any(ftypes == "batch")) {
      correct_batch <- TRUE
      # colBatch <- names(ftypes)[ftypes == "batch"]
    }

    object <- object[, , omicsToIntegrate]
    # omics_types <-  getOmicsTypes(object)

    # On each selected omics, according to its type, apply transformation if demanded.
    # Filter DE entities
    # TODO : add a possibility of choice for keeping every entity (small tables)
    for (SEname in omicsToIntegrate) {
      SEobject <- object[[SEname]]
      omicsType <- getOmicsTypes(SEobject)
      list_args <- list(
        object = object,
        SEname = SEname,
        correctBatch = correct_batch,
        contrasts_names = contrasts_names,
        type = type,
        choice = choice
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
      MOFAObject <- MOFA2::create_mofa(object,
        group = group,
        extract_metadata = TRUE
      )
      return(MOFAObject)
    } else if (method == "MixOmics") {
      # Common samples names:
      object <- object[, Reduce("intersect", colnames(object))]

      MixOmicsObject <- list(
        blocks = lapply(object@ExperimentList, FUN = function(SE) t(SummarizedExperiment::assay(SE))),
        metadata = object@colData
      )
      MixOmicsObject$blocks <- lapply(MixOmicsObject$blocks, FUN = function(mat) mat[order(rownames(mat)), ])

      return(MixOmicsObject)
    }
  }
)


######################## INTEGRATION USING MOFA ########################

#' @title run_MOFA_analysis
#' @description Runs a MOFA analysis based on an untrained MOFA object and user arguments.
#' @param object An untrained MOFA object
#' @param scale_views boolean. MOFA option to scale the views so they have the same variance. Default is FALSE.
#' @param maxiter integer. MOFA option, maximum number of iterations to be considered if there it does not converge. Default is 1000.
#' @param num_factors integer. MOFA option, maximum number of factor to consider. Default is 10.
#' @param ... Not in use.
#' @return A list with an untrained MOFA object (containing all options for the run) and a trained MOFA object
#' @export
#' @exportMethod run_MOFA_analysis
#' @importClassesFrom MOFA2 MOFA
#' 
methods::setMethod(
  f = "run_MOFA_analysis",
  signature = "MOFA",
  definition = function(object,
                        scale_views = FALSE,
                        maxiter = 1000,
                        num_factors = 10,
                        ...) {
    
    data_opts  <- MOFA2::get_default_data_options(object)
    model_opts <- MOFA2::get_default_model_options(object)
    train_opts <- MOFA2::get_default_training_options(object)

    data_opts$scale_views  <- scale_views
    train_opts$maxiter     <- maxiter
    train_opts$verbose     <- FALSE
    model_opts$num_factors <- num_factors

    MOFAObject.untrained <- MOFA2::prepare_mofa(
      object           = object,
      data_options     = data_opts,
      model_options    = model_opts,
      training_options = train_opts
    )

    MOFAObject.trained <- MOFA2::run_mofa(MOFAObject.untrained, use_basilisk = TRUE)
    # peut poser probleme au niveau python et mofapy.
    # Installer python, numpy et mofapy, ensuite reinstaller totalement package MOFA2 et restart R.

    return(list("MOFAObject.untrained" = MOFAObject.untrained, "MOFAObject.trained" = MOFAObject.trained))
  }
)

######################## INTEGRATION USING MixOMICS ########################

#' @title run_MixOmics_analysis
#' @description Run a mixOmics analysis. Given the specification of the user (type of response, if response there is, multi block or not)
#'  the function will determine which mixOmics function to use. Please see mixOmics manual or website for more information.
#' @param object list of blocks (matrices) with the same samples (rows)
#' @param scale_views Boolean. Do the matrices have to be scaled? Is used inside the mixOmics function.
#' @param selectedResponse Default is NULL (pls functions). The response, can be a matrix or a single factor (discriminant analysis is set in this case).
#' @param ncomp Number of component to be computed.
#' @param link_dataset Numeric between 0 and 1. Only used for multi block analysis. Indicates the correlation to be considered between the matrices.
#'        It impacts the weights of the features, hence the feature selection. Please see mixOmics user's guide for better explaination.
#' @param link_response Numeric between 0 and 1. Indicates the correlation to be considered between the matrices and the response matrix.
#'        It impacts the weights of the features, hence the feature selection. Please see mixOmics user's guide for better explaination.
#' @param sparsity Boolean. Indicates if there is a feature selection purpose. If TRUE, functions like spls(da), block.spls(da) will be used.
#' @param cases_to_try If sparsity is TRUE, indicates the number of cases to try for the feature selection. The best outcome, as computed by tuning function, will be displayed.
#' @return A mixOmics result.
#' @export
#' @exportMethod run_MixOmics_analysis
methods::setMethod(
  f = "run_MixOmics_analysis",
  signature = "list",
  definition = function(object,
                        scale_views = FALSE,
                        selectedResponse = NULL,
                        ncomp = 2,
                        link_datasets = 1,
                        link_response = 1,
                        sparsity = FALSE,
                        cases_to_try = 5,
                        ...) {
    list_res <- list()
    dis_anal <- FALSE # is this a discriminant analysis

    # TODO : check this, c'est un peu etrange d'avoir eu a reordonner les lignes dans prepareForIntegration
    # du coup ca demande de le faire aussi pour Y
    # TODO : faudrait le faire directement dans preparedList ?
    Y <- data.frame(object$metadata, stringsAsFactors = TRUE) %>%
      tibble::rownames_to_column(var = "rowNam") %>%
      dplyr::arrange(rowNam) %>%
      tibble::column_to_rownames(var = "rowNam") %>%
      dplyr::select(tidyselect::all_of(selectedResponse))

    # Is this a discriminant analysis?
    # TODO : is this used?
    if (ncol(Y) == 1 && is.factor(Y[, 1])) {
      dis_anal <- TRUE
      Y <- Y[, 1]
    } else {
      YrowNames <- rownames(Y)
      YFactors <- do.call("cbind", lapply(1:ncol(Y), FUN = function(j) {
        if (is.factor(Y[, j])) {
          mat_return <- mixOmics::unmap(Y[, j])
          colnames(mat_return) <- attr(mat_return, "levels")
          return(mat_return)
        }
      }))

      Y <- cbind(Y %>% dplyr::select_if(is.numeric), YFactors)
      Y <- apply(Y, 2, as.numeric)
      rownames(Y) <- YrowNames
    }

    # Design matrix
    Design_mat <- matrix(link_datasets,
      nrow = length(object$blocks) + 1,
      ncol = length(object$blocks) + 1
    )
    Design_mat[, ncol(Design_mat)] <-
      Design_mat[nrow(Design_mat), ] <- link_response
    diag(Design_mat) <- 0

    # What function to use for the analysis (can't be sparse if there is a continous response)
    functionName <- "pls"
    if (dis_anal) functionName <- paste0(functionName, "da")
    if (sparsity && !is.numeric(Y)) functionName <- paste0("s", functionName)
    if (length(object$blocks) > 1) functionName <- paste0("block.", functionName)

    # Model Tuning (if required, for sparsity)
    if (sparsity && dis_anal && functionName != "block.spls") {
      # no tune.block.spls so far, there is one for spls
      # It is a bit weird to consider tuning with folds when there is very few samples per condition
      # Add a warning or something when it's the case?
      # Also consider adding a warning when the number of feature is still very high even after differential analysis

      tune_function <- paste0("tune.", functionName)

      test_keepX <- lapply(object$blocks, FUN = function(dat) {
        ceiling(seq(from = ceiling(0.1 * ncol(dat)), to = ncol(dat), length.out = cases_to_try))
      })

      list_tuning_args <- list(
        X = object$blocks,
        Y = Y,
        ncomp = ncomp,
        scale = scale_views,
        test.keepX = test_keepX,
        folds = min(floor(nrow(Y) / 2), 10) # TODO ???
      )
      if (length(object$blocks) > 1) list_tuning_args$design <- Design_mat

      list_res$tuning_res <- do.call(getFromNamespace(tune_function, ns = "mixOmics"), list_tuning_args)
    }

    # Model fitting

    list_args <- list(
      X = object$blocks,
      Y = Y,
      ncomp = ncomp,
      scale = scale_views
    )
    if (length(object$blocks) > 1) list_args$design <- Design_mat
    if (sparsity && !functionName %in% c("block.spls", "block.pls")) list_args$keepX <- list_res$tuning_res$choice.keepX

    list_res$analysis_res <- do.call(getFromNamespace(functionName, ns = "mixOmics"), list_args)

    return(list_res)
  }
)
