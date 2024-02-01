
####################################### CO-EXPRESSION ##############################


#' @title try_rflomics
#' @details
#' This function come from https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#  The author indicated that he merged Martins solution (https://stackoverflow.com/a/4952908/2161065) and
#  the one from the R-help mailing list you get with demo(error.catching).
#' @param expr The expression that has been to be evaluated.
#' @return a named list
#' \itemize{
#' \item{\code{value:} }{The results of the expr evaluation or NULL if an error occured }
#' \item{\code{warning:} }{warning message or NULL}
#' \item{\code{error:} }{error message or NULL}
#' }
#' @export
#' @keywords internal
#' @noRd

try_rflomics <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr,
             error    = function(e){ err <- e
             NULL
             }),
    warning = function(w){ warn <- w
    invokeRestart("muffleWarning")}
  )
  
  return(list(value = value,
              warning = warn,
              error = err)
  )
  
}


#' @title coseq.error.manage
#' @param coseq.res.list list of coseq object
#' @param K number of group
#' @param replicates number of replication to run
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @export
#' @importFrom dplyr n group_by summarise mutate filter
#' @importFrom coseq ICL
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @noRd
#'
coseq.error.manage <- function(coseq.res.list, K, replicates, cmd = FALSE){
  
  # Create a table of jobs summary
  error.list <- unlist(lapply(coseq.res.list, function(x){
    ifelse(is.null(x$error), "success", as.character(x$error))
  }))
  
  # status of jobs
  nK_success.job <- table(error.list)["success"]
  
  if (is.na(nK_success.job)) { nK_success.job <- 0 }
  
  # if there is at least one failed job
  # => generate table with error summary
  K.list <- rep(paste0("K",min(K), "-", max(K)), each = replicates)
  
  jobs.tab <- data.frame(K = K.list, error.message = as.factor(error.list))
  
  jobs.tab.sum1 <- jobs.tab %>% 
    dplyr::group_by(K, error.message) %>%
    dplyr::summarise(n = dplyr::n()) %>%  
    dplyr::mutate(prop.failed = round((n/replicates)*100)) %>%
    dplyr::filter(error.message != "success")
  
  jobs.tab.sum <- jobs.tab.sum1
  
  if (nK_success.job != 0) {
    
    # Generate the list of results
    coseq.res.list[["value"]] <- list()
    
    for (x in names(coseq.res.list)) {
      
      if (!is.null(coseq.res.list[[x]]$value)) {
        coseq.res.list[["value"]][[x]] <- coseq.res.list[[x]]$value
      }
    }
    
    # if ICL == NA | loglike == 0
    
    # if (cmd) print("#     => error management: level 2 ")
    like.vec <- unlist(lapply(1:nK_success.job, function(x){ 
      coseq::likelihood(coseq.res.list[["value"]][[x]])
    })) %>%
      lapply(., function(x){ 
        ifelse(is.na(x) | (x==0), "failed", "success") }) %>% unlist()
    
    nK_success <- table(like.vec)["success"]
    
    replicates <- nK_success.job
    
    # expected list of cases
    K.list.ex <- rep(K, each = replicates)
    
    # observed list of cases
    K.list.ob <- stringr::str_replace(string = names(like.vec), pattern = "K=", replacement = "") %>% 
      as.numeric() %>% 
      sort()
    
    # missed cases
    if (length(K.list.ob) != length(K.list.ex)) {
      
      missed.K.vec <- names(table(K.list.ob)[table(K.list.ob) < nK_success.job])
      
      like.vec.bis <- rep("failed", length(missed.K.vec))
      names(like.vec.bis) <- paste0("K=", missed.K.vec)
      
      like.vec <- c(like.vec, like.vec.bis)
    }
    
    jobs.tab <- data.frame(K = names(like.vec), error.message = as.factor(like.vec))
    
    jobs.tab.sum2 <- jobs.tab %>% 
      dplyr::group_by(K, error.message) %>%
      dplyr::summarise(n = dplyr::n()) %>%  
      dplyr::mutate(prop.failed = round((n/replicates)*100)) %>%
      dplyr::filter(error.message != "success")
    
    jobs.tab.sum <- data.table::rbindlist(list(jobs.tab.sum1, jobs.tab.sum2), use.names = TRUE) %>% 
      tibble::tibble()
    
  } else{
    nK_success <- 0
  }
  
  return(list(jobs.tab.sum = jobs.tab.sum, 
              nK_success = nK_success, 
              coseq.res.list.values = coseq.res.list[["value"]]))
}



#' @title coseq.results.process
#' @param coseqObjectList list of coseq object
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @importFrom coseq ICL likelihood clusters
#' @export
#' @noRd
#'
coseq.results.process <- function(coseqObjectList, K, conds){
  
  # ICL plot
  ICL.list <- list()
  
  # get ICL as one vector for all replicates
  
  ICL.vec <- lapply(1:length(coseqObjectList), function(x){ coseq::ICL(coseqObjectList[[x]]) }) %>% unlist()
  ICL.list[["ICL.vec"]] <- ICL.vec
  
  # Find the ICL min by replicates
  
  ICL.min.per.rep <- lapply(1:length(coseqObjectList), function(x){
    ICL.vec.min <- coseq::ICL(coseqObjectList[[x]])
    ICL.vec.min[which(ICL.vec.min == min(ICL.vec.min))]
  })  %>% unlist()
  
  # Construct a table of results with ICL and K:
  
  ICL.tab <- data.frame(K = stringr::str_replace(names(ICL.vec), "K=", ""), ICL = ICL.vec) %>%
    dplyr::mutate(K = as.numeric(K))
  ICL.list[["ICL.tab"]] <- ICL.tab
  
  # Summarize the table: by K, compute the median of the replicate's ICL.
  
  ICL.n <- ICL.tab  %>% 
    dplyr::group_by(.,K) %>% 
    dplyr::filter(!is.na(ICL)) %>%
    dplyr::summarise(median = median(ICL, na.rm = TRUE), n = dplyr::n()) %>%
    dplyr::mutate(K = as.numeric(K))
  
  ICL.list[["ICL.n"]] <- ICL.n
  
  # Search for a replicate with a ICL min corresponding to the K with the min median 
  
  K.ICL.median.min <- ICL.n[which.min(ICL.n$median),]$K
  
  index  <- which(names(ICL.min.per.rep) == paste0("K=", K.ICL.median.min))
  
  # Case where the median.min has a K.min.rep
  if(length(index)>0){
    # Case where there is several rep with a min.rep, the min of them is taken
    index2 <- which(ICL.min.per.rep[index] == min(ICL.min.per.rep[index]))
    coseq.res <- coseqObjectList[index][index2][[1]]
    
    # logLike plot
    logLike.list <- list()
    
    logLike.vec <- lapply(1:length(coseqObjectList), function(x){ coseq::likelihood(coseqObjectList[[x]]) }) %>% unlist()
    logLike.list[["logLike.vec"]] <- logLike.vec
    
    logLike.tab <- data.frame(K = stringr::str_replace(names(logLike.vec), "K=", ""), logLike = logLike.vec) %>% 
      dplyr::mutate(K = as.numeric(K))
    logLike.list[["logLike.tab"]] <- logLike.tab
    
    logLike.n <- logLike.tab %>% 
      dplyr::group_by(.,K) %>% 
      dplyr::filter(!is.na(logLike)) %>%
      dplyr::summarise(median = median(logLike), n = dplyr::n()) %>%
      dplyr::mutate(K = as.numeric(K))
    logLike.list[["logLike.n"]] <- logLike.n
    
    # list of genes per cluster
    clusters <- lapply(1:length(table(coseq::clusters(coseq.res))), function(i){
      names(coseq::clusters(coseq.res)[coseq::clusters(coseq.res) == i])
    })
    names(clusters) <- paste("cluster", 1:length(table(coseq::clusters(coseq.res))), sep = ".")
    
    # nbr of cluster
    # Gestion des NA dans les ICLs
    ICLv <- na.omit(coseq.res@metadata$ICL)
    nb_cluster <- na.omit(coseq.res@metadata$nbCluster)[min(ICLv) == ICLv]
    
    #output
    CoExpAnal <- list()
    CoExpAnal[["results"]]      <- TRUE
    CoExpAnal[["coseqResults"]] <- coseq.res
    CoExpAnal[["clusters"]]     <- clusters
    CoExpAnal[["cluster.nb"]]   <- nb_cluster
    CoExpAnal[["plots"]]        <- list("ICL" = ICL.list, "logLike" = logLike.list)
  } else{
    # Pb of convergence: if there is no K.min.rep which correspond to the median.min, return an error
    CoExpAnal <- list()
    CoExpAnal[["results"]]      <- FALSE
    CoExpAnal[["error"]]  <- "No min.median correspond to min.ICL.rep"
  }
  
  return(CoExpAnal)
}


#' @title run Coseq for co-expression analysis on cluster
#' @param counts matrix
#' @param K number of groups
#' @param replicates number of replication to run
#' @param param.list list of coseq parameters
#' @return coseqResults
#' @importFrom coseq coseq
#' @importFrom clustermq Q_rows 
#' @export
#' @noRd
#'
runCoseq_clustermq <- function(counts, conds, K=2:20, replicates = 5, param.list, silent = TRUE, cmd = FALSE){
  
  # iter <-  rep(K, each = replicates)
  # seed_arg = rep(1:replicates, max(K) - 1)  
  
  df_args <- data.frame(x = rep(K, each = replicates),
                        seed_arg = rep(1:replicates, max(K) - 1) )
  # colnames(df_args) = c("x", "seed_arg")
  
  # nbr_iter <- length(iter)
  nbr_iter <- nrow(df_args)
  coseq.res.list <- list()
  set.seed(12345)
  
  # setting to run coseq on clustermq
  param.list[["object"]] <- counts
  param.list[["K"]] <- K
  
  fx <- function(x, seed_arg){
    
    try_rflomics <- function(expr) {
      warn <- err <- NULL
      value <- withCallingHandlers(
        tryCatch(expr,
                 error    = function(e){ err <- e
                 NULL
                 }),
        warning = function(w){ warn <- w
        invokeRestart("muffleWarning")}
      )
      list(value = value, warning = warn, error = err)
    }
    if (silent) { 
      co <- suppressMessages(capture.output(
        res <- try_rflomics(coseq::coseq(object = param.list[["object"]], 
                                         K = x,
                                         model = param.list$model,
                                         transformation = param.list$transformation,
                                         GaussianModel = param.list$GaussianModel,
                                         normFactors = param.list$normFactors,
                                         meanFilterCutoff = param.list$meanFilterCutoff,
                                         seed = seed_arg, 
                                         verbose = FALSE))
      ))
    }else{
      res <- try_rflomics(coseq::coseq(object = param.list[["object"]], 
                                       K = x,
                                       model = param.list$model,
                                       transformation = param.list$transformation,
                                       GaussianModel = param.list$GaussianModel,
                                       normFactors = param.list$normFactors,
                                       meanFilterCutoff = param.list$meanFilterCutoff,
                                       seed = seed_arg))
    }
    return(res)
  }
  
  # coseq.res.list <- clustermq::Q(fx,
  #                                x = iter, seed_arg = rep(1:replicates, replicates), 
  #                                export = param.list, n_jobs = nbr_iter, pkgs = "coseq")
  # 
  coseq.res.list <- clustermq::Q_rows(fun = fx,
                                      df = df_args, 
                                      export = param.list, n_jobs = nbr_iter, pkgs = "coseq")
  
  names(coseq.res.list) <- c(1:nbr_iter)
  
  CoExpAnal <- list()
  
  if (cmd) print("#     => error management ")
  
  # Create a table of jobs summary
  error.list <- unlist(lapply(coseq.res.list, function(x){
    ifelse(is.null(x$error), "success", as.character(x$error))
  }))
  
  nK_success <- table(error.list)["success"]
  if (cmd) print(paste0("#     => nbr of success jobs: ", nK_success))
  
  K.list <- rep(paste0("K=", K), each = replicates)
  
  jobs.tab <- data.frame(K = K.list, error.message = as.factor(error.list))
  
  jobs.tab.sum <- jobs.tab %>% dplyr::group_by(K, error.message) %>%
    dplyr::summarise(n=dplyr::n()) %>%  dplyr::mutate(prop.failed=round((n/replicates)*100)) %>%
    dplyr::filter(error.message != "success")
  
  
  # If they are at least the half of K which succeed, valid results
  if(nK_success !=0 ){
    
    if (cmd) print("#     => process results ")
    # Generate the list of results
    #coseq.res.list[["value"]] <- lapply(coseq.res.list,function(x){x$value})
    coseq.res.list[["value"]] <- list()
    for(x in names(coseq.res.list)){
      
      if(!is.null(coseq.res.list[[x]]$value)){
        coseq.res.list[["value"]][[x]] <- coseq.res.list[[x]]$value
      }
    }
    
    CoExpAnal <- coseq.results.process(coseq.res.list[["value"]], K = K, conds = conds)
    CoExpAnal[["warning"]] <- coseq.res.list$warning
    
    if(nK_success/length(iter) < 0.8){
      CoExpAnal[["error"]] <- TRUE
    }
    
  }
  # Réinitialisation de l'objet CoExpAnal
  else{
    CoExpAnal[["results"]] <- FALSE
    CoExpAnal[["error"]] <- TRUE
    
  }
  
  CoExpAnal[["stats"]] <- jobs.tab.sum
  
  return(CoExpAnal)
}




#' @title run Coseq for co-expression analysis on cluster
#' @param counts matrix
#' @param K number of groups
#' @param replicates number of replication to run
#' @param param.list list of coseq parameters
#' @return coseqResults
#' @importFrom coseq coseq clusters
#' @export
#' @noRd
#'
runCoseq_local <- function(counts, conds, K=2:20, replicates = 5, param.list, silent = TRUE, cmd = FALSE){
  
  iter <- rep(K, replicates)
  coseq.res.list <- list()
  
  # set.seed(12345)
  if (silent) {
    
    coseq.res.list <- lapply(1:replicates, function(x){
      
      co <- capture.output(suppressMessages(
        res <- try_rflomics(
          coseq::coseq(counts, K = K, parallel = TRUE,
                       model            = param.list[["model"]],
                       transformation   = param.list[["transformation"]],
                       meanFilterCutoff = param.list[["meanFilterCutoff"]],
                       normFactors      = param.list[["normFactors"]],
                       GaussianModel    = param.list[["GaussianModel"]],
                       seed = x, 
                       verbose = FALSE)
        )))
      return(res)
    })  
    
  }else{
    coseq.res.list <- lapply(1:replicates, function(x){
      
      try_rflomics(coseq::coseq(counts, K = K, parallel = TRUE,
                                model            = param.list[["model"]],
                                transformation   = param.list[["transformation"]],
                                meanFilterCutoff = param.list[["meanFilterCutoff"]],
                                normFactors      = param.list[["normFactors"]],
                                GaussianModel    = param.list[["GaussianModel"]],
                                seed = x))
    })  
  }
  
  names(coseq.res.list) <- c(1:replicates)
  
  CoExpAnal <- list()
  
  # error managment
  if (cmd) print("#     => error management: level 1 ")
  coseq.error.management <- coseq.error.manage(coseq.res.list = coseq.res.list, 
                                               K = K, 
                                               replicates = replicates,
                                               cmd = cmd)
  
  nK_success   <- coseq.error.management$nK_success
  
  # If they are more than 60 % of succeeded jobs, valid results, find min.median.ICL
  
  if (nK_success/length(iter) >=0.8) {
    
    CoExpAnal <-  coseq.results.process(coseqObjectList = coseq.error.management$coseq.res.list.values, 
                                        K = K,
                                        conds = conds)
    # If ICL.median has been found
    
    if(CoExpAnal[["results"]]==TRUE){
      
      CoExpAnal[["results"]] <- TRUE
      CoExpAnal[["warning"]] <- coseq.res.list$warning
      
      if (cmd) { 
        print(paste0("#     => Number of clusters: ", 
                     max(unique(coseq::clusters(CoExpAnal$coseqResults)))))
      }
    }
  } else{
    CoExpAnal[["results"]] <- FALSE
    CoExpAnal[["error"]] <- "nK_success/length(iter) < 0.8"
  }
  
  CoExpAnal[["stats"]] <- coseq.error.management$jobs.tab.sum 
  
  return(CoExpAnal)
}



#' @title profile boxplot per cluster
#' @param coseq.res coseq object
#' @param selectedCluster cluster num
#' @param conds condition matrix
#' @return boxplot profiles.
#' @export
#' @importFrom ggplot2 geom_boxplot facet_wrap theme element_blank
#' @importFrom dplyr arrange
#' @importFrom purrr reduce
#' @noRd
#'
coseq.y_profile.one.plot <- function(coseq.res, selectedCluster, conds){
  
  samples <- variable <- value <- cluster <- NULL
  
  nb_cluster <- coseq.res@metadata$nbCluster[min(coseq.res@metadata$ICL) == coseq.res@metadata$ICL]
  groups <- conds %>% dplyr::arrange(factor(samples, levels = names(coseq.res@y_profiles)))
  y_profiles <- list()
  
  for (i in 1:nb_cluster) {
    y_profiles[[i]] <- coseq.res@y_profiles[coseq.res@allResults[[paste0("K=", nb_cluster)]][,i] > 0.8,] %>%
      data.frame() %>% 
      reshape2::melt() %>%  
      dplyr::rename(samples = variable) %>%
      dplyr::full_join(conds , by = "samples") %>% 
      dplyr::mutate(cluster = i)
  }
  
  y_profiles.gg <-  y_profiles %>% purrr::reduce(rbind)
  y_profiles.gg$groups <- factor(y_profiles.gg$groups, levels = unique(conds$groups))
  y_profiles.gg$samples <- factor(y_profiles.gg$samples, levels = unique(conds$samples))
  
  
  p <- ggplot2::ggplot(data = dplyr::filter(y_profiles.gg, cluster == selectedCluster)) +
    ggplot2::geom_boxplot(ggplot2::aes(x = samples, y = value, fill = groups), outlier.size = 0.3) + 
    ggplot2::facet_wrap(~cluster) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  print(p)
}
