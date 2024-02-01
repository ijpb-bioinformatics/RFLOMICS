##### Global import ####

#' @importFrom ggplot2 ggplot geom_col theme_classic aes
#' theme element_text element_blank ylab xlab ggtitle
#' scale_fill_gradientn geom_tile theme_bw guides scale_fill_gradient2
#' guide_colourbar labs

# @export
#' @importFrom magrittr "%>%" 
magrittr::`%>%`




#' @title TMM.Normalization
#' Interface to the calcNormFactors function of the edgeR package  with the choosen TMM parameters as the normalization method
#' @param counts numeric matrix of read counts
#' @param groups vector or factor giving the experimental group/condition for each sample/library.
#' @return a data.frame with a row for each sample and columns group, lib.size and norm.factors containing the group labels, library sizes and normalization factors. Other columns can be optionally added to give more detailed sample information.
#' @export
#' @importFrom edgeR DGEList calcNormFactors
#' @noRd

TMM.Normalization <- function(counts, groups){
  dge <- edgeR::DGEList(counts=counts, group=groups)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  nf  <- dge$samples
  return(nf)
}


#' plotDistr
#'
#' @param abundances matrix or dataframe of feature/gene abundances/counts
#' @export
#' @importFrom ggplot2 geom_density xlab
#' @importFrom reshape2 melt
#' @noRd
plotDistr <- function(abundances, dataType, transform_method){
  
  
  value <- samples <- NULL
  
  switch(dataType,
         "RNAseq" = {
           pseudo_counts <- log2(abundances+1) %>% reshape2::melt()
           colnames(pseudo_counts) <- c("features", "samples", "value")
           x_lab <-"log2(feature abundances)"
         },
         "proteomics"={
           pseudo_counts <- abundances %>% reshape2::melt()
           colnames(pseudo_counts) <- c("features", "samples", "value")
           x_lab <- switch (transform_method,
                            "log2"= "log2(feature abundances)",
                            "none"="feature abundances")
         },
         "metabolomics"={
           pseudo_counts <- abundances %>% reshape2::melt()
           colnames(pseudo_counts) <- c("features", "samples", "value")
           x_lab <- switch (transform_method,
                            "log2"= "log2(feature abundances)",
                            "none"="feature abundances")
         }
  )
  
  p <- ggplot2::ggplot(pseudo_counts) + geom_density(aes(value, color=samples) ) + xlab(x_lab) +
    theme(legend.position='none')
  print(p)
  
}


# TODO do we use this?
# utils::globalVariables(names(data))




