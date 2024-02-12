


###### METHOD to check the completness of the ExpDesign





##### RESETTING ####

#' resetFlomicsMultiAssay
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param results vector of results names
#' @param dataset dataset name. If dataset == NULL, all datasets will be reset
#' @return An object of class \link{RflomicsMAE}
#' @export
#' @exportMethod resetFlomicsMultiAssay
#' @noRd
#'
methods::setMethod(f="resetFlomicsMultiAssay", signature="RflomicsMAE",
                   
                   definition <- function(object, results, datasets = NULL){
                     
                     # if dataset is null we take all datasets present in RflomicsMAE object
                     if(is.null(datasets)){
                       datasets <- unlist(object@metadata$omicList)
                     }
                     else{
                       # check if given dataset name include in datasets presente in RflomicsMAE object
                       if(!datasets %in% unlist(object@metadata$omicList)){
                         warning("The given dataset name is not present in RflomicsMAE object")
                         return(object)
                       }
                     }
                     
                     for(res in results){
                       
                       for(data in datasets){
                         if(!is.null(object[[data]])){
                           
                           if(!is.null(object[[data]]@metadata[[res]])){ object[[data]]@metadata[[res]] <- NULL }
                         }
                       }
                       
                       object@metadata[[res]] <- NULL
                     }
                     
                     
                     # for(data in datasets){
                     #   
                     #   if(!is.null(object[[data]])){
                     #     
                     #     dataset <- object[[data]]
                     #     
                     #     for(res in results){
                     #       if(!is.null(dataset@metadata[[res]])){ dataset@metadata[[res]] <- NULL }
                     #     }
                     #     
                     #     object[[data]] <- dataset
                     #   }
                     #   
                     # }
                     return(object)
                   })


