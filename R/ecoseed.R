### ============================================================================
### [data] function and internal function
### ----------------------------------------------------------------------------
# D. Charif

#' @title Ecoseed project data
#' @description 
#' The dataset ecoseed is a multi-omics dataset composed of three data matrices:
#' transcriptomic, metabolomic and proteomic.
#' It has been obtained from the model plant Arabidopsis thaliana in the context 
#' of the study of seed germination and vigor.  
#' In particular, the authors were interested in deciphering key entities 
#' involved in response to environmental stresses (on the mother plant): 
#' influences of temperature (high, medium and low) and imbibition stage 
#' (Dry: DI, early imbibition: EI and late imbibition: LI). 
#' @name ecoseed
#' @aliases ecoseed
#' @docType data
#' @usage data(ecoseed)
#' @format design 
#' @format set1 transcritptomic data: raw read count data matrix
#' @format set2 proteomic data : relative abundance matrix as XIC (area under the extracted ion chromatogram)
#' @format set3 metabolomic data : relative abundance matrix as XIC.
#' @keywords datasets
#' @references See (https://www.uibk.ac.at/botany/ecoseed/home/)
#' @source Loic Rajjou and Gwendal Cueff
#' @examples
#' data(ecoseed)
"ecoseed"
