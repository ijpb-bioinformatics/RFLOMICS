################################### ANNOTATION #############################

#' @title EnrichmentHyperG
#' @param alpha pvalue threshold
#' @param annotation gene annotation
#' @param geneList gene list
#' @return list
#' @noRd
#'
EnrichmentHyperG <- function(annotation, geneList, alpha = 0.01){
  #enrichment_analysis <- function(Reference,Gene_List,Alpha){
  
  Term <- Name <- Domain <- geneID <- NULL
  Pvalue_over <- Pvalue_under <- Decision <- NULL
  
  # ## success in the urn /  For each annotation term, number of annotated genes in the Reference file
  # m=table(Reference[,2])
  # ## failures in the urn / For each annotation term, number of not annotated genes in the Reference file
  # n=length(unique(Reference[,1]))-m
  #
  # trial<-merge(Gene_List,Reference)
  # ## trial effective / number of genes in the gene list
  # k=length(unique(trial[,1]))
  # ## trial success /  For each annotation term, number of annotated genes in the gene list file
  # x=table(factor(trial[,2],levels=rownames(m)))
  
  ## success in the urn /  For each annotation term, number of annotated genes in the Reference file
  
  Urn_Success <- annotation %>% dplyr::group_by(Term, Name, Domain) %>% dplyr::count(name = "Urn_Success")
  
  ## size of reference / nbr of genes in Ref file
  Urn_effective <- length(unique(annotation$geneID))
  
  ##
  #trial<-merge(geneList,annotation)
  trial <- dplyr::filter(annotation , geneID %in% geneList)
  
  ## trial effective / number of genes in the gene list
  Trial_effective <- length(unique(trial$geneID))
  
  ## trial success /  For each annotation term, number of annotated genes in the gene list file
  Trial_Success <- trial %>% dplyr::group_by(Term, Name, Domain) %>% dplyr::count(name = "Trial_Success")
  
  ## Result files
  res=NaN
  # Term=rownames(m)
  # m=as.numeric(m)
  # n=as.numeric(n)
  # x=as.numeric(x)
  # res=data.frame(Term,Urn_Success=m,Urn_Failures=n,Trial_Success=x,Trial_effective=k,
  #                Urn_percentage_Success=signif(100*m/(m+n),3),
  #                Trial_percentage_Success=signif(100*x/k,3),
  #                Pvalue_over=phyper(x-1,m,n,k,lower.tail=FALSE),
  #                Pvalue_under=phyper(x,m,n,k,lower.tail=TRUE))
  
  res= dplyr::full_join(Urn_Success, Trial_Success, by = c("Term", "Name", "Domain")) %>%
    dplyr::mutate(Urn_percentage_Success   = signif(100*Urn_Success/Urn_effective, 3), Urn_effective = Urn_effective,
                  
                  Trial_percentage_Success = signif(100*Trial_Success/Trial_effective, 3), Trial_effective = Trial_effective,
                  Pvalue_over  = stats::phyper(Trial_Success-1,Urn_Success, (Urn_effective-Urn_Success),Trial_effective,lower.tail=FALSE),
                  Pvalue_under = stats::phyper(Trial_Success,  Urn_Success, (Urn_effective-Urn_Success),Trial_effective,lower.tail=TRUE))
  
  # res_over_under <-NULL
  # index=which(res$Pvalue_over<Alpha)
  # if(length(index)!=0)
  # {
  #   res_over <- res[index,]
  #   res_over[,10] <- "overrepresented"
  #   colnames(res_over)[10] <- c("Decision")
  #   res_over_under <- res_over
  # }
  #
  # index=which(res$Pvalue_under<Alpha)
  # if(length(index)!=0)
  # {
  #   res_under <- res[index,]
  #   res_under[,10] <- "underrepresented"
  #   colnames(res_under)[10] <- c("Decision")
  #   res_over_under <- rbind(res_over_under,res_under)
  # }
  
  res_over_under <- NULL
  res_over_under <- res %>% dplyr::mutate(Decision = dplyr::if_else(Pvalue_over <alpha, "overrepresented",
                                                                    dplyr::if_else(Pvalue_under<alpha, "underrepresented", NULL))) %>%
    dplyr::filter(!is.na(Decision))
  
  Results <- list("All_results"   = res, "Over_Under_Results" = res_over_under,
                  "Urn_effective" = Urn_effective, "Trial_effective" = Trial_effective)
  
  return(Results)
  
}