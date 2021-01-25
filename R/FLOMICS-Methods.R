#' @title [\code{\link{ExpDesign-class}}] Class constructor
#' @description
#' @param ExpDesign data.frame with experimental design
#' @param refList vector with factor ref
#' @param typeList vector with type of factor
#' @return An object of class [\code{\link{ExpDesign-class}}]
#' @name ExpDesign-Constructor
#' @rdname ExpDesign-Constructor
#' @export
ExpDesign.constructor <- function(ExpDesign, refList, typeList){

  # List.Factors
  dF.List <- lapply(1:dim(ExpDesign)[2], function(i){
    
    relevel(as.factor(ExpDesign[[i]]), ref=refList[[i]])
  })
  names(dF.List) <- names(ExpDesign)

  # Factors.Type
  names(typeList) <- names(ExpDesign)
  
  
  Design = new(Class = "ExpDesign",
               ExpDesign=ExpDesign,
               List.Factors=dF.List,
               Factors.Type=typeList,
               Model.formula=vector(),
               Model.matrix=vector(),
               Contrasts.List=list(),
               Contrasts.Sel=data.frame(),
               Contrasts.Coeff=data.frame())
           
  return(Design)
  }

#' @title RunDiffAnalysis
#' @param An object of class [\code{\link{MultiAssayExperiment}]
#' @param data omic data type
#' @param DiffMethod Differential analysis method
#' @param FDR FDR threshold
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @return MultiAssayExperiment
#' @exportMethod RunDiffAnalysis
#' @examples
#'
setMethod(f="RunDiffAnalysis",
          signature="MultiAssayExperiment",
          definition <- function(object, data, FDR = 0.05, DiffAnalysisMethod, clustermq){

            # Run the Diff analysis and get the results as a list of object depending of the

            ListOfDiffResults <- switch(DiffAnalysisMethod,
                                     "edgeRglmfit"=edgeR.AnaDiff(object, data, clustermq)
                                )

            # The results of the methods have to be formated

            # Set an AnaDiff object to
            object@ExperimentList[[data]]@metadata[["AnaDiff"]] <- ListOfDiffResults

            #
            object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]] <- lapply(ListOfDiffResults, function(x){
              
              res<-topTags(x, n = dim(x)[1])
              DEGs<- res$table[res$table$FDR <= FDR,]
              #DEGs<-res$table
              return(DEGs)
            })
            names(object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]]) <- names(ListOfDiffResults)

            ## merge results in bin matrix
            DEG_list <- lapply(1:length(object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]]), function(x){
              
              res <- object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]][[x]]
              tmp <- data.frame(DEG = rownames(res), bin = rep(1,length(rownames(res))))
              colnames(tmp) <- c("DEG", paste("H", x, sep=""))
              return(tmp)
            })
            names(DEG_list) <- names(object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]])
            
            object@ExperimentList[[data]]@metadata[["AnaDiffDeg.mat"]] <- DEG_list %>% purrr::reduce(full_join, by="DEG") %>% 
              mutate_at(.vars = 2:(length(DEG_list)+1), .funs = function(x){if_else(is.na(x), 0, 1)}) %>% data.table()
            
            return(object)
          })


#' @title [\code{\link{DiffAnalysis-class}}] Class constructor
#'
#' @description
#'
#' @return An object of class [\code{\link{DiffAnalysis-class}}]
#' @name DiffAnalysis-Constructor
#' @rdname DiffAnalysis-Constructor
#' @export
DiffAnalysis <- function(){
  .DiffAnalysis(
    method="glmLRT",
    ListOfDEResults=list(),
    ListOfRawMethodResults=list()
  )}


#' @title Multivariate Quality Check
#'
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param data data name
#' @param axis The number of PCA axis
#' @param PCA pca axis to plot
#' @param pngFile plot file name to save

#'
#' @exportMethod mvQCdesign
#'
#' @rdname mvQCdesign
setMethod(f="mvQCdesign",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA=c("raw","norm"), axis=5, pngFile=NULL){

            resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <-  object@metadata$colDataStruc["n_dFac"]

            bigdf <- list()
            for(i in 1:n_dFac){
              Factor <-  object@colData[,i]
              nameFac <- names(object@colData)[i]
              df <- list()
              for(j in 1:axis){
                qc = as.vector(resPCA$ind$coord[,j])
                o = order(qc)[order(Factor[order(qc)])]
                col = cc[Factor][o]
                y=qc[o]
                Axis = rep(paste("PCA",j, "\n(Var=",round(resPCA$eig[j,2],1),")",sep=""),length(qc))
                Fac = Factor[o]
                df[[j]] <- data.frame(y,col,Axis,dfac=names(object@colData)[1:n_dFac][i],
                                      Levels=Fac,x=1:length(y))
              }
              bigdf[[i]] <- dplyr::bind_rows(df)
            }
            big <- dplyr::bind_rows(bigdf)
            out <- by(data = big, INDICES = big$dfac, FUN = function(m) {
              m <- droplevels(m)
              names(m) <- c("y", "col", "Axis", "dfac", unique(m$dfac), "x")
              m <- ggplot(m,aes(y=y,x=x),aes_string(color=unique(m$dfac)))+
                geom_bar(stat = "identity",position = position_dodge(),aes_string(fill=unique(m$dfac)))+
                facet_grid(as.factor(dfac)~Axis) +
                labs(x = "Samples", y="Coordinates on \n the PCA axis")+
                theme(axis.title.y = element_text(size = 5),
                      axis.title.x=element_text(size = 5),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
            })
            p <- do.call(grid.arrange, out)
            print(p)
            
            if(! is.null(pngFile)){
              ggsave(filename = pngFile,  plot = p)
            }
})


#' @title multivariate QC data
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param data data name
#' @param axis The number of PCA axis
#' @param PCA pca axis to plot
#' @param pngFile plot file name to save
#'
#' @exportMethod mvQCdata
#' @rdname mvQCdata
setMethod(f="mvQCdata",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA=c("raw","norm"),axis=3, pngFile=NULL){

            resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <- object@metadata$colDataStruc["n_dFac"]
            n_qcFac <- dim(object[[data]]@colData[,-c(1:2)])[2]
            var <- names(object[[data]]@colData[,-c(1:2)])

            corA=list()
            VarAxis = list()
            for(i in 1:axis){
            corA[[i]] <- cor(resPCA$ind$coord[,i],
                             as.data.frame(object[[data]]@colData[,-c(1:2)]),
                             method="spearman")
            VarAxis[[i]] <- paste("\n(Var=",round(resPCA$eig[i,2],1),")",sep="")
            }

            df = data.frame("QCparam"=rep(var,axis),
                            "Spearman"= unlist(corA),
                            "Axis"=rep(paste(rep("Axis",axis),1:axis,VarAxis,sep=""), each=n_qcFac))

            p <- ggplot(df,aes(x=Axis, y=abs(Spearman),fill=QCparam))+
              geom_bar(stat="identity",position=position_dodge(),width=0.7)+ylim(0,1)+
              labs(x = "Axis number", y="Cor(Coord_dFactor_PCA,QCparam)")
            
            print(p)
            
            if(! is.null(pngFile)){
              ggsave(filename = pngFile, plot = p)
            }
            
          })



#' @title abundanceBoxplot
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param dataType omic data type
#' @param pngFile
#' @exportMethod abundanceBoxplot
#' @rdname abundanceBoxplot
#'
setMethod(f= "abundanceBoxplot",
          signature = "MultiAssayExperiment",
          definition <- function(object, dataType, pngFile=NULL){

            # this function generate boxplot (abandance distribution) from raw data and normalized data

            
            sample_names <- row.names(object@colData)

            groups  <- object@metadata$design@List.Factors[object@metadata$design@Factors.Type == "Bio"] %>% as.data.frame() %>%
                       unite(col="groups", sep="_", remove = FALSE) %>% mutate(samples=sample_names)

            # normalized data
            pseudo  <- log2(scale(assay(object[[dataType]]), center=FALSE,
                                  scale=object[[dataType]]@metadata$Normalization$coefNorm$norm.factors)+1) %>% reshape2::melt()
            colnames(pseudo) <- c("feature", "samples", "value")
            pseudo_bis <- full_join(pseudo, groups, by="samples")

            pseudo_bis$samples <- factor(pseudo_bis$samples, levels = unique(pseudo_bis$samples))

            # boxplot
            p <- ggplot(pseudo_bis, aes(x=samples, y=value)) + geom_boxplot(aes(fill=groups)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
              xlab(paste0(dataType, " samples"))
              
              #scale_fill_manual(values=col)
            print(p)
            
            if(! is.null(pngFile)){
              ggsave(filename = pngFile, plot = p)
            }
            
          }
)


#' @title plotPCAnorm
#'
#' @param MultiAssayExperiment
#' @param condition
#' @param color color palette
#'
#' @return
#' @exportMethod plotPCAnorm
#'
#' @examples
setMethod(f= "plotPCAnorm",
          signature = "MultiAssayExperiment",
          definition <- function(object, data, PCA, PCs=c(1,2), condition="groups", pngFile){

            #
            PC1 <- paste("Dim.",PCs[1], sep="")
            PC2 <- paste("Dim.",PCs[2], sep="")

            groups    <- unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]),
                               col="groups", sep="_", remove = TRUE)$groups
            conditions<- object@colData %>% as.data.frame() %>% mutate(samples=row.names(.), groups=groups)

            score     <- object[[data]]@metadata$PCAlist[[PCA]]$ind$coord[, PCs] %>% as.data.frame() %>%
                         mutate(samples=row.names(.)) %>% full_join(., conditions, by="samples")

            var1 <- round(object[[data]]@metadata$PCAlist[[PCA]]$eig[PCs,2][1], digit=3)
            var2 <- round(object[[data]]@metadata$PCAlist[[PCA]]$eig[PCs,2][2], digit=3)

            
            p <- ggplot(score, aes_string(x=PC1, y=PC2, color=condition))  +
              geom_point(size=3) +
              geom_text(aes(label=samples), size=3, vjust = 0) +
              xlab(paste(PC1, " (",var1,"%)", sep="")) +
              ylab(paste(PC2, " (",var2,"%)", sep="")) +
              geom_hline(yintercept=0, linetype="dashed", color = "red") +
              geom_vline(xintercept=0, linetype="dashed", color = "red") +
              theme(strip.text.x = element_text(size=8, face="bold.italic"),
                    strip.text.y = element_text(size=8, face="bold.italic")) #+
              #scale_color_manual(values=col$colors)
            
            print(p)
            if(! is.null(pngFile)){
              ggsave(filename = pngFile, plot = p)
            }
            })



#' @title barplotPCAnorm
#'
#' @param MultiAssayExperiment
#' @param condition
#' @param colors color palette
#'
#' @return
#' @exportMethod barplotPCAnorm
#'
#'
#' @examples
setMethod(f= "barplotPCAnorm",
          signature = "MultiAssayExperiment",
          definition <- function(object, condition="samples"){

            col <- colorPlot(object@design, object@colData, condition=condition)

            score_raw  <- object@listPCA$raw$ind$coord  %>% data.table::melt %>%
                          mutate(tag="1.Unnormalised data")
            score_norm <- object@listPCA$norm$ind$coord %>% data.table::melt %>%
                          mutate(tag=paste("2.Normalised data : ", object@Normalization@Method,  sep=""))

            score <- rbind(score_raw, score_norm)
            colnames(score) <- c("samples", "PCs", "value", "tag")

            ggplot(data=score, aes(x=PCs, y=value, fill=samples)) +
              geom_bar(stat="identity", position=position_dodge(), color="black") +
              facet_grid(tag~PCs, scale ="free", space = "free") + scale_fill_manual(values=col)
          })




#' FilterLowAbundance
#'
#' @param MultiAssayExperiment
#' @param threshold
#'
#' @return MultiAssayExperiment
#' @exportMethod FilterLowAbundance
#'
#' @examples
setMethod(f= "FilterLowAbundance",
          signature = "MultiAssayExperiment",
          definition <- function(object, data, Filter_Strategy = "NbConditions", CPM_Cutoff = 5){
            
            objectFilt <- object[[data]]
            
            ## nbr of genes with 0 count
            genes_flt0  <- objectFilt[rowSums(assay(objectFilt)) <= 0, ]@NAMES
            
            ## remove 0 count 
            objectFilt <- objectFilt[rowSums(assay(objectFilt))  > 0, ]
            
           
            
            
            ## filter cpm
            NbReplicate <- table(unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]), col="groups", sep="_"))
            NbConditions <- unique(unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]), col="groups", sep="_"))$groups %>% length()
            
            switch(Filter_Strategy,
                   "NbConditions" = { keep <- rowSums(cpm(assay(objectFilt)) >= CPM_Cutoff) >=  NbConditions },
                   "NbReplicates" = { keep <- rowSums(cpm(assay(objectFilt)) >= CPM_Cutoff) >=  min(NbReplicate) },
                   "filterByExpr" = { dge  <- edgeR::DGEList(counts = assay(objectFilt), genes = rownames(assay(objectFilt)))
                                      #keep <- filterByExpr(dge, GLM_Model) 
                                      keep <- filterByExpr(dge)
                                      }
                   )
            
            ## nbr of genes filtered
            genes_flt1  <- objectFilt[!keep]@NAMES
            
            objectFilt@metadata[["FilteredFeature"]] <-  c(genes_flt0, genes_flt1)
            
            object@ExperimentList[[paste0(data, ".filtred")]] <- objectFilt[keep]
            
            return(object)

          })



#' @title RunNormalization
#' @param MultiAssayExperiment
#' @param data omic data type
#' @param NormMethod normalisation methode
#' @return MultiAssayExperiment
#' @exportMethod RunNormalization
#' @examples
#'
setMethod(f="RunNormalization",
          signature="MultiAssayExperiment",
          definition <- function(object, data, NormMethod){

            groups <- unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]), col="groups", sep="_")$groups

            coefNorm  = switch(NormMethod,
                               "TMM"=TMM.Normalization(assay(object[[data]]), groups)
            )
            object@ExperimentList[[data]]@metadata[["Normalization"]] <- list(methode = NormMethod, coefNorm = coefNorm)
            return(object)
          }
)



#' @title RunPCA
#' @param MultiAssayExperiment
#' @param data omic data type
#' @param PCA norm or raw
#' @return MultiAssayExperiment
#' @exportMethod RunPCA
#' @examples
#'
setMethod(f="RunPCA",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA){

            if(PCA=="raw"){
            pseudo <- log2(scale(assay(object[[data]]),
                                      center=FALSE)+1)
            }
            else if(PCA=="norm"){
            pseudo <- log2(scale(assay(object[[data]]),
                                        center=FALSE,
                                        scale=object[[data]]@metadata[["Normalization"]]$coefNorm$norm.factors)+1)
            }

            pca <- FactoMineR::PCA(t(pseudo),ncp = 5,graph=F)

            object@ExperimentList[[data]]@metadata[["PCAlist"]][[PCA]]<- pca

            return(object)
          }
)



#' @title SetModelMatrix
#'
#' @param ExpDesign 
#'
#' @return
#' @exportMethod SetModelMatrix
#'
#' @examples
setMethod(f="SetModelMatrix",
          signature="ExpDesign",
          definition <- function(object){
            
            #  => list of fact from model
            Fact.list <- strsplit(gsub("[ ~]", "",object@Model.formula), split = "+", fixed = TRUE)[[1]]
            Fact.list <- Fact.list[str_detect(Fact.list, pattern=":", negate = TRUE)]
            
            #  => list of fact bio
            Fact.Bio  <- names(object@Factors.Type[which(object@Factors.Type == "Bio")])
            
            #  => factor vectors to set model.matrix `
            #  * stock ref values
            ref.list <- vector()
            for (i in Fact.list) {
              assign(i, object@List.Factors[[i]])
              ref.list <- c(ref.list,levels(object@List.Factors[[i]])[1])
            }
            
            #  => set model matrix
            #  * change colnames
            model.design.matrix <- stats::model.matrix(as.formula(object@Model.formula))
            design.colnames <- colnames(model.design.matrix)
            #  colnames
            for (i in names(object@List.Factors)) {
              design.colnames <- gsub(i, "", design.colnames)
            }
            design.colnames <- gsub(":", "_", design.colnames)
            colnames(model.design.matrix) <- design.colnames
            object@Model.matrix <- model.design.matrix

            # => Get all the contrasts
            contrasts <- list()
            contrast.coef <- list()
            contrasts.nbr <- 0
            for(i in Fact.list[which(Fact.list %in% Fact.Bio)]){
              
              mat <- as.vector(unique(object@List.Factors[[i]])) %>% utils::combn(.,2)
              
              contrast.tmp <- list()
              
              contrast.tmp[["hypoth"]] <- apply (mat, 2, function(x) {
                paste(x, collapse=" - ")
              })
              
              contrast.tmp[["hypoth_tmp"]] <- contrast.tmp[["hypoth"]]
              
              for (ref in ref.list) {
                contrast.tmp[["hypoth_tmp"]] <- gsub(paste(ref, "-"), "", contrast.tmp[["hypoth_tmp"]])
                contrast.tmp[["hypoth_tmp"]] <- gsub(paste("-", ref), "", contrast.tmp[["hypoth_tmp"]])
              }
              
              contrast.tmp[["idContrast"]] <- paste("C", (contrasts.nbr+1):(contrasts.nbr+length(contrast.tmp[["hypoth"]])), sep = "")
              
              contrasts.nbr <- contrasts.nbr+length(contrast.tmp[["hypoth"]])
              contrasts[[i]] <- data.frame(contrast.tmp) %>% mutate(factors=i)
                   
              ## contrast coef
              #contrast.coef[[i]] <- makeContrasts(contrasts = contrasts[[i]]$hypoth_tmp, levels = model.design.matrix) %>% as.data.frame()
              #colnames(contrast.coef[[i]]) <- contrasts[[i]]$idContrast
            }
            object@Contrasts.List  <- contrasts %>% purrr::reduce(rbind)
            
            ## contrast coef
            object@Contrasts.Coeff <- makeContrasts(contrasts = object@Contrasts.List$hypoth_tmp, levels = model.design.matrix) %>% as.data.frame()
            colnames(object@Contrasts.Coeff) <- object@Contrasts.List$idContrast
            
            return(object)

})




#' @title CheckExpDesignCompleteness
#'
#' @param ExpDesign 
#' @return list
#' @exportMethod CheckExpDesignCompleteness
#'
#' @examples
setMethod(f="CheckExpDesignCompleteness",
            signature="ExpDesign",
            definition <- function(object){
  
  # output list
  output <- list()
  
  # check presence of bio factors
  if (! "Bio" %in% object@Factors.Type){
    
    message <- "noBio"
    
    group_count  <- object@List.Factors[object@Factors.Type == "batch"] %>% as.data.frame() %>% table() %>% as.data.frame()
    names(group_count)[names(group_count) == "Freq"] <- "Count"
    output[["count"]]   <- group_count
    
  }else{
    
    # count occurence of bio conditions
    group_count  <- object@List.Factors[object@Factors.Type == "Bio"] %>% as.data.frame() %>% table() %>% as.data.frame()
    names(group_count)[names(group_count) == "Freq"] <- "Count"
    
    output[["count"]]   <- group_count
    
    # check presence of relicat / batch
    # check if design is complete
    # check if design is balanced
    # check nbr of replicats
    message <- if_else(! "batch" %in% object@Factors.Type , "noBatch",
                       if_else(0 %in% group_count$Count ,   "noCompl", 
                               if_else(length(unique(group_count$Count)) != 1, "noBalan", 
                                       if_else(group_count$Count[1] < 3, "lowRep", "true"))))
  }
  
  
  # switch pour message complet 
  output[["message"]] <- switch(message ,
         "true"       = { c("true",    "The experimental design is complete and balanced.") },
         "lowRep"     = { c("warning", "WARNING : 3 biological replicates are needed.") },
         "noCompl"    = { c("false",   "ERROR : The experimental design is not complete.") },
         "noBalan"    = { c("warning", "WARNING : The experimental design is complete but not balanced.") },
         
         "noBio"      = { c("false",   "ERROR : no bio factor !") },
         "noBatch"    = { c("false",   "ERROR : no replicat") }
         )
  
  
  return(output)
})



######################################### Contrasts ############################################



#' @title getExpressionContrast
#' get simple, pairwise comparison and averaged expression contrast data frames, offer to the user the possibility to select for each type of contrast
#' the contrast he want to keep and bind the selected expression contrast data frames 
#'
#' @param ExpDesign 
#' @param model.formula 
#'
#' @return ExpDesign
#' @exportMethod getExpressionContrast
#'
#' @examples
#' @author Christine Paysant-Le Roux
#'
setMethod(f="getExpressionContrast",
          signature="ExpDesign",
          definition <- function(object, model.formula){

  # model formula
  modelFormula <- formula(model.formula)
  
  #Design@Model.formula <- formula(model.formula)
  object@Model.formula <- model.formula
  
  # bio factor list in formulat 
  labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
  
  FactorBioInDesign <- intersect(names(object@Factors.Type[object@Factors.Type == "Bio"]), labelsIntoDesign)
  
  #BioFactors <- object@List.Factors[FactorBioInDesign]
  
  treatmentFactorsList <- lapply(FactorBioInDesign, function(x){paste(x, unique(object@List.Factors[[x]]), sep="")})
  names(treatmentFactorsList) <- FactorBioInDesign
            
  interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
  # define all simple contrasts pairwise comparisons
  allSimpleContrast_df <- defineAllSimpleContrasts(treatmentFactorsList)
  listOfContrastsDF <- list(simple = allSimpleContrast_df)
  
  # define all simples contrast means
  # exists("allSimpleContrast_df", inherits = FALSE)
  if(length(treatmentFactorsList) != 1){
    allAveragedContrasts_df <- define_averaged_contrasts (allSimpleContrast_df)
    listOfContrastsDF[["averaged"]] <- allAveragedContrasts_df
  }
  
  # define all interaction contrasts
  if(length(treatmentFactorsList) != 1){
    if(interactionPresent){
      labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
      labelOrder <- attr(terms.formula(modelFormula), "order")
      twoWayInteractionInDesign <- labelsIntoDesign[which(labelOrder==2)]
      groupInteractionToKeep <- gsub(":", " vs ", twoWayInteractionInDesign)
      allInteractionsContrasts_df <- defineAllInteractionContrasts(treatmentFactorsList, groupInteractionToKeep)
      
      listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
    }
    #allInteractionsContrasts_df <- defineAllInteractionContrasts(treatmentFactorsList)
    #listOfContrastsDF[["interaction"]] <- allInteractionsContrasts_df
  }
  # choose the contrasts and rbind data frames of contrasts
  #selectedContrasts <- returnSelectedContrasts(listOfContrastsDF)
  
  # replace interactive selection of contrasts by return all contrasts -> shiny
  object@Contrasts.List  <- listOfContrastsDF
  object@Contrasts.Coeff <- data.frame()
  object@Contrasts.Sel   <- data.frame()
  
  return(object)
})


#' @title getContrastMatrix
#' define contrast matrix or contrast list with contrast name and contrast coefficients
#'
#' @param ExpDesign
#' @param contrastList
#' @return ExpDesign
#' @exportMethod getContrastMatrix
#'
#' @author Christine Paysant-Le Roux
setMethod(f="getContrastMatrix",
          signature="ExpDesign",
          definition <- function(object, contrastList){
            
  contrast.sel.list <- list()
  contrast.sel.list <- lapply(names(Design@Contrasts.List), function(contrastType) {
  
    tmp <- object@Contrasts.List[[contrastType]] %>% dplyr::filter(contrast %in% contrastList) %>%
                    dplyr::select(contrast, contrastName, type, groupComparison)
    return(tmp)
  })
  object@Contrasts.Sel <- contrast.sel.list %>% purrr::reduce(rbind) %>% dplyr::mutate(tag = paste("H", 1:dim(.)[1], sep=""))
            
            
  sampleData <-  object@ExpDesign         
  selectedContrasts <- object@Contrasts.Sel$contrast
  
  modelFormula <- formula(object@Model.formula)
  # bio factor list in formulat 
  labelsIntoDesign <- attr(terms.formula(modelFormula),"term.labels")
  FactorBioInDesign <- intersect(names(object@Factors.Type[object@Factors.Type == "Bio"]), labelsIntoDesign)
  
  #BioFactors <- object@List.Factors[FactorBioInDesign]
  
  treatmentFactorsList <- lapply(FactorBioInDesign, function(x){paste(x, unique(object@List.Factors[[x]]), sep="")})
  names(treatmentFactorsList) <- FactorBioInDesign
  
  treatmentCondenv <- new.env()
  
  interactionPresent <- any(attr(terms.formula(modelFormula),"order") > 1)
  isThreeOrderInteraction <- any(attr(terms.formula(modelFormula),"order") == 3)
  
  # get model matrix
  modelMatrix <- stats::model.matrix(modelFormula, data = object@List.Factors %>% as.data.frame())
  colnames(modelMatrix)[colnames(modelMatrix) == "(Intercept)"] <- "Intercept"
  # assign treatment conditions(group) to boolean vectors according to the design model matrix
  #treatmentCondenv <- new.env()
  assignVectorToGroups(treatmentFactorsList = treatmentFactorsList, modelMatrix = modelMatrix, interactionPresent = interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction, treatmentCondenv = treatmentCondenv)
  # get the coefficient vector associated with each selected contrast
  # contrast <- allSimpleContrast_df$contrast[1]
  colnamesGLMdesign <- colnames(modelMatrix)
  
  
  #coefficientsMatrix <- sapply(selectedContrasts$contrast, function(x) returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
  coefficientsMatrix <- sapply(selectedContrasts, function(x) returnContrastCoefficients(x, colnamesGLMdesign, treatmentCondenv = treatmentCondenv))
  
  #coefficientsMatrix <- MASS::as.fractions(coefficientsMatrix)
  colnames(coefficientsMatrix) <- selectedContrasts

  rownames(coefficientsMatrix) <- colnamesGLMdesign 
  contrastMatrix <- as.data.frame(t(coefficientsMatrix))
  #contrastMatrix <- as_tibble(t(coefficientsMatrix)) %>%
    #dplyr::mutate(contrast = selectedContrasts, .before = "Intercept") %>%
    #dplyr::mutate(type = selectedContrasts$type, .after = "contrast")
  #contrastMatrix <- MASS::as.fractions(contrastMatrix)
  #contrastMatrix
  # contrastList <- as.list(as.data.frame(coefficientsMatrix))
  
  object@Contrasts.Coeff <- contrastMatrix
  return(object)
})

#coefmatrices <- sapply(unique(names(coefvectors)),
#                       function(n) as.matrix(as.data.frame(coefvectors[names(coefvectors)==n])),
#                       simplify=FALSE, USE.NAMES=TRUE)

#Contrasts = list(D1vsD2          = c(1,  1, -1, -1,  0),
#                 C1vsC2          = c(1, -1,  1, -1,  0),
#                 InteractionDC   = c(1, -1, -1,  1,  0),
#                 C1vsC2forD1only = c(1, -1,  0,  0,  0),
#                 C1vsC2forD2only = c(0,  0,  1, -1,  0),
#                 TreatsvsControl = c(1,  1,  1,  1, -4),
#                 T1vsC           = c(1,  0,  0,  0, -1),
#                 T2vsC           = c(0,  1,  0,  0, -1),
#                 T3vsC           = c(0,  0,  1,  0, -1),
#                 T4vsC           = c(0,  0,  0,  1, -1))

# contrast From emmeans v1.3.5 by Russell Lenth 16th Percentile Contrasts and linear functions of EMMs
# coef returns a data.frame containing the object's grid, along with columns named c.1, c.2, ... containing the contrast coefficients. 




################################### CO-EXPRESSION #############################


#' @title runCoExpression
#' @param object MultiAssayExperiment
#' @param data dataset name
#' @param tools  
#' @param geneList 
#' @param K list of number of clusters
#' @param iter gene list 
#' @param model
#' @param transformation
#' @param normFactors
#' @return MultiAssayExperiment
#' @exportMethod runCoExpression
#'
setMethod(f="runCoExpression",
          signature="MultiAssayExperiment",
          definition <- function(object, data, tools = "coseq", geneList, K, iter=5 , model="normal", transformation="arcsin", normFactors="TMM"){
            
            counts = assay(object@ExperimentList[[data]])[geneList,] 
            
            switch (tools,
              "coseq" = {
                  coseq.res <- runCoseq(counts, K=K, iter=iter, model=model, transformation=transformation, normFactors=normFactors)
                  object@ExperimentList[[data]]@metadata[["CoExpResults"]][["coseqResults"]] <- coseq.res
                  
                  clusters <- lapply(1:length(table(clusters(coseq.res))), function(i){ 
                    names(clusters(coseq.res)[clusters(coseq.res) == i])
                    })
                  object@ExperimentList[[data]]@metadata[["CoExpResults"]][["clusters"]] <- clusters
                  names(object@ExperimentList[[data]]@metadata[["CoExpResults"]][["clusters"]]) <- paste("cluster", 1:length(table(clusters(coseq.res))), sep = ".")
                
                }
            )
              
      return(object)
})

################################### ANNOTATION #############################

#' @title runAnnotationEnrichment
#' @param object MultiAssayExperiment
#' @param data dataset name
#' @param alpha 
#' @param probaMethod 
#' @param annotation gene annotation
#' @param geneLists gene list 
#' @return MultiAssayExperiment
#' @exportMethod runAnnotationEnrichment
#'
setMethod(f="runAnnotationEnrichment",
          signature="MultiAssayExperiment",
          definition <- function(object, data, annotation, geneLists, alpha = 0.01, probaMethod = "hypergeometric"){
             
            Results <- list()
            
            for(geneList in names(geneLists)){
              
              Results[[geneList]] <- switch(probaMethod,
                   "hypergeometric"=EnrichmentHyperG(annotation, geneLists[[geneList]], alpha = 0.01)
                   )
            }
            object@ExperimentList[[data]]$metadata$AnnotEnrich <- Results
            
            return(object)

          })




