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
  
  # groups
  groups <- tidyr::unite(as.data.frame(ExpDesign[typeList == "Bio"]), col="groups", sep="_", remove = TRUE) %>% 
            dplyr::mutate(samples = rownames(.))
  
  Design = new(Class = "ExpDesign",
               ExpDesign=ExpDesign,
               List.Factors=dF.List,
               Factors.Type=typeList,
               Groups=groups,
               Model.formula=vector(),
               Model.matrix=vector(),
               Contrasts.List=list(),
               Contrasts.Sel=data.frame(),
               Contrasts.Coeff=data.frame())
           
  return(Design)
}



#' @title FlomicsMultiAssay.constructor Class constructor
#' @description
#' @param inputs list of input data
#' @param Design ExpDesign class
#' @return An object of class [\code{\link{MultiAssayExperiment}}]
#' @name FlomicsMultiAssay.constructor
#' @rdname FlomicsMultiAssay.constructor
#' @export
FlomicsMultiAssay.constructor <- function(inputs, Design){ 
  
 # if input == NULL
  
  SummarizedExperimentList <- list()
  listmap  <- list()
  omicList <- list()
  k <- 0
  for (dataName in names(inputs)){
  
      k <- k+1
    
      ## construct SummarizedExperiment for each data
      abundance <- read.table(inputs[[dataName]][["dataFile"]], header = TRUE, row.names = 1)
      
      if(!is.null(inputs[[dataName]][["qcFile"]])){
        print("# ... metadata QC...")
        QCmat <- read.table(inputs[[dataName]][["qcFile"]], header = TRUE)
      }
      else{
        QCmat <- data.frame(primary = colnames(abundance),
                            colname = colnames(abundance),
                            stringsAsFactors = FALSE)
      }
      
      SummarizedExperimentList[[dataName]] <- SummarizedExperiment::SummarizedExperiment(assays  = S4Vectors::SimpleList(abundance=as.matrix(abundance)),
                                                                   colData = QCmat)
      
      # metadata for sampleMap for MultiAssayExperiment
      listmap[[dataName]] <- data.frame(primary = as.vector(SummarizedExperimentList[[dataName]]@colData$primary),
                                        colname = as.vector(SummarizedExperimentList[[dataName]]@colData$colname),
                                        stringsAsFactors = FALSE)
      
      #
      omicType <- inputs[[dataName]][["omicType"]]
      
      colnames <- c(names(omicList[[omicType]]), k)
      omicList[[omicType]] <- c(omicList[[omicType]] ,dataName)
      names(omicList[[omicType]]) <- colnames
  
  }
  
  FlomicsMultiAssay <- MultiAssayExperiment::MultiAssayExperiment(experiments = SummarizedExperimentList,
                                             colData     = Design@ExpDesign,
                                             sampleMap   = MultiAssayExperiment::listToMap(listmap),
                                             metadata    = list(design = Design,
                                                                colDataStruc = c(n_dFac = dim(Design@ExpDesign)[2], n_qcFac = 0),
                                                                omicList = omicList))
  
  return(FlomicsMultiAssay)
}



#' @title RunDiffAnalysis
#' @param An object of class [\code{\link{MultiAssayExperiment}]
#' @param data omic data type
#' @param DiffMethod Differential analysis method
#' @param contrastList list of contrast to test
#' @param FDR FDR threshold
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @return MultiAssayExperiment
#' @exportMethod RunDiffAnalysis
#' @examples
#'
setMethod(f="RunDiffAnalysis",
          signature="MultiAssayExperiment",
          definition <- function(object, data, FDR = 0.05, contrastList, DiffAnalysisMethod, clustermq=FALSE){

            object@ExperimentList[[data]]@metadata$DiffExpAnal <- list()
            
            Contrasts.Sel <- dplyr::filter(object@metadata$design@Contrasts.Sel, contrastName %in% contrastList)
            object@ExperimentList[[data]]@metadata$DiffExpAnal[["contrasts"]] <- Contrasts.Sel
            object@ExperimentList[[data]]@metadata$DiffExpAnal[["method"]]    <- DiffAnalysisMethod
            object@ExperimentList[[data]]@metadata$DiffExpAnal[["FDR"]]       <- FDR
            
            # Run the Diff analysis and get the results as a list of object depending of the
            ListOfDiffResults <- switch(DiffAnalysisMethod,
                                     "edgeRglmfit"=edgeR.AnaDiff(object, data, clustermq)
                                     )

            # Set an AnaDiff object to
            object@ExperimentList[[data]]@metadata$DiffExpAnal[["DGELRT"]] <- ListOfDiffResults

            #
            object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]] <- lapply(ListOfDiffResults, function(x){
              
              res<-edgeR::topTags(x, n = dim(x)[1])
              DEGs<- res$table[res$table$FDR <= FDR,]
              #DEGs<-res$table
              return(DEGs)
            })
            names(object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]]) <- names(ListOfDiffResults)

            ## merge results in bin matrix
            DEG_list <- lapply(1:length(object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]]), function(x){
              
              res <- object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]][[x]]
              tmp <- data.frame(DEG = rownames(res), bin = rep(1,length(rownames(res))))
              colnames(tmp) <- c("DEG", paste("H", x, sep=""))
              return(tmp)
            })
            names(DEG_list) <- names(object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]])
            
            object@ExperimentList[[data]]@metadata$DiffExpAnal[["mergeDGE"]] <- DEG_list %>% purrr::reduce(dplyr::full_join, by="DEG") %>% 
              dplyr::mutate_at(.vars = 2:(length(DEG_list)+1), .funs = function(x){dplyr::if_else(is.na(x), 0, 1)}) %>% data.table::data.table()
            
            return(object)
          })




#' DiffAnal.plot
#'
#' @param data 
#' @param pngFile 
#' @param FDRcutoff 
#' @return plot
#' @exportMethod DiffAnal.plot
#' @export
#'
#' @examples
setMethod(f="DiffAnal.plot",
          signature="MultiAssayExperiment",
          definition <- function(object, data, hypothesis){
            
            plots <- list()
            
            res      <- object@ExperimentList[[data]]@metadata$DiffExpAnal[["DGELRT"]][[hypothesis]]
            resTable <- object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]][[hypothesis]]
            FDR      <- object@ExperimentList[[data]]@metadata$DiffExpAnal[["FDR"]]
            
            res.FDR <- edgeR::topTags(res, n = dim(res)[1])
            plots[["MA.plot"]] <- MA.plot(data = res.FDR$table, FDRcutoff = FDR, pngFile =NULL)
            
            
            plots[["Pvalue.hist"]] <- pvalue.plot(data =resTable[resTable$FDR <= FDR,], pngFile =NULL)
            
            return(plots)
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
            p <- do.call(gridExtra::grid.arrange, out)
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
                       tidyr::unite(col="groups", sep="_", remove = FALSE) %>% dplyr::mutate(samples=sample_names)

            # normalized data
            
            if(! is.null(object[[dataType]]@metadata$Normalization)){
              pseudo  <- log2(scale(MultiAssayExperiment::assay(object[[dataType]]), center=FALSE,
                                    scale=object[[dataType]]@metadata$Normalization$coefNorm$norm.factors)+1) %>% reshape2::melt()
            }
            else{
              pseudo  <- log2(scale(MultiAssayExperiment::assay(object[[dataType]]), center=FALSE)+1) %>% reshape2::melt()
            }
            
            colnames(pseudo) <- c("feature", "samples", "value")
            pseudo_bis <- dplyr::full_join(pseudo, groups, by="samples")

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

            groups    <- tidyr::unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]),
                               col="groups", sep="_", remove = TRUE)$groups
            conditions<- object@colData %>% as.data.frame() %>% dplyr::mutate(samples=row.names(.), groups=groups)

            score     <- object[[data]]@metadata$PCAlist[[PCA]]$ind$coord[, PCs] %>% as.data.frame() %>%
                         dplyr::mutate(samples=row.names(.)) %>% dplyr::full_join(., conditions, by="samples")

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
                          dplyr::mutate(tag="1.Unnormalised data")
            score_norm <- object@listPCA$norm$ind$coord %>% data.table::melt %>%
                          dplyr::mutate(tag=paste("2.Normalised data : ", object@Normalization@Method,  sep=""))

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
            assayFilt  <- MultiAssayExperiment::assay(objectFilt)
            
            ## nbr of genes with 0 count
            genes_flt0  <- objectFilt[rowSums(assayFilt) <= 0, ]@NAMES
            
            ## remove 0 count 
            objectFilt  <- objectFilt[rowSums(assayFilt)  > 0, ]
            assayFilt   <- MultiAssayExperiment::assay(objectFilt)
           
            
            
            ## filter cpm
            NbReplicate <- table(tidyr::unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]), col="groups", sep="_"))
            NbConditions <- unique(tidyr::unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]), col="groups", sep="_"))$groups %>% length()
            
            switch(Filter_Strategy,
                   "NbConditions" = { keep <- rowSums(edgeR::cpm(assayFilt) >= CPM_Cutoff) >=  NbConditions },
                   "NbReplicates" = { keep <- rowSums(edgeR::cpm(assayFilt) >= CPM_Cutoff) >=  min(NbReplicate) },
                   "filterByExpr" = { dge  <- edgeR::DGEList(counts = assayFilt, genes = rownames(assayFilt))
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

            groups <- tidyr::unite(as.data.frame(object@colData[object@metadata$design@Factors.Type == "Bio"]), col="groups", sep="_")$groups

            coefNorm  = switch(NormMethod,
                               "TMM"=TMM.Normalization(MultiAssayExperiment::assay(object[[data]]), groups)
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
            pseudo <- log2(scale(MultiAssayExperiment::assay(object[[data]]),
                                      center=FALSE)+1)
            }
            else if(PCA=="norm"){
            pseudo <- log2(scale(MultiAssayExperiment::assay(object[[data]]),
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
              contrasts[[i]] <- data.frame(contrast.tmp) %>% dplyr::mutate(factors=i)
                   
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
    message <- dplyr::if_else(! "batch" %in% object@Factors.Type , "noBatch",
                              dplyr::if_else(0 %in% group_count$Count ,   "noCompl", 
                                             dplyr::if_else(length(unique(group_count$Count)) != 1, "noBalan", 
                                                            dplyr::if_else(group_count$Count[1] < 3, "lowRep", "true"))))
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
          definition <- function(object, data, tools = "coseq", geneList, K, iter=5 , model="normal", 
                                 transformation="arcsin", normFactors="TMM", nameList, merge="union"){
            
            object@ExperimentList[[data]]@metadata$CoExpAnal <- list()
            object@ExperimentList[[data]]@metadata$CoExpAnal[["model"]]            <- model
            object@ExperimentList[[data]]@metadata$CoExpAnal[["transformation"]]   <- transformation
            object@ExperimentList[[data]]@metadata$CoExpAnal[["normFactors"]]      <- normFactors
            object@ExperimentList[[data]]@metadata$CoExpAnal[["meanFilterCutoff"]] <- 50
            object@ExperimentList[[data]]@metadata$CoExpAnal[["gene.list.names"]]  <- nameList
            object@ExperimentList[[data]]@metadata$CoExpAnal[["merge.type"]]       <- merge
            
            counts = MultiAssayExperiment::assay(object@ExperimentList[[data]])[geneList,] 
            
            switch (tools,
              "coseq" = {
                  coseq.res <- runCoseq(counts, K=K, iter=iter, model=model, transformation=transformation, normFactors=normFactors)
                  object@ExperimentList[[data]]@metadata$CoExpAnal[["coseqResults"]] <- coseq.res
                  
                  # list of genes per cluster
                  clusters <- lapply(1:length(table(coseq::clusters(coseq.res))), function(i){ 
                    names(coseq::clusters(coseq.res)[coseq::clusters(coseq.res) == i])
                    })
                  object@ExperimentList[[data]]@metadata$CoExpAnal[["clusters"]] <- clusters
                  names(object@ExperimentList[[data]]@metadata$CoExpAnal[["clusters"]]) <- paste("cluster", 1:length(table(coseq::clusters(coseq.res))), sep = ".")
                  
                  # nbr of cluster
                  nb_cluster <- coseq.res@metadata$nbCluster[min(coseq.res@metadata$ICL) == coseq.res@metadata$ICL]
                  object@ExperimentList[[data]]@metadata$CoExpAnal[["cluster.nb"]] <- nb_cluster
                
                  # plot
                  plot.coseq.res <- coseq::plot(coseq.res, conds = FlomicsMultiAssay@metadata$design@Groups$groups)
                  object@ExperimentList[[data]]@metadata$CoExpAnal[["plots"]] <- plot.coseq.res
                  
                }
            )
              
      return(object)
})

################################### ANNOTATION #############################

#' @title runAnnotationEnrichment
#' @param object MultiAssayExperiment
#' @param data dataset name
#' @param CoExpListNames
#' @param DiffListNames
#' @param alpha 
#' @param probaMethod 
#' @param annotation gene annotation
#' @return MultiAssayExperiment
#' @exportMethod runAnnotationEnrichment
#'
setMethod(f="runAnnotationEnrichment",
          signature="MultiAssayExperiment",
          definition <- function(object, data, DiffListNames = NULL, CoExpListNames = NULL, annotation, 
                                 alpha = 0.01, probaMethod = "hypergeometric"){
             
            EnrichAnal <- list()
            EnrichAnal[["list.names"]] <- c(CoExpListNames, DiffListNames)
            EnrichAnal[["alpha"]]      <- alpha
            EnrichAnal[["proba.test"]] <- probaMethod
            
            
            ## list of gene list to annotate
            geneLists <- list()
            geneLists.diff <- list()
            geneLists.diff <- lapply(DiffListNames, function(listname){
              
              row.names(object@ExperimentList[[data]]@metadata$DiffExpAnal[["TopDGE"]][[listname]])
            })
            names(geneLists.diff) <- DiffListNames
            
            geneLists.coseq <- list()
            #geneLists.coseq <- lapply(CoExpListNames, function(listname){
              
            geneLists.coseq <-  object@ExperimentList[[data]]@metadata[["CoExpAnal"]][["clusters"]][CoExpListNames]
            #})
            #names(geneLists.coseq) <- CoExpListNames
            
            geneLists <- c(geneLists.diff, geneLists.coseq)
            
            
            Results <- list()
            for(geneList in names(geneLists)){
              
              Results[[geneList]] <- switch(probaMethod,
                     "hypergeometric"=EnrichmentHyperG(annotation, geneLists[[geneList]], alpha = 0.01)
                     )
            }
            EnrichAnal[["results"]] <- Results
            
            object@ExperimentList[[data]]$metadata$EnrichAnal <- EnrichAnal
            
            return(object)
          })




