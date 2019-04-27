require(DESeq2)
require(pheatmap)
require(parallel)
require(genefilter)

DASeq = function(CountData,colData,formula = ~ Lane + condition,normRows=1,fitType=c("parametric","local","mean"),
     locType = c("shorth","median"),collapseBy=NULL){
     fitType = match.arg(fitType)
     locType = match.arg(locType)
     locfunc = switch(locType,
          shorth = genefilter::shorth,
          median = stats::median)
     if(sum(colnames(CountData[-1]) != colData$Names)){
          colData = colData[match(colnames(CountData[-1]),colData$Names),]
     }
     suppressMessages(das <- DESeqDataSetFromMatrix(CountData,colData,tidy = T,design = formula))
     if(!is.null(collapseBy)){
          das = collapseReplicates(das,collapseBy)
     }
     trial = try(das <- estimateSizeFactors(das,type="poscounts",locfunc = locfunc),silent = TRUE)
     if(inherits(trial,"try-error") & locType == "shorth"){
          message("Tie occured while caclculating size Factor location with shorth, switching to median")
          das = estimateSizeFactors(das,type="poscounts",locfunc = median)
     }
     das = estimateDispersions(das,fitType = fitType,quiet=T)
     das = nbinomWaldTest(das)
     das
}

PlotHeatmap = function(DAObj,formula = design(DAObj),normRows=1,row_annote=NULL,show_rownames = is.null(row_annote)){
     DAObj = DAObj[,do.call(order,colData(DAObj)[rev(all.vars(formula))])]
     df = as.data.frame(colData(DAObj)[all.vars(formula)])
     pheatmap(log2(counts(DAObj,normalized=T)+1)[1:(nrow(DAObj)-normRows),],annotation_row = row_annote,
              cluster_rows = T,show_rownames = show_rownames,cluster_cols = F,annotation_col = df,
              show_colnames = F,col= colorRampPalette(c("White","Blue","Red"))(100))
}

DAOrgList = function(DAObj,rNames = resultsNames(DAObj)[-1],threshold = 0.05){
     extractSig = function(r,DAObj){
          res = results(DAObj,name=r)
          res = res[ifelse(is.na(res$padj),FALSE,res$padj < threshold),]
          setNames(res$log2FoldChange,rownames(res))
     }
     sapply(rNames,extractSig,DAObj=DAObj)
}

SimDASeq = function(DAObj,normrows = 1,sim = 0,p = 1/(2+nrow(DAObj)-normrows),formula = design(DAObj),blank="None",
          fitType=c("parametric","local","mean"),locType = c("shorth","median")){
     fitType = match.arg(fitType)
     locType = match.arg(locType)
     BlankVar = all.vars(formula)[1]
     NoSim = counts(DAObj[,colData(DAObj)[,BlankVar] == blank])[-c(1:(nrow(DAObj)-normrows)),]
     CountData = as.data.frame(counts(DAObj))
     CountData$Organism=rownames(CountData)
     rownames(CountData) = NULL
     CountData = CountData[,c(ncol(CountData),1:(ncol(CountData)-1))]
     colData = as.data.frame(colData(DAObj))
     colData$Names = rownames(colData)
     rownames(colData) = NULL
     BlankCol = match(BlankVar,names(colData))
     if(sim > 0){
          for(i in 1:sim){
               Name = paste0("Pseudo_",i)
               CountData[,Name] = c(rnbinom(nrow(CountData)-normrows,1,1-p),NoSim)
               colData = rbind(colData,c(Name,rep(NA,ncol(colData)-1)))
               colData[nrow(colData),BlankCol] = blank
          }
     }
     DASeq(CountData,colData,normRows = normRows,formula = formula,fitType = fitType,locType = locType)
}

# SimDASeq = function(CountData,colData,condition="Spike",normrows=1,sim=0,p=1/(2+(nrow(CountData)-normrows)),collapseBy = NULL){
#      bunmapped = as.integer(rowMeans(CountData[-1][nrow(CountData),colData$condition == "Blank"]))
#      if(sim > 0){
#           colData$Lane = as.character(colData$Lane)
#           colData$condition = as.character(colData$condition)
#           for(i in 1:sim){
#                for(j in 1:2){
#                     Name = paste0("Pseudo_",i,"_L00",j)
#                     CountData[,Name] = c(rnbinom(nrow(CountData)-normrows,1,1-p),bunmapped)
#                     colData = rbind(colData,c(Name,paste0("L00",j),"Blank"))
#                }
#           }
#           colData$Lane = as.factor(colData$Lane)
#           colData$condition = as.factor(colData$condition)
#      }
#      obj = DASeq(CountData,colData)
#      if(!is.null(collapseBy)){
#           obj = collapseReplicates(das,collapseBy)
#      }
#      results(obj,resultsNames(obj)[2])
# }

DA_UnderSimulatedBlanks = function(DAObj,p,normrows=1,
          Threads = 8,Trials=100,MaxSim=5){
     cl = makeCluster(Threads)
     invisible(clusterEvalQ(cl,library(DESeq2)))
     invisible(clusterEvalQ(cl,library(pheatmap)))
     invisible(clusterExport(cl,c("SimDASeq","DASeq")))
     SimRes = function(i,DAObj,normrows,sim,p,thresh = 0.05){
          obj = SimDASeq(DAObj,normrows = normrows,sim,p)
          r = results(obj,name=resultsNames(obj)[2])
          result = (r$padj[1:(nrow(r)-normrows)] < thresh)
          result
     }
     SimrowSums = function(sim,DAObj,normrows){
          result = parSapply(cl,1:Trials,SimRes,DAObj = DAObj,normrows=normrows,sim=sim,p=p)
          apply(result,1,sum,na.rm=T)
     }
     Result = sapply(1:MaxSim,SimrowSums,DAObj=DAObj,normrows=normrows)
     stopCluster(cl)
     rownames(Result) = rownames(DAObj)[1:(nrow(DAObj)-normrows)]
     colnames(Result) = 1:MaxSim
     Result / Trials * 100
}

PlotSimulation = function(Data){
     barplot(t(Data),beside=T,col=colorRampPalette(c("white","blue"))(ncol(Data)),
          las=3,names.arg = sapply(strsplit(rownames(Data),"|",fixed=T),paste0,collapse="\n"),
          cex.names = 0.65,legend = TRUE, args.legend = list(x="topleft",horiz=TRUE,title="Simulated Blanks"),
          ylab = "Proportion of Simulations with Significant DA")
}


AddNormalizationRow = function(CountData,colData,NormValue = c("Total","Trimmed","NonHuman","Pseudo"),pseudo=1,nrowName="Normalization"){
     NormValue = match.arg(NormValue)
     NormRow = switch(NormValue,
          Total = colData$LibrarySize,
          Trimmed = colData$TrimmedSize,
          NonHuman = colData$TrimmedSize - colData$HumanReads,
          Pseudo = colData$PathogenReads + pseudo)
     NormRow = NormRow - colData$PathogenReads
     CountData[nrow(CountData)+1,] = c(0,NormRow)
     CountData[nrow(CountData),1] = nrowName
     CountData
}

SubSampleToTaxonLvl = function(CountData,TaxonLvl = c("Strain","Species","Genus"),delim="|",na.rm=TRUE,zero.rm=TRUE){
     TaxonLvl = match.arg(TaxonLvl)
     TaxonIndex = switch(TaxonLvl,
          Strain = 3,
          Species = 2,
          Genus = 1)
     tmp = strsplit(CountData[,1],delim,fixed=TRUE)
     tmp = lapply(tmp,"[",1:TaxonIndex)
     CountData[1] = sapply(tmp,paste0,collapse=delim)
     if(na.rm){
          TaxonhasNA = sapply(tmp,function(x){length(grep("^NA$",x)) > 0})
          CountData = CountData[!TaxonhasNA,]
     }
     CountData = aggregate(CountData[-1],by=list(Organism = CountData[,1]),FUN=sum)
     if(zero.rm){
          RowIsEmpty = apply(CountData[-1],1,sum,na.rm=T) == 0
          CountData = CountData[!RowIsEmpty,]
     }
     CountData
}
