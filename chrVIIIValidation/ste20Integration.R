library(beeswarm)
library(seqinr)
load("data/phosphoProt.RData")
phosphoProtResidualsQTL <- phosphoProtResiduals
load("data/diffExpressionReplacementStrains.RData")
load("data/diffProtReplacementStrainsNoCor.RData")
load("data/diffPhosphoReplacementStrainsNoCor.RData")
load("data/diffRnaProtReplacementStrains.RData")
load("data/diffPhosphoProtReplacementStrainsNoCor.RData")
load("data/phosphoLevelReplacementStrainsNoCor.RData")
load("data/hsInfo.RData")
load("data/pQTL_results170815.RData")
strainLabels <- gsub("_.*","",colnames(phosphoLevel))
load("data/phosphoLevelQtlMappingData170727.RData")
load("data/phosphoLevel.RData")
load("data/expressionLevel.RData")
expressionLevelsQTL <- eBatchGeneLengthCorrected
load("data/expressionLevelsReplacementStrains.RData")
load("data/binInfo.RData")
load("data/eQTL_results160831.RData")
load("data/proteinLevel.RData")
load("data/protLevelsReplacementStrainsNoCor.RData")
load("data/pepLevelsReplacementStrains.RData")
load("data/plotPack.RData")
gtf <- read.table(file = "Saccharomyces_cerevisiae/RM/genes_RM_CDSonly_refined.gtf",as.is = T)
rawCountsQTL <- read.table("data/rawRNAseqCount.tsv",sep="\t",as.is=T)
genoWithLocs <- read.table("data/genotype_for_mapping.tsv",as.is=T,header=T)
u2e <- read.table("Saccharomyces_cerevisiae/uniprot2ens.tab",sep="\t",header = T,quote = "",as.is = T)
nickNames <- tolower(gsub(x = u2e[,4],pattern = "_YEAST",replacement = ""))
names(nickNames) <- u2e[,1]
revNames <- names(nickNames)
names(revNames) <- nickNames
meta <- read.table("metadata/metadata.tsv",sep="\t",header=T)
rownames(meta) <- meta$culture
BYProts   <- read.fasta("Saccharomyces_cerevisiae/sacCer3/REFINED/proteins_S288C_R6411.fa", as.string=TRUE, forceDNAtolower=FALSE)
rownames(phospho2prot) <- phospho2prot[,1]
rownames(genotype)[!grepl("X",rownames(genotype))] <- paste0("X",rownames(genotype)[!grepl("X",rownames(genotype))])
library(RColorBrewer)
colTrait <-brewer.pal(6,name="Paired")[-1]
names(colTrait) <-c("e","pt","p","phosphoRegressed", "phospho")
colBY <- brewer.pal(10,name="Paired")[8]
colRM <- brewer.pal(10,name="Paired")[10]

nPep <- read.table("ste20Validation/proteinData/proteinNoCor/mapDIA/protein_level.txt",sep="\t",header=T,as.is=T)
nPep <- nPep[,c(1,18,19)]
rownames(nPep) <- nPep[,1]
rawPepsVal <- read.table("ste20Validation/proteinData/phosphoNoCor/mapDIA/peptide_level.txt",sep="\t",header=T,as.is=T)
nPepsPh <- rawPepsVal[,c(2,19)]
rownames(nPepsPh) <- nPepsPh[,1]
#compute foldchanges

phByStrain <- sapply(rownames(genotype),FUN=function(s){
  rowMeans(phosphoLevelBatchCorrected[,intersect(meta[paste0("X",meta$strain)==s,1],colnames(phosphoLevelBatchCorrected)),drop=F],na.rm=T)
})
hsFC <- rowMeans(phByStrain[,genotype[,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,1597]==1],na.rm=T)
hsFC2 <- rowMeans(phByStrain[,genotype[,234]==0&genotype[,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,234]==0&genotype[,1597]==1],na.rm=T)
hsFC3 <- rowMeans(phByStrain[,genotype[,396]==0&genotype[,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,396]==0&genotype[,1597]==1],na.rm=T)
hsFC4 <- rowMeans(phByStrain[,genotype[,396]==1&genotype[,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,396]==1&genotype[,1597]==1],na.rm=T)
hsFC_M <- rowMeans(phByStrain[,genotype[,396]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,396]==1],na.rm=T)
hsFC2565 <- rowMeans(phByStrain[,genotype[,2565]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,2565]==1],na.rm=T)
hsFC_2565_0 <- rowMeans(phByStrain[,genotype[,2565]==0&genotype[,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,2565]==0&genotype[,1597]==1],na.rm=T)
hsFC_2565_1 <- rowMeans(phByStrain[,genotype[,2565]==1&genotype[,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[,2565]==1&genotype[,1597]==1],na.rm=T)
phrByStrain <- sapply(rownames(genotype),FUN=function(s){
  rowMeans(phosphoProtResidualsQTL[,intersect(meta[paste0("X",meta$strain)==s,1],colnames(phosphoProtResidualsQTL)),drop=F],na.rm=T)
})
hsFCPhr <- rowMeans(phrByStrain[,genotype[,1597]==0],na.rm=T)-rowMeans(phrByStrain[,genotype[,1597]==1],na.rm=T)
hsFCPhr2 <- rowMeans(phrByStrain[,genotype[,234]==0&genotype[,1597]==0],na.rm=T)-rowMeans(phrByStrain[,genotype[,234]==0&genotype[,1597]==1],na.rm=T)
hsFCPhr3 <- rowMeans(phrByStrain[,genotype[,396]==0&genotype[,1597]==0],na.rm=T)-rowMeans(phrByStrain[,genotype[,396]==0&genotype[,1597]==1],na.rm=T)
hsFCPhr4 <- rowMeans(phrByStrain[,genotype[,396]==1&genotype[,1597]==0],na.rm=T)-rowMeans(phrByStrain[,genotype[,396]==1&genotype[,1597]==1],na.rm=T)
hsFCPhr_M <- rowMeans(phrByStrain[,genotype[,396]==0],na.rm=T)-rowMeans(phrByStrain[,genotype[,396]==1],na.rm=T)

eByStrain <- sapply(rownames(genotype),FUN=function(s){
  rowMeans(eBatchGeneLengthCorrected[,intersect(meta[paste0("X",meta$strain)==s,1],colnames(eBatchGeneLengthCorrected)),drop=F],na.rm=T)
})
hsFcE <- rowMeans(eByStrain[,genotype[,1597]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,1597]==1],na.rm=T)
hsFcE_ira2 <- rowMeans(eByStrain[,genotype[,3091]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,3091]==1],na.rm=T)
hsFcE2 <- rowMeans(eByStrain[,genotype[,234]==0&genotype[,1597]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,234]==0&genotype[,1597]==1],na.rm=T)
hsFcE3 <- rowMeans(eByStrain[,genotype[,396]==0&genotype[,1597]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,396]==0&genotype[,1597]==1],na.rm=T)
hsFcE4 <- rowMeans(eByStrain[,genotype[,396]==1&genotype[,1597]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,396]==1&genotype[,1597]==1],na.rm=T)
hsFcE_M <- rowMeans(eByStrain[,genotype[,396]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,396]==1],na.rm=T)
hsFcE1137 <- rowMeans(eByStrain[,genotype[,1137]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,1137]==1],na.rm=T)
hsFcE2565 <- rowMeans(eByStrain[,genotype[,2565]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,2565]==1],na.rm=T)
hsFcE694 <- rowMeans(eByStrain[,genotype[,694]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,694]==1],na.rm=T)
hsFcE228 <- rowMeans(eByStrain[,genotype[,228]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,228]==1],na.rm=T)
pByStrain <- sapply(rownames(genotype),FUN=function(s){
  rowMeans(proteinLevelBatchCorrected[,intersect(meta[paste0("X",meta$strain)==s,1],colnames(proteinLevelBatchCorrected)),drop=F],na.rm=T)
})
hsFcP <- rowMeans(pByStrain[,genotype[,1597]==0],na.rm=T)-rowMeans(pByStrain[,genotype[,1597]==1],na.rm=T)
fcByLocE <- apply(genotype,2,FUN=function(x){
  rowMeans(eByStrain[,which(x==0)],na.rm=T)-rowMeans(eByStrain[,which(x==1)],na.rm=T)
})
fcByLocPh <- apply(genotype,2,FUN=function(x){
  rowMeans(phByStrain[,which(x==0)],na.rm=T)-rowMeans(phByStrain[,which(x==1)],na.rm=T)
})

eCorrected <- eByStrain
#2952
loc <- 2952
eCorrected[,genotype[,loc]==1] <- eCorrected[,genotype[,loc]==1]-rowMeans(eCorrected[,genotype[,loc]==1])+rowMeans(eCorrected[,genotype[,loc]==0])
#1137
loc <- 1137
eCorrected[,genotype[,loc]==1] <- eCorrected[,genotype[,loc]==1]-rowMeans(eCorrected[,genotype[,loc]==1])+rowMeans(eCorrected[,genotype[,loc]==0])
hsFcECor <- rowMeans(eCorrected[,genotype[,1597]==0],na.rm=T)-rowMeans(eCorrected[,genotype[,1597]==1],na.rm=T)

isect <- intersect(rownames(eByStrain),rownames(res_1473_1477))
fcMatE <- cbind(hs1=hsFcE[isect],hs2=hsFcE2[isect],hs3=hsFcE3[isect],hs4=hsFcE4[isect],hsM=hsFcE_M[isect],double=res_1473_1477[isect,"log2FoldChange"])
fcMatE <- as.data.frame(fcMatE)
plot(fcMatE)

#map phosphopeptides
locByPep <- t(apply(phospho2ProtRep,1,FUN=function(x){
  pep <- x[1]
  prot <- x[2]
  rawPep <- gsub("\\(UniMod:[0-9]*\\)","",pep,fixed=F)
  pepLen <- nchar(rawPep)
  pepStart <- nchar(gsub(paste0(rawPep,".*"),"",BYProts[[prot]]))+1
  pepRange <- c(pepStart,pepStart+pepLen-1)
  return(c(prot,pepRange))
}))

locByPepQTL <- t(apply(phospho2prot,1,FUN=function(x){
  pep <- x[1]
  prot <- x[2]
  rawPep <- gsub("\\([0-9]P\\)_[0-9]$","",pep,fixed=F)
  rawPep <- gsub("\\(UniMod:[0-9]*\\)","",rawPep,fixed=F)
  pepLen <- nchar(rawPep)
  pepStart <- nchar(gsub(paste0(rawPep,".*"),"",BYProts[[prot]]))+1
  pepRange <- c(pepStart,pepStart+pepLen-1)
  return(c(prot,pepRange))
}))
colnames(locByPep) <- colnames(locByPepQTL) <- c("prot","start","end")
rownames(locByPepQTL) <- phospho2prot[,1]

exactMatches <- t(sapply(rownames(phospho2prot),FUN=function(pep){
  strippedPep <- gsub("_.*","",pep)
  if(strippedPep%in%phospho2ProtRep[,1]){
    return(c(pep,strippedPep))
  }else{
    return(c(pep,NA))
  }
}))
exactMatchesSub <- exactMatches[!is.na(exactMatches[,2]),]
# oldToNewMatch <- matrix(F,ncol=nrow(locByPep),nrow=nrow(locByPepQTL))
# colnames(oldToNewMatch) <- rownames(locByPep)
# rownames(oldToNewMatch) <- rownames(locByPepQTL)
# 
# for(p in colnames(oldToNewMatch)){
#   #hits <- which(locByPepQTL[,1]==locByPep[p,1]&as.numeric(locByPepQTL[,2])<=as.numeric(locByPep[p,3])&as.numeric(locByPepQTL[,3])>=as.numeric(locByPep[p,2]))
#   hits <- which(locByPepQTL[,1]==locByPep[p,1]&as.numeric(locByPepQTL[,2])==as.numeric(locByPep[p,2])&as.numeric(locByPepQTL[,3])==as.numeric(locByPep[p,3]))
#   if(length(hits)>0){
#     oldToNewMatch[hits,p] <- T
#   }
# }
# 
# oldToNewMatchSub <- oldToNewMatch[rownames(pQTL_results$phosphoProtResidualsQTL$qv),rownames(phosphoProtResiduals)]
# 
# minByMatchNewToOld <- apply(oldToNewMatch,2,FUN=function(x){
#   min(pQTL_results$phosphoLevelQTL$pv[which(x),1597])
# })
# minByMatchOldToNew <- apply(oldToNewMatch,1,FUN=function(x){
#   min(diffPhospho[[3]][which(x),"p"])
# })
# minFdrByMatchOldToNew <- apply(oldToNewMatch,1,FUN=function(x){
#   if(!any(x)){return(NA)}
#   min(diffPhospho[[3]][which(x),"fdr"])
# })
# # plot(-log10(pQTL_results$phosphoLevelQTL$pv[,1597]),-log10(minByMatchOldToNew))
# # plot(-log10(diffPhospho[[3]][,"p"]),-log10(minByMatchNewToOld))
# 
# minFdrByMatchOldToNew1477_res <- apply(oldToNewMatchSub,1,FUN=function(x){
#   if(!any(x)){return(NA)}
#   min(diffPhosphoProt[[3]][which(x),"fdr"])
# })
# 
# sigFcByMatchNewToOld <- apply(oldToNewMatch,2,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   hsFC[locs[which.min(pQTL_results$phosphoLevelQTL$pv[locs,1597])]]
# })
# sigFcByMatchOldToNew <- apply(oldToNewMatch,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhospho[[3]][locs[order(diffPhospho[[3]][locs,"p"],decreasing=F,na.last = T)][1],"l2FC"]
# })
# 
# sigFcByMatchOldToNew1475 <- apply(oldToNewMatch,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhospho[[1]][locs[order(diffPhospho[[1]][locs,"p"],decreasing=F,na.last = T)][1],"l2FC"]
# })
# sigFcByMatchOldToNew1476 <- apply(oldToNewMatch,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhospho[[2]][locs[order(diffPhospho[[2]][locs,"p"],decreasing=F,na.last = T)][1],"l2FC"]
# })
# sigFcByMatchOldToNew1477 <- apply(oldToNewMatch,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhospho[[3]][locs[order(diffPhospho[[3]][locs,"p"],decreasing=F,na.last = T)][1],"l2FC"]
# })
# 
# sigFcByMatchOldToNew1475_res <- apply(oldToNewMatchSub,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhosphoProt[[1]][locs[order(diffPhosphoProt[[1]][locs,"p"],decreasing=F,na.last = T)][1],"fc"]
# })
# sigFcByMatchOldToNew1476_res <- apply(oldToNewMatchSub,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhosphoProt[[2]][locs[order(diffPhosphoProt[[2]][locs,"p"],decreasing=F,na.last = T)][1],"fc"]
# })
# sigFcByMatchOldToNew1477_res <- apply(oldToNewMatchSub,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   diffPhosphoProt[[3]][locs[order(diffPhosphoProt[[3]][locs,"p"],decreasing=F,na.last = T)][1],"fc"]
# })
# 
# # plot(sigFcByMatchNewToOld,diffPhospho[[3]][,"l2FC"])
# # plot(sigFcByMatchOldToNew,hsFC)
# 
# meanFcByMatchNewToOld <- apply(oldToNewMatch,2,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   mean(hsFC[locs],na.rm=T)
# })
# meanFcByMatchOldToNew <- apply(oldToNewMatch,1,FUN=function(x){
#   locs <- which(x)
#   if(!any(x)){return(NA)}
#   mean(diffPhospho[[3]][locs,"l2FC"],na.rm=T)
# })

# plot(meanFcByMatchNewToOld,diffPhospho[[3]][,"l2FC"])
# plot(meanFcByMatchOldToNew,hsFC)

####plot changes specific for each omics level between by and mutated strains####
##transcripts
mainVec <- c("BY vs BY:STE20-RM","BY vs BY:GPA1-RM","BY vs BY:STE20-RM,GPA1-RM")
resList <- list(res_1473_1475,res_1473_1476,res_1473_1477)
#Volcano plots with all transcripts
pdf("graph/alleleReplacement/byVsRep_transcripts_volcano_all.pdf")
par(mfrow=c(2,2))
invisible(sapply(1:length(resList),FUN=function(i){
  mat <- resList[[i]]
  plot(mat[,"log2FoldChange"],-log10(mat[,"pvalue"]),pch=20,col=rgb(0,0,0,0.2),main=mainVec[i],xlab="l2FC",ylab="-log10(p)")
}))
dev.off()
pdf("graph/alleleReplacement/byVsRep_transcripts_volcano_all_color.pdf")
par(mfrow=c(2,2))
invisible(sapply(1:length(resList),FUN=function(i){
  mat <- resList[[i]]
  colvec <- rownames(mat)%in%rownames(eByStrain)
  plot(mat[,"log2FoldChange"],-log10(mat[,"pvalue"]),pch=20,col=rgb(as.numeric(!colvec),0,0,0.2),main=mainVec[i],xlab="l2FC",ylab="-log10(p)")
}))
dev.off()
#beeswarm of fold changes for transcripts in and out of the geneset used for eqtl mapping
par(mfrow=c(2,2))
invisible(sapply(1:length(resList),FUN=function(i){
  mat <- resList[[i]]
  colvec <- rownames(mat)%in%rownames(eByStrain)
  boxplot(abs(mat[,"log2FoldChange"])~as.factor(colvec),outline=T)
  #beeswarm(as.vector(mat[,"log2FoldChange"])~as.factor(colvec),add=T)
}))
#Volcano plots with shared transcripts
pdf("graph/alleleReplacement/byVsRep_transcripts_volcano_intersect.pdf")
par(mfrow=c(2,2))
invisible(sapply(1:length(resList),FUN=function(i){
  mat <- resList[[i]]
  mat <- mat[intersect(rownames(mat),rownames(eByStrain)),]
  plot(mat[,"log2FoldChange"],-log10(mat[,"pvalue"]),pch=20,col=rgb(0,0,0,0.2),main=mainVec[i],xlab="l2FC",ylab="-log10(p)")
}))
dev.off()
#Volcano plots with different transcripts
pdf("graph/alleleReplacement/byVsRep_transcripts_volcano_diff.pdf")
par(mfrow=c(2,2))
invisible(sapply(1:length(resList),FUN=function(i){
  mat <- resList[[i]]
  mat <- mat[setdiff(rownames(mat),rownames(eByStrain)),]
  plot(mat[,"log2FoldChange"],-log10(mat[,"pvalue"]),pch=20,col=rgb(0,0,0,0.2),main=mainVec[i],xlab="l2FC",ylab="-log10(p)")
}))
dev.off()

####plot changes on omics level between by and mutated strains in comparison to hs effects####
##transcripts
#all
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
pdf("graph/alleleReplacement/repVsHS_transcripts_all.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()


par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE4[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE4[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE4[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)

onlyRep <- isectGenes[which(res_1473_1477[isectGenes,"padj"]<0.1&pv[isectGenes,1597]>0.8)]

#sig for HS

isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- intersect(isectGenes,targetsByHS[[10]]$e)
pdf("graph/alleleReplacement/repVsHS_transcripts_sigHS.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()
round(cor(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes]),digits=2)
round(cor(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes]),digits=2)
cor(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes])

#sig for HS not defined by bin but by overlap

isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
ovHits <- sapply(QtlList$FDR10,FUN=function(qtl){
  qtlMarkers <- unlist(lapply(1:nrow(qtl$predictors),FUN=function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))
  if(1597%in%qtlMarkers){
    return(rownames(qv)[qtl$target])
  }else{
    return(NA)
  }
})
ovHits <- ovHits[!is.na(ovHits)]
isectGenes <- intersect(isectGenes,ovHits)
pdf("graph/alleleReplacement/repVsHS_transcripts_sigHSNarrow.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes])
cor(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes])
cor(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes])

#sig for 1597
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- isectGenes[which(qv[isectGenes,1597]<0.1)]
pdf("graph/alleleReplacement/repVsHS_transcripts_sig1597.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes])
cor(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes])
cor(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes])

par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE3[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE3[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE3[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)


#sig for double replacement
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- isectGenes[res_1473_1477[isectGenes,"padj"]<0.1]
pdf("graph/alleleReplacement/repVsHS_transcripts_sig1477.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1476[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(res_1473_1477[isectGenes,"log2FoldChange"],hsFcE[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()

#compare effects of 1475 to 1476
pdf("graph/alleleReplacement/gpa1_ste20_comparison_transcripts.pdf")
par(mfrow=c(2,2))
gs <- rownames(res_1473_1475)
plot(res_1473_1475[gs,"log2FoldChange"],res_1473_1476[gs,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(GPA1 replacement)",xlab="l2FC(STE20 replacement)",main="all transcripts")
abline(a=0,b=1,col="red")
gs <- intersect(rownames(res_1473_1475),targetsByHS[[10]]$e)
plot(res_1473_1475[gs,"log2FoldChange"],res_1473_1476[gs,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(GPA1 replacement)",xlab="l2FC(STE20 replacement)",main="HS targets")
abline(a=0,b=1,col="red")
gs <- intersect(rownames(res_1473_1475),rownames(qv)[qv[,1597]<0.1])
plot(res_1473_1475[gs,"log2FoldChange"],res_1473_1476[gs,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(GPA1 replacement)",xlab="l2FC(STE20 replacement)",main="HS targets (narrow definition)")
abline(a=0,b=1,col="red")
gs <- rownames(res_1473_1475)[which(res_1473_1477$padj<0.1)]
plot(res_1473_1475[gs,"log2FoldChange"],res_1473_1476[gs,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(GPA1 replacement)",xlab="l2FC(STE20 replacement)",main="affected by double replacement")
abline(a=0,b=1,col="red")
dev.off()

#compare the fold-changes caused by the double replacement to the single replacements
pdf("graph/alleleReplacement/single_vs_double_transcripts.pdf",height=3.5)
par(mfrow=c(1,2))
gs <- rownames(res_1473_1475)[which(res_1473_1477$padj<0.1)]
plot(res_1473_1475[gs,"log2FoldChange"],res_1473_1477[gs,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(double replacement)",xlab="l2FC(single replacement)",main="STE20 replacement")
abline(a=0,b=1,col="red")
plot(res_1473_1476[gs,"log2FoldChange"],res_1473_1477[gs,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(double replacement)",xlab="l2FC(single replacement)",main="GPA1 replacement")
abline(a=0,b=1,col="red")
dev.off()
#compare the double replacement to the sum of the single replacements
pdf("graph/alleleReplacement/expected_vs_observed_double_replacement_transcripts.pdf",height=3.5)
par(mfrow=c(1,2))
gs <- rownames(res_1473_1475)[which(res_1473_1477$padj<0.1)]
plot(res_1473_1475[gs,"log2FoldChange"]+res_1473_1476[gs,"log2FoldChange"],res_1473_1477[gs,"log2FoldChange"],ylab="observed l2FC",xlab="expected l2FC",pch=19,col=rgb(0,0,0,0.3),main="multiplicative model")
abline(a=0,b=1,col="red")
plot(log2(2^(res_1473_1475[gs,"log2FoldChange"])+2^(res_1473_1476[gs,"log2FoldChange"])-1),res_1473_1477[gs,"log2FoldChange"],ylab="observed l2FC",xlab="expected l2FC",pch=19,col=rgb(0,0,0,0.3),main="additive model")
abline(a=0,b=1,col="red")
dev.off()

#What genes are affected by the double replacement but not the HS?
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- isectGenes[which(res_1473_1477[isectGenes,"padj"]<0.1&!isectGenes%in%targetsByHS[[10]]$e)]
boxplot(-log10(pv[,1597])~as.factor(rownames(qv)%in%isectGenes))
wilcox.test(-log10(pv[,1597])~as.factor(rownames(qv)%in%isectGenes))
isectGenes <- isectGenes[which(res_1473_1477[isectGenes,"padj"]<0.1&!pv[isectGenes,1597]>0.5)]
cat(isectGenes,sep="\n")
plot(eByStrain[isectGenes,"XBY4716"],rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))]))
abline(a=0,b=1)

#is there a difference how well the rep by values fit with qtl by between significant and ns genes?
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- isectGenes[!is.na(res_1473_1477[isectGenes,"padj"])]
plot(rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T),eByStrain[isectGenes,"XBY4716"],pch=19,col=rgb(res_1473_1477[isectGenes,"padj"]<0.2,0,0,0.2+(res_1473_1477[isectGenes,"padj"]<0.2)/2))
plot(rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T)-eByStrain[isectGenes,"XBY4716"],res_1473_1477[isectGenes,"log2FoldChange"])
#answer: the differences between 1473 and 1477 are not due to measurement errors in 1473 compared to by

#comparison of variation between rep strains and in cross
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- isectGenes[!is.na(res_1473_1477[isectGenes,"padj"])]
plot(apply(eByStrain[isectGenes,],1,sd,na.rm=T),abs(res_1473_1477[isectGenes,"log2FoldChange"]),pch=19,col=rgb(res_1473_1477[isectGenes,"padj"]<0.1,0,0,(res_1473_1477[isectGenes,"padj"]<0.1)/2+0.2)) #some transcripts are stable in the cross but not the rep strains
#answer: transcripts for sterol synthesis are much more variable in the rep strains

##proteins
#all available
isectGenes <- intersect(rownames(prot),names(hsFcP))
pdf("graph/alleleReplacement/repVsHS_proteins_all.pdf")
par(mfrow=c(2,2))
plot(diffProt[[1]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[2]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[3]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()
#sig for HS
isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- intersect(isectGenes,targetsByHS[[10]]$p)
pdf("graph/alleleReplacement/repVsHS_proteins_sigHS.pdf")
par(mfrow=c(2,2))
plot(diffProt[[1]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[2]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[3]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()
#sig for 1597
isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- isectGenes[which(pQTL_results$pQTL$qv[isectGenes,1597]<0.1)]
pdf("graph/alleleReplacement/repVsHS_proteins_sig1597.pdf")
par(mfrow=c(2,2))
plot(diffProt[[1]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[2]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[3]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()

#best 100 for 1597
isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- isectGenes[order(pQTL_results$pQTL$pv[isectGenes,1597],na.last = T)[1:100]]
pdf("graph/alleleReplacement/repVsHS_proteins_best1597.pdf")
par(mfrow=c(2,2))
plot(diffProt[[1]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[2]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[3]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()


#sig for double replacement
isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- isectGenes[which(diffProt[[3]][isectGenes,"fdr"]<0.1)]
pdf("graph/alleleReplacement/repVsHS_proteins_sig1477.pdf")
par(mfrow=c(2,2))
plot(diffProt[[1]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[2]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffProt[[3]][isectGenes,"l2FC"],hsFcP[isectGenes],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()

cor.test(diffProt[[1]][isectGenes,"l2FC"],hsFcP[isectGenes])
cor.test(diffProt[[2]][isectGenes,"l2FC"],hsFcP[isectGenes])
cor.test(diffProt[[3]][isectGenes,"l2FC"],hsFcP[isectGenes]) 


isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- isectGenes[which(diffProt[[3]][isectGenes,"fdr"]<0.1)]

#look at changes of transcripts and proteins
isectGenes <- intersect(rownames(prot),rownames(expressionLevels))
pdf("graph/alleleReplacement/transcripts_vs_proteins_all.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],diffProt[[1]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="STE20 replacement")
plot(res_1473_1476[isectGenes,"log2FoldChange"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="GPA1 replacement")
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="double replacement")
dev.off()

isectGenes <- intersect(rownames(prot),rownames(expressionLevels))
isectGenes <- isectGenes[which(res_1473_1477[isectGenes,"padj"]<0.1)]
pdf("graph/alleleReplacement/transcripts_vs_proteins_sig_double.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],diffProt[[1]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="STE20 replacement")
plot(res_1473_1476[isectGenes,"log2FoldChange"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="GPA1 replacement")
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="double replacement")
dev.off()
cor(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"])

isectGenes <- intersect(rownames(prot),rownames(expressionLevels))
isectGenes <- isectGenes[which(abs(res_1473_1477[isectGenes,"log2FoldChange"])>.5)]
pdf("graph/alleleReplacement/transcripts_vs_proteins_strong_fc.pdf")
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],diffProt[[1]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="STE20 replacement")
plot(res_1473_1476[isectGenes,"log2FoldChange"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="GPA1 replacement")
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(protein)",xlab="l2FC(transcript)",main="double replacement")
dev.off()
cor.test(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"])



isectGenes <- intersect(rownames(prot),rownames(expressionLevels))
isectGenes <- intersect(isectGenes,rownames(expressionLevelsQTL))
isectGenes <- isectGenes[which(diffProt[[3]][isectGenes,"fdr"]<0.1&res_1473_1477[isectGenes,"padj"]<0.1)]
plot(hsFcE[isectGenes],diffProt[[3]][isectGenes,"l2FC"]) #no correlation

isectGenes <- intersect(rownames(prot),rownames(expressionLevels))
isectGenes <- intersect(isectGenes,rownames(pByStrain))
isectGenes <- isectGenes[which(diffProt[[3]][isectGenes,"fdr"]<0.1)]
plot(hsFcP[isectGenes],diffProt[[3]][isectGenes,"l2FC"]) #weak stuff

#What proteins change their levels independently of their transcripts?
# isectGenes <- rownames(diffRnaProt[[3]])[which(diffRnaProt[[3]][,"fdr"]<0.2)]
# par(mfrow=c(1,1))
# plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"])
# cat(isectGenes,sep="\n") #enriched for de groot publication as having pt effects
# 
# cat(isectGenes[which(diffProt[[3]][isectGenes,"l2FC"]>0)],sep="\n")  #nothing


isectGenes <- intersect(rownames(expressionLevels),rownames(prot))
isectGenes <- isectGenes[which(diffProt[[3]][isectGenes,"fdr"]<0.1&res_1473_1477[isectGenes,"padj"]>0.5)]
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"])

#Do significant proteins change similarly in each single replacment strain?
isectGenes <- rownames(prot)[which(diffProt[[3]][,"fdr"]<0.1)]
pdf("graph/alleleReplacement/prot_fc_sig1477.pdf")
par(mfrow=c(2,2))
plot(diffProt[[1]][isectGenes,"l2FC"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(double)",xlab="l2FC(single)",main="STE20 vs STE20+GPA1")
abline(a=0,b=1,col="red")
plot(diffProt[[2]][isectGenes,"l2FC"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(double)",xlab="l2FC(single)",main="GPA1 vs STE20+GPA1")
abline(a=0,b=1,col="red")
plot(diffProt[[1]][isectGenes,"l2FC"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC(GPA1)",xlab="l2FC(STE20)",main="GPA1 vs STE20")
abline(a=0,b=1,col="red")
dev.off()
#is there a difference how well the rep by values fit with qtl by between significant and ns genes?
isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- isectGenes[!is.na(diffProt[[3]][isectGenes,"fdr"])]
pvec <- diffProt[[3]][isectGenes,"fdr"]
plot(rowMeans(prot[isectGenes,grepl("1473",colnames(prot))],na.rm=T),pByStrain[isectGenes,"XBY4716"],pch=19,col=rgb(pvec<0.2,0,0,0.2+(pvec<0.2)/2))
plot(rowMeans(prot[isectGenes,grepl("1473",colnames(prot))])-pByStrain[isectGenes,"XBY4716"],diffProt[[3]][isectGenes,"l2FC"])

isectGenes <- intersect(rownames(prot),names(hsFcP))
isectGenes <- isectGenes[!is.na(diffProt[[3]][isectGenes,"fdr"])]
pvec <- diffProt[[3]][isectGenes,"fdr"]
plot(apply(pByStrain[isectGenes,],1,sd,na.rm=T),abs(diffProt[[3]][isectGenes,"l2FC"]),pch=19,col=rgb(pvec<0.2,0,0,(pvec<0.2)/2+0.2)) #some transcripts are stable in the cross but not the rep strains
#answer: some proteins (lipid metabolism) are more variable between the rep strains than in the cross
##phospho

#look at general distribution of fold changes
pdf("graph/alleleReplacement/phosphoLevelsVolcano.pdf")
par(mfrow=c(2,2))
plot(diffPhospho[[1]][,"l2FC"],-log10(diffPhospho[[1]][,"p"]),pch=19,col=rgb(0,0,0,0.3),xlab="l2FC",ylab="-log10(p)",main="STE20")
plot(diffPhospho[[2]][,"l2FC"],-log10(diffPhospho[[2]][,"p"]),pch=19,col=rgb(0,0,0,0.3),xlab="l2FC",ylab="-log10(p)",main="GPA1")
plot(diffPhospho[[3]][,"l2FC"],-log10(diffPhospho[[3]][,"p"]),pch=19,col=rgb(0,0,0,0.3),xlab="l2FC",ylab="-log10(p)",main="STE20+GPA1")
#plot(diffPhospho[[4]][,"l2FC"],-log10(diffPhospho[[4]][,"p"]),pch=19,col=rgb(0,0,0,0.3),xlab="l2FC",ylab="-log10(p)",main="GPA1 in STE20Delta background")
dev.off()

#compare foldchanges across comparisons
pdf("graph/alleleReplacement/replacementComparisonPhospho.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhospho[[1]][,"l2FC"],diffPhospho[[2]][,"l2FC"],xlab="l2FC in STE20 replacement",ylab="l2FC in GPA1 replacement",pch=19,col=rgb(0,0,0,0.3),main="STE20 vs GPA1")
abline(a=0,b=1,col="red")
plot(diffPhospho[[1]][,"l2FC"],diffPhospho[[3]][,"l2FC"],xlab="l2FC in STE20 replacement",ylab="l2FC in GPA1+STE20 replacement",pch=19,col=rgb(0,0,0,0.3),main="STE20 vs GPA1+STE20")
abline(a=0,b=1,col="red")
plot(diffPhospho[[2]][,"l2FC"],diffPhospho[[3]][,"l2FC"],xlab="l2FC in GPA1 replacement",ylab="l2FC in GPA1+STE20 replacement",pch=19,col=rgb(0,0,0,0.3),main="GPA1 vs GPA1+STE20")
abline(a=0,b=1,col="red")
plot(diffPhospho[[2]][,"l2FC"],diffPhospho[[4]][,"l2FC"],xlab="BY background",ylab="RM background",pch=19,col=rgb(0,0,0,0.3),main="GPA1 replacement in two STE20 backgrounds")
abline(a=0,b=1,col="red")
plot(diffPhospho[[1]][,"l2FC"],diffPhospho[[5]][,"l2FC"],xlab="BY background",ylab="RM background",pch=19,col=rgb(0,0,0,0.3),main="STE20 replacement in two GPA1 backgrounds")
abline(a=0,b=1,col="red")
dev.off()

cor.test(diffPhospho[[1]][,"l2FC"],diffPhospho[[2]][,"l2FC"],use='pair')
cor.test(diffPhospho[[1]][,"l2FC"],diffPhospho[[3]][,"l2FC"],use='pair')
cor.test(diffPhospho[[2]][,"l2FC"],diffPhospho[[3]][,"l2FC"],use='pair')
cor.test(diffPhospho[[2]][,"l2FC"],diffPhospho[[4]][,"l2FC"],use='pair')
cor.test(diffPhospho[[1]][,"l2FC"],diffPhospho[[5]][,"l2FC"],use='pair')

#look at which fold changes fit the best with the phospho data
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
fcMat <- cbind(hs1=hsFC[sites],hs2=hsFC2[sites],hs3=hsFC3[sites],hs4=hsFC4[sites],hsM=hsFC_M[sites],double=diffPhospho[[3]][exactMatches[sites,2],"l2FC"])
fcMat <- as.data.frame(fcMat)
plot(fcMat)
#all matched sites from the qtl data
# sites <- names(sigFcByMatchOldToNew)[!is.na(sigFcByMatchOldToNew)]
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
pdf("graph/alleleReplacement/repVsHS_phospho_exact_all.pdf")
# par(mfrow=c(2,2))
# plot(sigFcByMatchOldToNew1475[sites],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
# abline(h=0,v=0,col="grey",lty=2)
# plot(sigFcByMatchOldToNew1476[sites],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
# abline(h=0,v=0,col="grey",lty=2)
# plot(sigFcByMatchOldToNew1477[sites],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
# abline(h=0,v=0,col="grey",lty=2)

par(mfrow=c(2,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)

cor.test(diffPhospho[[1]][exactMatchesSub[sites,2],"l2FC"],hsFC[sites])
cor.test(diffPhospho[[2]][exactMatchesSub[sites,2],"l2FC"],hsFC[sites])
cor.test(diffPhospho[[3]][exactMatchesSub[sites,2],"l2FC"],hsFC[sites])

dev.off()
#sig by HS
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
round(corVec[1:2],digits=2)
nameVec <- c("STE20","GPA1","STE20 + GPA1")
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sigHS.pdf",height = 7)
par(mfrow=c(2,2))
invisible(sapply(1:3,FUN=function(i){
  mat <- diffPhospho[[i]]
  plot(mat[exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0(nameVec[i]," replacement,\n r=",round(corVec[i],digits=2)),xlim=range(mat[exactMatchesSub[sites,2],1],na.rm=T)*1.3,ylim=range(hsFC[sites],na.rm=T)*1.3)
  abline(h=0,v=0,col="grey",lty=2)
  text(x=mat[exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5,cex = 0.7)
}))
# 
# plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
# abline(h=0,v=0,col="grey",lty=2)
# text(x=diffPhospho[[1]][exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5)
# plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
# abline(h=0,v=0,col="grey",lty=2)
# text(x=diffPhospho[[2]][exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5)
# plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
# abline(h=0,v=0,col="grey",lty=2)
# text(x=diffPhospho[[3]][exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5)
# plot(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
# abline(h=0,v=0,col="grey",lty=2)
# plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
# abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites])

colVec <- ifelse(pQTL_results$phosphoLevelQTL$qv[sites,1597]<0.1,"green","blue")
pdf(file = "graph/alleleReplacement/phValidationWithLabels.pdf")
par(mfrow=c(1,1))
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=colVec,ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement",xlim=range(diffPhospho[[3]][exactMatchesSub[sites,2],1],na.rm=T)*1.1,ylim=range(hsFC[sites],na.rm=T)*1.1)
abline(h=0,v=0,col="grey",lty=2)
text(x=diffPhospho[[3]][exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5)
legend(x = "topleft",legend = c("significant for 1597", "significant for broader HS"),pch = 19,col=c("green","blue"))
dev.off()

plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],diffPhospho[[2]][exactMatchesSub[sites,2],1],pch=19,col=rgb(0,0,0,0.3))
abline(a=0,b=1,col="red")
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],diffPhospho[[3]][exactMatchesSub[sites,2],1],pch=19,col=rgb(0,0,0,0.3))
abline(a=0,b=1,col="red")
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],diffPhospho[[3]][exactMatchesSub[sites,2],1],pch=19,col=rgb(0,0,0,0.3))
abline(a=0,b=1,col="red")



#sig by HS, narrower definition
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
ovHits <- sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10,FUN=function(qtl){
  qtlMarkers <- unlist(lapply(1:nrow(qtl$predictors),FUN=function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))
  if(1597%in%qtlMarkers){
    return(rownames(pQTL_results$phosphoLevelQTL$qv)[qtl$target])
  }else{
    return(NA)
  }
})
ovHits <- ovHits[!is.na(ovHits)]
sites <- intersect(sites,ovHits)
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sigHSNarrow.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites])


#sig by HS and 1473 vs 1477
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
sites <- sites[which(diffPhospho[[3]][exactMatchesSub[sites,2],"fdr"]<0.1)]
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sigHS_sig1477.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
colVec <- ifelse(pQTL_results$phosphoLevelQTL$qv[sites,1597]<0.1,"green","blue")
pdf(file = "graph/alleleReplacement/phValidationWithLabels_HS_14733.pdf")
par(mfrow=c(1,1))
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=colVec,ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement",xlim=range(diffPhospho[[3]][exactMatchesSub[sites,2],1],na.rm=T)*1.1,ylim=range(hsFC[sites],na.rm=T)*1.1)
abline(h=0,v=0,col="grey",lty=2)
text(x=diffPhospho[[3]][exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5)
legend(x = "topleft",legend = c("significant for 1597", "significant for broader HS"),pch = 19,col=c("green","blue"))
dev.off()
#sig by HS as ph or phr, shown as ph
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1|apply(pQTL_results$phosphoProtResidualsQTL$qv[,binPerMarker%in%c(119,120)],1,min,na.rm=T)[sites]<0.1)]
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sigHS_ph_or_phr.pdf")
par(mfrow=c(2,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()

#sig by 1597
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(pQTL_results$phosphoLevelQTL$qv[sites,1597]<0.1)]
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sig1597.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites])

plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main="STE20 replacement in GPA1 background")
abline(h=0,v=0,col="grey",lty=2)

#sig by 1597 and many frags
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(pQTL_results$phosphoLevelQTL$qv[sites,1597]<0.1)]
sites <- sites[nPepsPh[exactMatchesSub[sites,2],2]>6]
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sig1597_fragfilter.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
#sig by 1473 vs 1477 but not the broader HS
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
hsSites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
sites <- sites[which(diffPhospho[[3]][exactMatches[sites,2],3]<0.1)]
sites <- setdiff(sites,hsSites)
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phospho_exact_sig1477.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS l2FC",xlab="replacement l2FC",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhospho[[1]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[2]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[3]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[4]][exactMatchesSub[sites,2],1],hsFC[sites])
cor.test(diffPhospho[[5]][exactMatchesSub[sites,2],1],hsFC[sites])

#is there a difference how well the rep by values fit with qtl by between significant and ns genes?
myPeps <- which(!is.na(diffPhospho[[3]][exactMatchesSub[,2],'fdr']))
pvec <- diffPhospho[[3]][exactMatchesSub[myPeps,2],"fdr"]
plot(rowMeans(phosphoLevel[exactMatchesSub[myPeps,2],grepl("1473",colnames(phosphoLevel))],na.rm=T),phByStrain[exactMatchesSub[myPeps,1],"XBY4716"],pch=19,col=rgb(pvec<0.2,0,0,0.2+(pvec<0.2)/2))
plot(rowMeans(phosphoLevel[exactMatchesSub[myPeps,2],grepl("1473",colnames(prot))])-phByStrain[exactMatchesSub[myPeps,1],"XBY4716"],diffPhospho[[3]][exactMatchesSub[myPeps,2],"l2FC"])


##phospho residuals
#compare foldchanges across comparisons
pdf("graph/alleleReplacement/replacementComparisonPhosphoRes.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[2]][,"fc"],xlab="fc in STE20 replacement",ylab="fc in GPA1 replacement",pch=19,col=rgb(0,0,0,0.3),main="STE20 vs GPA1")
abline(a=0,b=1,col="red")
plot(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[3]][,"fc"],xlab="fc in STE20 replacement",ylab="fc in GPA1+STE20 replacement",pch=19,col=rgb(0,0,0,0.3),main="STE20 vs GPA1+STE20")
abline(a=0,b=1,col="red")
plot(diffPhosphoProt[[2]][,"fc"],diffPhosphoProt[[3]][,"fc"],xlab="fc in GPA1 replacement",ylab="fc in GPA1+STE20 replacement",pch=19,col=rgb(0,0,0,0.3),main="GPA1 vs GPA1+STE20")
abline(a=0,b=1,col="red")
plot(diffPhosphoProt[[2]][,"fc"],diffPhosphoProt[[4]][,"fc"],xlab="fc in BY background",ylab="fc in RM background",pch=19,col=rgb(0,0,0,0.3),main="GPA1 replacement in two STE20 backgrounds")
abline(a=0,b=1,col="red")
plot(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[5]][,"fc"],xlab="fc in BY background",ylab="fc in RM background",pch=19,col=rgb(0,0,0,0.3),main="STE20 replacement in two GPA1 backgrounds")
abline(a=0,b=1,col="red")
dev.off()

cor.test(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[2]][,"fc"],use='pair')
cor.test(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[3]][,"fc"],use='pair')
cor.test(diffPhosphoProt[[2]][,"fc"],diffPhosphoProt[[3]][,"fc"],use='pair')
cor.test(diffPhosphoProt[[2]][,"fc"],diffPhosphoProt[[4]][,"fc"],use='pair')
cor.test(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[5]][,"fc"],use='pair')

cor.test(diffPhosphoProt[[1]][,"fc"],diffPhosphoProt[[4]][,"fc"],use='pair')

#all matched sites from the qtl data
sites <- exactMatchesSub[exactMatchesSub[,1]%in%rownames(phrByStrain),1]
sites <- sites[exactMatchesSub[sites,2]%in%rownames(diffPhosphoProt[[1]])]
corVec <- sapply(diffPhosphoProt,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phosphoRes_all.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[4]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[5]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[4]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[5]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])

#sig by HS as phres
sites <- exactMatchesSub[exactMatchesSub[,1]%in%rownames(phrByStrain),1]
sites <- sites[exactMatchesSub[sites,2]%in%rownames(diffPhosphoProt[[1]])]
sites <- sites[which(apply(pQTL_results$phosphoProtResidualsQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
corVec <- sapply(diffPhosphoProt,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phosphoRes_sigHS.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[4]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[5]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[4]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[5]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])

nickNames[phospho2prot[sites,2]]

#sig by HS as phospho
sites <- exactMatchesSub[exactMatchesSub[,1]%in%rownames(phrByStrain),1]
sites <- sites[exactMatchesSub[sites,2]%in%rownames(diffPhosphoProt[[1]])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
pdf("graph/alleleReplacement/repVsHS_phosphoRes_sigHSByPh.pdf")
par(mfrow=c(2,2))
plot(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS FC",xlab="replacement FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS FC",xlab="replacement FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS FC",xlab="replacement FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()

#sig by 1473 vs 1477 but not HS
sites <- exactMatchesSub[exactMatchesSub[,1]%in%rownames(phrByStrain),1]
sites <- sites[exactMatchesSub[sites,2]%in%rownames(diffPhosphoProt[[1]])]
hsSites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
sites <- sites[which(diffPhosphoProt[[3]][exactMatches[sites,2],"fdr"]<0.2)]
sites <- setdiff(sites,hsSites)
corVec <- sapply(diffPhosphoProt,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],use="pair")})
pdf("graph/alleleReplacement/repVsHS_phosphoRes_sig1477.pdf",height = 10)
par(mfrow=c(3,2))
plot(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("STE20 replacement,\n r=",round(corVec[1],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1 replacement,\n r=",round(corVec[2],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1+STE20 replacement,\n r=",round(corVec[3],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[4]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("GPA1 replacement in RM-STE20 background,\n r=",round(corVec[4],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[5]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS fc",xlab="replacement fc",main=paste0("STE20 replacement in RM-GPA1 background,\n r=",round(corVec[5],digits=2)))
abline(h=0,v=0,col="grey",lty=2)
dev.off()
cor.test(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[4]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])
cor.test(diffPhosphoProt[[5]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites])

#sig by 1597
sites <- exactMatchesSub[exactMatchesSub[,1]%in%rownames(phrByStrain),1]
sites <- sites[exactMatchesSub[sites,2]%in%rownames(diffPhosphoProt[[1]])]
sites <- sites[which(pQTL_results$phosphoProtResidualsQTL$qv[sites,1597]<0.1)]
pdf("graph/alleleReplacement/repVsHS_phosphoRes_sig1597.pdf")
par(mfrow=c(2,2))
plot(diffPhosphoProt[[1]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS FC",xlab="replacement FC",main="STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[2]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS FC",xlab="replacement FC",main="GPA1 replacement")
abline(h=0,v=0,col="grey",lty=2)
plot(diffPhosphoProt[[3]][exactMatchesSub[sites,2],"fc"],hsFCPhr[sites],pch=19,col=rgb(0,0,0,0.3),ylab="HS FC",xlab="replacement FC",main="GPA1+STE20 replacement")
abline(h=0,v=0,col="grey",lty=2)
dev.off()


####compare omics differences between rep strains on multiple layers####
##e and p
#all available genes
isectGenes <- intersect(rownames(prot),rownames(res_1473_1475))
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],diffProt[[1]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="STE20 replacement")
plot(res_1473_1476[isectGenes,"log2FoldChange"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="GPA1 replacement")
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="GPA1+STE20 replacement")
#transcripts that change between 1473 and 1477
isectGenes <- intersect(rownames(prot),rownames(res_1473_1475))
isectGenes <- isectGenes[which(res_1473_1477[isectGenes,"padj"]<0.1)]
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],diffProt[[1]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="STE20 replacement")
plot(res_1473_1476[isectGenes,"log2FoldChange"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="GPA1 replacement")
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="GPA1+STE20 replacement")
#proteins that change between 1473 and 1477
isectGenes <- intersect(rownames(prot),rownames(res_1473_1475))
isectGenes <- isectGenes[which(diffProt[[3]][isectGenes,"fdr"]<0.1)]
par(mfrow=c(2,2))
plot(res_1473_1475[isectGenes,"log2FoldChange"],diffProt[[1]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="STE20 replacement")
plot(res_1473_1476[isectGenes,"log2FoldChange"],diffProt[[2]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="GPA1 replacement")
plot(res_1473_1477[isectGenes,"log2FoldChange"],diffProt[[3]][isectGenes,"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC transcripts",main="GPA1+STE20 replacement")

#compare variation on each layer between strains
strains <- c("1473","1475","1476","1477")
sdByGeneAndStrainE <- sapply(strains,FUN=function(s){
  mat <- expressionLevels
  samples2use <- which(grepl(s,colnames(mat)))[1:3]
  varByGene <- apply(mat[,samples2use],1,sd,na.rm=T)
  return(varByGene)
})
boxplot(sdByGeneAndStrainE-rowMeans(sdByGeneAndStrainE,na.rm=T),outline=F)
cor(sdByGeneAndStrainE)
cor(sdByGeneAndStrainE-rowMeans(sdByGeneAndStrainE))
sdByGeneAndStrainP <- sapply(strains,FUN=function(s){
  mat <- prot
  samples2use <- which(grepl(s,colnames(mat)))[1:3]
  varByGene <- apply(mat[,samples2use],1,sd,na.rm=T)
  return(varByGene)
})
boxplot(sdByGeneAndStrainP-rowMeans(sdByGeneAndStrainP,na.rm=T),outline=F)

##p and ph
#all available phosphopeptides
phWithProt <- rownames(diffPhosphoProt[[1]])
par(mfrow=c(2,2))
plot(diffPhospho[[1]][phWithProt,"l2FC"],diffProt[[1]][phospho2ProtRep[phWithProt,2],"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC phosphopeptides",main="STE20 replacement")
plot(diffPhospho[[2]][phWithProt,"l2FC"],diffProt[[2]][phospho2ProtRep[phWithProt,2],"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC phosphopeptides",main="GPA1 replacement")
plot(diffPhospho[[3]][phWithProt,"l2FC"],diffProt[[3]][phospho2ProtRep[phWithProt,2],"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC phosphopeptides",main="STE20+GPA1 replacement")
cor.test(diffPhospho[[1]][phWithProt,"l2FC"],diffProt[[1]][phospho2ProtRep[phWithProt,2],"l2FC"])
cor.test(diffPhospho[[2]][phWithProt,"l2FC"],diffProt[[2]][phospho2ProtRep[phWithProt,2],"l2FC"])
cor.test(diffPhospho[[3]][phWithProt,"l2FC"],diffProt[[3]][phospho2ProtRep[phWithProt,2],"l2FC"])
#phosphopeptides from proteins that change abundance in the double
phWithProt <- rownames(diffPhosphoProt[[1]])
phWithProt <- phWithProt[which(diffProt[[3]][phospho2ProtRep[phWithProt,2],"p"]<0.1)]
par(mfrow=c(2,2))
plot(diffPhospho[[1]][phWithProt,"l2FC"],diffProt[[1]][phospho2ProtRep[phWithProt,2],"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC phosphopeptides",main="STE20 replacement")
plot(diffPhospho[[2]][phWithProt,"l2FC"],diffProt[[2]][phospho2ProtRep[phWithProt,2],"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC phosphopeptides",main="GPA1 replacement")
plot(diffPhospho[[3]][phWithProt,"l2FC"],diffProt[[3]][phospho2ProtRep[phWithProt,2],"l2FC"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC proteins",xlab="l2FC phosphopeptides",main="STE20+GPA1 replacement")
cor.test(diffPhospho[[1]][phWithProt,"l2FC"],diffProt[[1]][phospho2ProtRep[phWithProt,2],"l2FC"])
cor.test(diffPhospho[[2]][phWithProt,"l2FC"],diffProt[[2]][phospho2ProtRep[phWithProt,2],"l2FC"])
cor.test(diffPhospho[[3]][phWithProt,"l2FC"],diffProt[[3]][phospho2ProtRep[phWithProt,2],"l2FC"])

#by protein as a box plot
phWithProt <- rownames(diffPhosphoProt[[1]])
phWithProt <- phWithProt[which(diffProt[[3]][phospho2ProtRep[phWithProt,2],"p"]<0.1)]
intProts <- unique(phospho2ProtRep[phWithProt,2])
intProts <- intProts[order(diffProt[[3]][intProts,"l2FC"])]
phFCbyProt <- lapply(intProts,FUN=function(p){
  diffPhospho[[3]][phospho2ProtRep[phospho2ProtRep[,2]==p,1],"l2FC"]
}) 
par(mfrow=c(1,1))
boxplot(phFCbyProt)
abline(h=0)

#correlation between peptides and proteins
phWithProt <- rownames(diffPhosphoProt[[1]])
pepProtCor <- sapply(phWithProt,FUN=function(pep){
  cor(prot[phospho2ProtRep[pep,2],],phosphoLevel[pep,colnames(prot)],use='pair')
})
pdf("graph/alleleReplacement/corPepProt.pdf",height = 4)
hist(pepProtCor,xlab="correlation",main="correlation between phosphopeptides and host proteins",xlim=c(-1,1))
dev.off()
boxplot(pepProtCor~as.factor(diffPhospho[[3]][phWithProt,"p"]<0.1))
boxplot(pepProtCor~as.factor(diffPhosphoProt[[3]][phWithProt,"p"]<0.1))

plot(abs(diffProt[[3]][phospho2ProtRep[phWithProt,2],"l2FC"]),abs(diffPhospho[[3]][phWithProt,"l2FC"]),pch=19,col=rgb(0,0,0,0.3))
abline(a=0,b=1,col="red")
hist(abs(diffProt[[3]][phospho2ProtRep[phWithProt,2],"l2FC"])-abs(diffPhospho[[3]][phWithProt,"l2FC"]))
#do rna changes also agree with phospho changes?
phPep <- rownames(diffPhospho[[1]])
par(mfrow=c(2,2))
plot(diffPhospho[[1]][phPep,"l2FC"],res_1473_1475[phospho2ProtRep[phPep,2],"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC transcripts",xlab="l2FC phosphopeptides",main="STE20 replacement")
plot(diffPhospho[[2]][phPep,"l2FC"],res_1473_1476[phospho2ProtRep[phPep,2],"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC transcripts",xlab="l2FC phosphopeptides",main="GPA1 replacement")
plot(diffPhospho[[3]][phPep,"l2FC"],res_1473_1477[phospho2ProtRep[phPep,2],"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC transcripts",xlab="l2FC phosphopeptides",main="STE20+GPA1 replacement")
cor.test(diffPhospho[[1]][phPep,"l2FC"],res_1473_1475[phospho2ProtRep[phPep,2],"log2FoldChange"])
cor.test(diffPhospho[[2]][phPep,"l2FC"],res_1473_1476[phospho2ProtRep[phPep,2],"log2FoldChange"])
cor.test(diffPhospho[[3]][phPep,"l2FC"],res_1473_1477[phospho2ProtRep[phPep,2],"log2FoldChange"])

#phosphopeptides from transcripts that change abundance in the double
phPep <- rownames(diffPhospho[[1]])
phPep <- phPep[which(res_1473_1477[phospho2ProtRep[phPep,2],"padj"]<0.2)]
par(mfrow=c(2,2))
plot(diffPhospho[[1]][phPep,"l2FC"],res_1473_1475[phospho2ProtRep[phPep,2],"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC transcripts",xlab="l2FC phosphopeptides",main="STE20 replacement")
plot(diffPhospho[[2]][phPep,"l2FC"],res_1473_1476[phospho2ProtRep[phPep,2],"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC transcripts",xlab="l2FC phosphopeptides",main="GPA1 replacement")
plot(diffPhospho[[3]][phPep,"l2FC"],res_1473_1477[phospho2ProtRep[phPep,2],"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3),ylab="l2FC transcripts",xlab="l2FC phosphopeptides",main="STE20+GPA1 replacement")

phSub <- phosphoLevel
colnames(phSub) <- gsub(pattern = "_r",replacement = "_",colnames(phSub))
colnames(phSub) <- gsub(pattern = "X",replacement = "s",colnames(phSub))

phPep <- rownames(diffPhospho[[1]])
pepRnaCor <- sapply(phPep,FUN=function(pep){
  cor(expressionLevels[phospho2ProtRep[pep,2],],phSub[pep,colnames(expressionLevels)],use='pair')
})
hist(pepRnaCor)
pdf("graph/alleleReplacement/protCorVsRnaCor.pdf")
plot(pepRnaCor[names(pepProtCor)],pepProtCor,xlab="correlation with RNA",ylab="correlation with protein",xlim=c(-1,1),ylim=c(-1,1),pch=19,col=rgb(0,0,0,0.5))
abline(h=0,v=0)
abline(a=0,b=1,lty=2,lwd=2,col="red")
dev.off()

plot(pepRnaCor[names(pepProtCor)],pepProtCor,xlab="correlation with RNA",ylab="correlation with protein",xlim=c(-1,1),ylim=c(-1,1),pch=19,col=rgb(phospho2ProtRep[names(pepProtCor),2]%in%targetsByHS[[10]]$ph,0,0,0.3))
abline(h=0,v=0)



boxplot(pepProtCor,pepRnaCor[names(pepProtCor)],pepRnaCor[setdiff(names(pepRnaCor),names(pepProtCor))])

pdf("graph/alleleReplacement/corPepRNA.pdf",height = 4)
hist(pepRnaCor,xlab="correlation",main="correlation between phosphopeptides and transcripts of host proteins",xlim=c(-1,1))
dev.off()

usePh <- intersect(rownames(phosphoProtResiduals),exactMatchesSub[,2])
plot(rowMeans(phosphoLevel[exactMatchesSub[,2],],na.rm=T),rowMeans(phByStrain[exactMatchesSub[,1],],na.rm=T),pch=19,col=rgb(phospho2ProtRep[names(pepProtCor),2]%in%targetsByHS[[10]]$ph,0,0,0.3))
abline(a=0,b=1)
usePh <- exactMatchesSub[exactMatchesSub[,1]%in%rownames(phosphoProtResidualsQTL)&exactMatchesSub[,2]%in%rownames(phosphoProtResiduals),]
phProtRatioMeans <- rowMeans(prot[phospho2ProtRep[usePh[,2],2],],na.rm=T)-rowMeans(phosphoLevel[usePh[,2],],na.rm=T)
names(phProtRatioMeans) <- usePh[,2]
phProtRatioMeansQTL <- rowMeans(pByStrain[phospho2prot[usePh[,1],2],],na.rm=T)-rowMeans(phByStrain[usePh[,1],],na.rm=T)
names(phProtRatioMeansQTL) <- usePh[,1]

plot(phProtRatioMeansQTL,phProtRatioMeans,pch=19,col=rgb(0,0,0,0.3))
abline(a=0,b=1)

pdf("graph/alleleReplacement/comparisonsOfMeans.pdf")
par(mfrow=c(2,2))
plot(rowMeans(prot[unique(phospho2ProtRep[usePh[,2],2]),],na.rm=T),rowMeans(pByStrain[unique(phospho2prot[usePh[,1],2]),],na.rm=T),pch=19,col=rgb(0,0,0,0.3),xlab="validation data",ylab="QTL data",main="average protein abundances")
plot(rowMeans(phosphoLevel[usePh[,2],],na.rm=T),rowMeans(phByStrain[usePh[,1],],na.rm=T),pch=19,col=rgb(0,0,0,0.3),xlab="validation data",ylab="QTL data",main="average phosphopeptide abundances")
plot(phProtRatioMeans,phProtRatioMeansQTL,pch=19,col=rgb(0,0,0,0.3),xlab="validation data",ylab="QTL data",main="average ratios of phosphopeptides\n and proteins")
dev.off()

cor(rowMeans(prot[unique(phospho2ProtRep[usePh[,2],2]),],na.rm=T),rowMeans(pByStrain[unique(phospho2prot[usePh[,1],2]),],na.rm=T))
cor(rowMeans(phosphoLevel[usePh[,2],],na.rm=T),rowMeans(phByStrain[usePh[,1],],na.rm=T))
cor(phProtRatioMeans,phProtRatioMeansQTL)

# ####check quality of phospho results####
# highPepProt <- names(which(table(phospho2ProtRep[,2])>10))
# phosphoProts <- intersect(rownames(diffProt[[3]]),highPepProt)
# plotMat <- t(sapply(phosphoProts,FUN=function(p){
#   c(diffProt[[3]][p,"l2FC"],mean(diffPhospho[[3]][phospho2ProtRep[,2]==p,"l2FC"]))
# }))
# plot(plotMat)
# cor(plotMat,use="pair")
# 
# epIsect <- intersect(rownames(diffProt[[3]]),rownames(res_1473_1477))
# plot(res_1473_1477[epIsect,"log2FoldChange"],diffProt[[3]][epIsect,"l2FC"])
# 
# sharedSamples <- setdiff(intersect(colnames(phosphoLevel),colnames(prot)),"X1473_r3")
# protPhosphoCor <- sapply(phospho2ProtRep[phospho2ProtRep[,2]%in%rownames(prot),1],FUN=function(pep){
#   cor(phosphoLevel[pep,sharedSamples],prot[phospho2ProtRep[pep,2],sharedSamples],use="pair")
# })
# strainLabelsSub <- gsub("_.*","",sharedSamples)
# names(strainLabelsSub) <- sharedSamples
# phCorrectedForStrain <- sapply(sharedSamples,FUN=function(s){
#   phosphoLevel[,s]-rowMeans(phosphoLevel[,sharedSamples[strainLabelsSub==strainLabelsSub[s]]],na.rm=T)
# })
# pCorrectedForStrain <- sapply(sharedSamples,FUN=function(s){
#   prot[,s]-rowMeans(prot[,sharedSamples[strainLabelsSub==strainLabelsSub[s]]],na.rm=T)
# })
# protPhosphoCor2 <- sapply(phospho2ProtRep[phospho2ProtRep[,2]%in%rownames(prot),1],FUN=function(pep){
#   cor(phCorrectedForStrain[pep,sharedSamples],pCorrectedForStrain[phospho2ProtRep[pep,2],sharedSamples],use="pair")
# })
# 
# protPhosphoCor3 <- sapply(phospho2ProtRep[phospho2ProtRep[,2]%in%rownames(prot),1],FUN=function(pep){
#   if(any(is.na(phCorrectedForStrain[pep,]))|any(is.na(pCorrectedForStrain[phospho2ProtRep[pep,2],]))){return(NA)}
#   cor(phCorrectedForStrain[pep,sharedSamples],pCorrectedForStrain[phospho2ProtRep[pep,2],sharedSamples],use="pair")
# })
# 
# usePeps <- pep2Prot[pep2Prot[,1]%in%rownames(prot),2]
# pepProtCor <- sapply(usePeps,FUN=function(pep){
#   cor(as.numeric(prot[pep2Prot[pep2Prot[,2]==pep,1],]),as.numeric(pepLevels[pep,]),use="pair")
# })
# 
# singlePepProts <- names(which(table(pep2Prot[,1])==1))
# nSnonPh <- sapply(rownames(pepLevels),nchar)-sapply(gsub("S","",rownames(pepLevels)),nchar)
# nYnonPh <- sapply(rownames(pepLevels),nchar)-sapply(gsub("Y","",rownames(pepLevels)),nchar)
# 
# 
# pep2 <- pepLevels
# colnames(pep2) <- gsub("X","s",gsub("r","",colnames(pep2)))
# sharedSamples <- intersect(colnames(pep2),colnames(expressionLevels))
# pepRnaCor <- sapply(pep2Prot[,2],FUN=function(pep){
#   if(!pep2Prot[pep2Prot[,2]==pep,1]%in%rownames(expressionLevels)){return(NA)}
#   cor(expressionLevels[pep2Prot[pep2Prot[,2]==pep,1],sharedSamples],as.numeric(pep2[pep,sharedSamples]),use="pair")
# })
# 
# rawPeps <- gsub("(","X",phospho2ProtRep[,1],fixed=T)
# rawPeps <- gsub(")","X",rawPeps,fixed=T)
# rawPeps <- gsub("UniMod:","",rawPeps,fixed=T)
# rawPeps <- gsub("X[0-9]*X","",rawPeps,fixed=F)
# 
# 
# rawPeps2 <- gsub("(","X",rownames(pepLevels),fixed=T)
# rawPeps2 <- gsub(")","X",rawPeps2,fixed=T)
# rawPeps2 <- gsub("UniMod:","",rawPeps2,fixed=T)
# rawPeps2 <- gsub("X[0-9]*X","",rawPeps2,fixed=F)
# 
# boxplot(pepRnaCor~as.factor(rawPeps2%in%rawPeps))
# 
# pepCor <- sapply(which(rawPeps%in%rawPeps2),FUN=function(i){
#   cor(phosphoLevel[i,],as.numeric(pepLevels[which(rawPeps2==rawPeps[i])[1],]),use="pair")
# })
# names(pepCor) <- rownames(phosphoLevel)[which(rawPeps%in%rawPeps2)]
# 
# subProts <- names(which(table(phospho2ProtRep[phospho2ProtRep[,2]%in%rownames(prot),2])>20))
# 
# kInTheMiddle <- sapply(rawPeps,FUN=function(pep){
#   noEnd <- substr(pep,1,nchar(pep)-1)
#   sum(strsplit(x = noEnd,split = "")[[1]]=="K")
# })
# rInTheMiddle <- sapply(rawPeps,FUN=function(pep){
#   noEnd <- substr(pep,1,nchar(pep)-1)
#   sum(strsplit(x = noEnd,split = "")[[1]]=="R")
# })
# 
# corWithSte20 <- apply(phosphoLevel[,-3],1,cor,y=expressionLevels["YHL007C",],use="pair")
# corWithGpa1 <- apply(phosphoLevel[,-3],1,cor,y=expressionLevels["YHR005C",],use="pair")
# 
# myRawPeps <- rawPeps[rownames(diffPhosphoProt[[3]])]
# corByRawPep <- sapply(unique(rawPeps),FUN=function(pep){
#   
# })

#Bodenmiller data
kinaseInteractionsFull <- read.table(file = "data/Bodenmiller2010.tsv",header = T,sep = "\t",quote = "",as.is = T)
kinaseInteractions <- kinaseInteractionsFull[grepl("[166.9984]",kinaseInteractionsFull[,3],fixed = T),]
uniRowsKinase <- kinaseInteractionsFull[,c(2,4)]
uniRowsKinase <- unique(uniRowsKinase)
biogridDataFull <- read.table("/data/user/jgrossb1/Downloads/BIOGRID-ALL-3.4.158.tab2.txt",sep="\t",header=T,quote="",comment.char = "",as.is=T)
biogridDataBK <- biogridDataFull[biogridDataFull$Author=="Breitkreutz A (2010)"&biogridDataFull$Experimental.System.Type=="physical",]
#biogridDataBK <- unname(biogridDataBK)
#uniRowsKinase <- unname(uniRowsKinase)
uniRowsKinase <- rbind(cbind(Regulator=uniRowsKinase[,1],Target=uniRowsKinase[,2],Author="Bodenmiller"),cbind(Regulator=biogridDataBK[,"Systematic.Name.Interactor.A"],Target=biogridDataBK[,"Systematic.Name.Interactor.B"],Author="Breitkreuz"))
uniRowsKinase <- uniRowsKinase[uniRowsKinase[,1]!=uniRowsKinase[,2],]
uniRowsKinase <- unique(uniRowsKinase)

ste20SpecTr <- c("crh1", "vhs3", "kar5", "cdc20", "rgc1")
revNames[ste20SpecTr]
res_1473_1475[revNames[ste20SpecTr],]
res_1473_1476[revNames[ste20SpecTr],]
res_1473_1477[revNames[ste20SpecTr],]

ste20SpePhProts <- tolower(strsplit(x = "Rts1, Aut1, Nup2, Rgc1, Rck2, Rsc9, Ssk1",split = ", ")[[1]])
diffPhospho[[1]][phospho2ProtRep[,2]%in%revNames[ste20SpePhProts],]
diffPhospho[[2]][phospho2ProtRep[,2]%in%revNames[ste20SpePhProts],]

###see if the variation in trait levels is even plausible/could be missed in the qtl data###
isectProts <- intersect(rownames(diffProt[[3]]),rownames(pByStrain))
par(mfrow=c(1,1))
plot(abs(diffProt[[3]][isectProts,"l2FC"]),abs(pByStrain[isectProts,"XBY4716"]-pByStrain[isectProts,"XRM11-1a"]))
abline(a=0,b=1)
hist(abs(diffProt[[3]][isectProts,"l2FC"])-abs(pByStrain[isectProts,"XBY4716"]-pByStrain[isectProts,"XRM11-1a"]),breaks=50)

par(mfrow=c(1,1))
plot(abs(diffProt[[3]][isectProts,"l2FC"]),abs(apply(pByStrain[isectProts,],1,sd,na.rm=T)))
hist(abs(diffProt[[3]][isectProts,"l2FC"])-abs(apply(pByStrain[isectProts,],1,sd,na.rm=T)),breaks=50)

plot(rowMeans(prot[isectProts,1:4],na.rm=T),pByStrain[isectProts,"XBY4716"])
experimentBYFC <- lm(rowMeans(prot[isectProts,1:4],na.rm=T)~pByStrain[isectProts,"XBY4716"])$residuals
plot(experimentBYFC,diffProt[[3]][names(experimentBYFC),"l2FC"])

pepsWithProts <- phospho2ProtRep[phospho2ProtRep[,2]%in%rownames(prot),1]
plot(diffPhospho[[3]][pepsWithProts,"l2FC"],diffProt[[3]][phospho2ProtRep[pepsWithProts,2],"l2FC"])

protNoNA <- t(apply(prot,1,FUN=function(x){
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}))
repProtPCA <- princomp(protNoNA)

sigProts <- rownames(prot)[which(diffProt[[3]][,"fdr"]<0.2)]
plot(diffProt[[3]][sigProts,"l2FC"],res_1473_1477[sigProts,"log2FoldChange"],pch=19,col=rgb(0,0,0,0.3))
sigProtsOnly <- setdiff(sigProts,rownames(res_1473_1477)[which(res_1473_1477$padj<0.5)])
plot(diffProt[[3]][sigProts,"l2FC"],res_1473_1477[sigProts,"log2FoldChange"],pch=19,col=rgb(sigProts%in%sigProtsOnly,0,0,0.3))
boxplot(rnaProtCor[sigProts]~as.factor(sigProts%in%sigProtsOnly))

singlePepFCs <- rowMeans(singlePepLevel[,1:4],na.rm=T)-rowMeans(singlePepLevel[,13:16],na.rm=T)
names(singlePepFCs) <- pepProtMat[,2]
boxplot(as.numeric(table(pepProtMat)[sigProts])~as.factor(sigProts%in%sigProtsOnly))

pepProtMat2 <- pepProtMat
rownames(pepProtMat2) <- pepProtMat2[,2]
pepsOfInterest <- pepProtMat[pepProtMat[,1]%in%sigProts,2]
plot(singlePepFCs[pepsOfInterest],diffProt[[3]][pepProtMat2[pepsOfInterest,1],"l2FC"],pch=19,col=rgb(pepProtMat2[pepsOfInterest,1]%in%sigProtsOnly,0,0,0.3))
abline(a=0,b=1,col="green")
meanPepFC <- sapply(unique(pepProtMat[,1]),FUN=function(prot){mean(singlePepFCs[pepProtMat[,1]==prot],na.rm=T)})
plot(diffProt[[3]][,"l2FC"],meanPepFC[rownames(prot)],col=rgb(0,0,0,0.3),pch=19)
plot(diffProt[[3]][sigProts,"l2FC"],meanPepFC[sigProts],col=rgb(table(pepProtMat[,1])[sigProts]>1,0,0,0.3),pch=19)
abline(a=0,b=1,col="green")

protsWithLevels <- intersect(rownames(prot),pepProtMat[,1])
protsWith10Peps <- intersect(protsWithLevels,names(which(table(pepProtMat[,1])>9)))
plot(meanPepFC[protsWithLevels],res_1473_1477[protsWithLevels,"log2FoldChange"],col=rgb(table(pepProtMat[,1])[protsWithLevels]>3,0,0,0.5),pch=19)
plot(meanPepFC[protsWith10Peps],res_1473_1477[protsWith10Peps,"log2FoldChange"],col=rgb(0,0,0,0.5),pch=19)
plot(meanPepFC[protsWith10Peps],diffProt[[3]][protsWith10Peps,"l2FC"])
abline(a=0,b=1)
nPep <- read.table("ste20Validation/proteinData/protein/mapDIA/protein_level.txt",sep="\t",header=T)
nPep <- nPep[,c(1,18,19)]
rownames(nPep) <- nPep[,1]

rawPepsVal <- read.table("ste20Validation/proteinData/phosphoNoCor/mapDIA/peptide_level.txt",sep="\t",header=T,as.is=T)
nPepsPh <- rawPepsVal[,c(2,19)]
rownames(nPepsPh) <- nPepsPh[,1]


myProts <- rownames(nPep)[nPep[,3]>4]
myPeps <- rownames(nPepsPh)[nPepsPh[,2]>14&phospho2ProtRep[nPepsPh[,1],2]%in%myProts]
myPeps <- intersect(myPeps,rownames(diffPhosphoProtRatio[[3]]))
protsToPlot <- intersect(intersect(rownames(nPep)[nPep[,3]>4],rownames(diffProt[[3]])),rownames(sv[[3]]))
# 
# sv <- diffProt
# load("data/diffProtReplacementStrainsNoCor.RData")
# par(mfrow=c(1,3))
# plot(diffProt[[3]][protsToPlot,"l2FC"],res_1473_1477[protsToPlot,"log2FoldChange"],xlim=c(-2,2))
# abline(a=0,b=1)
# plot(sv[[3]][protsToPlot,"l2FC"],res_1473_1477[protsToPlot,"log2FoldChange"],xlim=c(-2,2))
# abline(a=0,b=1)
# plot(meanPepFC[protsToPlot],res_1473_1477[protsToPlot,"log2FoldChange"],xlim=c(-2,2))
# abline(a=0,b=1)
# 
# 
# plot(meanPepFC[protsToPlot],diffProt[[3]][protsToPlot,"l2FC"])


###look at the levels of gpa1 and ste20 themselves
#gpa1
gpa1 <- "YHR005C"
boxplot(expressionLevels[gpa1,]~as.factor(gsub("_.*","",colnames(expressionLevels))))
boxplot(prot[gpa1,]~as.factor(gsub("_.*","",colnames(prot))))
hist(eByStrain[gpa1,]-eByStrain[gpa1,"XBY4716"])
abline(v=res_1473_1477[gpa1,"log2FoldChange"],col="red")#in range
hist(pByStrain[gpa1,]-pByStrain[gpa1,"XBY4716"])
abline(v=diffProt[[3]][gpa1,"l2FC"],col="red") #out of range but noisy
#ste20
ste20 <- "YHL007C"
boxplot(expressionLevels[ste20,]~as.factor(gsub("_.*","",colnames(expressionLevels))))
boxplot(prot[ste20,]~as.factor(gsub("_.*","",colnames(prot))))
hist(eByStrain[ste20,]-eByStrain[ste20,"XBY4716"])
abline(v=res_1473_1477[ste20,"log2FoldChange"],col="red")#in range
hist(pByStrain[ste20,]-pByStrain[ste20,"XBY4716"],xlim=c(-1,1))
abline(v=diffProt[[3]][ste20,"l2FC"],col="red") #out of range but noisy

hsFcE[c(gpa1,ste20)]
res_1473_1477[c(gpa1,ste20),"log2FoldChange"]
hsFcP[c(gpa1,ste20)]
diffProt[[3]][c(gpa1,ste20),"l2FC"]

###look at dataset specific differences
isectGenes <- intersect(intersect(rownames(eByStrain),rownames(pByStrain)),rownames(prot))

plot(rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T)-eByStrain[isectGenes,"XBY4716"],rowMeans(prot[isectGenes,grepl("1473",colnames(prot))],na.rm=T)-pByStrain[isectGenes,"XBY4716"])
cor(rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T)-eByStrain[isectGenes,"XBY4716"],rowMeans(prot[isectGenes,grepl("1473",colnames(prot))],na.rm=T)-pByStrain[isectGenes,"XBY4716"])
plot(rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T)-eByStrain[isectGenes,"XBY4716"],hsFcE228[isectGenes])
plot(rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T)-eByStrain[isectGenes,"XBY4716"],hsFcE_ira2[isectGenes])

isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
expDiff <- rowMeans(expressionLevels[isectGenes,grepl("1473",colnames(expressionLevels))],na.rm=T)-eByStrain[isectGenes,"XBY4716"]
corWithDiff <- apply(fcByLocE[isectGenes,1:3593],2,cor,y=expDiff,use="pair")
tail(order(corWithDiff))
plot(fcByLocE[isectGenes,1056],expDiff)
plot(fcByLocE[isectGenes,712],expDiff)
###look at other loci with effects on mating###
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
plot(hsFC[sites],hsFC2565[sites])
plot(hsFC_2565_0[sites],hsFC_2565_1[sites])
abline(a=0,b=1)
plot(diffPhospho[[3]][exactMatchesSub[sites,2],"l2FC"],hsFC_2565_0[sites])

isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- intersect(isectGenes,targetsByHS[[10]]$e)
pdf("graph/hsTargets_three_loci.pdf")
invisible(sapply(isectGenes,FUN=function(g){
  boxplot(eByStrain[g,]~genotype[,396]*genotype[,1597]*genotype[,2565],main=g)
}))
dev.off()
pdf("graph/hsTargets_four_loci.pdf")
invisible(sapply(isectGenes,FUN=function(g){
  boxplot(eByStrain[g,]~genotype[,396]*genotype[,694]*genotype[,1597]*genotype[,2565],main=g)
}))
dev.off()
strains0 <- rownames(genotype)[which(genotype[,694]==0&genotype[,396]==0&genotype[,2565]==0)]
repFCE <- rowMeans(eByStrain[,genotype[strains0,1597]==0],na.rm=T)-rowMeans(eByStrain[,genotype[strains0,1597]==1],na.rm=T)
repFC <- rowMeans(phByStrain[,genotype[strains0,1597]==0],na.rm=T)-rowMeans(phByStrain[,genotype[strains0,1597]==1],na.rm=T)
plot(repFCE[isectGenes],res_1473_1477[isectGenes,"log2FoldChange"])
plot(repFC[sites],diffPhospho[[3]][exactMatchesSub[sites,2],"l2FC"])

isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- intersect(isectGenes,targetsByHS[[10]]$e)
corByFCE <- apply(fcByLocE[,1:3593],2,cor,y=fcByLocE[,1597],use="pair")
m <- genotype[,1597]
ct <- 0
suppressionOf1597 <- t(apply(genotype[,1:3593],2,FUN=function(x){
  ct<<-ct+1
  if(fisher.test(m,x)$p.value<0.01){return(c(NA,NA,NA))}
  fc0 <- rowMeans(eByStrain[,m==0&x==0],na.rm=T)-rowMeans(eByStrain[,m==1&x==0],na.rm=T)
  fc1 <- rowMeans(eByStrain[,m==0&x==1],na.rm=T)-rowMeans(eByStrain[,m==1&x==1],na.rm=T)
  p <- wilcox.test(abs(fc0[isectGenes]),abs(fc1[isectGenes]),paired = T)$p.value
  s <- lm(fc0[isectGenes]~fc1[isectGenes])$coefficients[2]
  s2 <- lm(fc1[isectGenes]~fc0[isectGenes])$coefficients[2]
  return(c(p,s,s2))
}))
corWith1597 <- apply(genotype[,1:3593],2,FUN=function(x){
  cor.test(x,m)$p.value
})


par(mfrow=c(2,2))
hist(rowMeans(expressionLevels[,grepl("s1473",colnames(expressionLevels))]),breaks=100)
abline(v=mean(expressionLevels["YPL187W",grepl("s1473",colnames(expressionLevels))]),col="red")

hist(eByStrain[,"XBY4716"],breaks=100)
abline(v=eByStrain["YPL187W","XBY4716"],col="red")

hist(eByStrain[,"XRM11-1a"],breaks=100)
abline(v=eByStrain["YDR461W","XRM11-1a"],col="red")

par(mfrow=c(1,1))
beeswarm(list(eByStrain["YPL187W",],eByStrain["YDR461W",]),pwcol=rep(rgb(genotype[,1597]==0,0,0,0.5),2),pch=19)

par(mfrow=c(2,1))
barplot(-log10(qv["YPL187W",]))
barplot(-log10(qv["YDR461W",]))

localHsGenes <- unique(gtf[gtf[,1]=="chrVIII"&gtf[,4]<=140000&gtf[,5]>=70000,9])
table(localHsGenes%in%rownames(res_1473_1477))
as.matrix(res_1473_1477[localHsGenes,])
plot(res_1473_1477[localHsGenes,"log2FoldChange"],hsFcE[localHsGenes])

hist(log2(as.numeric(rawCountsQTL["YPL187W",])),breaks=30)
hist(log2(as.numeric(rawCountsQTL["YDR461W",])),breaks=30)

plot(log2(as.numeric(rawCountsQTL["YPL187W",])),log2(as.numeric(rawCountsQTL["YDR461W",])))

plot(eByStrain["YHR005C",],eByStrain["YPL187W",],pch=19,col=rgb(genotype[,396]==0,0,0,0.5))
plot(eByStrain["YHR005C",],eByStrain["YDR461W",],pch=19,col=rgb(genotype[,396]==0,0,0,0.5))

corWithEByAllele <- sapply(c(0,1),FUN=function(a){
  apply(eByStrain[,genotype[,1597]==a],1,cor,y=p)
})

###are the effects in the rep strains within the physiological range of the cross?###
eRange <- eByStrain-eByStrain[,"XBY4716"]
intGenes <- intersect(rownames(res_1473_1477)[which(res_1473_1477$pvalue<0.001)],rownames(eByStrain))
boxplot(t(eRange[intGenes,]),horizontal=T)
abline(v=0)
points(y = 1:length(intGenes),x=(-1)*res_1473_1477[intGenes,"log2FoldChange"],col="red",pch=19)

plot(apply(eRange[intGenes,],1,median),(-1)*res_1473_1477[intGenes,"log2FoldChange"])
abline(h=0,v=0,a=0,b=1)

pdf("graph/alleleReplacement/expressionRanges.pdf")
invisible(sapply(intGenes,FUN=function(g){
  hist(eRange[g,],xlim=range(c(eRange[g,],(-1)*res_1473_1477[g,"log2FoldChange"]))*1.2,main=paste0(g,", ",nickNames[g]))
  abline(v=res_1473_1477[g,"log2FoldChange"]*(-1),col="red")
}))
dev.off()

pRange <- pByStrain-pByStrain[,"XBY4716"]
intGenes <- intersect(rownames(diffProt[[3]])[which(diffProt[[3]][,"p"]<0.005)],rownames(pByStrain))
intGenes <- intGenes[!is.na(pByStrain[intGenes,"XBY4716"])]
pdf("graph/alleleReplacement/expressionRangesP.pdf")
invisible(sapply(intGenes,FUN=function(g){
  hist(pRange[g,],xlim=range(c(pRange[g,],(-1)*diffProt[[3]][g,"l2FC"]),na.rm = T)*1.2,main=paste0(g,", ",nickNames[g]))
  abline(v=(-1)*diffProt[[3]][g,"l2FC"],col="red")
}))
dev.off()

strainCenteredPhospho <- lapply(unique(strainLabels),FUN=function(s){
  phosphoLevel[,strainLabels==s]-rowMeans(phosphoLevel[,strainLabels==s],na.rm=T)
})
strainCenteredPhospho <- do.call("cbind",strainCenteredPhospho)
strainCenteredPhospho[is.na(strainCenteredPhospho)] <- 0
pca_cen_ph <- princomp(strainCenteredPhospho)
cor_cen_ph <- apply(pca_cen_ph$loadings,2,FUN=function(pc){
  apply(strainCenteredPhospho,1,cor,y=pc,use="pair")
})

strainLabelsE <- sapply(colnames(expressionLevels),FUN=function(s){
  gsub("_.*","",s)
})
strainCenteredRNA <- lapply(unique(strainLabelsE),FUN=function(s){
  expressionLevels[,strainLabelsE==s]-rowMeans(expressionLevels[,strainLabelsE==s],na.rm=T)
})
strainCenteredRNA <- do.call("cbind",strainCenteredRNA)
strainCenteredRNA[is.na(strainCenteredRNA)] <- 0
pca_cen_e <- princomp(strainCenteredRNA)
cor_cen_e <- apply(pca_cen_e$loadings,2,FUN=function(pc){
  apply(strainCenteredRNA,1,cor,y=pc,use="pair")
})

###plot whole mutant effects###
#e
isectGenes <- isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- intersect(isectGenes,targetsByHS[[10]]$e)
validationMatE <- data.frame(s=res_1473_1475[isectGenes,"log2FoldChange"],g=res_1473_1476[isectGenes,"log2FoldChange"],sg=res_1473_1477[isectGenes,"log2FoldChange"])
rownames(validationMatE) <- isectGenes
plotPointsE <- t(apply(validationMatE,2,FUN=function(x){
  fcLm <- summary(lm(x~hsFcE[isectGenes]))$coefficients
  slope <- fcLm[2,1]
  slopeError <- fcLm[2,2]
  r <- cor(x,hsFcE[isectGenes],use="pair")
  p <- fcLm[2,4]
  return(c(slope,r,slopeError,p))
}))
colnames(plotPointsE) <- c("slope","r","slopeError","p")
plotPointsE <- data.frame(plotPointsE,type="e",stringsAsFactors = F)

#ph
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
validationMatPh <- data.frame(s=diffPhospho[[1]][exactMatches[sites,2],"l2FC"],g=diffPhospho[[2]][exactMatches[sites,2],"l2FC"],sg=diffPhospho[[3]][exactMatches[sites,2],"l2FC"])
plotPointsPh <- t(apply(validationMatPh,2,FUN=function(x){
  fcLm <- summary(lm(x~hsFC[sites]))$coefficients
  slope <- fcLm[2,1]
  slopeError <- fcLm[2,2]
  r <- cor(x,hsFC[sites],use="pair")
  p <- fcLm[2,4]
  return(c(slope,r,slopeError,p))
}))
colnames(plotPointsPh) <- c("slope","r","slopeError","p")
plotPointsPh <- data.frame(plotPointsPh,type="phospho",stringsAsFactors = F)

plotPointsCombined <- rbind(plotPointsE,plotPointsPh)
labelVec <- rep(c("STE20","GPA1","STE20+GPA1"),2)
pdf("graph/alleleReplacement/summary_e_ph.pdf")
par(cex=1.4,cex.lab=1.4)
plot(plotPointsCombined[,1:2],pch=19,col=colTrait[plotPointsCombined[,"type"]],xlim=c(0.3,1.6),ylim=c(0.2,0.8),cex=2.5)
arrows(x0 = plotPointsCombined[,1]-plotPointsCombined[,3]/2,y0 = plotPointsCombined[,2],x1 = plotPointsCombined[,1]+plotPointsCombined[,3]/2,y1 = plotPointsCombined[,2],code = 3,angle = 90,length = .05,lwd=2)
points(plotPointsCombined[,1:2],pch=19,col=colTrait[plotPointsCombined[,"type"]],cex=2.5)
points(plotPointsCombined[,1:2],cex=2.5,lwd=2)
text(x = plotPointsCombined[,1]+plotPointsCombined[,3]/2,y=plotPointsCombined[,2],labels=labelVec,pos=4)
legend(x = "topleft",legend = c("transcriptome","phosphoproteome"),pch=19,col = colTrait[c("e","phospho")])
dev.off()

####complete lower panel for pathway figure####
pdf("graph/alleleReplacement/pathway_lower_part.pdf",width=12,height=4)
par(cex=1,cex.axis=1.3,cex.lab=1.7,cex.main=1.7,mfrow=c(1,3),las=1,mar=c(5,5,4,1)+0.1)
#b
isectGenes <- intersect(rownames(expressionLevels),names(hsFcE))
isectGenes <- intersect(isectGenes,targetsByHS[[10]]$e)
plot(y=res_1473_1477[isectGenes,"log2FoldChange"],x=hsFcE[isectGenes],pch=19,col=colTrait["e"],xlab="hotspot log2(FC)",ylab="STE20+GPA1 allele replacement log2(FC)",main="transcriptome: r=0.54, m=1.07",cex=1.5)
abline(h=0,v=0,col="grey",lty=2)
abline(a=0,b=plotPointsCombined[3,1],col=colTrait["e"])
#c
sites <- names(hsFC)[!is.na(exactMatches[names(hsFC),2])]
sites <- sites[which(apply(pQTL_results$phosphoLevelQTL$qv[sites,binPerMarker%in%c(119,120)],1,min,na.rm=T)<0.1)]
corVec <- sapply(diffPhospho,FUN=function(mat){cor(mat[exactMatchesSub[sites,2],1],hsFC[sites],use="pair")})
round(corVec[1:2],digits=2)
nameVec <- c("STE20","GPA1","STE20 + GPA1")
invisible(sapply(3,FUN=function(i){
  mat <- diffPhospho[[i]]
  plot(y=mat[exactMatchesSub[sites,2],1],x=hsFC[sites],pch=19,col=colTrait["phospho"],xlab="hotspot log2(FC)",ylab="STE20+GPA1 allele replacement log2(FC)",main="phosphoproteome: r=0.45, m=0.7",xlim=range(mat[exactMatchesSub[sites,2],1],na.rm=T)*1.3,ylim=range(hsFC[sites],na.rm=T)*1.3,cex=1.5)
  abline(h=0,v=0,col="grey",lty=2)
  abline(a=0,b=plotPointsCombined[6,1],col=colTrait["phospho"])
  #text(x=mat[exactMatchesSub[sites,2],1],y=hsFC[sites],labels = nickNames[phospho2prot[sites,2]],adj=c(-1,-1)*0.5,cex = 0.7)
}))
#d
plot(plotPointsCombined[,1:2],pch=19,col=colTrait[plotPointsCombined[,"type"]],xlim=c(0.3,1.6),ylim=c(0.2,0.8),cex=2.5,main="allele replacement effects",ylab="Pearson's r")
arrows(x0 = plotPointsCombined[,1]-plotPointsCombined[,3]/2,y0 = plotPointsCombined[,2],x1 = plotPointsCombined[,1]+plotPointsCombined[,3]/2,y1 = plotPointsCombined[,2],code = 3,angle = 90,length = .05,lwd=2)
points(plotPointsCombined[,1:2],pch=19,col=colTrait[plotPointsCombined[,"type"]],cex=2.5)
points(plotPointsCombined[,1:2],cex=2.5,lwd=2)
text(x = plotPointsCombined[,1]+plotPointsCombined[,3]/2,y=plotPointsCombined[,2],labels=labelVec,pos=4)
legend(x = "topleft",legend = c("transcriptome","phosphoproteome"),pch=19,col = colTrait[c("e","phospho")])
dev.off()