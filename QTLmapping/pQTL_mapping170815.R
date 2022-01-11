###generate a vector so I can loop through the datasets I want to map###
homePath <- "/cellnet/phosphoQTL/data/phosphoProtein20170814/" #pQTL was copied from the last mapping to have everything in the same folder
mapVec <- list.files(homePath)

###get mapping stats###
sapply(mapVec,FUN=function(dir){
  permObjects <- length(list.files(paste0(homePath,dir,"/perms/")))
  realObjects <- length(list.files(paste0(homePath,dir,"/real/")))
  return(c(permObjects,realObjects))
})

#compare to the input
#phosphoLevelQTL
load("data/phosphoLevelQtlMappingData170727.RData")
nrow(phenotype4perm)
nrow(phenotype)
#phosphoProtResidualsQTL
load("data/phosphoProtResidualsQtlMappingData170727.RData")
nrow(phenotype4perm)
nrow(phenotype)
#phosphoRnaResidualsQTL
load("data/phosphoRnaResidualsQtlMappingData170727.RData")
nrow(phenotype4perm)
nrow(phenotype)
#phosphoRnaResidualsQTL
load("data/pQtlMappingData161117.RData")
nrow(phenotype4perm)
nrow(phenotype)
#ptQTL
load("data/ptQtlMappingData170801.RData")
nrow(phenotype4perm)
nrow(phenotype)
#everything is complete

###get p-values###
source("lib/pvalues.R")
source("lib/QTL_draft.R")
fdrs <- c(1,5,10,15,20,25) 
geno <- read.table("data/genotype_for_mapping.tsv",header=T)
chrVec <- geno[,1]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)
pList <- lapply(mapVec,FUN=function(dir){
  realFiles <- list.files(paste0(homePath,dir,"/real/"))
  realFiles <- paste0("trait",1:length(realFiles),".RData")
  scoreMat <- t(sapply(realFiles,FUN=function(file){
    load(paste0(homePath,dir,"/real/",file))
    scores
  }))
  pv <- pEstByFile(path = paste0(homePath,dir,"/perms/"),scores = scoreMat,fileBatch = 100)
  qv <- pv
  qv[1:length(qv)] <- p.adjust(as.vector(qv),method="fdr")
  QtlList <- lapply(fdrs,FUN=function(fdr){
    if(sum(qv<=(fdr/100))==0){return("no significant associations")}
    QTLgrouper(pmat = qv,sigThreshold = fdr/100,corThreshold = 0.8,distThreshold = 9,genotype = geno,chrVec = chrVec)
  })
  names(QtlList) <- paste0("FDR",fdrs)
  return(list(pv=pv,qv=qv,QtlList=QtlList))
})

metaStats <- lapply(pList,FUN=function(pl){
  out <- sapply(fdrs,FUN=function(fdr){
    if(sum(pl$qv<=(fdr/100))==0){return(c(0,0,0))}
    nTraits <- sum(apply(pl$qv<=(fdr/100),1,any))
    nAssoc <- sum(pl$qv<=(fdr/100))
    nQTL <- length(pl$QtlList[[paste0("FDR",fdr)]])
    return(c(nTraits,nAssoc,nQTL))
  })
  rownames(out) <- c("nTraits","nAssociations","nQTL")
  colnames(out) <- paste0("FDR",fdrs)
  return(out)
})
source("lib/redmineTbls.R")
sapply(metaStats,printRedmine)

names(pList) <- names(metaStats) <- mapVec
#name traits
nameList <- list()
load("data/phosphoLevelQtlMappingData170727.RData")
nameList[[1]] <- rownames(phenotype)
load("data/phosphoProtResidualsQtlMappingData170727.RData")
nameList[[2]] <- rownames(phenotype)
load("data/phosphoRnaResidualsQtlMappingData170727.RData")
nameList[[3]] <- rownames(phenotype)
load("data/pQtlMappingData161117.RData")
nameList[[4]] <- rownames(phenotype)
load("data/ptQtlMappingData170801.RData")
nameList[[5]] <- rownames(phenotype)
for(i in 1:length(pList)){
  rownames(pList[[i]]$pv) <- nameList[[i]]
  rownames(pList[[i]]$qv) <- nameList[[i]]
}
pQTL_results <- pList
save(pQTL_results,file="/cellnet/phosphoQTL/data/pQTL_results170815.RData")

#plot qtlmaps
load("data/plotPack.RData")
load("data/phosphoLevel.RData")
pho2prot <- phospho2prot[,2]
names(pho2prot) <- phospho2prot[,1]

#colors
library(RColorBrewer)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

traitNames <- c("phQTL","phResQTL","pQTL","ptQTL")
names(chr_length) <- gsub("chr","",names(chr_length))
invisible(sapply(1:4,FUN=function(i){
  sapply(1:length(pQTL_results[c(1,2,4,5)][[i]]$QtlList),FUN=function(j){
    QtlList <- pQTL_results[c(1,2,4,5)][[i]]$QtlList[[j]]
    if(length(QtlList)==1){return(NA)}
    png(paste0("/cellnet/phosphoQTL/graph/pQTL/",traitNames[i],"_FDR",fdrs[j],".png"),width = 960,height = 960)
    par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
    if(any(nameList[c(1,2,4,5)][[i]]%in%names(geneLocs))){
      glocs <- geneLocs[nameList[c(1,2,4,5)][[i]]]
    }else{
      glocs <- geneLocs[pho2prot[nameList[c(1,2,4,5)][[i]]]]
    }
    QTLplotter(QTLlist = QtlList,targetLocs = glocs,predictorLocs = markerLocs,chrLen = chr_length,xlab = paste0(traitNames[i]," position"),ylab = "target-gene position",main = "",qtl.lwd = 5,chrLty=1,col=adjustcolor(col = colors[traitNames[i]],alpha.f = 0.8),chrCol = rgb(0,0,0,0.3),labcex = 1)
    dev.off()
  })
}))
invisible(sapply(1:4,FUN=function(i){
  sapply(1:length(pQTL_results[c(1,2,4,5)][[i]]$QtlList),FUN=function(j){
    QtlList <- pQTL_results[c(1,2,4,5)][[i]]$QtlList[[j]]
    if(length(QtlList)==1){return(NA)}
    pdf(paste0("/cellnet/phosphoQTL/graph/pQTL/",traitNames[i],"_FDR",fdrs[j],".pdf"),width = 11,height = 11)
    par(cex=1.5,cex.lab=1.7,mar=c(5,5,4,2)+0.1)
    if(any(nameList[c(1,2,4,5)][[i]]%in%names(geneLocs))){
      glocs <- geneLocs[nameList[c(1,2,4,5)][[i]]]
    }else{
      glocs <- geneLocs[pho2prot[nameList[c(1,2,4,5)][[i]]]]
    }
    QTLplotter(QTLlist = QtlList,targetLocs = glocs,predictorLocs = markerLocs,chrLen = chr_length,xlab = paste0(traitNames[i]," position"),ylab = "target-gene position",main = "",qtl.lwd = 5,chrLty=1,col=adjustcolor(col = colors[traitNames[i]],alpha.f = 0.8),chrCol = rgb(0,0,0,0.3),labcex = .8)
    dev.off()
  })
}))
