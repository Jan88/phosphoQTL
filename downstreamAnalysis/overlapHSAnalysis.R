####load data####
load("data/pQTL_results170815.RData")
load("data/eQTL_results160831.RData")
load("data/phosphoLevel.RData")
load("data/hsInfo.RData")
load("data/binInfo.RData")
rownames(phospho2prot) <- phospho2prot[,1]

####What genes are available on all molecular layers?####
useGenes <- intersect(rownames(qv),phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$pv),2])

####compute null distributions for each hotspot and each pair of layers####
phTargetsByHS <- apply(hsOvMat,1,FUN=function(hsVec){
  pmat <- pQTL_results$phosphoLevelQTL$qv
  hsMarkers <- which(binPerMarker%in%as.numeric(hsVec["startBin"]):as.numeric(hsVec["endBin"]))
  targets <- rownames(pmat)[rowSums(pmat[,hsMarkers]<.1)>0]
  return(targets)
})
phrTargetsByHS <- apply(hsOvMat,1,FUN=function(hsVec){
  pmat <- pQTL_results$phosphoProtResidualsQTL$qv
  hsMarkers <- which(binPerMarker%in%as.numeric(hsVec["startBin"]):as.numeric(hsVec["endBin"]))
  targets <- rownames(pmat)[rowSums(pmat[,hsMarkers]<.1)>0]
  return(targets)
})

proteinTargets <- unique(unlist(lapply(targetsByHS,FUN=function(hs){return(hs$p)})))
proteinTargetsSub <- intersect(proteinTargets,useGenes)
eToPSubNull <- sapply(1:nrow(hsOvMat),FUN=function(nHS){
  eHits <- intersect(useGenes,targetsByHS[[nHS]]$e)
  npHits <- length(intersect(useGenes,targetsByHS[[nHS]]$p))
  nullVec <- sapply(1:1000,FUN=function(i){
    length(intersect(eHits,sample(proteinTargetsSub,size = npHits)))
  })
  return(nullVec)
})
eToPSub <- t(sapply(1:nrow(hsOvMat),FUN=function(nHS){
  eHits <- intersect(useGenes,targetsByHS[[nHS]]$e)
  pHits <- intersect(useGenes,targetsByHS[[nHS]]$p)
  nShared <- length(intersect(eHits,pHits))
  nullVec <- eToPSubNull[,nHS]
  nExpected <- mean(nullVec)
  pBigger <- max(sum(nullVec>=nShared)/length(nullVec),1/length(nullVec))
  pSmaller <- max(sum(nullVec<=nShared)/length(nullVec),1/length(nullVec))
  return(c(nShared,nExpected,pBigger,pSmaller))
}))


phrTargets <- unique(unlist(phrTargetsByHS))
phrTargetsSub <- phrTargets[phospho2prot[phrTargets,2]%in%useGenes]
eToPhrSubNull <- sapply(1:nrow(hsOvMat),FUN=function(nHS){
  eHits <- intersect(useGenes,targetsByHS[[nHS]]$e)
  nPhrHits <- length(intersect(phrTargetsSub,phrTargetsByHS[[nHS]]))
  nullVec <- sapply(1:1000,FUN=function(i){
    length(intersect(eHits,phospho2prot[sample(phrTargetsSub,size = nPhrHits),2]))
  })
  return(nullVec)
})
eToPhrSub <- t(sapply(1:nrow(hsOvMat),FUN=function(nHS){
  eHits <- intersect(useGenes,targetsByHS[[nHS]]$e)
  phrHits <- unique(phospho2prot[phrTargetsByHS[[nHS]],2])
  nShared <- length(intersect(eHits,phrHits))
  nullVec <- eToPhrSubNull[,nHS]
  nExpected <- mean(nullVec)
  pBigger <- max(sum(nullVec>=nShared)/length(nullVec),1/length(nullVec))
  pSmaller <- max(sum(nullVec<=nShared)/length(nullVec),1/length(nullVec))
  return(c(nShared,nExpected,pBigger,pSmaller))
}))

####plot null distributions and observed values####
pdf(file = "graph/overlapByHSSub.pdf",width=9)
par(mfrow=c(2,1),las=2)
boxplot(eToPSubNull,ylim=c(0,max(c(eToPSubNull,eToPSub[,1]))),main="overlap of eQTL and pQTL",ylab="overlap of targets",names=hsOvMat[,1])
points(x = 1:ncol(eToPSubNull),y=eToPSub[,1],pch=19,col="red")
boxplot(eToPhrSubNull,ylim=c(0,max(c(eToPhrSubNull,eToPhrSub[,1]))),main="overlap of eQTL and phResQTL",ylab="overlap of targets",names=hsOvMat[,1])
points(x = 1:ncol(eToPhrSubNull),y=eToPhrSub[,1],pch=19,col="red")
dev.off()
####compute global null distributions for each pair of layers (future)####




