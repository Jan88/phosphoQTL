scores <- lapply(1:5429,FUN=function(x){
  load(paste0("/cellnet/phosphoQTL/data/raw_mapping_data/eQTL_160829/real/trait",x,".RData"))
  scores
})
scores <- do.call("rbind",scores)

source("lib/pvalues.R")
source("lib/QTL_draft.R")
#pv <- pEst(path = "/mnt/cellnet/phosphoQTL/data/raw_mapping_data/eQTL_160829/perms/",markersPerIteration = 300,scores = scores)
pv <- pEstByFile(path = "/mnt/cellnet/phosphoQTL/data/raw_mapping_data/eQTL_160829/perms/",scores = scores,fileBatch = 200)

qv <- pv
qv[1:length(qv)] <- p.adjust(as.vector(qv),method="fdr")
load("data/eQtlMappingData160817.RData")
rownames(pv) <- rownames(qv) <- rownames(phenotype)

geno <- read.table("data/genotype_for_mapping.tsv",header=T)
chrVec <- geno[,1]
geno <- geno[,4:ncol(geno)]
geno <- t(geno)
fdrs <- c(1,5,10,15,20,25) 
QtlList <- lapply(fdrs,FUN=function(fdr){
  QTLgrouper(pmat = qv,sigThreshold = fdr/100,corThreshold = 0.8,distThreshold = 9,genotype = geno,chrVec = chrVec)
})
names(QtlList) <- paste0("FDR",fdrs)

#colors
library(RColorBrewer)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")



load("data/plotPack.RData")
names(chr_length) <- gsub(pattern = "chr",replacement = "",x = names(chr_length))
png(filename = "graph/eQTL/eQTL_FDR10.png",width = 960,height = 960)
par(cex=1.5,cex.lab=2,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = QtlList$FDR10,targetLocs = geneLocs[rownames(qv)],main="",predictorLocs = markerLocs,chrLen = chr_length,xlab = "eQTL position",ylab = "target-gene position",qtl.lwd=5,labcex = 1,col=adjustcolor(col = colors["eQTL"],alpha.f = 0.8),chrCol = rgb(0,0,0,0.3),chrLty=1)
dev.off()

pdf("graph/eQTL/eQTL_FDR10.pdf",width = 11,height = 11)
par(cex=1.5,cex.lab=1.7,mar=c(5,5,4,2)+0.1)
QTLplotter(QTLlist = QtlList$FDR10,targetLocs = geneLocs[rownames(qv)],main="",predictorLocs = markerLocs,chrLen = chr_length,xlab = "eQTL position",ylab = "target-gene position",qtl.lwd=5,labcex = .8,col=adjustcolor(col = colors["eQTL"],alpha.f = 0.8),chrCol = rgb(0,0,0,0.3),chrLty=1)
dev.off()

save(pv,qv,QtlList,file="data/eQTL_results160831.RData")

nTargets <- sapply(fdrs/100,FUN=function(fdr){
  sum(apply(qv<=fdr,1,any))
})
propTargets <- round(nTargets/nrow(qv),digits=3)
nQTL <- sapply(QtlList,length)
nLoci <- sapply(fdrs/100,FUN=function(fdr){
  sum(apply(qv<=fdr,2,any))
})
propLoci <- round(nLoci/ncol(qv),digits=3)

metaTbl <- rbind(nTargets,propTargets,nQTL,nLoci,propLoci)
printRedmine(metaTbl)