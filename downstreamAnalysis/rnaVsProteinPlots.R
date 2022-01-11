library(seqinr)
library(topGO)
source("lib/redmineTbls.R")
source("lib/general_function.R")

load("data/expressionLevel.RData")
load("data/proteinLevel.RData")
load("data/eQTL_results160831.RData")
load("data/pQTL_results170815.RData")
load("data/binInfo.RData")
load("data/hsInfo.RData")
load("data/phosphoLevel.RData")
rownames(phospho2prot) <- phospho2prot[,1]
meta <- read.table("metadata/metadata.tsv", sep="\t", header=T)
rownames(meta) <- meta$culture
gff <- read.table("/cellnet/phosphoQTL/Saccharomyces_cerevisiae/sacCer3/S288C_R6411_CDSonly.gff",sep="\t",comment.char = "#",quote = "",as.is=T)

genotype <- read.table("data/genotype_for_mapping.tsv",header=T,check.names = F)
colnames(genotype)[1]<-"chr"
geno <- t(genotype[,-(1:3)])
ePheno <- sapply(rownames(geno),FUN=function(s){
  rowMeans(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]==s,drop=F],na.rm=T)
})

eEffect <- lapply(1:ncol(geno),FUN=function(m){
  rowMeans(ePheno[,geno[,m]==0,drop=F],na.rm=T)-rowMeans(ePheno[,geno[,m]==1,drop=F],na.rm=T)
})
eEffect <- do.call("cbind",eEffect)
rownames(eEffect) <- rownames(ePheno)

pPheno <- sapply(rownames(geno),FUN=function(s){
  rowMeans(proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]==s,drop=F],na.rm=T)
})
pEffect <- lapply(1:ncol(geno),FUN=function(m){
  rowMeans(pPheno[,geno[,m]==0,drop=F],na.rm=T)-rowMeans(pPheno[,geno[,m]==1,drop=F],na.rm=T)
})
pEffect <- do.call("cbind",pEffect)
rownames(pEffect) <- rownames(pPheno)

phPheno <- sapply(rownames(geno),FUN=function(s){
  rowMeans(phosphoLevelBatchCorrected[,meta[colnames(phosphoLevelBatchCorrected),"strain"]==s,drop=F],na.rm=T)
})
phEffect <- lapply(1:ncol(geno),FUN=function(m){
  rowMeans(phPheno[,geno[,m]==0,drop=F],na.rm=T)-rowMeans(phPheno[,geno[,m]==1,drop=F],na.rm=T)
})
phEffect <- do.call("cbind",phEffect)
rownames(phEffect) <- rownames(phPheno)

#term lists
goTbl <- read.table("data/gene_association.sgd",comment.char = "!",sep="\t",quote="",as.is=T)
goTbl[,11] <- gsub("\\|{1}.*","",goTbl[,11],fixed = F)

uniGenes <- unique(gff[,9])
BP_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="P"),]
  return(as.character(sub_go[,5]))
})
names(BP_annotation) <- uniGenes

MF_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="F"),]
  return(as.character(sub_go[,5]))
})
names(MF_annotation) <- uniGenes

CC_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="C"),]
  return(as.character(sub_go[,5]))
})
names(CC_annotation) <- uniGenes

allGenes <- rep(0,length(uniGenes))
names(allGenes) <- uniGenes
allGenes[1:(length(allGenes)-1)] <- 1
allGenes <- as.factor(allGenes)

GOdataBP <- new("topGOdata",
                ontology = "BP",
                allGenes = allGenes,
                annot = annFUN.gene2GO,
                gene2GO = BP_annotation,
                nodeSize=1)
resultFisherBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")

GOdataMF <- new("topGOdata",
                ontology = "MF",
                allGenes = allGenes,
                annot = annFUN.gene2GO,
                gene2GO = MF_annotation,
                nodeSize=1)
resultFisherMF <- runTest(GOdataMF, algorithm = "classic", statistic = "fisher")

GOdataCC <- new("topGOdata",
                ontology = "CC",
                allGenes = allGenes,
                annot = annFUN.gene2GO,
                gene2GO = CC_annotation,
                nodeSize=1)
resultFisherCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")

BPterms <- names(resultFisherBP@score)
BPtermList <- lapply(BPterms,FUN=function(term){
  unlist(genesInTerm(object = GOdataBP,whichGO = term))
})
names(BPtermList) <- BPterms

MFterms <- names(resultFisherMF@score)
MFtermList <- lapply(MFterms,FUN=function(term){
  unlist(genesInTerm(object = GOdataMF,whichGO = term))
})
names(MFtermList) <- MFterms

CCterms <- names(resultFisherCC@score)
CCtermList <- lapply(CCterms,FUN=function(term){
  unlist(genesInTerm(object = GOdataCC,whichGO = term))
})
names(CCtermList) <- CCterms

termList <- c(BPtermList,MFtermList,CCtermList)
#effect sizes
eqtlEffectSizes <- t(sapply(QtlList$FDR25,FUN=function(qtl){
  target <- qtl$target
  marker <- qtl$mostSignificantPredictor
  pv <- qv[target,marker]
  pvlevel <- which(c(0.01,0.05,0.1,0.15,0.2,0.25)>pv)[1]
  eEffect <- eEffect[target,marker]
  targetName <- rownames(ePheno)[target]
  if(targetName%in%rownames(pQTL_results$ptQTL$qv)){
    pEffect <- pEffect[targetName,marker]
    inPT <- 1
  }else{
    pEffect <- NA
    inPT <- 0
  }
  return(c(target,pv,pvlevel,marker,eEffect,pEffect,inPT))
}))
colnames(eqtlEffectSizes) <- c("target","pv","pvlevel","marker","eEffect","pEffect","inPT")
#eqtlEffectSizes <- cbind(eqtlEffectSizes,local=as.numeric(localeQTLlists$FDR25))
fcDiffThreshold <- 0.15
pThreshold <- 0.1
eqtlEffectSizes <- cbind(eqtlEffectSizes,effectClass=apply(eqtlEffectSizes[,5:6],1,FUN=function(x){
  if(any(is.na(x))){return(NA)}
  r <- (x[2]-x[1])*ifelse(x[1]>0,1,-1)
  if(r>fcDiffThreshold){return(2)}
  if(r<(fcDiffThreshold*(-1))){return(3)}
  return(1)
}))

#Are QTL from the same gene more likely to be of the same class?
subEqtl <- eqtlEffectSizes[rowSums(is.na(eqtlEffectSizes))==0,]
allPairs <- expand.grid(1:nrow(subEqtl),1:nrow(subEqtl))
allPairs <- allPairs[allPairs[,1]<allPairs[,2],]
sameGene <- subEqtl[allPairs[,1],"target"]==subEqtl[allPairs[,2],"target"]
sameClass <- subEqtl[allPairs[,1],"effectClass"]==subEqtl[allPairs[,2],"effectClass"]

qtlClassContMat <- matrix(c(sum(sameGene&sameClass),sum(!sameGene&sameClass),sum(sameGene&!sameClass),sum(!sameGene&!sameClass)),nrow=2)
colnames(qtlClassContMat) <- c("sameGene_SameClass","sameGene_DiffClass")
rownames(qtlClassContMat) <- c("diffGene_SameClass","diffGene_DiffClass")
fisher.test(qtlClassContMat)

par(mfrow=c(1,1))
plot(eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold,5:6],pch=20,col=sapply(c("darkgrey","darkgreen","purple")[eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold,8]],col2alpha,0.5),xlab="rna log2FC",ylab="protein log2FC")
abline(a=fcDiffThreshold,1)
abline(a=(-1)*fcDiffThreshold,1)
legend(x = "topleft",legend = c("as expected","enhanced","buffered"),pch=20,col=c("darkgrey","darkgreen","purple"))

pdfAndPng(file = "graph/eQTL_effects",8,8, expression({
  par(cex=1.5)
  plot(eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold,5:6],pch=20,col=sapply(c("darkgrey","darkgreen","purple")[eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold,8]],col2alpha,0.5),xlab="rna log2FC",ylab="protein log2FC")
  abline(a=fcDiffThreshold,1)
  abline(a=(-1)*fcDiffThreshold,1)
  legend(x = "topleft",legend = c("as expected","enhanced","buffered"),pch=20,col=c("darkgrey","darkgreen","purple"))
  abline(h=0,v=0,lty=2,col="darkgrey")}))

cor.test(eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold,5],eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold,6],use="pair")

classTargets <- lapply(1:3,FUN=function(i){
  unique(rownames(qv)[eqtlEffectSizes[which(eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,"effectClass"]==i),1]])
})
backGround <- intersect(rownames(qv)[apply(qv<=pThreshold,1,any)],rownames(pPheno))

pQtlProp <- sapply(classTargets,FUN=function(gvec){
  sum(rowSums(pQTL_results$pQTL$qv[gvec,]<=0.1)>0)/length(gvec)
})

#GO

uniGenes <- unique(goTbl[,11])
BP_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="P"),]
  return(as.character(sub_go[,5]))
})
names(BP_annotation) <- uniGenes

MF_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="F"),]
  return(as.character(sub_go[,5]))
})
names(MF_annotation) <- uniGenes

CC_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="C"),]
  return(as.character(sub_go[,5]))
})
names(CC_annotation) <- uniGenes

#plot by HS

par(mfrow=c(3,5))
invisible(sapply(which(hsOvMat[,"eQtlHotspot"]=="TRUE"),FUN=function(i){
  markers <- which(binPerMarker%in%(as.numeric(hsOvMat[i,"startBin"]):as.numeric(hsOvMat[i,"endBin"])))
  mat2plot <- eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,4]%in%markers&eqtlEffectSizes[,"inPT"]==1,]
  plot(mat2plot[,5:6],pch=20,col=sapply(c("darkgrey","darkgreen","purple")[mat2plot[,8]],col2alpha,0.5),main=hsOvMat[i,1],xlab="rna FC",ylab="protein FC")
  abline(a=fcDiffThreshold,1)
  abline(a=(-1)*fcDiffThreshold,1)
}))


i <- 20
pdfAndPng(file = "graph/ira2_effects_colors",8,8, expression({
  par(mfrow=c(1,1),cex=1.5)
  markers <- which(binPerMarker%in%(as.numeric(hsOvMat[i,"startBin"]):as.numeric(hsOvMat[i,"endBin"])))
  mat2plot <- eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,4]%in%markers&eqtlEffectSizes[,"inPT"]==1,]
  plot(mat2plot[,5:6],pch=20,col=sapply(c("darkgrey","darkgreen","purple")[mat2plot[,8]],col2alpha,0.5),main=hsOvMat[i,1],xlab="RNA log2FC",ylab="protein log2FC")
  abline(a=fcDiffThreshold,1)
  abline(a=(-1)*fcDiffThreshold,1)
}))


i <- 20
markers <- which(binPerMarker%in%(as.numeric(hsOvMat[i,"startBin"]):as.numeric(hsOvMat[i,"endBin"])))
mat2plot <- eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,4]%in%markers&eqtlEffectSizes[,"inPT"]==1,]
ira2enhanced <- rownames(qv)[mat2plot[mat2plot[,"effectClass"]==2,1]]
ira2buffered <- rownames(qv)[mat2plot[mat2plot[,"effectClass"]==3,1]]

i <- 14
markers <- which(binPerMarker%in%(as.numeric(hsOvMat[i,"startBin"]):as.numeric(hsOvMat[i,"endBin"])))
mat2plot <- eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,4]%in%markers&eqtlEffectSizes[,"inPT"]==1,]

i <- 19
markers <- which(binPerMarker%in%(as.numeric(hsOvMat[i,"startBin"]):as.numeric(hsOvMat[i,"endBin"])))
mat2plot <- eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,4]%in%markers&eqtlEffectSizes[,"inPT"]==1,]
#mat2plot2 <- eqtlEffectSizes[eqtlEffectSizes[,2]<=pThreshold&eqtlEffectSizes[,4]%in%markers,]
mkt1enhanced <- rownames(qv)[mat2plot[mat2plot[,"effectClass"]==2,1]]
mkt1buffered <- rownames(qv)[mat2plot[mat2plot[,"effectClass"]==3,1]]
mkt1up <- rownames(qv)[mat2plot[mat2plot[,5]>0,1]]
mkt1down <- rownames(qv)[mat2plot[mat2plot[,5]<0,1]]
mkt1Targtes <- intersect(targetsByHS[[i]]$e,rownames(pQTL_results$ptQTL$qv))
GOanalysisMKT1 <- lapply(list(mkt1up,mkt1down,mkt1enhanced,mkt1buffered),FUN=function(hits){
  hits <-  as.factor(as.numeric(mkt1Targtes%in%hits))
  names(hits) <- mkt1Targtes
  
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = BP_annotation,
                  nodeSize=3)
  resultFisherBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
  allResBP <- GenTable(object=GOdataBP, elimFisher = resultFisherBP,topNodes=geneData(resultFisherBP)[[4]])
  allResBP <- subset(allResBP,resultFisherBP@score[allResBP[,1]]<=0.01)
  
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = MF_annotation,
                  nodeSize=3)
  resultFisherMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
  allResMF <- GenTable(object=GOdataMF, elimFisher = resultFisherMF,topNodes=geneData(resultFisherMF)[[4]])
  allResMF <- subset(allResMF,resultFisherMF@score[allResMF[,1]]<=0.01)
  
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = CC_annotation,
                  nodeSize=3)
  resultFisherCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "fisher")
  allResCC <- GenTable(object=GOdataCC, elimFisher = resultFisherCC,topNodes=geneData(resultFisherCC)[[4]])
  allResCC <- subset(allResCC,resultFisherCC@score[allResCC[,1]]<=0.01)
  
  return(list(allResBP,allResMF,allResCC))
})
names(GOanalysisMKT1) <- c("up","down","enhanced","buffered")



invisible(sapply(GOanalysisMKT1[[1]],printRedmine))
invisible(sapply(GOanalysisMKT1[[2]],printRedmine))
invisible(sapply(GOanalysisMKT1[[3]],printRedmine))
invisible(sapply(GOanalysisMKT1[[4]],printRedmine))

hasPT <- as.factor(as.numeric(apply(pQTL_results$ptQTL$qv<=0.1,1,any)))
names(hasPT) <- rownames(pQTL_results$ptQTL$qv)
hasE <- as.factor(as.numeric(apply(qv<=0.1,1,any)))
names(hasE) <- rownames(qv)
hasP <- as.factor(as.numeric(apply(pQTL_results$pQTL$qv<=0.1,1,any)))
names(hasP) <- rownames(pQTL_results$pQTL$qv)
onlyPT <- as.factor(as.numeric(apply(pQTL_results$ptQTL$qv<=0.1,1,any)&!apply(qv[rownames(pQTL_results$ptQTL$qv),]<=0.1,1,any)))
names(onlyPT) <- rownames(pQTL_results$ptQTL$qv)
targetEnrichment <- lapply(list(hasE,hasPT,hasP,onlyPT),FUN=function(hits){
  
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = BP_annotation,
                  nodeSize=3)
  resultFisherBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
  allResBP <- GenTable(object=GOdataBP, elimFisher = resultFisherBP,topNodes=geneData(resultFisherBP)[[4]])
  allResBP <- subset(allResBP,resultFisherBP@score[allResBP[,1]]<=0.01)
  
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = MF_annotation,
                  nodeSize=3)
  resultFisherMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
  allResMF <- GenTable(object=GOdataMF, elimFisher = resultFisherMF,topNodes=geneData(resultFisherMF)[[4]])
  allResMF <- subset(allResMF,resultFisherMF@score[allResMF[,1]]<=0.01)
  
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = CC_annotation,
                  nodeSize=3)
  resultFisherCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "fisher")
  allResCC <- GenTable(object=GOdataCC, elimFisher = resultFisherCC,topNodes=geneData(resultFisherCC)[[4]])
  allResCC <- subset(allResCC,resultFisherCC@score[allResCC[,1]]<=0.01)
  
  return(list(allResBP,allResMF,allResCC))
})
names(targetEnrichment) <- c("e","pt","p","onlyPt")

eTargets <- sapply(QtlList$FDR10,FUN=function(qtl){qtl$target})
pTargets <- sapply(pQTL_results$pQTL$QtlList$FDR10,FUN=function(qtl){qtl$target})
ptTargets <- sapply(pQTL_results$ptQTL$QtlList$FDR10,FUN=function(qtl){qtl$target})
phTargets <- sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10,FUN=function(qtl){qtl$target})
# qtlUnionEP <- lapply(rownames(pQTL_results$ptQTL$qv),FUN=function(g){
#   gn <<- g
#   eI <- which(rownames(qv)==g)
#   pI <- which(rownames(pQTL_results$pQTL$qv)==g)
#   eqtl <- QtlList$FDR10[which(eTargets==eI)]
#   eqtl <- do.call("rbind",lapply(eqtl,FUN=function(qtl){return(cbind(qtl$predictors,qtl$mostSignificantPredictor))}))
#   pqtl <- pQTL_results$pQTL$QtlList$FDR10[which(pTargets==pI)]
#   pqtl <- do.call("rbind",lapply(pqtl,FUN=function(qtl){return(cbind(qtl$predictors,qtl$mostSignificantPredictor))}))
#   out <- NULL
#   if(!is.null(pqtl)){
#     inE <- apply(pqtl,1,FUN=function(x){
#       any(eqtl[,1]<=x[2]&eqtl[,2]>=x[1])
#     })
#     pqtl <- pqtl[!inE,,drop=F]
#     if(nrow(pqtl)==0){
#       pqtl <- NULL
#     }
#   }
#   if(!is.null(pqtl)){
#     out <- rbind(out,cbind(pqtl,1))
#   }
#   if(!is.null(eqtl)){
#     out <- rbind(out,cbind(eqtl,0))
#   }
#   return(out)
#   # eLead <- unique(eqtl[,3])
#   # eLead <- cbind(eLead,T)
#   # pLead <- cbind(pLead,F)
#   # if(ncol(pLead)==1){pLead <- NULL}
#   # if(!is.null(pLead)){
#   #   allLead <- rbind(eLead,pLead)
#   # }else{
#   #   allLead <- eLead
#   # }
#   # return(allLead)
# })
# names(qtlUnionEP) <- rownames(pQTL_results$ptQTL$qv)

corThreshold <- 0.5
qtlUnionEP <- lapply(rownames(pQTL_results$ptQTL$qv),FUN=function(g){
  gn <<- g
  eI <- which(rownames(qv)==g)
  ptI <- which(rownames(pQTL_results$ptQTL$qv)==g)
  eqtl <- QtlList$FDR10[which(eTargets==eI)]
  eqtl <- unique(sapply(eqtl,FUN=function(qtl){return(qtl$mostSignificantPredictor)}))
  ptqtl <- pQTL_results$ptQTL$QtlList$FDR10[which(ptTargets==ptI)]
  ptqtl <- unique(sapply(ptqtl,FUN=function(qtl){return(qtl$mostSignificantPredictor)}))
  
  out <- NULL
  if(length(ptqtl)>0&length(eqtl)>0){
    inE <- sapply(ptqtl,FUN=function(x){
      max(abs(apply(geno[,eqtl,drop=F],2,cor,y=geno[,x],use="pair")))>=corThreshold
    })
    ptqtl <- ptqtl[!inE]
    
  }
  if(length(ptqtl)>0){
    out <- rbind(out,cbind(ptqtl,1))
  }
  if(length(eqtl)>0){
    out <- rbind(out,cbind(eqtl,0))
  }
  return(out)
  # eLead <- unique(eqtl[,3])
  # eLead <- cbind(eLead,T)
  # pLead <- cbind(pLead,F)
  # if(ncol(pLead)==1){pLead <- NULL}
  # if(!is.null(pLead)){
  #   allLead <- rbind(eLead,pLead)
  # }else{
  #   allLead <- eLead
  # }
  # return(allLead)
})
names(qtlUnionEP) <- rownames(pQTL_results$ptQTL$qv)
# 
# plotMatListEP <- lapply(names(qtlUnionEP),FUN=function(g){
#   cbind(eEffect[g,qtlUnionEP[[g]][,3]],pEffect[g,qtlUnionEP[[g]][,3]],qtlUnionEP[[g]][,3],qtlUnionEP[[g]][,4])
# })
# names(plotMatListEP) <- names(qtlUnionEP)


plotMatListEP <- lapply(names(qtlUnionEP),FUN=function(g){
  cbind(eEffect[g,qtlUnionEP[[g]][,1]],pEffect[g,qtlUnionEP[[g]][,1]],qtlUnionEP[[g]][,1],qtlUnionEP[[g]][,2])
})
names(plotMatListEP) <- names(qtlUnionEP)

qtlUnionPPH <- lapply(rownames(pQTL_results$phosphoProtResidualsQTL$qv),FUN=function(ph){
  phI <- which(rownames(pQTL_results$phosphoLevelQTL$qv)==ph)
  g <- phospho2prot[ph,2]
  pI <- which(rownames(pQTL_results$pQTL$qv)==g)
  pqtl <- pQTL_results$pQTL$QtlList$FDR10[which(pTargets==pI)]
  pqtl <- do.call("rbind",lapply(pqtl,FUN=function(qtl){return(cbind(qtl$predictors,qtl$mostSignificantPredictor))}))
  phqtl <- pQTL_results$phosphoLevelQTL$QtlList$FDR10[which(phTargets==phI)]
  phqtl <- do.call("rbind",lapply(phqtl,FUN=function(qtl){return(cbind(qtl$predictors,qtl$mostSignificantPredictor))}))
  if(!is.null(phqtl)){
    inP <- apply(phqtl,1,FUN=function(x){
      any(pqtl[,1]<=x[2]&pqtl[,2]>=x[1])
    })
    phLead <- unique(phqtl[!inP,3])
  }else{
    phLead <- NULL
  }
  pLead <- unique(pqtl[,3])
  allLead <- c(phLead,pLead)
  return(allLead)
})
names(qtlUnionPPH) <- rownames(pQTL_results$phosphoProtResidualsQTL$qv)
plotMatListPPH <- lapply(names(qtlUnionPPH),FUN=function(ph){
  cbind(pEffect[phospho2prot[ph,2],qtlUnionPPH[[ph]]],phEffect[ph,qtlUnionPPH[[ph]]],qtlUnionPPH[[ph]])
})
names(plotMatListPPH) <- names(qtlUnionPPH)
phProts <- unique(phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv),2])
pToPh <- lapply(phProts,FUN=function(g){
  ph <- phospho2prot[phospho2prot[,2]==g,1]
  pI <- which(rownames(pQTL_results$pQTL$qv)==g)
  if(pI%in%pTargets){
    leadMarkers <- sapply(pQTL_results$pQTL$QtlList$FDR10[which(pTargets==pI)],FUN=function(qtl){qtl$mostSignificantPredictor})
    out <- do.call("rbind",lapply(ph,FUN=function(pep){
      cbind(pEffect[g,leadMarkers],phEffect[pep,leadMarkers])
    }))
    return(out)
  }else{
    return(NULL)
  }
  
})
pToPh <- do.call("rbind",pToPh)
cor(pToPh)

subEP <- do.call("rbind",plotMatListEP)
subPPH <- do.call("rbind",plotMatListPPH)
fcThreshold <- 0.15
p1Range <- range(as.vector(subEP[,1:2]))
colVec <- apply(subEP,1,FUN=function(x){
  if(x[4]==0){
    fcDiff <- abs(x[1]-x[2])
    if(fcDiff<fcThreshold){
      out <- "lightgrey"
    }else{
      r <- (x[2]-x[1])*ifelse(x[1]>0,1,-1)
      if(r>fcThreshold){out <- "darkgreen"}
      if(r<(fcThreshold*(-1))){out <- "purple"}
    }
  }else{
    out <- "gold"
  }
  return(out)
})
pdfAndPng(file = "graph/full_effects_R_P_190129",8,8, expression({
  par(mfrow=c(1,1),cex=1.5)
  plot(subEP[,1:2],xlab="RNA log2FC",ylab="protein log2FC",pch=20,xlim=p1Range,ylim=p1Range,col=sapply(colVec,col2alpha,0.5))
  abline(h=0,v=0,lty=2,col="darkgrey")
  abline(a=fcThreshold,b=1,lty=1,lwd=2,col="black")
  abline(a=fcThreshold*(-1),b=1,lty=1,lwd=2,col="black")
  legend(x = "topleft",legend = c("similar","enhanced","buffered","protein only"),pch=20,col=c("darkgrey","darkgreen","purple","gold"))
}))

pdfAndPng(file = "graph/eqtl_effects_R_P_with_color_190107",8,8, expression({
  par(mfrow=c(1,1),cex=1.5)
  plot(subEP[show,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,xlim=p1Range,ylim=p1Range,col=sapply(colVec[show],col2alpha,0.5))
  abline(h=0,v=0,lty=2,col="darkgrey")
  abline(a=fcThreshold,b=1,lty=1,lwd=2,col="black")
  abline(a=fcThreshold*(-1),b=1,lty=1,lwd=2,col="black")
  legend(x = "topleft",legend = c("as expected","enhanced","buffered"),pch=20,col=c("darkgrey","darkgreen","purple"))
}))

show <- colVec!="gold"
pdfAndPng(file = "graph/eqtl_effects_R_P_no_color_190107",8,8, expression({
  par(mfrow=c(1,1),cex=1.5)
  plot(subEP[show,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,xlim=p1Range,ylim=p1Range,col=rgb(0,0,0,0.5))
  abline(h=0,v=0,lty=2,col="darkgrey")
  
}))
pdfAndPng(file = "graph/eqtl_effects_R_P_no_color_190107",8,8, expression({
  par(mfrow=c(1,1),cex=1.5)
  plot(subEP[show,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,xlim=p1Range,ylim=p1Range,col=rgb(0,0,0,0.5))
  abline(h=0,v=0,lty=2,col="darkgrey")
  
}))
###complex vs rest
uniProt2geneSymbol <- read.table(file = "Saccharomyces_cerevisiae/uniprot2ens.tab",header = T,sep = "\t",quote="",as.is=T)
rownames(uniProt2geneSymbol) <- uniProt2geneSymbol[,3]
complexAnnotation <- read.table("Saccharomyces_cerevisiae/saccharomyces_cerevisiae_complexes.tsv",sep="\t",as.is=T,quote="",comment.char = "",header = T)
complexMembers <- lapply(complexAnnotation[,5],FUN=function(idString){
  ids <- strsplit(x = idString,split = "|",fixed=T)[[1]]
  ids <- gsub(pattern = "\\(.*",replacement = "",x = ids)
  geneSymbols <- uniProt2geneSymbol[ids,1]
  geneSymbols <- geneSymbols[!is.na(geneSymbols)]
  return(geneSymbols)
})
names(complexMembers) <- complexAnnotation[,1]
#complexMembers <- complexMembers[sapply(complexMembers,length)>1]
complexProteins <- intersect(unique(unlist(complexMembers)),intersect(rownames(pQTL_results$pQTL$qv),rownames(qv)))
riboGenes <- unique(c(termList[["GO:0005840"]],termList[["GO:0005730"]]))
plotList <- cbind(do.call("rbind",plotMatListEP),unlist(lapply(names(plotMatListEP),FUN=function(g){rep(g%in%complexProteins,nrow(plotMatListEP[[g]]))})),unlist(lapply(names(plotMatListEP),FUN=function(g){rep(ifelse(g%in%riboGenes,1,0),nrow(plotMatListEP[[g]]))})))

#are weak effects transmitted?
plotListSub <- plotList[abs(plotList[,1])<=0.2&plotList[,4]==0&plotList[,6]==0,]
cor(plotListSub[plotListSub[,5]==0,1:2])
cor(plotListSub[plotListSub[,5]==1,1:2])
pdfAndPng(file = "graph/eqtl_effects_complex_weakeqtl_190220",16,8, expression({
  ymax <- max(abs(plotListSub[,2]))
  xmax <- 0.2
  par(mfrow=c(1,2),cex=1.5)
  plot(plotListSub[plotListSub[,5]==0,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.5),xlim=c(xmax*(-1),xmax),ylim=c(ymax*(-1),ymax),main="not in complex")
  abline(h=0,v=0,lty=2,col="red")
  plot(plotListSub[plotListSub[,5]==1,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.5),xlim=c(xmax*(-1),xmax),ylim=c(ymax*(-1),ymax),main="in complex")
  abline(h=0,v=0,lty=2,col="red")
}))

complexLM <- lm(plotListSub[,2]~plotListSub[,1]+plotListSub[,5]+plotListSub[,1]:plotListSub[,5])
anova(complexLM)
plotListSub2 <- plotList[plotList[,4]==0&plotList[,6]==0,]
pdfAndPng(file = "graph/eqtl_effects_complex_alleqtl_190220",16,8, expression({
  ymax <- max(abs(plotListSub2[,2]))
  xmax <- max(abs(plotListSub2[,1]))
  par(mfrow=c(1,2),cex=1.5)
  plot(plotListSub2[plotListSub2[,5]==0,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.5),xlim=c(xmax*(-1),xmax),ylim=c(ymax*(-1),ymax),main="not in complex")
  abline(h=0,v=0,lty=2,col="red")
  plot(plotListSub2[plotListSub2[,5]==1,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.5),xlim=c(xmax*(-1),xmax),ylim=c(ymax*(-1),ymax),main="in complex")
  abline(h=0,v=0,lty=2,col="red")
}))
complexLM2 <- lm(plotListSub2[,2]~plotListSub2[,1]+plotListSub2[,5]+plotListSub2[,1]:plotListSub2[,5])
anova(complexLM2)
#do complex proteins have lower variance compared to their rna variance?
complexSet <- setdiff(complexProteins,union(termList[["GO:0005730"]],riboGenes))
nonComplexSet <- setdiff(intersect(rownames(eBatchGeneLengthCorrected),rownames(proteinLevelBatchCorrected)),c(termList[["GO:0005730"]],riboGenes,complexProteins))
load("data/expressionLevel.RData")
load("data/proteinLevel.RData")
load("data/phosphoNew/phosphoProt.RData")
load("data/protRna.RData")
strains2use <- intersect(meta[colnames(eBatchGeneLengthCorrected),"strain"],intersect(meta[colnames(proteinLevelBatchCorrected),"strain"],meta[colnames(phosphoProtResiduals),"strain"]))
eByStrain <- sapply(strains2use,FUN=function(s){
  rowMeans(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]==s,drop=F],na.rm=T)
})
pByStrain <- sapply(strains2use,FUN=function(s){
  rowMeans(proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]==s,drop=F],na.rm=T)
})
rnaSD <- apply(eByStrain,1,sd,na.rm=T)
protSD <- apply(pByStrain,1,sd,na.rm=T)
maxCoord <- 0.5
plot(rnaSD[nonComplexSet],protSD[nonComplexSet],ylim=c(0,maxCoord),xlim=c(0,maxCoord))
plot(rnaSD[complexSet],protSD[complexSet],ylim=c(0,maxCoord),xlim=c(0,maxCoord))
boxplot(rnaSD[nonComplexSet],rnaSD[complexSet])
boxplot(protSD[nonComplexSet],protSD[complexSet])
byRnaSD <- apply(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]=="BY4716"],1,sd,na.rm=T)
byProtSD <- apply(proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]=="BY4716"],1,sd,na.rm=T)
rmRnaSD <- apply(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]=="RM11-1a"],1,sd,na.rm=T)
rmProtSD <- apply(proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]=="RM11-1a"],1,sd,na.rm=T)
rnaWithinSampleSD <- lapply(names(which(table(meta[colnames(eBatchGeneLengthCorrected),"strain"])>1)),FUN=function(s){
  mat <- eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]==s]
  return(mat-rowMeans(mat,na.rm=T))
})
rnaWithinSampleSD <- apply(do.call("cbind",rnaWithinSampleSD),1,sd,na.rm=T)
remRnaSD <- rnaSD-rnaWithinSampleSD
protWithinSampleSD <- lapply(names(which(table(meta[colnames(proteinLevelBatchCorrected),"strain"])>1)),FUN=function(s){
  mat <- proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]==s]
  return(mat-rowMeans(mat,na.rm=T))
})
protWithinSampleSD <- apply(do.call("cbind",protWithinSampleSD),1,sd,na.rm=T)
remProtSD <- protSD-protWithinSampleSD

#together
plot(protSD[c(complexSet,nonComplexSet)],byProtSD[c(complexSet,nonComplexSet)],col=ifelse(c(complexSet,nonComplexSet)%in%complexSet,"red","blue"),pch=20)
#separate
plot(protSD[nonComplexSet],byProtSD[nonComplexSet],col="blue",pch=20,xlim=c(0,2),ylim=c(0,2))
abline(a=0,b=1,lwd=3)
plot(protSD[complexSet],byProtSD[complexSet],col="red",pch=20,xlim=c(0,2),ylim=c(0,2))
abline(a=0,b=1,lwd=3)
plot(protSD[nonComplexSet],rmProtSD[nonComplexSet],col="blue",pch=20,xlim=c(0,2),ylim=c(0,2))
abline(a=0,b=1,lwd=3)
plot(protSD[complexSet],rmProtSD[complexSet],col="red",pch=20,xlim=c(0,2),ylim=c(0,2))
abline(a=0,b=1,lwd=3)
boxplot((protSD-(byProtSD+rmProtSD)/2)[nonComplexSet],(protSD-(byProtSD+rmProtSD)/2)[complexSet],outline=F)
boxplot((rnaSD-(byRnaSD+rmRnaSD)/2)[nonComplexSet],(rnaSD-(byRnaSD+rmRnaSD)/2)[complexSet],outline=F)

###local vs distant
gff <- read.table(file = "Saccharomyces_cerevisiae/RM/genes_RM_CDSonly_refined.gtf",stringsAsFactors = F)
CDS_position<-t(sapply(unique(gff[,9]),function(gene){
  sel<-gff[,9]==gene
  pos<-c(gff[sel,4],gff[sel,5])
  pos<-range(pos)
  chr<-gff[sel,1][1]
  return(list(chr,pos[1],pos[2]))
}))
gene<-CDS_position
gene<-data.frame(chr=unlist(gene[,1]),
                 start=unlist(gene[,2]),
                 end=unlist(gene[,3]),
                 stringsAsFactors = F)
gene$med <- rowMeans(cbind(gene$start,gene$end))
gene <- gene[rownames(pQTL_results$ptQTL$qv),]
geneMarker<-sapply(1:nrow(gene),function(i){
  before <- max(which(as.character(genotype$chr)==as.character(gene$chr[i]) &
                        genotype$end<gene$med[i]))
  if(is.infinite(before)) {before =min(which(as.character(genotype$chr)==as.character(gene$chr[i])))}
  after <- min(which(as.character(genotype$chr)==as.character(gene$chr[i]) &
                       genotype$start>gene$med[i]))
  if(is.infinite(after)) {after =max(which(as.character(genotype$chr)==as.character(gene$chr[i])))}
  return(c(before,after))
})
gene <- cbind(gene,marker_1=geneMarker[1,],marker_2=geneMarker[2,],stringsAsFactors=F)
local <- unlist(lapply(1:length(qtlUnionEP),FUN=function(i){
  i <<- i
  if(is.null(qtlUnionEP)){return(NULL)}
  unlist(lapply(qtlUnionEP[[i]][,1],FUN=function(j){
    any(abs(apply(geno[,unlist(gene[i,5:6]),drop=F],2,FUN=cor,y=geno[,j],use="pair"))>0.8)
  }))
  
}))
localLM <- lm(subEP[show&local,2]~subEP[show&local,1]+0)
distantLM <- lm(subEP[show&!local,2]~subEP[show&!local,1]+0)
pdfAndPng(file = "graph/eqtl_effects_R_P_local_dist_190107",8,8, expression({
  par(mfrow=c(1,1),cex=1.5)
  plot(subEP[show,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,xlim=p1Range,ylim=p1Range,col=sapply(ifelse(local,"red","blue")[show],col2alpha,0.5))
  abline(h=0,v=0,lty=2,col="darkgrey")
  legend(x = "topleft",legend = paste(c("local:","distant:"),paste0("m = ",round(c(localLM$coefficients[1],distantLM$coefficients[1]),digits=2))),lty=1,col=c("red","blue"))
  abline(a=0,b=localLM$coefficients[1],col="red")
  abline(a=0,b=distantLM$coefficients[1],col="blue")
}))

hasPT <- names(plotMatListEP)[sapply(plotMatListEP,FUN=function(x){any(x[,4]==1)})]
backGround <- rownames(pQTL_results$ptQTL$qv)
GObyClass <- lapply(c(classTargets,list(hasPT)),FUN=function(x){
  
  hits <- as.factor(as.numeric(backGround%in%x))
  names(hits) <- backGround
  
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = BP_annotation,
                  nodeSize=3)
  resultFisherBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
  allResBP <- GenTable(object=GOdataBP, pvalue = resultFisherBP,topNodes=geneData(resultFisherBP)[[4]])
  allResBP <- subset(allResBP,resultFisherBP@score[allResBP[,1]]<=0.01)
  
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = MF_annotation,
                  nodeSize=3)
  resultFisherMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
  allResMF <- GenTable(object=GOdataMF, pvalue = resultFisherMF,topNodes=geneData(resultFisherMF)[[4]])
  allResMF <- subset(allResMF,resultFisherMF@score[allResMF[,1]]<=0.01)
  
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = hits,
                  annot = annFUN.gene2GO,
                  gene2GO = CC_annotation,
                  nodeSize=3)
  resultFisherCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "fisher")
  allResCC <- GenTable(object=GOdataCC, pvalue = resultFisherCC,topNodes=geneData(resultFisherCC)[[4]])
  allResCC <- subset(allResCC,resultFisherCC@score[allResCC[,1]]<=0.01)
  
  return(list(allResBP,allResMF,allResCC))
})
names(GObyClass) <- c("similar","enhanced","buffered","proteinOnly")
invisible(sapply(GObyClass[[1]],printRedmine))
invisible(sapply(GObyClass[[2]],printRedmine))
invisible(sapply(GObyClass[[3]],printRedmine))
invisible(sapply(GObyClass[[4]],printRedmine))

ontologies <- c("BP","MF","CC")
mat2write <- lapply(GObyClass,FUN=function(goList){
  matlist <- lapply(1:3,FUN=function(i){
    if(nrow(goList[[i]])==0){return(NULL)}
    out <- cbind(ontology=ontologies[i],goList[[i]])
    return(out)
  })
  combinedMat <- do.call("rbind",matlist)
  return(combinedMat)
})
# 
# write.table(x = mat2write[[1]],file = "manuscript/suppData/GO_similar.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
# write.table(x = mat2write[[2]],file = "manuscript/suppData/GO_enhanced.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
# write.table(x = mat2write[[3]],file = "manuscript/suppData/GO_buffered.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
# write.table(x = mat2write[[4]],file = "manuscript/suppData/GO_protein_only.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

pRangeGO <- 21

# sapply(1:length(GObyClass),FUN=function(i){
#   
# })

fisher.test(rownames(pQTL_results$ptQTL$qv)%in%hasPT,rownames(pQTL_results$ptQTL$qv)%in%classTargets[[1]],alternative="greater")$p.value*3
fisher.test(rownames(pQTL_results$ptQTL$qv)%in%hasPT,rownames(pQTL_results$ptQTL$qv)%in%classTargets[[2]],alternative="greater")$p.value*3
fisher.test(rownames(pQTL_results$ptQTL$qv)%in%hasPT,rownames(pQTL_results$ptQTL$qv)%in%classTargets[[3]],alternative="greater")$p.value*3

p2Range <- range(subEP[,2])
pdfAndPng(file = "graph/full_effects_R_P_PH",16,8, expression({
  par(mfrow=c(1,2),cex=1.5)
  plot(subEP[,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p1Range,ylim=p1Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
  plot(subPPH[,1:2],xlab="protein l2FC",ylab="phospho-peptide l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p2Range,ylim=p2Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
}))

#IRA2 effects
p1Range <- c(-2.5,1.5)
p2Range <- c(-2,0.8)
fcThreshold <- 0.15
i <- 20
hsRange <- head(which(binPerMarker==as.numeric(hsOvMat[i,"startBin"])),n=1):tail(which(binPerMarker==as.numeric(hsOvMat[i,"endBin"])),n=1)
subEP <- do.call("rbind",plotMatListEP)
subEP <- subEP[subEP[,3]%in%hsRange,]
subPPH <- do.call("rbind",plotMatListPPH)
subPPH <- subPPH[subPPH[,3]%in%hsRange,]

pdfAndPng(file = "graph/ira2_effects",16,8, expression({
  par(mfrow=c(1,2))
  plot(subEP[,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p1Range,ylim=p1Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
  plot(subPPH[,1:2],xlab="protein l2FC",ylab="phospho-peptide l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p2Range,ylim=p2Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
}))


#HAP1 effects
p1Range <- c(-2.5,1.5)
p2Range <- c(-2,0.8)
fcDiffThreshold <- 0.15
i <- 14
hsRange <- head(which(binPerMarker==as.numeric(hsOvMat[i,"startBin"])),n=1):tail(which(binPerMarker==as.numeric(hsOvMat[i,"endBin"])),n=1)
subEP <- do.call("rbind",plotMatListEP)
subEP <- subEP[subEP[,3]%in%hsRange,]
subPPH <- do.call("rbind",plotMatListPPH)
subPPH <- subPPH[subPPH[,3]%in%hsRange,]

pdfAndPng(file = "graph/hap1_effects",16,8, expression({
  par(mfrow=c(1,2),cex=2)
  plot(subEP[,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p1Range,ylim=p1Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
  plot(subPPH[,1:2],xlab="protein l2FC",ylab="phospho-peptide l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p2Range,ylim=p2Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
}))
pdfAndPng(file = "graph/hap1_effects_colors",16,8, expression({
  par(mfrow=c(1,2))
  colEP <- apply(subEP[,1:2],1,FUN=function(x){
    absDiff <- abs(x[1]-x[2])
    if(absDiff<fcDiffThreshold){
      return("darkgrey")
    }else{
      
    }
    
  })
  plot(subEP[,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p1Range,ylim=p1Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
  plot(subPPH[,1:2],xlab="protein l2FC",ylab="phospho-peptide l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p2Range,ylim=p2Range)
  abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
}))
# 
# #MKT1
# i <- 19
# hsRange <- head(which(binPerMarker==as.numeric(hsOvMat[i,"startBin"])),n=1):tail(which(binPerMarker==as.numeric(hsOvMat[i,"endBin"])),n=1)
# subEP <- do.call("rbind",plotMatListEP)
# subEP <- subEP[subEP[,3]%in%hsRange,]
# subPPH <- do.call("rbind",plotMatListPPH)
# subPPH <- subPPH[subPPH[,3]%in%hsRange,]
# 
# pdfAndPng(file = "graph/MKT1_effects",16,8, expression({
#   par(mfrow=c(1,2))
#   plot(subEP[,1:2],xlab="RNA l2FC",ylab="protein l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p1Range,ylim=p1Range)
#   abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
#   plot(subPPH[,1:2],xlab="protein l2FC",ylab="phospho-peptide l2FC",pch=20,col=rgb(0,0,0,0.4),xlim=p2Range,ylim=p2Range)
#   abline(h=0,v=0,a=0,b=1,lty=2,col="darkgrey")
# }))