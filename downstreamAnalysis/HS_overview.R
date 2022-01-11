source("lib/general_function.R")
load("data/binInfo.RData")
load("data/pQTL_results170815.RData")
load("data/eQTL_results160831.RData")
load("data/phosphoLevel.RData")
rownames(phospho2prot) <- phospho2prot[,1]
gff <- read.table("/cellnet/phosphoQTL/Saccharomyces_cerevisiae/sacCer3/S288C_R6411_CDSonly.gff",sep="\t",comment.char = "#",quote = "",as.is=T)
library(topGO)
goTbl <- read.table("data/gene_association.sgd",comment.char = "!",sep="\t",quote="",as.is=T)
goTbl[,11] <- gsub("\\|{1}.*","",goTbl[,11],fixed = F)

###How many hotspots do we find?
#eQTL
sigBinsE <- which(bins[,"eQTL"]>=30)
incVec <- sigBinsE[2:length(sigBinsE)]-sigBinsE[1:(length(sigBinsE)-1)]
sameChrom <- as.character(bins[sigBinsE[-1],"chr"])==as.character(bins[sigBinsE[-length(sigBinsE)],"chr"])
eHS <- sum(incVec>1|!sameChrom)+1
#pQTL
sigBinsP <- which(bins[,"pQTL"]>=13)
incVec <- sigBinsP[2:length(sigBinsP)]-sigBinsP[1:(length(sigBinsP)-1)]
sameChrom <- as.character(bins[sigBinsP[-1],"chr"])==as.character(bins[sigBinsP[-length(sigBinsP)],"chr"])
pHS <- sum(incVec>1|!sameChrom)+1
#ptQTL
sigBinsPt <- which(bins[,"ptQTL"]>=9)
incVec <- sigBinsPt[2:length(sigBinsPt)]-sigBinsPt[1:(length(sigBinsPt)-1)]
sameChrom <- as.character(bins[sigBinsPt[-1],"chr"])==as.character(bins[sigBinsPt[-length(sigBinsPt)],"chr"])
ptHS <- sum(incVec>1|!sameChrom)+1
#phosphoLevelQTL
sigBinsPh <- which(bins[,"phosphoLevelQTL"]>=11)
incVec <- sigBinsPh[2:length(sigBinsPh)]-sigBinsPh[1:(length(sigBinsPh)-1)]
sameChrom <- as.character(bins[sigBinsPh[-1],"chr"])==as.character(bins[sigBinsPh[-length(sigBinsPh)],"chr"])
phHS <- sum(incVec>1|!sameChrom)+1
#phosphoResidualQTL
sigBinsPhr <- which(bins[,"phosphoProtResidualsQTL"]>=4)
incVec <- sigBinsPhr[2:length(sigBinsPhr)]-sigBinsPhr[1:(length(sigBinsPhr)-1)]
sameChrom <- as.character(bins[sigBinsPhr[-1],"chr"])==as.character(bins[sigBinsPhr[-length(sigBinsPhr)],"chr"])
phrHS <- sum(incVec>1|!sameChrom)+1

#collect target genes for each molecular layer and bin
targetsByBin <- lapply(1:nrow(bins),FUN=function(b){
  markers <- which(binPerMarker==b)
  eTargets <- rownames(qv)[which(apply(qv[,markers,drop=F]<0.1,1,any))]
  pTargets <- rownames(pQTL_results$pQTL$qv)[which(apply(pQTL_results$pQTL$qv[,markers,drop=F]<0.1,1,any))]
  ptTargets <- rownames(pQTL_results$ptQTL$qv)[which(apply(pQTL_results$ptQTL$qv[,markers,drop=F]<0.1,1,any))]
  phTargets <- unique(phospho2prot[rownames(pQTL_results$phosphoLevelQTL$qv)[which(apply(pQTL_results$phosphoLevelQTL$qv[,markers,drop=F]<0.1,1,any))],2])
  phrTargets <- unique(phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv)[which(apply(pQTL_results$phosphoProtResidualsQTL$qv[,markers,drop=F]<0.1,1,any))],2])
  out <- list(e=eTargets,p=pTargets,pt=ptTargets,ph=phTargets,phr=phrTargets)
  return(out)
})

#HS unification

logMat <- cbind(e=bins[,"eQTL"]>=30,pt=bins[,"ptQTL"]>=9,p=bins[,"pQTL"]>=13,phr=bins[,"phosphoProtResidualsQTL"]>=4,ph=bins[,"phosphoLevelQTL"]>=11)
hsVec <- which(apply(logMat,1,any))
incVec <- hsVec[2:length(hsVec)]-hsVec[1:(length(hsVec)-1)]
sameChrom <- as.character(bins[hsVec[-1],"chr"])==as.character(bins[hsVec[-length(hsVec)],"chr"])
hsMat <- lapply(c(1,which(incVec>1|!sameChrom)+1),FUN=function(i){
  hsEnd <- which((1:length(hsVec))>i&(c(0,incVec)>1|c(F,!sameChrom)))[1]
  if(is.na(hsEnd)){
    hsEnd <- length(hsVec)
  }
  binVec <- hsVec[i:(hsEnd-1)]
  out <- cbind(bins[binVec,],logMat[binVec,,drop=F])
  return(out)
})
chrVec <- as.character(sapply(hsMat,FUN=function(mat){mat[1,1]}))
hsNames <- sapply(1:length(chrVec),FUN=function(i){
  paste0(chrVec[i],":",sum(chrVec[(1:length(chrVec))<i]==chrVec[i])+1)
})

targetsByHS <- lapply(hsMat,FUN=function(mat){
  startBin <- as.numeric(rownames(mat)[1])
  endBin <- as.numeric(rownames(mat)[nrow(mat)])
  eTargets <- unique(unlist(lapply(startBin:endBin,FUN=function(bin){
    targetsByBin[[bin]]$e
  })))
  ptTargets <- unique(unlist(lapply(startBin:endBin,FUN=function(bin){
    targetsByBin[[bin]]$pt
  })))
  pTargets <- unique(unlist(lapply(startBin:endBin,FUN=function(bin){
    targetsByBin[[bin]]$p
  })))
  phrTargets <- unique(unlist(lapply(startBin:endBin,FUN=function(bin){
    targetsByBin[[bin]]$phr
  })))
  phTargets <- unique(unlist(lapply(startBin:endBin,FUN=function(bin){
    targetsByBin[[bin]]$ph
  })))
  out <- list(e=eTargets,p=pTargets,pt=ptTargets,ph=phTargets,phr=phrTargets)
  return(out)
})

hsOvMat <- t(sapply(1:length(hsNames),FUN=function(i){
  chr <- as.character(hsMat[[i]][1,1])
  startBin <- as.numeric(rownames(hsMat[[i]])[1])
  endBin <- as.numeric(rownames(hsMat[[i]])[nrow(hsMat[[i]])])
  startPos <- min(hsMat[[i]][,2])
  endPos <- max(hsMat[[i]][,3])
  nTargets <- sapply(targetsByHS[[i]],length)[c("e","pt","p","phr","ph")]
  hsStatus <- apply(hsMat[[i]][,c("e","pt","p","phr","ph")],2,any)
  out <- c(hsNames[i],chr,startPos,endPos,startBin,endBin,nTargets,hsStatus)
  return(out)
}))
colnames(hsOvMat) <- c("name","chr","startPos","endPos","startBin","endBin","eQtlTargets","ptQtlTargets","pQtlTargets","phResQtlTargets","phQtlTargets","eQtlHotspot","ptQtlHotspot","pQtlHotspot","phResQtlHotspot","phQtlHotspot")

#three biggest phrQTL HS: HAP1:192, 203, 237
ovMat <- lapply(c(192,203,237),FUN=function(b){
  ov_p_e <- length(intersect(targetsByBin[[b]]$p,targetsByBin[[b]]$e))/length(targetsByBin[[b]]$e)
  ov_phr_p <- length(intersect(targetsByBin[[b]]$phr,targetsByBin[[b]]$p))/length(targetsByBin[[b]]$phr)
  ov_phr_e <- length(intersect(targetsByBin[[b]]$phr,targetsByBin[[b]]$e))/length(targetsByBin[[b]]$phr)
  ov_phr_ep <- length(intersect(targetsByBin[[b]]$phr,union(targetsByBin[[b]]$e,targetsByBin[[b]]$p)))/length(targetsByBin[[b]]$phr)
  return(c(ov_p_e,ov_phr_p,ov_phr_e,ov_phr_ep))
})
ovMat <- do.call("rbind",ovMat)
round(1-ovMat[,4],digits=2)



#look at GO-enrichments per HS

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

ct <- 0
GOanalysis <- lapply(targetsByHS,FUN=function(bin){
  ct <<- ct+1
  print(ct)
  eTargets <- as.factor(as.numeric(rownames(qv)%in%bin$e))
  names(eTargets) <- rownames(qv)
  pTargets <- as.factor(as.numeric(rownames(pQTL_results$pQTL$qv)%in%bin$p))
  names(pTargets) <- rownames(pQTL_results$pQTL$qv)
  ptTargets <- as.factor(as.numeric(rownames(pQTL_results$ptQTL$qv)%in%bin$pt))
  names(ptTargets) <- rownames(pQTL_results$ptQTL$qv)
  phResTargets <- as.factor(as.numeric(unique(phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv),2])%in%bin$phr))
  names(phResTargets) <- unique(phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv),2])
  phTargets <- as.factor(as.numeric(unique(phospho2prot[rownames(pQTL_results$phosphoLevelQTL$qv),2])%in%bin$ph))
  names(phTargets) <- unique(phospho2prot[rownames(pQTL_results$phosphoLevelQTL$qv),2])
  
  byLevel <- lapply(list(eTargets,ptTargets,pTargets,phResTargets,phTargets),FUN=function(hits){
    if(all(hits==0)){return(NULL)}
    GOdataBP <- new("topGOdata",
                    ontology = "BP",
                    allGenes = hits,
                    annot = annFUN.gene2GO,
                    gene2GO = BP_annotation,
                    nodeSize=10)
    resultFisherBP <- runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
    allResBP <- GenTable(object=GOdataBP, pvalue = resultFisherBP,topNodes=geneData(resultFisherBP)[[4]])
    allResBP <- subset(allResBP,resultFisherBP@score[allResBP[,1]]<=0.01)
    
    GOdataMF <- new("topGOdata",
                    ontology = "MF",
                    allGenes = hits,
                    annot = annFUN.gene2GO,
                    gene2GO = MF_annotation,
                    nodeSize=10)
    resultFisherMF <- runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
    allResMF <- GenTable(object=GOdataMF, pvalue = resultFisherMF,topNodes=geneData(resultFisherMF)[[4]])
    allResMF <- subset(allResMF,resultFisherMF@score[allResMF[,1]]<=0.01)
    
    GOdataCC <- new("topGOdata",
                    ontology = "CC",
                    allGenes = hits,
                    annot = annFUN.gene2GO,
                    gene2GO = CC_annotation,
                    nodeSize=10)
    resultFisherCC <- runTest(GOdataCC, algorithm = "weight01", statistic = "fisher")
    allResCC <- GenTable(object=GOdataCC, pvalue = resultFisherCC,topNodes=geneData(resultFisherCC)[[4]])
    allResCC <- subset(allResCC,resultFisherCC@score[allResCC[,1]]<=0.01)
    
    out <- c(list(allResBP,allResMF,allResCC))
    return(out)
  })
  names(byLevel) <- c("e","pt","p","phr","ph")
  return(byLevel)
})
names(GOanalysis) <- hsOvMat[,"name"]
save(hsOvMat,targetsByHS,GOanalysis,file="data/hsInfo.RData")

#find lead markers for each HS
load("data/binInfo.RData")
load("data/hsInfo.RData")
load("data/pQTL_results170815.RData")
load("data/eQTL_results160831.RData")
load("data/expressionLevel.RData")
load("data/proteinLevel.RData")
load("data/phosphoNew/phosphoProt.RData")
load("data/eQtlMappingData160817.RData")
meta <- read.table("metadata/metadata.tsv",sep="\t",as.is=T,header=T)
rownames(meta) <- meta[,1]
meta$strain <- sapply(meta$strain,FUN=function(s){
  if(grepl(pattern = ".",x = s,fixed = T)){
    s <- paste0("X",s)
  }
  return(s)
})
strains2use <- intersect(meta[colnames(eBatchGeneLengthCorrected),"strain"],intersect(meta[colnames(proteinLevelBatchCorrected),"strain"],meta[colnames(phosphoProtResiduals),"strain"]))
eByStrain <- sapply(strains2use,FUN=function(s){
  rowMeans(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]==s,drop=F],na.rm=T)
})
pByStrain <- sapply(strains2use,FUN=function(s){
  rowMeans(proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]==s,drop=F],na.rm=T)
})
eEffects <- sapply(1:ncol(genotype),FUN=function(i){
  rowMeans(eByStrain[,genotype[,i]==0,drop=F],na.rm=T)-rowMeans(eByStrain[,genotype[,i]==1,drop=F],na.rm=T)
})
pEffects <- sapply(1:ncol(genotype),FUN=function(i){
  rowMeans(pByStrain[,genotype[,i]==0,drop=F],na.rm=T)-rowMeans(pByStrain[,genotype[,i]==1,drop=F],na.rm=T)
})
pTargetsByMarker <- colSums(pQTL_results$ptQTL$qv<=0.1)
pLead <- apply(hsOvMat,1,FUN=function(x){
  binVec <- as.numeric(x["startBin"]):as.numeric(x["endBin"])
  hsMarkers <- which(binPerMarker%in%binVec)
  pLeadMarker <- hsMarkers[which.max(pTargetsByMarker[hsMarkers])]
  return(pLeadMarker)
})
corByPHs <- sapply(which(as.logical(hsOvMat[,"pQtlHotspot"])),FUN=function(i){
  targets <- intersect(targetsByHS[[i]]$p,rownames(eByStrain))
  cor(eEffects[targets,pLead[i]],pEffects[targets,pLead[i]],use="pair")
})
mean(corByPHs)

###Venn diagrams for IRA2 and HAP1
library(eulerr)
library(RColorBrewer)
source("lib/general_function.R")
library(grid)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]

names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

load("data/phosphoLevel.RData")
rownames(phospho2prot) <- phospho2prot[,1]
phrProts <- unique(phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv),2])

#HAP1
i <- 14

fit = euler(list("pQTL" = intersect(targetsByHS[[i]]$p,phrProts),
                 "eQTL" = intersect(targetsByHS[[i]]$e,phrProts), 
                 "phQTL" = intersect(targetsByHS[[i]]$ph,phrProts)))
colvec <- colors[c(3,1,5)]
eg = plot(fit, fills = colvec, 
          alpha = 0.5, edges = F,
          quantities = list(font = 4, cex = 2),
          labels = list(font = 4, cex = 3))
eg
pdfAndPng(file = "graph/Venn_HAP1",4,4, expression({plot(eg)}))

#HAP1 with e pt and phres for paper
asRGB <- function(vec){
  rgb(vec[1],vec[2],vec[3])
}
colorMat <- col2rgb(colors)/255
fit = euler(list("ptQTL" = intersect(targetsByHS[[i]]$pt,phrProts),
                 "eQTL" = intersect(targetsByHS[[i]]$e,phrProts), 
                 "phResQTL" = intersect(targetsByHS[[i]]$phr,phrProts)))
#colvec <- c(colors[c(3,1,5)],rep("grey",20))
colvec <- c(asRGB(rowMeans(colorMat[,c(1,2)])), #e pt
            asRGB(rowMeans(colorMat[,c(4,2)])), #pt phr
            asRGB(rowMeans(colorMat[,c(4,1)])), #e phr
            asRGB(rowMeans(colorMat[,c(4,1,2)])), #e pt phr
            asRGB(colorMat[,c(2)]), #pt
            asRGB(colorMat[,c(1)]), #e
            asRGB(colorMat[,c(4)])) #pt
            

eg = plot(fit, fills = colors[c(2,1,4)], 
          alpha = 0.5, edges = F,
          quantities = list(col = colvec, font = 2, cex = 2),
          labels = list(col = colors[c(2,1,4)], font = 2, cex = 2))
#eg$gp <- gpar(mar=c(10,0,0,0))
eg$vp <- viewport(w = .7, h = .9, gp = gpar(col="blue"))
eg
pdfAndPng(file = "graph/Venn_HAP1_2",4,4, expression({plot(eg)}))
#IRA2
i <- 15

fit = euler(list("pQTL" = intersect(targetsByHS[[i]]$p,phrProts),
                 "eQTL" = intersect(targetsByHS[[i]]$e,phrProts), 
                 "phQTL" = intersect(targetsByHS[[i]]$ph,phrProts)))
eg = plot(fit, fills = colvec, 
          alpha = 0.5, edges = F,
          quantities = list(font = 4, cex = 2),
          labels = list(font = 4, cex = 3))
eg
pdfAndPng(file = "graph/Venn_IRA2",4,4, expression({plot(eg)}))