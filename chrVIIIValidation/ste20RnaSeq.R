library(sva)
library(DESeq2)
load("data/hsInfo.RData")
load("data/binInfo.RData")
hsMarkers <- which(binPerMarker%in%(119:120))

countMat <- lapply(1:15,FUN=function(i){
  numAsChar <- as.character(i)
  fileName <- paste0("ste20Validation/mappingAgainstSpec/counts/",paste(rep("0",3-nchar(numAsChar)),collapse=""),numAsChar,".txt")
  mat <- read.table(fileName,as.is=T,header=F)
  return(mat)
})
allGenes <- unique(unlist(lapply(countMat,FUN=function(mat){mat[,2]})))
countMat <- sapply(countMat,FUN=function(mat){
  countVec <- rep(0,length(allGenes))
  names(countVec) <- allGenes
  countVec[mat[,2]] <- mat[,1]
  return(countVec)
})
strainLabels <- c("s1473",
                  "s1473",
                  "s1475",
                  "s1475",
                  "s1475",
                  "s1476",
                  "s1476",
                  "s1476",
                  "s1477",
                  "s1477",
                  "s1477",
                  "s1473",
                  "s1475",
                  "s1476",
                  "s1477")

batch <- c(1,4,1,3,4,1,3,4,1,3,4,2,2,2,2)

#Confirm genotypes
ste20Loc <- c(95118,97937)
ste20Len <- ste20Loc[2]-ste20Loc[1]+1
gpa1Loc <- c(113499,114917)
gpa1Len <- gpa1Loc[2]-gpa1Loc[1]+1
voc <- c(1,2,3,4)
names(voc) <- c("A","C","G","T")
baseModusSte20 <- sapply(1:15,FUN=function(i){
  mat <- read.table(paste0("ste20Validation/mappingAgainstSpec/mappedReads/",paste(rep(0,3-nchar(i)),collapse=""),i,".bam.ste20Reads.txt"),sep="\t",as.is=T)
  mat <- mat[order(mat[,1]),]
  out <- matrix(0,ncol=4,nrow=ste20Loc[2]-ste20Loc[1]+1)
  for(j in 1:nrow(mat)){
    rl <- nchar(mat[j,2])
    sP <- as.numeric(mat[j,1])-ste20Loc[1]+1
    r <- strsplit(mat[j,2],split = "")[[1]]
    targetPos <- sP:(sP+rl-1)
    use <- targetPos>0&targetPos<=ste20Len&r%in%names(voc)
    out[cbind(targetPos[use],voc[r[use]])] <- out[cbind(targetPos[use],voc[r[use]])]+1
  }
  maxProp <- apply(out,1,max)/rowSums(out)
  baseCall <- rep(NA,nrow(out))
  baseCall[maxProp>=0.9] <- names(voc)[apply(out,1,which.max)][maxProp>=0.9]
  return(baseCall)
})
varBasesSte20 <- which(apply(baseModusSte20,1,FUN=function(x){max(table(x))+sum(is.na(x))})<ncol(baseModusSte20))
varCallsSte20 <- apply(baseModusSte20[varBasesSte20,],2,paste,collapse="")
baseModusGpa1 <- sapply(1:15,FUN=function(i){
  mat <- read.table(paste0("ste20Validation/mappingAgainstSpec/mappedReads/",paste(rep(0,3-nchar(i)),collapse=""),i,".bam.gpa1Reads.txt"),sep="\t",as.is=T)
  mat <- mat[order(mat[,1]),]
  out <- matrix(0,ncol=4,nrow=gpa1Loc[2]-gpa1Loc[1]+1)
  for(j in 1:nrow(mat)){
    rl <- nchar(mat[j,2])
    sP <- as.numeric(mat[j,1])-gpa1Loc[1]+1
    r <- strsplit(mat[j,2],split = "")[[1]]
    targetPos <- sP:(sP+rl-1)
    use <- targetPos>0&targetPos<=gpa1Len&r%in%names(voc)
    out[cbind(targetPos[use],voc[r[use]])] <- out[cbind(targetPos[use],voc[r[use]])]+1
  }
  maxProp <- apply(out,1,max)/rowSums(out)
  baseCall <- rep(NA,nrow(out))
  baseCall[maxProp>=0.9] <- names(voc)[apply(out,1,which.max)][maxProp>=0.9]
  return(baseCall)
})
varBasesGpa1 <- which(apply(baseModusGpa1,1,FUN=function(x){max(table(x))+sum(is.na(x))})<ncol(baseModusGpa1))
varCallsGpa1 <- apply(baseModusGpa1[varBasesGpa1,],2,paste,collapse="")

cbind(as.numeric(as.factor(varCallsSte20)),as.numeric(as.factor(varCallsGpa1))) #1477 mutated back at the last variable base of gpa1

# randos <- sample(rownames(res_1473_1475)[res_1473_1477$pvalue>0.6],100)
# countMat2 <- countMat
# countMat2[randos,batch==4] <- countMat2[randos,batch==4]*3
####for all transcripts####
dds <- DESeqDataSetFromMatrix(countData = countMat,colData = data.frame(strain=strainLabels,batch=as.factor(batch)),design=~batch+strain)
dds <- DESeq(dds)
res_1473_1475 <- results(dds,contrast=c("strain","s1473","s1475"))
res_1473_1476 <- results(dds,contrast=c("strain","s1473","s1476"))
res_1473_1477 <- results(dds,contrast=c("strain","s1473","s1477"))
res_1476_1477 <- results(dds,contrast=c("strain","s1476","s1477"))
res_1475_1477 <- results(dds,contrast=c("strain","s1475","s1477"))
res_1475_1476 <- results(dds,contrast=c("strain","s1475","s1476"))
res_b <- results(dds,contrast=c("batch","1","2"))

save(res_1473_1475,res_1475_1477,res_1473_1476,res_1473_1477,res_1476_1477,res_1475_1476,file="data/diffExpressionReplacementStrains.RData")


expressionLevels <- counts(dds,normalized=T)
expressionLevels[expressionLevels==0] <- 0.1
expressionLevels <- log2(expressionLevels)
colnames(expressionLevels) <- c("s1473_1",
                                "s1473_4",
                                "s1475_1",
                                "s1475_3",
                                "s1475_4",
                                "s1476_1",
                                "s1476_3",
                                "s1476_4",
                                "s1477_1",
                                "s1477_3",
                                "s1477_4",
                                "s1473_2",
                                "s1475_2",
                                "s1476_2",
                                "s1477_2")
mod <- model.matrix(~strainLabels)
expressionLevels <- ComBat(dat=expressionLevels, batch=batch, par.prior=F, prior.plots=F,mod = mod)
save(expressionLevels,file="data/expressionLevelsReplacementStrains.RData")

####for the transcripts from eQTL mapping####
load("data/eQtlMappingData160817.RData")
dds <- DESeqDataSetFromMatrix(countData = countMat[rownames(countMat)%in%rownames(phenotype),],colData = data.frame(strain=strainLabels,batch=as.factor(batch)),design=~batch+strain)
dds <- DESeq(dds)
res_1473_1475 <- results(dds,contrast=c("strain","s1473","s1475"))
res_1473_1476 <- results(dds,contrast=c("strain","s1473","s1476"))
res_1473_1477 <- results(dds,contrast=c("strain","s1473","s1477"))
res_1476_1477 <- results(dds,contrast=c("strain","s1476","s1477"))
res_1475_1477 <- results(dds,contrast=c("strain","s1475","s1477"))
res_1475_1476 <- results(dds,contrast=c("strain","s1475","s1476"))
res_b <- results(dds,contrast=c("batch","1","2"))

save(res_1473_1475,res_1475_1477,res_1473_1476,res_1473_1477,res_1476_1477,res_1475_1476,file="data/diffExpressionReplacementStrainsSubset.RData")


expressionLevels <- counts(dds,normalized=T)
expressionLevels[expressionLevels==0] <- 0.1
expressionLevels <- log2(expressionLevels)
colnames(expressionLevels) <- c("s1473_1",
                                "s1473_4",
                                "s1475_1",
                                "s1475_3",
                                "s1475_4",
                                "s1476_1",
                                "s1476_3",
                                "s1476_4",
                                "s1477_1",
                                "s1477_3",
                                "s1477_4",
                                "s1473_2",
                                "s1475_2",
                                "s1476_2",
                                "s1477_2")
mod <- model.matrix(~strainLabels)
expressionLevels <- ComBat(dat=expressionLevels, batch=batch, par.prior=F, prior.plots=F,mod = mod)
save(expressionLevels,file="data/expressionLevelsReplacementStrainsSub.RData")
load("data/expressionLevel.RData")
load("data/eQtlMappingData160728.RData")
load("data/eQTL_results160831.RData")
meta <- read.table("metadata/metadata.tsv",sep="\t",as.is=T,header=T)
rownames(meta) <- meta[,1]
meta$strain <- sapply(meta$strain,FUN=function(s){
  if(grepl(pattern = ".",x = s,fixed = T)){
    s <- paste0("X",s)
  }
  return(s)
})

strains2use <- unique(meta[colnames(eBatchGeneLengthCorrected),"strain"])
eByStrain <- sapply(strains2use,FUN=function(s){
  rowMeans(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]==s,drop=F],na.rm=T)
})
eByStrain <- eByStrain[,rownames(genotype)]
hsEffect <- rowMeans(eByStrain[,genotype[,1597]==0],na.rm=T)-rowMeans(eByStrain[,genotype[,1597]==1],na.rm=T)

sharedGenes <- intersect(rownames(eByStrain),rownames(res_1473_1475))
sigGenes <- sharedGenes[qv[sharedGenes,1597]<0.2]
plot(res_1473_1477[sigGenes,"log2FoldChange"],hsEffect[sigGenes])
plot(res_1473_1475[sigGenes,"log2FoldChange"],hsEffect[sigGenes])
plot(res_1473_1476[sigGenes,"log2FoldChange"],hsEffect[sigGenes])
plot(res_1475_1477[sigGenes,"log2FoldChange"],hsEffect[sigGenes])
plot(res_1476_1477[sigGenes,"log2FoldChange"],hsEffect[sigGenes])
plot(res_1475_1476[sigGenes,"log2FoldChange"],hsEffect[sigGenes])
myCounts <- log2(counts(dds,normalized=T))
myCounts[is.infinite(myCounts)] <- NA
countsCen <- t(scale(t(myCounts),center=T,scale=F))

sigGenes <- sharedGenes[apply(qv[sharedGenes,hsMarkers],1,min)<0.1]
pdf(file = "graph/ste20FirstLook.pdf")
par(mfrow=c(2,2))
plot(hsEffect[sigGenes],res_1473_1475[sigGenes,"log2FoldChange"],xlab="log2FC(HS)",ylab="log2FC STE20-RM",pch=20)
plot(hsEffect[sigGenes],res_1473_1476[sigGenes,"log2FoldChange"],xlab="log2FC(HS)",ylab="log2FC GPA1-RM",pch=20)
plot(hsEffect[sigGenes],res_1473_1477[sigGenes,"log2FoldChange"],xlab="log2FC(HS)",ylab="log2FC GPA1-RM, STE20-RM",pch=20)
plot(hsEffect[sigGenes],res_1476_1477[sigGenes,"log2FoldChange"],xlab="log2FC(HS)",ylab="log2FC GPA1-RM (STE20-RM background)",pch=20)
dev.off()

#PCAs
library(ggfortify)
autoplot(prcomp(t(myCounts[rowSums(is.na(myCounts))==0,])),data=cbind(t(as.data.frame(myCounts[rowSums(is.na(myCounts))==0,])),strain=strainLabels),col="strain")
pc <- prcomp(t(myCounts[rowSums(is.na(myCounts))==0,]))
pcCors <- t(apply(myCounts,1,FUN=function(x){return(c(cor(x,pc$rotation[1,],use='pair'),cor(x,pc$rotation[2,],use='pair')))}))

sigGenes <- sharedGenes[qv[sharedGenes,1597]<0.3]
autoplot(prcomp(t(myCounts[sigGenes,])),data=cbind(t(as.data.frame(myCounts[sigGenes,])),strain=strainLabels),col="strain")
pc <- prcomp(t(myCounts[sigGenes,]))
pcCors <- t(apply(myCounts,1,FUN=function(x){return(c(cor(x,pc$rotation[1,],use='pair'),cor(x,pc$rotation[2,],use='pair')))}))

#What genes differ between replicates of the same strain?
sdBetweenReps <- sapply(unique(strainLabels),FUN=function(s){
  apply(myCounts[,strainLabels==s],1,sd,na.rm=T)
})
meanByStrain <- sapply(unique(strainLabels),FUN=function(s){
  apply(myCounts[,strainLabels==s],1,mean,na.rm=T)
})
centeredByStrain <- sapply(1:15,FUN=function(i){myCounts[,i]-meanByStrain[,strainLabels[i]]})
meanByBatch <- sapply(1:4,FUN=function(b){
  apply(myCounts[,batch==b&strainLabels!="s1473"],1,mean,na.rm=T)
})
sigGenes <- sharedGenes[apply(qv[sharedGenes,hsMarkers],1,min)<0.1]
countsNormCorr <- myCounts-meanByBatch[,batch]
autoplot(prcomp(t(countsNormCorr[rowSums(is.na(myCounts))==0,])),data=cbind(t(as.data.frame(countsNormCorr[rowSums(is.na(myCounts))==0,])),strain=strainLabels),col="strain")
autoplot(prcomp(t(countsNormCorr[sigGenes,])),data=cbind(t(as.data.frame(countsNormCorr[sigGenes,])),strain=strainLabels),col="strain")


autoplot(prcomp(t(centeredByStrain[rowSums(is.na(myCounts))==0,])),data=cbind(t(as.data.frame(centeredByStrain[rowSums(is.na(myCounts))==0,])),strain=strainLabels),col="strain")
autoplot(prcomp(t(centeredByStrain[sigGenes,])),data=cbind(t(as.data.frame(centeredByStrain[sigGenes,])),strain=strainLabels),col="strain")
pc <- prcomp(t(myCounts[rowSums(is.na(myCounts))==0,]))
pcCors <- t(apply(centeredByStrain[rowSums(is.na(myCounts))==0,],1,FUN=function(x){return(c(cor(x,pc$rotation[1,],use='pair'),cor(x,pc$rotation[2,],use='pair')))}))

#why are the pvalues for 1473 vs 1475 so high?
sigGenes <- sharedGenes[apply(qv[sharedGenes,hsMarkers],1,min)<0.1]
pdf("graph/sigGenesBoxplots.pdf")
sapply(sigGenes,FUN=function(g){
  boxplot(myCounts[g,]~strainLabels,main=g)
})
dev.off()

library(beeswarm)
sigGenes <- sharedGenes[apply(qv[sharedGenes,hsMarkers],1,min)<0.1]
pdf("graph/sigGenesBeeswarms.pdf")
par(mfrow=c(2,2))
sapply(sigGenes,FUN=function(g){
  beeswarm(countsNormCorr[g,]~strainLabels,main=g)
})
dev.off()

#model by allele
gpa1 <- as.factor(c(0,0,0,0,0,1,1,1,1,1,1,0,0,1,1))
ste20 <- as.factor(c(0,0,1,1,1,0,0,0,1,1,1,0,1,0,1))
dds <- DESeqDataSetFromMatrix(countData = countMat,colData = data.frame(gpa1=gpa1,ste20=ste20,batch=as.factor(batch)),design=~batch+gpa1+ste20)
dds <- DESeq(dds)
res_gpa1 <- results(dds,contrast=c("gpa1","0","1"))
res_ste20 <- results(dds,contrast=c("ste20","0","1"))
res_int <- results(dds,name="gpa11.ste201")

#compare each sample to its strain average
sampleStrainCor <- sapply(1:15,FUN=function(i){
  mean(apply(countsNormCorr[sigGenes,setdiff(which(strainLabels==strainLabels[i]),i)],2,cor,y=countsNormCorr[sigGenes,i],use="pair"))
})