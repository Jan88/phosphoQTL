# ------- Library and constant ---------------------------------------
source("lib/general_function.R")
library(RColorBrewer)
library(xlsx)
library(dplyr)
chrLength <- c(230218,813184,316620,1531933,576874,270161,1090940,
               562643,439888,745751,666816,1078177,924431,784333,
               1091291,948066)
names(chrLength) <- c("chrI","chrII","chrIII","chrIV","chrV",
                      "chrVI","chrVII","chrVIII","chrIX","chrX",
                      "chrXI","chrXII","chrXIII","chrXIV","chrXV",
                      "chrXVI")


meta    <- read.xlsx("metadata/metadata.xlsx", sheetName="Main", 
                     row.names=1, header=TRUE,as.data.frame=TRUE,
                     stringsAsFactors=F,na.strings = "NA")
meta$id <- rownames(meta) 
meta[meta=="NA"]<-NA

meta <- meta[!is.na(meta$RNAseqRun) & meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]

# ----- Read the raw counts ----------

# sense
rawSenseCounts <- sapply(meta$id, function(strain) {
  rc <- read.table(paste("RNAseq/bam_ssg/",strain,"_ReadsPerGene.out.tab", sep=""),
                   stringsAsFactors = F, 
                   header = F, 
                   sep="\t")
  rc <- rc[!(grepl("N_",rc[,1])|grepl("Q[0..9]+",rc[,1])),] # remove first lines of STAR output and QNNNN annotations
  rc <- rc[order(rc[,1]),]
  counts <- rc[,4]   
  names(counts) <- rc[,1]
  return(counts)
})

write.table(x=rawSenseCounts, file="data/rawRNAseqCount.tsv", sep="\t", quote = F,row.names = T)

# anti-sense
rawASenseCounts <- sapply(meta$id, function(strain) {
  rc <- read.table(paste("RNAseq/bam_ssg/",strain,"_ReadsPerGene.out.tab", sep=""),
                   stringsAsFactors = F, 
                   header = F, 
                   sep="\t")
  rc <- rc[!(grepl("N_",rc[,1])|grepl("Q[0..9]+",rc[,1])),] # remove first lines of STAR output and QNNNN annotations
  rc <- rc[order(rc[,1]),]
  counts <- rc[,3]   
  names(counts) <- rc[,1]
  return(counts)
})


#-------- Normalisation ------------
lowExpressionCount<-apply(rawSenseCounts,1,function(x)sum(x<10)/length(x)) 
sum(lowExpressionCount<.1)  # 5429/6664
filteredSenseCounts <- rawSenseCounts[lowExpressionCount<0.1,]

#estimate scaling vs DESeq2
library("DESeq2")
dds<- DESeqDataSetFromMatrix(countData=filteredSenseCounts,
                             colData=as.data.frame(meta),
                             design=~1)

normDataDESeq2<-assay(rlog(dds,blind = T))

tmm<-apply(log2(filteredSenseCounts+1),2,mean,trim=.15)
normDataTMM <- sapply(1:length(tmm),function(i)log2(filteredSenseCounts+1)[,i]-tmm[i])

pdfAndPng("graph/RNA_normalization_comparison",8,10,expression({
layout(matrix(c(1,2,3,3,4,4),byrow=TRUE, nrow=3, ncol=2))
boxplot(normDataDESeq2,main="DESeq2",outline = F)
boxplot(normDataTMM, main="TMM",outline = F)
par(mar=c(2,0,3,0))
colBatchCul <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(meta$cultureBatch)))[as.numeric(as.factor(meta$cultureBatch))]
temp<- normDataDESeq2
colnames(temp) <- meta$strain
PlotHclustWithColClasses(temp, rbind(colBatchCul),main="DESeq2")
temp<- normDataTMM
colnames(temp) <- meta$strain
PlotHclustWithColClasses(temp, rbind(colBatchCul),main="TMM")
}))

eDataNorm <- normDataDESeq2
colnames(eDataNorm) <- colnames(filteredSenseCounts)
save(eDataNorm, file="data/eDataNorm.RData")

#----- batch structure
load("data/eDataNorm.RData")

pvalue<-apply(normDataDESeq2, 1, function(x){
  summary(aov(x~meta$cultureBatch))[[1]][["Pr(>F)"]][1]
})
pvalue2<-apply(normDataDESeq2, 1, function(x){
  summary(aov(sample(x)~meta$cultureBatch))[[1]][["Pr(>F)"]][1]
})

pvalueAmp<-apply(normDataDESeq2, 1, function(x){
  summary(aov(x~meta$riboAmpBatch))[[1]][["Pr(>F)"]][1]
})
hist(pvalue)

library(sva)
mod <- model.matrix(~1, data=meta)
batch <- meta$cultureBatch
eComBat = ComBat(dat=eDataNorm, batch=batch, mod=mod, par.prior=F, prior.plots=F)


pdfAndPng("graph/RNA_normalization_combat",9,7,expression({
  par(mfrow=c(2,1), mar=c(4,0,4,0))
  temp<- eDataNorm
  colnames(temp) <- meta$strain
  PlotHclustWithColClasses(temp, rbind(colBatchCul),main="DESseq")
  temp<- eComBat
  colnames(temp) <- meta$strain
  PlotHclustWithColClasses(temp, rbind(colBatchCul),main="DESeq + Combat")
}))


eLevelNonBatch <- eDataNorm
eLevelBatchCorrected <- eComBat
dimnames(eLevelBatchCorrected) <- dimnames(eLevelNonBatch)


# correct for gene length

gff <- read.table("genomes/sacCer3/sacCer3.gtf",stringsAsFactors = F, sep="\t")
gff<- gff[gff[,3]=="exon",]

geneLength<-sapply(rownames(eLevelBatchCorrected),function(g){
  sel<-grepl(g,gff[,9])
  sum(gff[sel,][,5]-gff[sel,][,4])  
})

eBatchGeneLengthCorrected <- t(sapply(1:nrow(eLevelBatchCorrected), function(i){
  log2(2^(eLevelBatchCorrected[i,])*1000/geneLength[i])
}))

dimnames(eBatchGeneLengthCorrected) <- dimnames(eLevelBatchCorrected) 

save(eLevelBatchCorrected,eBatchGeneLengthCorrected,eLevelNonBatch, file="data/expressionLevel.RData")

write.table(x=eLevelBatchCorrected, file="data/normalizedGeneExpression.tsv", sep="\t", quote = F,row.names = T)


library(sva)
mod <- model.matrix(~1, data=meta)
batch <- meta$cultureBatch
eComBat = ComBat(dat=eDataNorm, batch=batch, mod=mod, par.prior=T, prior.plots=F)


mod <- model.matrix(~1, data=meta)
batch <- meta$riboAmpBatch
eComBat = ComBat(dat=eComBat, batch=batch, mod=mod, par.prior=F, prior.plots=F)

eLevelNonBatch <- eDataNorm
eLevelBatchCorrected <- eComBat
dimnames(eLevelBatchCorrected) <- dimnames(eLevelNonBatch)

save(eLevelBatchCorrected,eBatchGeneLengthCorrected,eLevelNonBatch, file="data/expressionLevel_test.RData")
###########


#-------- Compute coverage######
countMatrix <- read.table(file = "data/rawRNAseqCount.tsv",
                          header = T,sep = "\t")
gff <- read.table(file = "Saccharomyces_cerevisiae/RM/genes_RM_CDSonly_refined.gtf",
                  stringsAsFactors = F)
transcriptLengths <- sapply(rownames(countMatrix),FUN=function(g){
  subgff <- gff[gff[,9]==g,,drop=F]
  if(nrow(subgff)==0){
    return(NA)
  }
  return(sum(apply(subgff[,4:5,drop=F],1,FUN=function(x){abs(x[1]-x[2])})))
})
any(is.na(transcriptLengths))
transcriptomeLength <- sum(transcriptLengths)
readLengthSumBySample <- colSums(countMatrix)*100
coverage = readLengthSumBySample/transcriptomeLength
range(coverage)


nReads = t(sapply(colnames(countMatrix), function(x){
  input = read.table(paste0("RNAseq/bowtieOutput/",x,
                            "/align_summary.txt"),
                     skip=1,nrows=1, as.is = T)[,3]
  mapped = read.table(paste0("RNAseq/bowtieOutput/",x,
                             "/align_summary.txt"),
                      skip=2,nrows=1, as.is = T)[,3]
  c(input = input, mapped = mapped, ratio = mapped/input)
}))
nReads = cbind(meta[rownames(nReads),], 
               nReads, coverage = coverage[rownames(nReads)])
cult = do.call(rbind,strsplit(nReads$platePosition, "-"))
nReads$plate = cult[,1]
nReads$plateRow = substr(cult[,2],1,1)
nReads$plateCol = substr(cult[,2],2,10)
nReads$cultureBatch = substr(nReads$cultureBatch,13,15)


png("graph/nReads.png", width = 1200, height = 400, pointsize=15)
par(oma = c(0,0,0,0), mar = c(5,4,1,1))
barplot(rbind(nReads[,"mapped"],nReads[,"input"]-nReads[,"mapped"]), 
        las = 2, col = c("orange3", "grey50"), ylab = "reads", 
        names.arg=rownames(nReads), 
        legend.text=c("mapped", "unmapped"), 
        args.legend=list(bty = "n", inset=0), cex.names=0.5)
dev.off()

png("graph/coverageByFeatures.png", 
    width = 1000, height = 1000, pointsize=20)
{
  featureY = nReads$coverage
  featureYname = "coverage"
  par(mfrow = c(3,2), mar = c(5,7,1,1))
  for(i in c("strain", "OD", "cultureBatch", 
             "plate", "plateRow", "plateCol")){
    if(is.numeric(nReads[,i])){
      plot(featureY ~ nReads[,i], xlab = "", ylab = "", las = 1)
    } else{
      boxplot(featureY ~ nReads[,i],
              xlab = "", ylab = "", 
              las = 2)
    }
    mtext(1,text=i, line = 4)
    mtext(2,text = featureYname, line = 5)
    p = summary(aov(featureY ~ nReads[,i]))[[1]]$`Pr(>F)`[1]
    legend("bottomright", legend=paste0("p-value = ",format(p,digits=3)))
  }
}
dev.off()

# plot n reads for different features
png("graph/nReadsByFeatures.png", 
    width = 1000, height = 1000, pointsize=20)
{
  featureY = nReads$input
  featureYname = "n Reads"
  par(mfrow = c(3,2), mar = c(5,7,1,1))
  for(i in c("strain", "OD", "cultureBatch", 
             "plate", "plateRow", "plateCol")){
    if(is.numeric(nReads[,i])){
      plot(featureY ~ nReads[,i], xlab = "", ylab = "", las = 1)
    } else{
      boxplot(featureY ~ nReads[,i],
              xlab = "", ylab = "", 
              las = 2)
    }
    mtext(1,text=i, line = 4)
    mtext(2,text = featureYname, line = 5)
    p = summary(aov(featureY ~ nReads[,i]))[[1]]$`Pr(>F)`[1]
    legend("bottomright", legend=paste0("p-value = ",format(p,digits=3)))
  }
}
dev.off()

# plot %mapped for different features
png("graph/MappedByFeatures.png", 
    width = 1000, height = 1000, pointsize=20)
{
  featureY = nReads$ratio
  featureYname = "% mapped"
  par(mfrow = c(3,2), mar = c(5,7,1,1))
  for(i in c("strain", "OD", "cultureBatch", 
             "plate", "plateRow", "plateCol")){
    if(is.numeric(nReads[,i])){
      plot(featureY ~ nReads[,i], xlab = "", ylab = "", las = 1)
    } else{
      boxplot(featureY ~ nReads[,i],
              xlab = "", ylab = "", 
              las = 2)
    }
    mtext(1,text=i, line = 4)
    mtext(2,text = featureYname, line = 5)
    p = summary(aov(featureY ~ nReads[,i]))[[1]]$`Pr(>F)`[1]
    legend("bottomright", legend=paste0("p-value = ",format(p,digits=3)))
  }
}
dev.off()

png("graph/nReadsByRiboAmpBatch.png")
par(mfrow = c(2,1))
boxplot(nReads$input ~ nReads[,"riboAmpBatch"],
        xlab = "riboAmpBatch", ylab = "n Reads", 
        las = 1, names = 1:15)
p = summary(aov(nReads$input ~ nReads[,"riboAmpBatch"]))[[1]]$`Pr(>F)`[1]
legend("bottomright", legend=paste0("p-value = ",format(p,digits=3)))
boxplot(nReads$ratio ~ nReads[,"riboAmpBatch"],
        xlab = "riboAmpBatch", ylab = "%mapped", 
        las = 1, names = 1:15)
p = summary(aov(nReads$ratio ~ nReads[,"riboAmpBatch"]))[[1]]$`Pr(>F)`[1]
legend("bottomright", legend=paste0("p-value = ",format(p,digits=3)))
dev.off()

