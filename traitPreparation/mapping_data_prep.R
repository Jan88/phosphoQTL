###general info
setwd("/mnt/cellnet/phosphoQTL/data")
meta <- read.table("/mnt/cellnet/phosphoQTL/metadata/metadata.tsv",header=T,sep="\t")
meta$id <- rownames(meta)
meta[meta=="NA"]<-NA

dontuse <- c("RM11-1-1", "RM11-1b", "BY4724")
genotype4mapping <- read.table("genotype_for_mapping.tsv",header=T)
colnames(genotype4mapping)[115] <- "RM11-1a"
load("populationStructureCovariates.Rdata")
rownames(populationStructureCovariates) <- colnames(genotype4mapping)[-(1:3)]

###prepare eQTL mapping data

load("expressionLevel.RData")

#use eBatchGeneLengthCorrected
subMeta <- subset(meta,culture%in%colnames(eBatchGeneLengthCorrected))
rnaStrains <- unique(subMeta$strain)
rnaStrains <- setdiff(rnaStrains,dontuse)

rna <- sapply(rnaStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- eBatchGeneLengthCorrected[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})

save(rna, file="buddingYeastRNA.RData") #for Maria, 7.6.2017

#make a minimal genotype for the rnaStrains
rnaStrainsX <- sapply(rnaStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})

rnaGeno <- genotype4mapping[,rnaStrainsX]
rnaGeno <- t(rnaGeno)
# id <- sapply(1:ncol(rnaGeno),FUN=function(i){
#   sapply(1:ncol(rnaGeno),FUN=function(j){
#     identical(rnaGeno[,i],rnaGeno[,j])&i!=j
#   })
# }) #no changes necessary
colnames(rna) <- rnaStrainsX


genotype <- cbind(rnaGeno,populationStructureCovariates)
phenotype <- rna
phenotype <- t(scale(t(phenotype),scale = T,center = T))
save(genotype,phenotype,file="eQtlMappingData160817.RData")

###prepare pQTL mapping data
load("protPreprocessed170829/proteinLevel.RData")
subMeta <- subset(meta,culture%in%colnames(proteinLevelBatchCorrected))
protStrains <- unique(subMeta$strain)
protStrains <- setdiff(protStrains,dontuse)

prot <- sapply(protStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- proteinLevelBatchCorrected[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})

save(prot, file="buddingYeastProtein.RData") #for Maria, 7.6.2017

#there are strains for which we dont seem to have genotypes

protStrainsX <- sapply(protStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
colnames(prot) <- protStrainsX
protStrainsX <- intersect(protStrainsX,colnames(genotype4mapping))
prot <- prot[,protStrainsX]
all(protStrainsX%in%rnaStrainsX)
length(protStrainsX)==length(rnaStrainsX)
prot <- prot[,rnaStrainsX]
phenotype <- prot
phenotype <- t(scale(t(phenotype),scale = T,center = T))
phenotype4perm <- phenotype[!apply(is.na(phenotype),1,any),]
identical(colnames(phenotype),rownames(genotype))
save(genotype,phenotype,phenotype4perm,file="pQtlMappingData170831.RData")

###prepare phosphoLevel
load("protPreprocessed170829/phosphoLevel.RData")
subMeta <- subset(meta,culture%in%colnames(phosphoLevelBatchCorrected))
phosphoLevelStrains <- unique(subMeta$strain)
phosphoLevelStrains <- setdiff(phosphoLevelStrains,dontuse)
phosphoLevel <- sapply(phosphoLevelStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- phosphoLevelBatchCorrected[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})
phosphoLevelStrainsX <- sapply(phosphoLevelStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
colnames(phosphoLevel) <- phosphoLevelStrainsX
phosphoLevelStrainsX <- intersect(phosphoLevelStrainsX,colnames(genotype4mapping))
phosphoLevel <- phosphoLevel[,phosphoLevelStrainsX]
phenotype <- phosphoLevel
phenotype <- t(scale(t(phenotype),scale = T,center = T))
identical(colnames(phenotype),rownames(genotype))
phenotype <- phenotype[,rownames(genotype)]
phenotype4perm <- phenotype[!apply(is.na(phenotype),1,any),]
save(genotype,phenotype,phenotype4perm,file="phosphoLevelQtlMappingData170831.RData")

###prepare phosphoRatio
load("protPreprocessed170829/phosphoProt.RData")
subMeta <- subset(meta,culture%in%colnames(phosphoProtRatio))
phosphoRatioStrains <- unique(subMeta$strain)
phosphoRatioStrains <- setdiff(phosphoRatioStrains,dontuse)
phosphoRatio <- sapply(phosphoRatioStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- phosphoProtRatio[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})
phosphoRatioStrainsX <- sapply(phosphoRatioStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
colnames(phosphoRatio) <- phosphoRatioStrainsX
phosphoRatioStrainsX <- intersect(phosphoRatioStrainsX,colnames(genotype4mapping))
phosphoRatio <- phosphoRatio[,phosphoRatioStrainsX]
phenotype <- phosphoRatio
phenotype <- t(scale(t(phenotype),scale = T,center = T))
identical(colnames(phenotype),rownames(genotype))
phenotype <- phenotype[,rownames(genotype)]
phenotype4perm <- phenotype[!apply(is.na(phenotype),1,any),]
save(genotype,phenotype,phenotype4perm,file="phosphoRatioQtlMappingData170831.RData")

###prepare phosphoProtResiduals
subMeta <- subset(meta,culture%in%colnames(phosphoProtResiduals))
phosphoResidualsStrains <- unique(subMeta$strain)
phosphoResidualsStrains <- setdiff(phosphoResidualsStrains,dontuse)
phosphoResiduals <- sapply(phosphoResidualsStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- phosphoProtResiduals[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})
phosphoResidualsStrainsX <- sapply(phosphoResidualsStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
colnames(phosphoResiduals) <- phosphoResidualsStrainsX
phosphoResidualsStrainsX <- intersect(phosphoResidualsStrainsX,colnames(genotype4mapping))
phosphoResiduals <- phosphoResiduals[,phosphoResidualsStrainsX]
phenotype <- phosphoResiduals
phenotype <- t(scale(t(phenotype),scale = T,center = T))
identical(colnames(phenotype),rownames(genotype))
phenotype <- phenotype[,rownames(genotype)]
phenotype4perm <- phenotype[!apply(is.na(phenotype),1,any),]
save(genotype,phenotype,phenotype4perm,file="phosphoProtResidualsQtlMappingData170831.RData")

###prepare phosphoRnaResiduals
load("protPreprocessed170829/phosphoRna.RData")
subMeta <- subset(meta,culture%in%colnames(phosphoRnaResiduals))
phosphoResidualsStrains <- unique(subMeta$strain)
phosphoResidualsStrains <- setdiff(phosphoResidualsStrains,dontuse)
phosphoResiduals <- sapply(phosphoResidualsStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- phosphoRnaResiduals[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})
phosphoResidualsStrainsX <- sapply(phosphoResidualsStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
colnames(phosphoResiduals) <- phosphoResidualsStrainsX
phosphoResidualsStrainsX <- intersect(phosphoResidualsStrainsX,colnames(genotype4mapping))
phosphoResiduals <- phosphoResiduals[,phosphoResidualsStrainsX]
phenotype <- phosphoResiduals
phenotype <- t(scale(t(phenotype),scale = T,center = T))
identical(colnames(phenotype),rownames(genotype))
phenotype <- phenotype[,rownames(genotype)]
phenotype4perm <- phenotype[!apply(is.na(phenotype),1,any),]
save(genotype,phenotype,phenotype4perm,file="phosphoRnaResidualsQtlMappingData170831.RData")


###prepare pt traits
load("protPreprocessed170829/protRna.RData")
subMeta <- subset(meta,culture%in%colnames(protRnaResiduals))
ptStrains <- unique(subMeta$strain)
ptStrains <- setdiff(ptStrains,dontuse)
ptTraits <- sapply(ptStrains,FUN=function(s){
  samples <- as.character(subMeta[subMeta$strain==s,1])
  mat <- protRnaResiduals[,samples,drop=F]
  out <- apply(mat,1,mean,na.rm=T)
  return(out)
})
ptStrainsX <- sapply(ptStrains,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
colnames(ptTraits) <- ptStrainsX
ptStrainsX <- intersect(ptStrainsX,colnames(genotype4mapping))
ptTraits <- ptTraits[,ptStrainsX]
phenotype <- ptTraits
phenotype <- t(scale(t(phenotype),scale = T,center = T))
identical(colnames(phenotype),rownames(genotype))
phenotype <- phenotype[,rownames(genotype)]
phenotype4perm <- phenotype[!apply(is.na(phenotype),1,any),]
save(genotype,phenotype,phenotype4perm,file="ptQtlMappingData170831.RData")


###prepare bQTL mapping data

load("bufferingProteinLevel.RData")

#buffering trait computed with all available genes and the s3-protein data
s3samples <- intersect(colnames(bufferingSum3prot),colnames(eBatchGeneLengthCorrected))
s3prot <- bufferingSum3prot[,s3samples]
s3rna <- eBatchGeneLengthCorrected[rownames(s3prot),s3samples]
buffering_s3_ag <- sapply(1:ncol(s3prot),FUN=function(i){
  cor(s3prot[,i],s3rna[,i],use="pair",method="spearman")
})
names(buffering_s3_ag) <- colnames(s3prot)
#s3strains <- meta[meta[,1]%in%s3samples,2]
subMeta <- subset(meta,meta[,1]%in%names(buffering_s3_ag))
strainInfo <- subMeta$strain
names(strainInfo) <- subMeta$culture
strainInfo <- strainInfo[names(buffering_s3_ag)]

#library(beeswarm)
#beeswarm(buffering_s3_ag~strainInfo,pch=20)
boxplot(buffering_s3_ag~strainInfo,main="buffering trait based on s3-protein levels and 548 genes",ylab="spearman's rho between prot and rna",xlab="strain")

#abline(v=seq(1,length(unique(strainInfo)),1))

buffering_s3_ag <- sapply(unique(strainInfo),FUN=function(s){
  vec <- buffering_s3_ag[which(strainInfo==s)]
  mean(vec)
})
names(buffering_s3_ag) <- unique(strainInfo)
hist(buffering_s3_ag,main="distribution of the buffering trait based on s3-protein levels and 548 genes",xlab="spearman's rho between prot and rna")

s3_prot_buddingYeast <- sapply(names(buffering_s3_ag),FUN=function(s){
  rowMeans(s3prot[,strainInfo==s,drop=F],na.rm=T)
})

save(s3_prot_buddingYeast, file="buddingYeastS3Protein.RData") #for Maria, 7.6.2017

#buffering trait computed with the third of the genes with the highest MAD based on the s3 prot
rnaByStrain <- sapply(names(buffering_s3_ag),FUN=function(s){
  mat <- s3rna[,as.character(subMeta[subMeta[,2]==s,1]),drop=F]
  apply(mat,1,mean)
})
rnaMAD <- apply(s3rna,1,mad)
hm_genes <- names(rnaMAD)[rnaMAD>=quantile(rnaMAD,2/3)]

buffering_s3_hm <- sapply(1:ncol(s3prot),FUN=function(i){
  cor(s3prot[hm_genes,i],s3rna[hm_genes,i],use="pair",method="spearman")
})
names(buffering_s3_hm) <- colnames(s3prot)
buffering_s3_hm <- buffering_s3_hm[names(strainInfo)]
boxplot(buffering_s3_hm~strainInfo,main="buffering trait based on s3-protein levels and 183 high MAD genes",ylab="spearman's rho between prot and rna",xlab="strain")
buffering_s3_hm <- sapply(unique(strainInfo),FUN=function(s){
  vec <- buffering_s3_hm[which(strainInfo==s)]
  mean(vec)
})
names(buffering_s3_hm) <- unique(strainInfo)
hist(buffering_s3_hm,main="distribution of the buffering trait based on s3-protein levels and 183 high MAD genes",xlab="spearman's rho between prot and rna",cex.main=1)

plot(buffering_s3_hm,buffering_s3_ag,main="buffering traits based on s3 protein levels and different gene sets",xlab="correlation in high MAD genes",ylab="correlation in all genes (548)")
abline(a=0,b=1,col="red",lty=2)


#buffering trait computed with all available genes and the mapdia-protein data
MDsamples <- intersect(colnames(bufferingMapDia),colnames(eBatchGeneLengthCorrected))
MDprot <- bufferingMapDia[,MDsamples]
MDrna <- eBatchGeneLengthCorrected[rownames(MDprot),MDsamples]
buffering_MD_ag <- sapply(1:ncol(MDprot),FUN=function(i){
  cor(MDprot[,i],MDrna[,i],use="pair",method="spearman")
})
names(buffering_MD_ag) <- colnames(MDprot)
#MDstrains <- meta[meta[,1]%in%MDsamples,2]
subMeta <- subset(meta,meta[,1]%in%names(buffering_MD_ag))
strainInfo <- subMeta$strain
names(strainInfo) <- subMeta$culture
strainInfo <- strainInfo[names(buffering_MD_ag)]

boxplot(buffering_MD_ag~strainInfo,main="buffering trait based on MD-protein levels and 548 genes",ylab="spearman's rho between prot and rna",xlab="strain")

#abline(v=seq(1,length(unique(strainInfo)),1))

buffering_MD_ag <- sapply(unique(strainInfo),FUN=function(s){
  vec <- buffering_MD_ag[which(strainInfo==s)]
  mean(vec)
})
names(buffering_MD_ag) <- unique(strainInfo)
hist(buffering_MD_ag,main="distribution of the buffering trait based on MD-protein levels and 548 genes",xlab="spearman's rho between prot and rna")

#buffering trait computed with the third of the genes with the highest MAD based on the mapdia prot

buffering_MD_hm <- sapply(1:ncol(MDprot),FUN=function(i){
  cor(MDprot[hm_genes,i],MDrna[hm_genes,i],use="pair",method="spearman")
})
names(buffering_MD_hm) <- colnames(MDprot)
buffering_MD_hm <- buffering_MD_hm[names(strainInfo)]
boxplot(buffering_MD_hm~strainInfo,main="buffering trait based on MD-protein levels and 183 high MAD genes",ylab="spearman's rho between prot and rna",xlab="strain")
buffering_MD_hm <- sapply(unique(strainInfo),FUN=function(s){
  vec <- buffering_MD_hm[which(strainInfo==s)]
  mean(vec)
})
names(buffering_MD_hm) <- unique(strainInfo)
hist(buffering_MD_hm,main="distribution of the buffering trait based on MD-protein levels and 183 high MAD genes",xlab="spearman's rho between prot and rna",cex.main=1)

plot(buffering_MD_hm,buffering_MD_ag,main="buffering traits based on MD protein levels and different gene sets",xlab="correlation in high MAD genes",ylab="correlation in all genes (548)")
abline(a=0,b=1,col="red",lty=2)

#plot md vs s3
plot(buffering_MD_ag,buffering_s3_ag,xlab="MD-derived buffering trait",ylab="S3-derived buffering trait",main="buffering traits based on all available genes")
abline(a=0,b=1,col="red",lty=2)
plot(buffering_MD_hm,buffering_s3_hm,xlab="MD-derived buffering trait",ylab="S3-derived buffering trait",main="buffering traits based on high MAD genes")
abline(a=0,b=1,col="red",lty=2)

#how similar are the prot matrices?
identical(colnames(bufferingMapDia),colnames(bufferingSum3prot))
identical(rownames(bufferingMapDia),rownames(bufferingSum3prot))
corByGene <- sapply(1:nrow(bufferingSum3prot),FUN=function(i){
  cor(bufferingSum3prot[i,],bufferingMapDia[i,],use="pair",method="spearman")
})
corBySample <- sapply(1:ncol(bufferingSum3prot),FUN=function(i){
  cor(bufferingSum3prot[,i],bufferingMapDia[,i],use="pair",method="spearman")
})

all(apply(cbind(names(buffering_s3_ag),names(buffering_s3_hm),names(buffering_MD_ag),names(buffering_MD_hm)),1,FUN=function(x){length(unique(x))})==1) #all name vectors are identical
strainsX <- sapply(names(buffering_MD_hm),FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
names(buffering_s3_ag) <- names(buffering_s3_hm) <- names(buffering_MD_ag) <- names(buffering_MD_hm) <- strainsX
#save everything
s3_ag <- list(genotype=genotype[strainsX,],phenotype=as.vector(scale(buffering_s3_ag)))
s3_hm <- list(genotype=genotype[strainsX,],phenotype=as.vector(scale(buffering_s3_hm)))
MD_ag <- list(genotype=genotype[strainsX,],phenotype=as.vector(scale(buffering_MD_ag)))
MD_hm <- list(genotype=genotype[strainsX,],phenotype=as.vector(scale(buffering_MD_hm)))
traitPack <- list(s3_ag=s3_ag,s3_hm=s3_hm,MD_ag=MD_ag,MD_hm=MD_hm)
save(traitPack,file="bQtlMappingData160729.RData")
