# ------- Library and constant ---------------------------------------
source("lib/general_function.R")
source("lib/genotyping_function.R")
library(dplyr)
chrLength <- c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
names(chrLength) <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

################################
### RNAseq based genotyping ####
################################

vcfPATH = "RNAseq/vcf/"
#vcfPATH = "/Volumes/TOSHIBA EXT/temp/vcf"

library(xlsx)
meta    <- read.xlsx("metadata/metadata.xlsx", sheetName="Main", row.names=1, header=TRUE,as.data.frame=TRUE, stringsAsFactors=F,na.strings = "NA")
meta$id <- rownames(meta) 
meta[meta=="NA"]<-NA

# single sample vcf
id<-lapply(meta$id[!is.na(meta$RNAseqRun) & meta$matching_with_brem_genotype], function(s){
  vcf <- read.vcf(paste(vcfPATH,"/",s,".vcf",sep=""))
  paste(vcf[,1],vcf[,2],sep="")
})

id<-Reduce(intersect, id)
#id<-unique(unlist(id))

vcfPerSample <- sapply(meta$id[!is.na(meta$RNAseqRun)], function(s){
  vcf <- read.vcf(paste(vcfPATH,"/",s,".vcf",sep=""))
  out<-vcf[,10]
  names(out) <- paste(vcf[,1],vcf[,2],sep="")
  out[id]
})

vcf <- read.vcf("RNAseq/unfiltered.vcf")
rownames(vcf)<-paste(vcf[,1],vcf[,2],sep="")
vcfPerSample<-cbind(vcf[id,1:9],vcfPerSample)

# -------------- check that genotype is always the same in all replicates ------
geno <- vcfPerSample[,10:ncol(vcfPerSample)]
geno[geno=="."] <-NA
geno[geno=="./."] <-NA
GT <- apply(geno,2,function(x){  # GT field of the genotypes
  as.numeric(sapply(as.character(x),function(x)strsplit(x,':',fixed=T)[[1]][1]))
})
genotypePerSample <- cbind(vcfPerSample[,1:4],GT)

corGT<-cor(GT,use="pair")
pdf("graph/hclust.pdf",5,15)
d<-hclust(as.dist(1-abs(corGT)))
par(cex=.5)
plot(as.dendrogram(d),horiz = T)
dev.off()


load("data/gdata.Rdata")
gdata$chromosome <- paste("chr",as.character(as.roman(gdata$chromosome)),sep="")
gdata <- gdata[,-c(1,2)]
colnames(gdata)[colnames(gdata)=="BY"] <- "BY4716"
colnames(gdata)[colnames(gdata)=="RM"] <- "RM11-1a"
gdata <- cbind(gdata[,1:3], useless=0, gdata[4:ncol(gdata)])

M<-meta[!is.na(meta$RNAseqRun),c("id","strain")]


pdf("graph/checkReplicate.pdf",width = 10,height = 4)
sapply(unique(M$strain), function(strain){
  par(mfrow=c(2,1),mar=c(0,5,3,2))
  s<-M$id[M$strain==strain]
  plot_genotype(genotypePerSample[,c("CHROM","POS","ID","REF",s)],chr_length = chrLength,p1 = 0, p2=1,axis=F)
  title(main=strain)
  par(mar=c(8,5,0,2))
  if(length(which(colnames(gdata)==strain))!=0){
    plot_genotype(labels = "brem",gdata[,c(1:4,which(colnames(gdata)==strain))],chr_length = chrLength,p1 = 1, p2=0,)
  } else {plot.new()}
})
dev.off()


png("graph/checkReplicate.png",width = 1200,height = 20000)
par(mfcol=c(124,2))
sapply(unique(M$strain), function(strain){
  par(mar=c(0,5,3,2))
  s<-M$id[M$strain==strain]
  plot_genotype(genotypePerSample[,c("CHROM","POS","ID","REF",s)],chr_length = chrLength,p1 = 0, p2=1,axis=F)
  title(main=strain)
  par(mar=c(8,5,0,2))
  if(length(which(colnames(gdata)==strain))!=0){
    plot_genotype(labels = "brem",gdata[,c(1:4,which(colnames(gdata)==strain))],chr_length = chrLength,p1 = 1, p2=0,)
  } else {plot.new()}
})
dev.off()

png("graph/problematicSample.png",width = 1200,height = 2000)
par(mfcol=c(16,2))
sapply(c("4.2.a","9.7.d","12.2.b","14.4.a","16.1.d","18.6.d","20.4.c","1.2.d","2.7.d","4.4.d","8.7.b","11.2.d","9.6.d","23.4.d"), function(strain){
  par(mar=c(0,5,3,2))
  s<-M$id[M$strain==strain]
  plot_genotype(genotypePerSample[,c("CHROM","POS","ID","REF",s)],chr_length = chrLength,p1 = 0, p2=1,axis=F)
  title(main=strain)
  par(mar=c(8,5,0,2))
  if(length(which(colnames(gdata)==strain))!=0){
    plot_genotype(labels = "brem",gdata[,c(1:4,which(colnames(gdata)==strain))],chr_length = chrLength,p1 = 1, p2=0,)
  } else {plot.new()}
})
dev.off()



# Infer the genotypes --------------------------------
# The script takes as input the VCF file for the GATK pipeline
vcf <- read.vcf("RNAseq/unfiltered.vcf")
rownames(vcf) <- paste(vcf[,1],vcf[,2],sep="")
vcfBloom<-read.vcf("genomes/RM_Bloom2013/RM3.vcf")
# only keep potstion where MQ in bloom > 30 as in the paper
m <- regexpr("MQ=([0-9]+)", vcfBloom$INFO,perl = T)
MQ<-as.numeric(gsub("MQ=","",regmatches(vcfBloom$INFO, m)))
vcfBloom<-vcfBloom[MQ>=30,]
#remove double annotaitons
temp<- paste(vcfBloom[,1],vcfBloom[,2],sep="")
sel<-table(temp)[temp]==1
vcfBloom<-vcfBloom[sel,]
rownames(vcfBloom)<-rownames(sel)[sel]

sel<-intersect(rownames(vcfBloom),rownames(vcf))
vcf<-vcf[sel,]
vcfBloom<-vcfBloom[sel,]


# Check genotypes raw calls globally to detect potential issues ---------
geno <- vcf[,10:ncol(vcf)]
geno[geno=="."] <-NA
geno[geno=="./."] <-NA
GT <- apply(geno,2,function(x){  # GT field of the genotypes
  as.numeric(sapply(as.character(x),function(x)strsplit(x,':',fixed=T)[[1]][1]))
})
table(GT)/46893/115*100
sum(is.na(GT))/46893/115*100
a<-cbind(vcf[,c("CHROM","POS","ID","REF")],GT)

chr.lim <-sapply(1:length(chrLength),function(chr){
  c(getGenomicPosition(chr,1,chromosome.size=chrLength),getGenomicPosition(chr,chrLength[chr],chromosome.size=chrLength))
})

pdfAndPng("graph/rawGenotypeAll",18,10,expression({
plot_genotype(a,chr_length = chrLength,p1 = 0, p2=1,axis=F,cex.y.axis=0.5, col=c('#377EB855','#E41A1C55'))
axis(1,at = colMeans(chr.lim),labels = as.roman(1:16), tick = F, axis = F)
}))

b<-cbind(vcf[,c("CHROM","POS","ID","REF")],GT[,c("18.1.d","2.5.d", "3.1.d")])
pdfAndPng("graph/rawGenotypeProblem",18,10,expression({
plot_genotype(b,chr_length = chrLength,p1 = 0, p2=1,axis=F)
axis(1,at = colMeans(chr.lim),labels = as.roman(1:16), tick = F)
}))

# romove strain with to many issue 18.1.d  2.5.d  3.1.d 
vcf <- vcf[,!colnames(vcf)%in%c("18.1.d","2.5.d", "3.1.d")]
# remove rows with more than 2 variants
# vcf <- vcf[-grep(",",vcf$ALT,fixed = T),] 

# --- filter the genotpye calls -----------------------------------------------

geno <- vcf[,10:ncol(vcf)]
geno[geno=="."] <-NA
geno[geno=="./."] <-NA

PL  <- apply(geno,2,function(x){ # DP field of the genotypes
  sapply(strsplit(x,':'),function(x){min(as.numeric(unlist(strsplit(x[2],","))))})
})

GT <- apply(geno,2,function(x){  # GT field of the genotypes
  as.numeric(sapply(as.character(x),function(x)strsplit(x,':',fixed=T)[[1]][1]))
})

GQ <- apply(geno,2,function(x){ # GQ field of the genotypes
  as.numeric(sapply(strsplit(x,':'),function(x){x[length(x)-1]}))
})

DP  <- apply(geno,2,function(x){ # DP field of the genotypes
  as.numeric(sapply(strsplit(x,':'),function(x){x[3]}))
})


which(PL>5,arr.ind=T)
sum(colMeans(PL,na.rm=T)>5)
which(colMeans(PL,na.rm=T)>5)


pdfAndPng("graph/genotyping_qual_overview", 11,5,expression({
  par(cex.main = 1.5, mar = c(5, 6, 3, 4) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las=1, mfrow=c(1,2))
  histGQ<-hist(GQ,plot = F)
  histGQ$counts <- log(histGQ$counts)
  plot(histGQ, main = "", col = "grey", axes=F, xlab='Genotype Quality (GQ)',ylab="Frequency\n")
  axis(2, log(c(1,10,100,1000,10000,100000,1000000,4000000)), c(1,10,100,1000,10000,100000,1000000,4000000))
  axis(1, seq(0,100,25))
  
  hist(DP,freq = F, ylab="Density\n", col="grey", xlab="depth at marker", main='')
}))

filteredGT <- GT
filteredGT[GQ<=60 | DP<=7] <- NA

#46893
# remove markers problematic in the parent: no different in RNAseq
#diffInParent <- GT[,"RM11-1a"] != GT[,"BY4716"]
#i.e: non 0 in BY, non 1 in RM
sel <- (filteredGT[,"RM11-1a"]==1 & filteredGT[,"BY4716"]==0)

# remove badly detected in parents
sel <- (filteredGT[,"RM11-1a"]==1 & filteredGT[,"BY4716"]==0)
sel[is.na(sel)]<-F
# 29896

filteredGT <- GT
filteredGT[GQ<=40 | DP<=4] <- NA
# remove to many NAs
sel<- sel & rowSums(is.na(filteredGT)) < (ncol(filteredGT)*.3) #25704
temp <- rowSums(filteredGT==1,na.rm=T)/(rowSums(filteredGT==1,na.rm=T)+rowSums(filteredGT==0,na.rm=T))
sel<-sel & temp<.8 & temp>.2 #25684
sel<- sel & rowSums(GT>1,na.rm = T)==0  #25683

filteredGT<-cbind(vcf[,c(1,2,4,5)],filteredGT)[sel,]



# detect suspect call (the one with call different from flamking markers) ----
a<-sapply(filteredGT[5:ncol(filteredGT)], find.suspect, filteredGT$CHR)
summary(sapply(a,length))
sum(is.na(filteredGT))/nrow(filteredGT)/112
plot(sapply(a,length))
length(unlist(a)) # 1008
sum(table(unlist(a))>1) #93
filteredGT <- filteredGT[-as.numeric(names(which(table(unlist(a))>1))),]

# infer the missing calls from identical flanking markers ---------
filledGenotype <- sapply(colnames(filteredGT)[5:ncol(filteredGT)], fillTheGap,geno=filteredGT, max_distance_flanking_markers=20000)

genotype <- cbind(filteredGT[,1:4],filledGenotype)
rownames(genotype) <- NULL

write.table(genotype, file="data/strain_genotype.tsv", sep="\t", quote=T,row.names = F)
# genotype 25,657 position

###############################
## ---------- group gentoype
###############################
genotype<-read.table("data/strain_genotype.tsv", check.names=F,header=T,sep='\t',stringsAsFactors = F)

#plot
chr.lim <-sapply(1:length(chrLength),function(chr){
  c(getGenomicPosition(chr,1,chromosome.size=chrLength),getGenomicPosition(chr,chrLength[chr],chromosome.size=chrLength))
})
pdfAndPng("graph/FilteredGenotype",18,10,expression({
  plot_genotype(genotype,chr_length = chrLength,p1 = 0, p2=1,axis=F,cex.y.axis=0.5, col=c('#377EB855','#E41A1C55'),lines=-1)
  axis(1,at = colMeans(chr.lim),labels = as.roman(1:16), tick = F)
  title(ylab = "strains")
}))

g<-genotype[,-c(1:4)]

allCor <- cor(t(g),use="pair")^2
d <- as.dist(1-allCor)
h <- hclust(d)
max(cutree(h,h = 0.0))

genotype2group <- cutree(h,h=0)

genotypeGroup <-   t(sapply(1:max(genotype2group), function(i) {
  (colMeans(g[genotype2group==i,],na.rm=T))
}))

# anotation of genotype group
chr<-sapply(1:max(genotype2group), function(i) {
  unique(genotype$CHROM[genotype2group==i])
})
b <- sapply(1:max(genotype2group), function(i) {
  min(genotype$POS[genotype2group==i])
})
e <- sapply(1:max(genotype2group), function(i) {
  max(genotype$POS[genotype2group==i])
})

genotypeGroup <- data.frame(CHR=chr, start=b, end=e,genotypeGroup,stringsAsFactors = F,check.names = F)


write.table(genotypeGroup, file="data/genotype_for_mapping.tsv", sep="\t", quote=T,row.names = F)

genotypeGroup<-read.table("data/genotype_for_mapping.tsv", check.names=F,header=T,sep='\t',stringsAsFactors = F)
sum(is.na(genotypeGroup)) #179
sum(rowSums(is.na(genotypeGroup))>0) #86


pdfAndPng("graph/mappingGenotype",10,10,expression({
  g<-data.frame(CHR="chr1", POS=1:nrow(genotypeGroup), REF=NA, ALT=NA, genotypeGroup[,-(1:3)],check.names = F,stringsAsFactors = F)
  cl<-nrow(genotypeGroup)
  names(cl)<-"chr1"
  chr.lim <-sapply(1:length(chrLength),function(chr){
    c(getGenomicPosition(chr,1,chromosome.size=chrLength),getGenomicPosition(chr,chrLength[chr],chromosome.size=chrLength))
  })
  plot_genotype(g,chr_length=cl,p1 = 0, p2=1,axis=F,cex.y.axis=0.5, col=c('#377EB855','#E41A1C55'))
}))


#############################################
# Population structure ------------------ ##
#############################################

genotype<-read.table("data/strain_genotype.tsv", check.names=F,header=T,sep='\t',stringsAsFactors = F)
genotypeGroup<-read.table("data/genotype_for_mapping.tsv", check.names=F,header=T,sep='\t',stringsAsFactors = F)

g<-(as.matrix(genotype[,-(1:4)]))
g<-g[rowSums(is.na(g))==0,]
gg<-(apply(g,1,function(x){(x-mean(x))/(mean(x)/2*sqrt(1-mean(x)/2))}))
K<-  gg%*%t(gg)
KPCA <- eigen(K)


pdfAndPng("graph/kinship",13,7, expression({ 
  library(RColorBrewer)
  layout(matrix(c(1,2,3),byrow=TRUE, nrow=1, ncol=3), widths=c(0.8,5,4))
  par(cex.lab=2)
  image.scale(K,col=colorRampPalette(c("white","darkblue"))(48), axis.pos=2,mar=c(5,2,4,1))
  image(K,col=colorRampPalette(c("white","darkblue"))(48),xaxt="n",yaxt="n", main="kinship")
  axis(1,labels=colnames(genotype)[-(1:4)],at=seq(0,1,1/(length(colnames(genotype[,-(1:4)]))-1)),las=2, cex.lab=0.2)
  axis(2,labels=colnames(genotype)[-(1:4)],at=seq(0,1,1/(length(colnames(genotype[,-(1:4)]))-1)),las=2, cex.lab=0.2)
  barplot(KPCA$values[1:10],main="PCA",ylab="eigen value",names.arg = 1:10, las=2)
}))


populationStructureCovariates <- KPCA$vectors[,1:7]
save(populationStructureCovariates, file="data/populationStructureCovariates.Rdata")



g<-(as.matrix(genotype[,-(1:4)]))
g<-g[rowSums(is.na(g))==0,]
gg<-(apply(g,1,function(x){(x-mean(x))/(mean(x)/2*sqrt(1-mean(x)/2))}))
K<-  gg%*%t(gg)
KPCA <- eigen(K)
g<-(as.matrix(genotype[,-(1:4)]))
g<-g[rowSums(is.na(g))==0,]
gg<-(apply(g,1,function(x){(x-mean(x))}))
K2<-  gg%*%t(gg)
KPCA2 <- eigen(K2)
pdfAndPng("graph/PCA_comparison", 10,10,expression({
  par(mfrow=c(3,3))
  sapply(1:9,function(i){
    plot(x=KPCA$vectors[,i],KPCA2$vectors[,i], 
         xlab="patterson et al.", ylab="centering", 
         main=paste("Eigen vector #", i, sep=""))
  })
}))
pdfAndPng("graph/PCA_comparison2", 4,4,expression({
  plot(x=KPCA$vectors[,5],KPCA2$vectors[,6], 
       xlab="patterson (vector 5)", ylab="centering (vector 6)")
}))


barplot(KPCA$values)

K<-  g%*%t(g)
v<-sum(diag(K))
cumsum(KPCA$values[1:30]/v) 


library(emma)
g<-(as.matrix(genotype[,-(1:4)]))
K<-emma.kinship(g,use="pairwise.complete.obs" )
K<-cov(g,use="pair")
KPCA <- eigen(K)
v<-sum(diag(K))
cumsum(KPCA$values[1:30]/v) 

pdfAndPng("graph/kinship",13,7, expression({ 
  library(RColorBrewer)
  layout(matrix(c(1,2,3),byrow=TRUE, nrow=1, ncol=3), widths=c(0.8,5,4))
  par(cex.lab=2)
  image.scale(K,col=colorRampPalette(c("white","darkblue"))(48), axis.pos=2,mar=c(5,2,4,1))
  image(K,col=colorRampPalette(c("white","darkblue"))(48),xaxt="n",yaxt="n", main="kinship")
  axis(1,labels=colnames(genotype)[-(1:4)],at=seq(0,1,1/(length(colnames(genotype[,-(1:4)]))-1)),las=2, cex.lab=0.2)
  axis(2,labels=colnames(genotype)[-(1:4)],at=seq(0,1,1/(length(colnames(genotype[,-(1:4)]))-1)),las=2, cex.lab=0.2)
  KPCA <- eigen(K)
  barplot(KPCA$values[1:10],main="PCA",ylab="eigen value",names.arg = 1:10, las=2)
}))

populationStructureCovariates <- KPCA$vectors[,1:10]

save(populationStructureCovariates, file="data/populationStructureCovariates.Rdata")

##########################################################
# --- generate data for strain specific genome  -----------
###########################################################

# for this we want to use markers that were determined at >50 
# fold coverage
# by bloom and al.

vcfBloom<-read.vcf("genomes/RM_Bloom2013/RM3.vcf")


#filter call at Q>30
m <- regexpr("MQ=([0-9]+)", vcfBloom$INFO,perl = T)
MQ<-as.numeric(gsub("MQ=","",regmatches(vcfBloom$INFO, m)))
sum(MQ>=30)
vcfBloom<-vcfBloom[MQ>=30,] # as in the paper
vcf <- read.vcf("RNAseq/unfiltered.vcf")
rownames(vcf) <- paste(vcf[,1],vcf[,2],sep="")
temp<- paste(vcfBloom[,1],vcfBloom[,2],sep="")

sel<-table(temp)[temp]==1
vcfBloom<-vcfBloom[sel,]
rownames(vcfBloom)<-rownames(sel)[sel]

sel<-intersect(rownames(vcfBloom),rownames(vcf))
vcf<-vcf[sel,]
vcfBloom<-vcfBloom[sel,]
vcf <- vcf[,!colnames(vcf)%in%c("18.1.d","2.5.d", "3.1.d")]

geno <- vcf[,10:ncol(vcf)]
geno[geno=="."] <-NA
geno[geno=="./."] <-NA
GT <- apply(geno,2,function(x){  # GT field of the genotypes
  as.numeric(sapply(as.character(x),function(x)strsplit(x,':',fixed=T)[[1]][1]))
})
sel <- GT[,"BY4716"]==1 | GT[,"RM11-1a"]==0
sel <- sel & !is.na(sel)

vcf<- vcf[-sel,]
vcfBloom<-vcfBloom[-sel,]

genotype<-read.table("data/strain_genotype.tsv", check.names=F,header=T,sep='\t',stringsAsFactors = F)
rownames(genotype) <-  paste(genotype[,1],genotype[,2],sep="")

geno <- vcf[,c(1,2,4,5,10:ncol(vcf))]
geno[,5:ncol(geno)]<-NA
geno[rownames(genotype),]<-genotype

filledGeno <- sapply(colnames(geno)[5:ncol(geno)], fillTheGap,geno=geno, max_distance_flanking_markers=20000000)
filledGeno[is.na(filledGeno)]<-0


library(xlsx)
meta    <- read.xlsx("metadata/metadata.xlsx", sheetName="Main", row.names=1, header=TRUE,as.data.frame=TRUE, stringsAsFactors=F,na.strings = "NA")
meta$id <- rownames(meta) 
meta[meta=="NA"]<-NA

library(xlsx)
meta    <- read.xlsx("metadata/metadata.xlsx", sheetName="Main", row.names=1, header=TRUE,as.data.frame=TRUE, stringsAsFactors=F,na.strings = "NA")
meta$id <- rownames(meta) 
meta[meta=="NA"]<-NA

sapply(colanmes(filledGeno),function())



# single sample vcf
s<-sapply(meta$id[!is.na(meta$RNAseqRun) & meta$matching_with_brem_genotype & !meta$suspect_multiple_strain], function(s){
  strain <- meta[s,"strain"]
  ssvcf <- vcfBloom[filledGeno[,strain]==1,]
  ssvcf<-ssvcf[!grepl(",",ssvcf$ALT,),]
  ssvcf[order(ssvcf$CHROM,ssvcf$POS),]
  f<-paste("RNAseq/generated_vcf/",s,".vcf",sep="")
  cat('##fileformat=VCFv4.0', 
      '##INFO=<ID=AB,Number=1,Type=Float,Description="Allele Balance for hets (ref/(ref+alt))">',
      '#############################################',
      '#####GENERTATED FOR INTERNAL USE ONLY #######',
      '#############################################',
      paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO',sep='\t'),
      sep='\n',
      file=f)
  write.table(ssvcf,file=f,append=TRUE,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  return(s)
})
cat(paste(s,".fa",sep=""),sep=" ")

