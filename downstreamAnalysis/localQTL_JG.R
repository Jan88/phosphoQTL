library(seqinr)
library(gtools)
library(abind)
source("lib/genotyping_function.R")
source("lib/general_function.R")


# prot sequence
RMprotsFasta   <- "Saccharomyces_cerevisiae/RM/proteins_RM.fa"
BYProtsFasta   <- "Saccharomyces_cerevisiae/sacCer3/REFINED/proteins_S288C_R6411.fa"
RMProts   <- read.fasta(RMprotsFasta, as.string=TRUE, forceDNAtolower=FALSE)
BYProts   <- read.fasta(BYProtsFasta, as.string=TRUE, forceDNAtolower=FALSE)

load("data/expressionLevel.RData")
load(file = "data/phosphoProt.RData")
load("data/protRna.RData")
load("data/phosphoLevel.RData")
load("data/proteinLevel.RData")
load("data/pQTL_results170815.RData")
phosphoLevelQTL <- pQTL_results$phosphoLevelQTL
phosphoProtResidualsQTL <- pQTL_results$phosphoProtResidualsQTL
phosphoRnaResidualsQTL <- pQTL_results$phosphoRnaResidualsQTL
ptQTL <- pQTL_results$ptQTL
pQTL <- pQTL_results$pQTL
load("data/eQTL_results160831.RData")
eQTL <- list(pv, qv, QtlList)
names(eQTL) <- names(ptQTL)

rownames(phospho2prot)<-phospho2prot$peptide


genotype <- read.table("data/genotype_for_mapping.tsv",header=T,check.names = F)
colnames(genotype)[1]<-"chr"
meta <- read.table("metadata/metadata.tsv", sep="\t", header=T)
rownames(meta) <- meta$culture
meta    <- meta[mixedsort(rownames(meta)), ]
#meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]


#genotypeChrLimit<-as.data.frame(t(sapply(unique(genotype[,1]),function(chr){range(which(genotype[,1]==chr))})))



# get gene prosition (chr coordinates)
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

# get gene position (relative to genetic markers)
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





# identify local QTL
g<-genotype[,4:ncol(genotype)]
g<-t(as.matrix(g))
cor.all <- cor(g,use="pair")

# eQTL
localeQTL<-sapply(eQTL$QtlList$FDR10,function(qtl){
  target <- rownames(eLevelBatchCorrected)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localeQTLTarget<-unique(sapply(eQTL$QtlList$FDR10[localeQTL],function(x)x$target))
haseQTLTarget <- unique(sapply(eQTL$QtlList$FDR10,function(x)x$target))

round(sum(localeQTL)/length(localeQTL) *100,2) #17.30%
round(length(localeQTLTarget)/(nrow(eLevelBatchCorrected))*100,2) #18.27%
length(localeQTLTarget)/length(haseQTLTarget)

# pQTL
localpQTL<-sapply(pQTL$QtlList$FDR10,function(qtl){
  target <- rownames(proteinLevelBatchCorrected)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localpQTLTarget<-unique(sapply(pQTL$QtlList$FDR10[localpQTL],function(x)x$target))
haspQTLTarget <- unique(sapply(pQTL$QtlList$FDR10,function(x)x$target))

round(sum(localpQTL)/length(localpQTL) *100,2) #11.12%
round(length(localpQTLTarget)/(nrow(proteinLevelBatchCorrected))*100,2) #12.41%


# ptQTL
localptQTL<-sapply(ptQTL$QtlList$FDR10,function(qtl){
  target <- rownames(ptQTL$pv)[qtl$target]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(FALSE)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localptQTLTarget<-unique(sapply(ptQTL$QtlList$FDR10[localptQTL],function(x)x$target))
hasptQTLTarget <- unique(sapply(ptQTL$QtlList$FDR10,function(x)x$target))

round(sum(localptQTL)/length(localptQTL) *100,2) # 9.72%
round(length(localptQTLTarget)/(nrow(ptQTL$pv))*100,2) # 6.95%


#phosphoQTL
localphosphoQTL<-sapply(phosphoLevelQTL$QtlList$FDR10,function(qtl){
  target <- phospho2prot[rownames(phosphoLevelBatchCorrected)[qtl$target],"protein"]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(F)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localphosphoQTLTarget<-sapply(phosphoLevelQTL$QtlList$FDR10[localphosphoQTL],function(x)x$target)
hasphosphoQTLTarget <- sapply(phosphoLevelQTL$QtlList$FDR10,function(x)x$target)

round(sum(localphosphoQTL)/length(localphosphoQTL) *100,2) #11.41%
round(length(localphosphoQTLTarget)/(nrow(phosphoLevelBatchCorrected))*100,2) # 8.6%



# phospho Prot
localphosphoProtQTL<-sapply(phosphoProtResidualsQTL$QtlList$FDR10,function(qtl){
  target <- phospho2prot[rownames(phosphoProtResiduals)[qtl$target],"protein"]
  M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
  if(any(is.infinite(M)))return(F)
  COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
  COR>0.8
})

localphosphoProtQTLTarget<-sapply(phosphoProtResidualsQTL$QtlList$FDR10[localphosphoProtQTL],function(x)x$target)
hasphosphoProtQTLTarget <- sapply(phosphoProtResidualsQTL$QtlList$FDR10,function(x)x$target)

round(sum(localphosphoProtQTL)/length(localphosphoProtQTL) *100,2) #10.94%
round(length(localphosphoProtQTLTarget)/(nrow(phosphoProtResiduals))*100,2) # 5.8%

# 
# # phospho RNA
# localphosphoRnaQTL<-sapply(phosphoRnaResidualsQTL$QtlList$FDR10,function(qtl){
#   target <- phospho2prot[rownames(phosphoRnaResiduals)[qtl$target],"protein"]
#   M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
#   if(any(is.infinite(M)))return(F)
#   COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
#   COR>0.8
# })
# 
# localphosphoRnaQTLTarget<-sapply(phosphoRnaResidualsQTL$QtlList$FDR10[localphosphoRnaQTL],function(x)x$target)
# hasphosphoRnaQTLTarget <- sapply(phosphoRnaResidualsQTL$QtlList$FDR10,function(x)x$target)
# 
# round(sum(localphosphoRnaQTL)/length(localphosphoRnaQTL) *100,2) #11.66%
# round(length(localphosphoRnaQTLTarget)/(nrow(phosphoRnaResiduals))*100,2) # 7.28%%

save(localeQTL,localphosphoQTL,localptQTL,localphosphoProtQTL,localpQTL,file="data/localQTL.RData")
# plots portion of local QLT
pdfAndPng("graph/proportion_of_local_QTL", 8,6,expression({
  par(mfrow=c(1,2))
  barplot(rev(c(sum(localeQTL)/length(localeQTL)*100, #e
                sum(localpQTL)/length(localpQTL)*100, #p
                sum(localptQTL)/length(localptQTL)*100, #pt
                sum(localphosphoQTL)/length(localphosphoQTL)*100, #phospho
                sum(localphosphoProtQTL)/length(localphosphoProtQTL)*100)),  #phospho Prot
          names.arg = rev(c("eQTL", "pQTL", "ptQTL", "phQTL", "phpQTL")), horiz = T, las=1, xaxt="n",xlab="proportion of local Linkage") 
  axis(1,at = c(0,5,10,15),labels = paste(c(0,5,10,15),"%",sep="" ))
  
  barplot(rev(c(length(localeQTLTarget)/(nrow(eLevelBatchCorrected))*100, #e
                length(localpQTLTarget)/(nrow(proteinLevelBatchCorrected))*100, #p
                length(localptQTLTarget)/(nrow(ptQTL$pv))*100, #pt
                length(localphosphoQTLTarget)/(nrow(phosphoLevelBatchCorrected))*100, #phospho
                length(localphosphoProtQTLTarget)/(nrow(phosphoProtResiduals))*100)),  #phospho Prot
          names.arg = rev(c("eQTL", "pQTL", "ptQTL", "phQTL", "phpQTL")), horiz = T, las=1, xaxt="n",xlab="proportion of traits with local QTL") 
  axis(1,at = c(0,5,10,15),labels = paste(c(0,5,10,15),"%",sep="" ))
}))

pdf("graph/jan/proportion_of_local_QTL")
par(mfrow=c(1,2))
barplot(rev(c(sum(localeQTL)/length(localeQTL)*100, #e
              sum(localpQTL)/length(localpQTL)*100, #p
              sum(localptQTL)/length(localptQTL)*100, #pt
              sum(localphosphoQTL)/length(localphosphoQTL)*100, #phospho
              sum(localphosphoProtQTL)/length(localphosphoProtQTL)*100)),  #phospho Prot
        names.arg = rev(c("eQTL", "pQTL", "ptQTL", "phQTL", "phpQTL")), horiz = T, las=1, xaxt="n",xlab="proportion of local Linkage") 
axis(1,at = c(0,5,10,15),labels = paste(c(0,5,10,15),"%",sep="" ))

barplot(rev(c(length(localeQTLTarget)/(nrow(eLevelBatchCorrected))*100, #e
              length(localpQTLTarget)/(nrow(proteinLevelBatchCorrected))*100, #p
              length(localptQTLTarget)/(nrow(ptQTL$pv))*100, #pt
              length(localphosphoQTLTarget)/(nrow(phosphoLevelBatchCorrected))*100, #phospho
              length(localphosphoProtQTLTarget)/(nrow(phosphoProtResiduals))*100)),  #phospho Prot
        names.arg = rev(c("eQTL", "pQTL", "ptQTL", "phQTL", "phpQTL")), horiz = T, las=1, xaxt="n",xlab="proportion of traits with local QTL") 
axis(1,at = c(0,5,10,15),labels = paste(c(0,5,10,15),"%",sep="" ))
dev.off()


library(RColorBrewer)
colTrait <-brewer.pal(6,name="Paired")[-1]
names(colTrait) <-c("e","pt","p","phosphoRegressed", "phospho")
colBY <- brewer.pal(10,name="Paired")[8]
colRM <- brewer.pal(10,name="Paired")[10]


d<- cbind(c(sum(localeQTL)/nrow(eQTL$qv), #e
            sum(localptQTL)/nrow(ptQTL$qv), #pt
            sum(localpQTL)/nrow(pQTL$qv), #p
            sum(localphosphoProtQTL)/nrow(phosphoProtResidualsQTL$qv),
            sum(localphosphoQTL)/nrow(phosphoLevelQTL$qv)),
          c(sum(rowSums(eQTL$qv<.1)>0)/nrow(eQTL$qv),
            sum(rowSums(ptQTL$qv<.1)>0)/nrow(ptQTL$qv),
            sum(rowSums(pQTL$qv<.1)>0)/nrow(pQTL$qv),
            sum(rowSums(phosphoProtResidualsQTL$qv<.1)>0)/nrow(phosphoProtResidualsQTL$qv),
            sum(rowSums(phosphoLevelQTL$qv<.1)>0)/nrow(phosphoLevelQTL$qv))
)

pdfAndPng("graph/proportion_trait",5,5,expression({
  par(mar=c(8,4,0.5,1))
  barplot(d[,2],col=colTrait,space = .3,yaxt="n", ylab="proportion of traits",
          names.arg = "",las=3,ylim=c(0,0.9))
  barplot(d[,1],add = T,density=20,angle=30, col="white", space=.3, yaxt="n")
  axis(2,at=c(0,.1,.2,0.3,.4,.5,.6,.7,.8),labels = c("0",paste(c(10,20,30,40,50,60,70,80),"%",sep="")),las=2)
  legend("topright", col = c("white","black"), density=c(0,20), angle=c(NA,30),
         legend=c("with QTL", "with local QTL"), cex=.9,bty="n")
  text((1:5)*1.3-0.5, par("usr")[3] - 0.05, labels = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), srt = 45, xpd = TRUE,cex=1.5,adj=1)
}))
pdfAndPng("graph/local_qtl_proportion",5,5,expression({
  par(mar=c(6,4,0.5,1),cex=1.5)
  barplot(sapply(list(localeQTL,localptQTL,localpQTL,localphosphoProtQTL,localphosphoQTL),FUN=function(x){sum(x)/length(x)}),col=colTrait,space = .3,yaxt="n", ylab="% local QTL",
          names.arg = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"),las=3,ylim=c(0,0.2))
  axis(2,at=c(0,.05,.1,.15,.2),labels = c("0",paste(c(5,10,15,20),"%",sep="")),las=2)
}))

pdf("graph/jan/proportion_trait.pdf")
par(cex.lab=1.5,mar=par()$mar*1.4)
barplot(d[,2],col=colTrait,space = .3,yaxt="n", ylab="proportion of traits",
        names.arg = "",las=3,ylim=c(0,0.9),mgp=par()$mgp*1.4)
barplot(d[,1],add = T,density=20,angle=30, col="white", space=.3, yaxt="n")
axis(2,at=c(0,.1,.2,0.3,.4,.5,.6,.7,.8),labels = c("0",paste(c(10,20,30,40,50,60,70,80),"%",sep="")),las=2,cex.axis=1.5)
legend("topright", col = c("white","black"), density=c(0,20), angle=c(NA,30),
       legend=c("with QTL", "with local QTL"), cex=1.5,bty="n")
text((1:5)*1.3-0.5, par("usr")[3] - 0.05, labels = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), srt = 45, xpd = TRUE,cex=1.5,adj=1)
dev.off()

###consider if genes with local eQTL have more mutations in their promoters
#compute the distance to the closest upstream and downstream mutation
fullGenotype <- read.table("data/strain_genotype.tsv",header=T,check.names = F)
bloomSnps <- read.table(file = "Saccharomyces_cerevisiae/RM_Bloom2013/RM2.vcf",comment.char = "#",header = F,sep="\t",as.is=T)[,1:2] #contains all snps from fullGenotype
byGff <- read.table(file = "Saccharomyces_cerevisiae/sacCer3/REFINED/genes_S288C_R6411_CDSonly_refined.gff",stringsAsFactors = F)
byGff[,9] <- gsub(pattern = ";.*",replacement = "",x = byGff[,9])
byFasta <- read.fasta(file = "Saccharomyces_cerevisiae/sacCer3/S228C_sorted.fa")
utr3p <- read.table(file="Saccharomyces_cerevisiae/SGD_all_ORFs_3prime_UTRs.fsa",sep="?",as.is = T)
utr5p <- read.table(file="Saccharomyces_cerevisiae/SGD_all_ORFs_5prime_UTRs.fsa",sep="?",as.is = T)
uniGenes <- unique(byGff[,9])
utrLengths <- lapply(list(utr3p,utr5p),FUN=function(utr){
  utr <- utr[grepl(">",utr[,1]),1]
  coordsByHeader <- t(sapply(utr,FUN=function(string){
    coords <- strsplit(x = string,split = " ")[[1]][2]
    coords <- gsub(pattern = ".*:",replacement = "",x = coords)
    coords <- as.numeric(strsplit(x = coords,split = "-")[[1]])
    return(c(min(coords),max(coords)))
  }))
  maxLenByGene <- sapply(uniGenes,FUN=function(g){
    hits <- grepl(pattern = paste0("_",g,"_"),x = utr)
    if(sum(hits)==0){return(0)}
    return(max(coordsByHeader[hits,1:2])-min(coordsByHeader[hits,1:2])+1)
  })
  return(maxLenByGene)
})

par(mfrow=c(1,1))
boxplot(rowMeans(eBatchGeneLengthCorrected)~as.factor(utrLengths[[2]][rownames(eBatchGeneLengthCorrected)]<1),names=c("5'UTR","no 5'UTR"),ylab="log2(RNA)")

closestSnpUp <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  startPos <- ifelse(dir=="+",min(subgff[,4:5]),max(subgff[,4:5]))
  snps <- bloomSnps[bloomSnps[,1]==chr,2]
  if(dir=="+"){
    startPos <- startPos-utrLengths[[2]][g]
    snps <- snps[snps<=startPos]
  }else{
    startPos <- startPos+utrLengths[[2]][g]
    snps <- snps[snps>=startPos]
  }
  if(length(snps)==0){return(NA)}
  minDist <- min(abs(snps-startPos))
  return(minDist)
})

closestSnpDown <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  startPos <- ifelse(dir=="-",min(subgff[,4:5]),max(subgff[,4:5]))
  snps <- bloomSnps[bloomSnps[,1]==chr,2]
  if(dir=="+"){
    startPos <- startPos+utrLengths[[1]][g]
    snps <- snps[snps>=startPos]
  }else{
    startPos <- startPos-utrLengths[[1]][g]
    snps <- snps[snps<=startPos]
  }
  if(length(snps)==0){return(NA)}
  minDist <- min(abs(snps-startPos))
  return(minDist)
})
windowSize <- 2000
nSnpUp <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  startPos <- ifelse(dir=="+",min(subgff[,4:5]),max(subgff[,4:5]))
  snps <- bloomSnps[bloomSnps[,1]==chr,2]
  if(dir=="+"){
    startPos <- startPos-utrLengths[[2]][g]
    snps <- snps[snps<startPos]
  }else{
    startPos <- startPos+utrLengths[[2]][g]
    snps <- snps[snps>startPos]
  }
  if(length(snps)==0){return(0)}
  snpCt <- sum(abs(snps-startPos)<=windowSize)
  return(snpCt)
})

nSnpDown <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  startPos <- ifelse(dir=="-",min(subgff[,4:5]),max(subgff[,4:5]))
  snps <- bloomSnps[bloomSnps[,1]==chr,2]
  if(dir=="+"){
    startPos <- startPos+utrLengths[[1]][g]
    snps <- snps[snps>startPos]
  }else{
    startPos <- startPos-utrLengths[[1]][g]
    snps <- snps[snps<startPos]
  }
  if(length(snps)==0){return(0)}
  snpCt <- sum(abs(snps-startPos)<=windowSize)
  return(snpCt)
})


#do genes with local eQTL have more snps?
nSnp <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  cds <- unlist(lapply(1:nrow(subgff),FUN=function(i){
    subgff[i,4]:subgff[i,5]
  }))
  n <- sum(bloomSnps[,1]==chr&bloomSnps[,2]%in%cds)
  return(n)
})

utr5Snps <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  if(dir=="+"){
    utrrange <- (min(subgff[,4:5])-utrLengths[[2]][g]):(min(subgff[,4:5])-1)
  }else{
    utrrange <- (max(subgff[,4:5])+1):(utrLengths[[2]][g]+max(subgff[,4:5]))
  }
  
  n <- sum(bloomSnps[,1]==chr&bloomSnps[,2]%in%utrrange)
  return(n)
})
utr5Snps[utrLengths[[2]]==0] <- 0
utr3Snps <- sapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NA)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  if(dir=="+"){
    utrrange <- (max(subgff[,4:5])+1):(utrLengths[[1]][g]+max(subgff[,4:5]))
  }else{
    utrrange <- (min(subgff[,4:5])-utrLengths[[1]][g]):(min(subgff[,4:5])-1)
  }
  
  n <- sum(bloomSnps[,1]==chr&bloomSnps[,2]%in%utrrange)
  return(n)
})
utr3Snps[utrLengths[[1]]==0] <- 0

#assemble profiles where the snps are for local and distant only (in the 2kb up and downstream of the utrs)
upSnps <- lapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NULL)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  startPos <- ifelse(dir=="+",min(subgff[,4:5]),max(subgff[,4:5]))
  snps <- bloomSnps[bloomSnps[,1]==chr,2]
  if(dir=="+"){
    startPos <- startPos-utrLengths[[2]][g]
    snps <- snps[snps<=startPos]
  }else{
    startPos <- startPos+utrLengths[[2]][g]
    snps <- snps[snps>=startPos]
  }
  if(length(snps)==0){return(NULL)}
  snpDistance <- abs(snps-startPos)
  snpDistance <- snpDistance[snpDistance<=windowSize]
  return(snpDistance)
})
downSnps <- lapply(uniGenes,FUN=function(g){
  if(!g%in%byGff[,9]){return(NULL)}
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  dir <- subgff[1,7]
  startPos <- ifelse(dir=="-",min(subgff[,4:5]),max(subgff[,4:5]))
  snps <- bloomSnps[bloomSnps[,1]==chr,2]
  if(dir=="+"){
    startPos <- startPos+utrLengths[[1]][g]
    snps <- snps[snps>=startPos]
  }else{
    startPos <- startPos-utrLengths[[1]][g]
    snps <- snps[snps<=startPos]
  }
  if(length(snps)==0){return(NULL)}
  snpDistance <- abs(snps-startPos)
  snpDistance <- snpDistance[snpDistance<=windowSize]
  return(snpDistance)
})
names(upSnps) <- names(downSnps) <- uniGenes
allSnps <- nSnpUp+utr5Snps+nSnp+utr3Snps+nSnpDown

geneLength <- sapply(uniGenes,FUN=function(g){
  sum(apply(byGff[byGff[,9]==g,4:5,drop=F],1,FUN=function(x){abs(diff(x))}))
})
nIntrons <- sapply(uniGenes,FUN=function(g){sum(byGff[,9]==g)-1})
#number of missense mutations

RMprotsFasta   <- "Saccharomyces_cerevisiae/RM/proteins_RM.fa"
BYProtsFasta   <- "Saccharomyces_cerevisiae/sacCer3/REFINED/proteins_S288C_R6411.fa"

RMProts   <- read.fasta(RMprotsFasta, as.string=TRUE, forceDNAtolower=FALSE)
BYProts   <- read.fasta(BYProtsFasta, as.string=TRUE, forceDNAtolower=FALSE)
polymorphicProt<-!(unlist(RMProts)==unlist(BYProts))

protPolymorphismPos<-sapply(1:length(BYProts),function(p){
  if(!polymorphicProt[p]){return(NA)}
  if(nchar(BYProts[[p]])!=nchar(RMProts[[p]])){return(NA)}
  which(do.call("!=",strsplit(c(BYProts[[p]], RMProts[[p]]), split = "")))
})
names(protPolymorphismPos)<-names(BYProts)
nPolyByProt <- sapply(protPolymorphismPos,FUN=function(p){ifelse(any(is.na(p)),NA,length(p))})
nPolyByProt[!polymorphicProt] <- 0
nPolyByProt <- nPolyByProt[uniGenes]

#plot for eqtl
par(mfrow=c(2,2))
boxplot(closestSnpUp[rownames(eQTL$qv)[haseQTLTarget]]~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,names=c("only distant","local"),ylab="distance (bp)",main="upstream")
boxplot(closestSnpDown[rownames(eQTL$qv)[haseQTLTarget]]~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,names=c("only distant","local"),ylab="distance (bp)",main="downstream")
boxplot(nSnp[rownames(eQTL$qv)[haseQTLTarget]]~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="snps per gene")
#boxplot(nSnpOnlyUtr[rownames(eQTL$qv)[haseQTLTarget]]~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="snps per gene in the UTRs")
boxplot(nSnpUp[rownames(eQTL$qv)[haseQTLTarget]]~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="upstream snps per gene")
#boxplot(nSnpDown[rownames(eQTL$qv)[haseQTLTarget]]~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="downstream snps per gene")
#plot for ptqtl
par(mfrow=c(2,2))
boxplot(closestSnpUp[setdiff(rownames(ptQTL$qv)[hasptQTLTarget],rownames(eQTL$qv)[localeQTLTarget])]~as.factor(setdiff(rownames(ptQTL$qv)[hasptQTLTarget],rownames(eQTL$qv)[localeQTLTarget])%in%rownames(ptQTL$qv)[localptQTLTarget]),outline=F,names=c("only distant","local"),ylab="distance (bp)",main="upstream")
boxplot(closestSnpDown[setdiff(rownames(ptQTL$qv)[hasptQTLTarget],rownames(eQTL$qv)[localeQTLTarget])]~as.factor(setdiff(rownames(ptQTL$qv)[hasptQTLTarget],rownames(eQTL$qv)[localeQTLTarget])%in%rownames(ptQTL$qv)[localptQTLTarget]),outline=F,names=c("only distant","local"),ylab="distance (bp)",main="downstream")
boxplot(nSnp[rownames(ptQTL$qv)[hasptQTLTarget]]~as.factor(rownames(ptQTL$qv)[hasptQTLTarget]%in%rownames(ptQTL$qv)[localptQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="snps per gene")
#boxplot(nSnpOnlyUtr[rownames(ptQTL$qv)[hasptQTLTarget]]~as.factor(rownames(ptQTL$qv)[hasptQTLTarget]%in%rownames(ptQTL$qv)[localptQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="snps per gene in the UTRs")
boxplot(nSnpUp[rownames(ptQTL$qv)[hasptQTLTarget]]~as.factor(rownames(ptQTL$qv)[hasptQTLTarget]%in%rownames(ptQTL$qv)[localptQTLTarget]),outline=F,names=c("only distant","local"),ylab="#snps",main="snps per gene")
#do the genes with local eqtl have lower or higher expression than others?
par(mfrow=c(1,2))
boxplot(rowMeans(eBatchGeneLengthCorrected[rownames(eQTL$qv)[haseQTLTarget],])~as.factor(rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]),outline=F,ylab="log2(RNA)",names=c("only distant","local"),main="eQTL")
boxplot(rowMeans(proteinLevelBatchCorrected[rownames(pQTL$qv)[haspQTLTarget],],na.rm=T)~as.factor(rownames(pQTL$qv)[haspQTLTarget]%in%rownames(pQTL$qv)[localpQTLTarget]),outline=F,ylab="log2(protein)",names=c("only distant","local"),main="pQTL")

#prepare proteins with local and distant qtl

eLog <- rownames(eQTL$qv)[haseQTLTarget]%in%rownames(eQTL$qv)[localeQTLTarget]
names(eLog) <- rownames(eQTL$qv)[haseQTLTarget]

wilcox.test(rowMeans(eBatchGeneLengthCorrected[names(which(eLog)),],na.rm=T),rowMeans(eBatchGeneLengthCorrected[names(which(!eLog)),],na.rm=T))

pLog <- rownames(pQTL$qv)[haspQTLTarget]%in%rownames(pQTL$qv)[localpQTLTarget]
names(pLog) <- rownames(pQTL$qv)[haspQTLTarget]
ptLog <- rownames(ptQTL$qv)[hasptQTLTarget]%in%rownames(ptQTL$qv)[localptQTLTarget]
names(ptLog) <- rownames(ptQTL$qv)[hasptQTLTarget]
phLog <- unique(phospho2prot[rownames(phosphoLevelQTL$qv)[hasphosphoQTLTarget],2])%in%unique(phospho2prot[rownames(phosphoLevelQTL$qv)[localphosphoQTLTarget],2])
names(phLog) <- unique(phospho2prot[rownames(phosphoLevelQTL$qv)[hasphosphoQTLTarget],2])
phrLog <- unique(phospho2prot[rownames(phosphoProtResidualsQTL$qv)[hasphosphoProtQTLTarget],2])%in%unique(phospho2prot[rownames(phosphoProtResidualsQTL$qv)[localphosphoProtQTLTarget],2])
names(phrLog) <- unique(phospho2prot[rownames(phosphoProtResidualsQTL$qv)[hasphosphoProtQTLTarget],2])

#do proteins with local phResQTL have more missense mutations than proteins with only distant phResQTL
par(mfrow=c(1,1),cex=1.5)
boxplot(nPolyByProt[names(phrLog)]~as.factor(phrLog),outline=F,ylab="number of missense mutations",names=c("distant phResQTL","local phResQTL"))

#show all kinds of mutations for all kinds of local qtl
plotList <- list(eLog,ptLog,pLog,phrLog,phLog)
plotNames <- c("e","pt","p","phRes","ph")
utrPresent <- utrLengths[[1]]>0&utrLengths[[2]]>0
sapply(plotList,FUN=function(gvec){
  gvec <- gvec[utrPresent[names(gvec)]]
  wilcox.test(allSnps[names(which(gvec))],allSnps[names(which(!gvec))],alternative="greater")
})
par(mfrow=c(1,5))
invisible(sapply(1:length(plotList),FUN=function(i){
  gvec <- plotList[[i]]
  gvec <- gvec[utrPresent[names(gvec)]]
  boxplot(allSnps[names(gvec)]~as.factor(gvec),outline=F)
}))

dev.off()
par(mfrow=c(6,6),mar=c(2,2,0,0),oma=c(0,0,1,1))

invisible(sapply(1:5,FUN=function(i){
  plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
  text(labels = plotNames[i],x = 0.5,y=0.5,cex=3)
  gvec <- plotList[[i]]
  boxplot(nSnpUp[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot(utr5Snps[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot(nSnp[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot(utr3Snps[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot(nSnpDown[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
}))
par(cex=1)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "2kb upstream",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "5' UTR",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "CDS",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "3' UTR",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "2kb downstream",x = 0.5,y=0.5)

dev.off()
par(mfrow=c(6,6),mar=c(2,2,0,0),oma=c(0,0,1,1))
plotList <- list(eLog,ptLog,pLog,phrLog,phLog)
plotNames <- c("e","pt","p","phRes","ph")
invisible(sapply(1:5,FUN=function(i){
  plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
  text(labels = plotNames[i],x = 0.5,y=0.5,cex=3)
  gvec <- plotList[[i]]
  boxplot((nSnpUp/2000)[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot((utr5Snps/utrLengths[[2]])[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot((nSnp/geneLength)[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot((utr3Snps/utrLengths[[1]])[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
  boxplot((nSnpDown/2000)[names(gvec)]~as.factor(gvec),outline=F,names=c("only distant","local"),col=colTrait[i])
}))
par(cex=1)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "2kb upstream",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "5' UTR",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "CDS",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "3' UTR",x = 0.5,y=0.5)
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
text(labels = "2kb downstream",x = 0.5,y=0.5)

dev.off()
par(mfrow=c(5,1))
boxplot(lapply(plotList,FUN=function(x){nSnpUp[names(which(x))]}),outline=F,col=colTrait,main="2kb upstream")
boxplot(lapply(plotList,FUN=function(x){utr5Snps[names(which(x))]}),outline=F,col=colTrait,main="5' UTR")
boxplot(lapply(plotList,FUN=function(x){nSnp[names(which(x))]}),outline=F,col=colTrait,main="CDS")
boxplot(lapply(plotList,FUN=function(x){utr3Snps[names(which(x))]}),outline=F,col=colTrait,main="3' UTR")
boxplot(lapply(plotList,FUN=function(x){nSnpDown[names(which(x))]}),outline=F,col=colTrait,main="2kb downstream")

#only genes with utr annotations are considered

par(mfrow=c(5,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){mean(nSnpUp[names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="2kb upstream",names.arg = plotNames,space = 0)
abline(h=mean(nSnpUp[utrPresent[names(nSnpUp)]]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utr5Snps[names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="5' UTR",names.arg = plotNames,space = 0)
abline(h=mean(utr5Snps[utrPresent[names(utr5Snps)]]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nSnp[names(which(x&utrPresent[x]))])}),col=colTrait,main="CDS",names.arg = plotNames,space = 0)
abline(h=mean(nSnp[utrPresent[names(nSnp)]]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utr3Snps[names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="3' UTR",names.arg = plotNames,space = 0)
abline(h=mean(utr3Snps[utrPresent[names(utr3Snps)]]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nSnpDown[names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="2kb downstream",names.arg = plotNames,space = 0)
abline(h=mean(nSnpDown[utrPresent[names(nSnpDown)]]),lty=2,lwd=2)

#correct for target size
par(mfrow=c(5,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){mean(nSnpUp[names(which(x&utrPresent[names(x)]))]/windowSize)}),col=colTrait,main="2kb upstream",names.arg = plotNames,space = 0)
abline(h=mean(nSnpUp[utrPresent[names(nSnpUp)]]/windowSize),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utr5Snps[names(which(x&utrPresent[names(x)]))]/utrLengths[[2]][names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="5' UTR",names.arg = plotNames,space = 0)
abline(h=mean(utr5Snps[utrPresent[names(utr5Snps)]]/utrLengths[[2]][names(which(utrPresent[names(utr5Snps)]))]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nSnp[names(which(x&utrPresent[names(x)]))]/geneLength[names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="CDS",names.arg = plotNames,space = 0)
abline(h=mean(nSnp[utrPresent[names(nSnp)]]/geneLength[names(nSnp)[utrPresent[names(nSnp)]]]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utr3Snps[names(which(x&utrPresent[names(x)]))]/utrLengths[[1]][names(which(x&utrPresent[names(x)]))])}),col=colTrait,main="3' UTR",names.arg = plotNames,space = 0)
abline(h=mean(utr3Snps[utrPresent[names(utr3Snps)]]/utrLengths[[1]][names(which(utrPresent[names(utr3Snps)]))]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nSnpDown[names(which(x&utrPresent[names(x)]))])/windowSize}),col=colTrait,main="2kb downstream",names.arg = plotNames,space = 0)
abline(h=mean(nSnpDown[utrPresent[names(nSnpDown)]]/windowSize),lty=2,lwd=2)

#only genes with protein levels
allDataPresent <- utrPresent[uniGenes]&uniGenes%in%rownames(pQTL$qv)
par(mfrow=c(6,1),mar=c(2,3.5,3,0),cex=.8)
barplot(sapply(plotList,FUN=function(x){mean(nSnpUp[names(which(x&allDataPresent[names(x)]))]/windowSize)}),col=colTrait,main="2kb upstream",names.arg = plotNames,space = 0,las=1)
abline(h=mean(nSnpUp[allDataPresent[names(nSnpUp)]]/windowSize),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utr5Snps[names(which(x&allDataPresent[names(x)]))]/utrLengths[[2]][names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="5' UTR",names.arg = plotNames,space = 0,las=1)
abline(h=mean(utr5Snps[allDataPresent[names(utr5Snps)]]/utrLengths[[2]][names(which(allDataPresent[names(utr5Snps)]))]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nSnp[names(which(x&allDataPresent[names(x)]))]/geneLength[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="CDS",names.arg = plotNames,space = 0,las=1)
abline(h=mean(nSnp[allDataPresent[names(nSnp)]]/geneLength[names(nSnp)[allDataPresent[names(nSnp)]]]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nPolyByProt[names(which(x&allDataPresent[names(x)]))]/(geneLength[names(which(x&allDataPresent[names(x)]))]/3),na.rm=T)}),col=colTrait,main="missense mutations",names.arg = plotNames,space = 0,las=1)
abline(h=mean(nPolyByProt[allDataPresent[names(nPolyByProt)]]/(geneLength[names(nPolyByProt)[allDataPresent[names(nPolyByProt)]]]/3),na.rm=T),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utr3Snps[names(which(x&allDataPresent[names(x)]))]/utrLengths[[1]][names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="3' UTR",names.arg = plotNames,space = 0,las=1)
abline(h=mean(utr3Snps[allDataPresent[names(utr3Snps)]]/utrLengths[[1]][names(which(allDataPresent[names(utr3Snps)]))]),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(nSnpDown[names(which(x&allDataPresent[names(x)]))])/windowSize}),col=colTrait,main="2kb downstream",names.arg = plotNames,space = 0,las=1)
abline(h=mean(nSnpDown[allDataPresent[names(nSnpDown)]]/windowSize),lty=2,lwd=2)

#only genes with protein levels not corrected for length
allDataPresent <- utrPresent[uniGenes]&uniGenes%in%rownames(pQTL$qv)
# par(mfrow=c(6,1),mar=c(2,3.5,3,0),cex=.8)
# barplot(sapply(plotList,FUN=function(x){mean(nSnpUp[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="2kb upstream",names.arg = plotNames,space = 0,las=1)
# abline(h=mean(nSnpUp[allDataPresent[names(nSnpUp)]]),lty=2,lwd=2)
# barplot(sapply(plotList,FUN=function(x){mean(utr5Snps[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="5' UTR",names.arg = plotNames,space = 0,las=1)
# abline(h=mean(utr5Snps[allDataPresent[names(utr5Snps)]]),lty=2,lwd=2)
# barplot(sapply(plotList,FUN=function(x){mean(nSnp[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="CDS",names.arg = plotNames,space = 0,las=1)
# abline(h=mean(nSnp[allDataPresent[names(nSnp)]]),lty=2,lwd=2)
# barplot(sapply(plotList,FUN=function(x){mean(nPolyByProt[names(which(x&allDataPresent[names(x)]))],na.rm=T)}),col=colTrait,main="missense mutations",names.arg = plotNames,space = 0,las=1)
# abline(h=mean(nPolyByProt[allDataPresent[names(nPolyByProt)]],na.rm=T),lty=2,lwd=2)
# barplot(sapply(plotList,FUN=function(x){mean(utr3Snps[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="3' UTR",names.arg = plotNames,space = 0,las=1)
# abline(h=mean(utr3Snps[allDataPresent[names(utr3Snps)]]),lty=2,lwd=2)
# barplot(sapply(plotList,FUN=function(x){mean(nSnpDown[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="2kb downstream",names.arg = plotNames,space = 0,las=1)
# abline(h=mean(nSnpDown[allDataPresent[names(nSnpDown)]]),lty=2,lwd=2)

# pdfAndPng("graph/localQTLmutations/avByLayerRaw", 6,10,expression({
#   lasPar <- 1
#   par(mfrow=c(7,1),mar=c(1,6,1.6,0),cex=1,oma=c(6,1,1,1))
#   barplot(sapply(plotList,FUN=function(x){mean(nSnpUp[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="2kb upstream",names.arg = "",space = 0,las=lasPar)
#   abline(h=mean(nSnpUp[allDataPresent[names(nSnpUp)]]),lty=2,lwd=2)
#   barplot(sapply(plotList,FUN=function(x){mean(utr5Snps[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="5' UTR",names.arg = "",space = 0,las=lasPar)
#   abline(h=mean(utr5Snps[allDataPresent[names(utr5Snps)]]),lty=2,lwd=2)
#   barplot(sapply(plotList,FUN=function(x){mean(nSnp[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="CDS",names.arg = "",space = 0,las=lasPar)
#   abline(h=mean(nSnp[allDataPresent[names(nSnp)]]),lty=2,lwd=2)
#   barplot(sapply(plotList,FUN=function(x){mean(nPolyByProt[names(which(x&allDataPresent[names(x)]))],na.rm=T)}),col=colTrait,main="missense mutations",names.arg = "",space = 0,las=lasPar)
#   abline(h=mean(nPolyByProt[allDataPresent[names(nPolyByProt)]],na.rm=T),lty=2,lwd=2)
#   barplot(sapply(plotList,FUN=function(x){mean(utr3Snps[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="3' UTR",names.arg = "",space = 0,las=lasPar)
#   abline(h=mean(utr3Snps[allDataPresent[names(utr3Snps)]]),lty=2,lwd=2)
#   barplot(sapply(plotList,FUN=function(x){mean(nSnpDown[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="2kb downstream",names.arg = "",space = 0,cex.names=1.4,las=2)
#   abline(h=mean(nSnpDown[allDataPresent[names(nSnpDown)]]),lty=2,lwd=2)
#   mtext(text="mean number of mutations",side = 2,line = 4,outer = F,at = 60,cex = 2.5)
#   par(mar=c(0,6,1.6,0))
#   plot(NULL,xlim=c(0,5),ylim=c(-1,1),axes=F,xlab="",ylab="")
#   text((1:5)-0.5, 1.8, labels = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), srt = 45, xpd = TRUE,cex=1.5,adj=1)
#   
# }))
pdf(file = "graph/localQTLmutations/avByLayerRaw.pdf",width =  6,height = 15)
  lasPar <- 1
  par(mfrow=c(7,1),mar=c(1,6,1.6,0),cex=1,oma=c(6,1,1,1))
  barplot(sapply(plotList,FUN=function(x){mean(nSnpUp[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="2kb upstream",names.arg = "",space = 0,las=lasPar)
  abline(h=mean(nSnpUp[allDataPresent[names(nSnpUp)]]),lty=2,lwd=2)
  barplot(sapply(plotList,FUN=function(x){mean(utr5Snps[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="5' UTR",names.arg = "",space = 0,las=lasPar)
  abline(h=mean(utr5Snps[allDataPresent[names(utr5Snps)]]),lty=2,lwd=2)
  barplot(sapply(plotList,FUN=function(x){mean(nSnp[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="CDS",names.arg = "",space = 0,las=lasPar)
  abline(h=mean(nSnp[allDataPresent[names(nSnp)]]),lty=2,lwd=2)
  barplot(sapply(plotList,FUN=function(x){mean(nPolyByProt[names(which(x&allDataPresent[names(x)]))],na.rm=T)}),col=colTrait,main="missense mutations",names.arg = "",space = 0,las=lasPar)
  abline(h=mean(nPolyByProt[allDataPresent[names(nPolyByProt)]],na.rm=T),lty=2,lwd=2)
  barplot(sapply(plotList,FUN=function(x){mean(utr3Snps[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="3' UTR",names.arg = "",space = 0,las=lasPar)
  abline(h=mean(utr3Snps[allDataPresent[names(utr3Snps)]]),lty=2,lwd=2)
  barplot(sapply(plotList,FUN=function(x){mean(nSnpDown[names(which(x&allDataPresent[names(x)]))])}),col=colTrait,main="2kb downstream",names.arg = "",space = 0,cex.names=1.4,las=2)
  abline(h=mean(nSnpDown[allDataPresent[names(nSnpDown)]]),lty=2,lwd=2)
  mtext(text="mean number of mutations",side = 2,line = 4,outer = F,at = 60,cex = 2.5)
  par(mar=c(0,6,1.6,0))
  plot(NULL,xlim=c(0,5),ylim=c(-1,1),axes=F,xlab="",ylab="")
  text((1:5)-0.5, 1.5, labels = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), srt = 45, xpd = TRUE,cex=1.5,adj=1)
dev.off()


#show intervals of expectected ranges of mutation frequencies through permutations of local targets per layer
snpMat <- cbind(nSnpUp=nSnpUp,utr5Snps=utr5Snps,nSnp=nSnp,nPolyByProt=nPolyByProt,utr3Snps=utr3Snps,nSnpDown=nSnpDown)
set.seed(2021)
permPerLayer <- lapply(plotList,FUN=function(gvec){
  gvec <- gvec[names(gvec)[allDataPresent[names(gvec)]]]
  nPos <- sum(gvec)
  resPerKindOfSNP <- t(sapply(1:10000,FUN=function(i){
    randomGenes <- sample(x = names(gvec),size = nPos,replace = F)
    out <- colMeans(snpMat[randomGenes,],na.rm=T)
    return(out)
  }))
  return(resPerKindOfSNP)
})
names(permPerLayer) <- c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL")
permPerLayer <- abind(permPerLayer,along=3)
realPerLayer <- t(sapply(plotList,FUN=function(gvec){
  gvec <- gvec[names(gvec)[allDataPresent[names(gvec)]]]
  out <- colMeans(snpMat[names(which(gvec)),],na.rm=T)
  return(out)
}))
mutationNames <- c("2kb upstream","5' UTR","CDS","missense mutations","3' UTR","2kb downstream")
names(mutationNames) <- dimnames(permPerLayer)[[2]]
pdf(file = "graph/localQTLmutations/avByLayerPerms.pdf",width =  6,height = 15)
  lasPar <- 1
  par(mfrow=c(7,1),mar=c(1,6,1.6,0),cex=1,oma=c(6,1,1,1))
  invisible(sapply(dimnames(permPerLayer)[[2]],FUN=function(kindOfSnp){
    boxplot(permPerLayer[,kindOfSnp,],xlab="",names=rep("",5),xaxt='n',main=mutationNames[kindOfSnp],ylim=c(min(c(as.vector(permPerLayer[,kindOfSnp,]),as.vector(realPerLayer[,kindOfSnp]))),max(c(as.vector(permPerLayer[,kindOfSnp,]),as.vector(realPerLayer[,kindOfSnp])))))
    points(x=1:5,y=realPerLayer[,kindOfSnp],col=colTrait,pch=19,lwd=3)
  }))
  mtext(text="mean number of mutations",side = 2,line = 4,outer = F,at = 60,cex = 2.5)
  par(mar=c(0,6,1.6,0))
  plot(NULL,xlim=c(0,5),ylim=c(-1,1),axes=F,xlab="",ylab="")
  text((1:5)-0.5, 1.5, labels = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), srt = 45, xpd = TRUE,cex=1.5,adj=1)
dev.off()
n <- nrow(permPerLayer)
pPerLayer <- t(sapply(1:5,FUN=function(i){
  sapply(1:6,function(j){
    max(sum(permPerLayer[,j,i]>=realPerLayer[i,j]),1)/n
  })
}))
colnames(pPerLayer) <- names(mutationNames)
rownames(pPerLayer) <- c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL")

par(mfrow=c(5,1))
barplot(sapply(plotList,FUN=function(x){mean((nSnpUp/2000)[names(which(x))])}),col=colTrait,main="2kb upstream",names.arg = plotNames,space = 0)
barplot(sapply(plotList,FUN=function(x){mean((utr5Snps/utrLengths[[2]])[names(which(x))],na.rm=T)}),col=colTrait,main="5' UTR",names.arg = plotNames,space = 0)
barplot(sapply(plotList,FUN=function(x){mean((nSnp/geneLength)[names(which(x))])}),col=colTrait,main="CDS",names.arg = plotNames,space = 0)
barplot(sapply(plotList,FUN=function(x){mean((utr3Snps/utrLengths[[1]])[names(which(x))],na.rm=T)}),col=colTrait,main="3' UTR",names.arg = plotNames,space = 0)
barplot(sapply(plotList,FUN=function(x){mean((nSnpDown/2000)[names(which(x))])}),col=colTrait,main="2kb downstream",names.arg = plotNames,space = 0)

#Are genes with local QTL enriched for introns?
sapply(plotList,FUN=function(gvec){
  ivec <- nIntrons[names(gvec)]>0
  fisher.test(gvec,ivec)$p.value
})
boxplot(rowMeans(eBatchGeneLengthCorrected)~as.factor(nIntrons[rownames(eBatchGeneLengthCorrected)]>0))
localNoIntron <- names(which(nIntrons[rownames(eQTL$qv)[localeQTLTarget]]==0))
distantNoIntron <- setdiff(names(which(nIntrons[rownames(eQTL$qv)[haseQTLTarget]]==0)),localNoIntron)
boxplot(rowMeans(eBatchGeneLengthCorrected[c(localNoIntron,distantNoIntron),])~as.factor(c(localNoIntron,distantNoIntron)%in%localNoIntron))
localOneIntron <- names(which(nIntrons[rownames(eQTL$qv)[localeQTLTarget]]==1))
distantOneIntron <- setdiff(names(which(nIntrons[rownames(eQTL$qv)[haseQTLTarget]]==1)),localOneIntron)
boxplot(rowMeans(eBatchGeneLengthCorrected[c(localOneIntron,distantOneIntron),])~as.factor(c(localOneIntron,distantOneIntron)%in%localOneIntron))
boxplot(rowMeans(eBatchGeneLengthCorrected[names(which(eLog)),])~as.factor(nIntrons[names(which(eLog))]>0))
boxplot(rowMeans(eBatchGeneLengthCorrected[names(which(!eLog)),])~as.factor(nIntrons[names(which(!eLog))]>0))
par(mfrow=c(1,1))
boxplot(rowMeans(eBatchGeneLengthCorrected[names(eLog),])~as.factor(nIntrons[names(eLog)]>0)*as.factor(eLog),ylab="log2(RNA mean)",names=c("","","",""))
boxplotNames <- c("no intron \nonly distant eQTL \nn=3036",">0 introns \nonly distant eQTL \nn=174","no intron \nlocal eQTL \nn=961",">0 introns \nlocal eQTL \nn=31")
invisible(sapply(1:4,FUN=function(n){
  name <- boxplotNames[n]
  mtext(text = name,side = 1,line = 3,at = n)
}))
boxplot(rowMeans(proteinLevelBatchCorrected[names(pLog),])~as.factor(nIntrons[names(pLog)]>0)*as.factor(pLog))
par(mfrow=c(1,2))
boxplot(rowMeans(eBatchGeneLengthCorrected[rownames(ptQTL$qv),])~as.factor(nIntrons[rownames(eBatchGeneLengthCorrected[rownames(ptQTL$qv),])]>0))
boxplot(rowMeans(proteinLevelBatchCorrected[rownames(ptQTL$qv),])~as.factor(nIntrons[rownames(proteinLevelBatchCorrected[rownames(ptQTL$qv),])]>0))

#do genes with any kind of local QTL have more missense mutations?

par(mfrow=c(5,1))
invisible(sapply(plotList,FUN=function(gvec){
  gvec <- gvec[!is.na(nPolyByProt[names(gvec)])]
  boxplot(nPolyByProt[names(gvec)]~as.factor(gvec),outline=F)
}))
par(mfrow=c(1,5))
invisible(sapply(1:length(plotList),FUN=function(i){
  gvec <- plotList[[i]]
  gvec <- gvec[!is.na(nPolyByProt[names(gvec)])]
  boxplot(nPolyByProt[names(gvec)]/geneLength[names(gvec)]/3~as.factor(gvec),outline=F,col=colTrait[i],main=plotNames[i],names=c("only distant","local"))
}))
sapply(plotList,FUN=function(gvec){
  gvec <- gvec[!is.na(nPolyByProt[names(gvec)])]
  fisher.test(nPolyByProt[names(gvec)]>0,gvec)
})
sapply(plotList,FUN=function(gvec){
  gvec <- gvec[!is.na(nPolyByProt[names(gvec)])]
  wilcox.test(nPolyByProt[names(which(gvec))],nPolyByProt[names(which(!gvec))],alternative="greater")
})
par(mfrow=c(1,1))
barplot(sapply(1:length(plotList),FUN=function(i){
  gvec <- plotList[[i]]
  gvec <- gvec[!is.na(nPolyByProt[names(gvec)])]
  mean(nPolyByProt[names(which(gvec))]/(geneLength[names(which(gvec))]/3))
}),space=0,names.arg = plotNames,col=colTrait)
abline(h=mean(nPolyByProt/(geneLength/3),na.rm=T),lty=2,lwd=2)
protPresent <- uniGenes%in%rownames(pQTL$qv)
names(protPresent) <- uniGenes
par(mfrow=c(1,1))
barplot(sapply(1:length(plotList),FUN=function(i){
  gvec <- plotList[[i]]
  gvec <- gvec[!is.na(nPolyByProt[names(gvec)])&names(gvec)%in%rownames(pQTL$qv)]
  mean(nPolyByProt[names(which(gvec&protPresent[names(gvec)]))]/(geneLength[names(which(gvec&protPresent[names(gvec)]))]/3))
}),space=0,names.arg = plotNames,col=colTrait)
abline(h=mean(nPolyByProt[protPresent]/(geneLength[protPresent]/3),na.rm=T),lty=2,lwd=2)
barplot(c(mean(nPolyByProt[!protPresent[names(nPolyByProt)]],na.rm=T),mean(nPolyByProt[protPresent[names(nPolyByProt)]],na.rm=T)),main="average number of missense mutations",names.arg = c("without protein level","with protein level"))
#do genes with any kind of local QTL have longer UTRs?

utrLengths2 <- lapply(utrLengths,FUN=function(vec){vec[vec==0] <- NA;vec})
par(mfrow=c(2,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths2[[2]][names(which(x))],na.rm=T)}),col=colTrait,main="5' UTR length",names.arg = plotNames,space = 0)
abline(h=mean(utrLengths2[[2]],na.rm=T),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths2[[1]][names(which(x))],na.rm=T)}),col=colTrait,main="3' UTR length",names.arg = plotNames,space = 0)
abline(h=mean(utrLengths2[[1]],na.rm=T),lty=2,lwd=2)

par(mfrow=c(2,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths[[2]][names(which(x))],na.rm=T)}),col=colTrait,main="5' UTR length",names.arg = plotNames,space = 0)
abline(h=mean(utrLengths[[2]],na.rm=T),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths[[1]][names(which(x))],na.rm=T)}),col=colTrait,main="3' UTR length",names.arg = plotNames,space = 0)
abline(h=mean(utrLengths[[1]],na.rm=T),lty=2,lwd=2)

#Are we missing UTR annotations in a biased way?
par(mfrow=c(5,1))
invisible(sapply(plotList,FUN=function(x){
  barplot(c(mean(utrPresent[names(x)[x]]),mean(utrPresent[names(x)[!x]])),space = 0,ylim=c(0,1))
}))
par(mfrow=c(2,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths[[2]][names(which(x))]==0,na.rm=T)}),col=colTrait,main="5' UTR missing",names.arg = plotNames,space = 0)
abline(h=mean(utrLengths[[2]]==0,na.rm=T),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths[[1]][names(which(x))]==0,na.rm=T)}),col=colTrait,main="3' UTR missing",names.arg = plotNames,space = 0)
abline(h=mean(utrLengths[[1]]==0,na.rm=T),lty=2,lwd=2)

par(mfrow=c(2,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths[[2]][names(x)]==0,na.rm=T)}),col=colTrait,main="5' UTR missing",names.arg = plotNames,space = 0,ylim=c(0,1))
abline(h=mean(utrLengths[[2]]==0,na.rm=T),lty=2,lwd=2)
barplot(sapply(plotList,FUN=function(x){mean(utrLengths[[1]][names(x)]==0,na.rm=T)}),col=colTrait,main="3' UTR missing",names.arg = plotNames,space = 0,ylim=c(0,1))
abline(h=mean(utrLengths[[1]]==0,na.rm=T),lty=2,lwd=2)

#Do proteins with a local phResQTL have more phResTraits?

boxplot(table(phospho2prot[rownames(phosphoProtResidualsQTL$qv),2])[names(phrLog)]~as.factor(phrLog))
rownames(phospho2prot) <- phospho2prot[,1]
#pepPerProt <- table(phospho2prot[rownames(phosphoProtResidualsQTL$qv),2])
pepPerProt <- sapply(names(phrLog),FUN=function(g){
  sum(phospho2prot[rownames(phosphoProtResidualsQTL$qv),2]==g)
})
boxplot(pepPerProt[names(phrLog)]~as.factor(phrLog))
wilcox.test(pepPerProt[names(which(phrLog))],pepPerProt[names(which(!phrLog))],alternative="greater")
snpCountsByNPep <- sapply(1:max(pepPerProt),FUN=function(i){
  x <- phrLog
  distantProts <- names(which(!x&pepPerProt[names(x)]==i))
  if(length(distantProts)==0){return(NA)}
  return(mean(nSnp[distantProts]))
})
sum(nSnp[names(which(phrLog))]>snpCountsByNPep[pepPerProt[names(which(phrLog))]],na.rm=T)
sum(nSnp[names(which(phrLog))]<snpCountsByNPep[pepPerProt[names(which(phrLog))]],na.rm=T)
hitsByNPep <- t(sapply(1:max(pepPerProt[names(phrLog)]),FUN=function(i){
  nHits <- sum(pepPerProt[names(phrLog)]==i&phrLog)
  nNonHits <- sum(pepPerProt[names(phrLog)]==i&!phrLog)
  return(c(nHits,nNonHits,sum(nHits,nNonHits)))
}))

sizeByNPep <- t(sapply(1:max(pepPerProt[names(phrLog)]),FUN=function(i){
  c(mean(geneLength[names(phrLog)[pepPerProt[names(phrLog)]==i&phrLog]]),mean(geneLength[names(phrLog)[pepPerProt[names(phrLog)]==i&!phrLog]]))
}))

#Is the size difference in proteins with a local QTL driven by the number of traits?
boxplot(geneLength[names(phrLog)]~as.factor(phrLog),outline=F)
wilcox.test(geneLength[names(which(phrLog))],geneLength[names(which(!phrLog))],alternative="greater")
subPhospo <- phospho2prot[rownames(phosphoProtResidualsQTL$qv),]
subPhospo <- cbind(subPhospo,local=1:nrow(subPhospo)%in%localphosphoProtQTLTarget)
subPhospo <- subPhospo[unique(hasphosphoProtQTLTarget),]
boxplot(geneLength[unique(subPhospo[,2])]~as.factor(unique(subPhospo[,2])%in%unique(subPhospo[subPhospo[,3],2])),outline=F)
pLength <- wilcox.test(geneLength[names(which(phrLog))],geneLength[names(which(!phrLog))],alternative = "greater")$p.value
pDistLength <- sapply(1:10000,FUN=function(i){
  subPhospo2 <- subPhospo
  subPhospo2[,3] <- subPhospo2[sample(1:nrow(subPhospo2)),3]
  phrLog2 <- unique(subPhospo2[,2])%in%subPhospo2[which(subPhospo2[,3]),2]
  names(phrLog2) <- unique(subPhospo2[,2])
  wilcox.test(geneLength[names(which(phrLog2))],geneLength[names(which(!phrLog2))],alternative = "greater")$p.value
})
boxplot(nPolyByProt[unique(subPhospo[,2])]~as.factor(unique(subPhospo[,2])%in%unique(subPhospo[subPhospo[,3],2])),outline=F)
pMut <- wilcox.test(nPolyByProt[names(which(phrLog))],nPolyByProt[names(which(!phrLog))],alternative = "greater")$p.value
pDistMut <- sapply(1:10000,FUN=function(i){
  subPhospo2 <- subPhospo
  subPhospo2[,3] <- subPhospo2[sample(1:nrow(subPhospo2)),3]
  phrLog2 <- unique(subPhospo2[,2])%in%subPhospo2[which(subPhospo2[,3]),2]
  names(phrLog2) <- unique(subPhospo2[,2])
  wilcox.test(nPolyByProt[names(which(phrLog2))],nPolyByProt[names(which(!phrLog2))],alternative = "greater")$p.value
})

sizeBins <- sapply(names(phrLog),FUN=function(x){
  min(sum(quantile(geneLength[names(phrLog)])<=geneLength[x]),4)
})
propPerSizeBin <- sapply(1:4,FUN=function(i){
  mean(subPhospo[sizeBins[subPhospo[,2]]==i,3])
})
mutBins <- sapply(names(phrLog),FUN=function(x){
  min(sum(quantile(nPolyByProt[names(phrLog)],na.rm=T)<=nPolyByProt[x]),4)
})
propPerMutBin <- sapply(1:4,FUN=function(i){
  mean(subPhospo[mutBins[subPhospo[,2]]==i,3],na.rm=T)
})
#Do phosphopeptides from genes with many snps have a higher chance of being regulated locally than others?
localPep <- rownames(phosphoProtResidualsQTL$qv)[localphosphoProtQTLTarget]
distantPep <- setdiff(rownames(phosphoProtResidualsQTL$qv)[hasphosphoProtQTLTarget],localPep)
nonRegPep <- rownames(phosphoProtResidualsQTL$qv)[setdiff(1:nrow(phosphoProtResidualsQTL$qv),hasphosphoProtQTLTarget)]
boxplot(list(nSnp[phospho2prot[nonRegPep,2]],nSnp[phospho2prot[distantPep,2]],nSnp[phospho2prot[localPep,2]]),outline=F,names=c("no QTL","only distant QTL","local QTL"),ylab="number of SNPs in the host protein",col=colTrait[4])
par(mfrow=c(2,2))
boxplot(list(nPolyByProt[phospho2prot[nonRegPep,2]],nPolyByProt[phospho2prot[distantPep,2]],nPolyByProt[phospho2prot[localPep,2]]),outline=F,names=c("no QTL","only distant QTL","local QTL"),ylab="number of missense mutations in the host protein",col=colTrait[4])
boxplot(list(geneLength[phospho2prot[nonRegPep,2]],geneLength[phospho2prot[distantPep,2]],geneLength[phospho2prot[localPep,2]]),outline=F,names=c("no QTL","only distant QTL","local QTL"),ylab="size of the host protein",col=colTrait[4])
boxplot(list(pepPerProt[phospho2prot[nonRegPep,2]],pepPerProt[phospho2prot[distantPep,2]],pepPerProt[phospho2prot[localPep,2]]),outline=F,names=c("no QTL","only distant QTL","local QTL"),ylab="number of phosphopeptides for the host protein",col=colTrait[4])

wilcox.test(nPolyByProt[phospho2prot[distantPep,2]],nPolyByProt[phospho2prot[localPep,2]])
wilcox.test(geneLength[phospho2prot[distantPep,2]],geneLength[phospho2prot[localPep,2]])
wilcox.test(pepPerProt[phospho2prot[distantPep,2]],pepPerProt[phospho2prot[localPep,2]])

localPep2 <- setdiff(rownames(phosphoLevelQTL$qv)[localphosphoQTLTarget],union(localPep,distantPep))
distantPep2 <- setdiff(setdiff(rownames(phosphoLevelQTL$qv)[hasphosphoQTLTarget],localPep2),union(localPep,distantPep))
nonRegPep2 <- setdiff(rownames(phosphoLevelQTL$qv)[setdiff(1:nrow(phosphoLevelQTL$qv),hasphosphoQTLTarget)],union(localPep,distantPep))
boxplot(list(nSnp[phospho2prot[nonRegPep2,2]],nSnp[phospho2prot[distantPep2,2]],nSnp[phospho2prot[localPep2,2]]),outline=F,names=c("no QTL","only distant QTL","local QTL"),ylab="number of SNPs in the host protein",col=colTrait[5])


#Are levels of proteins with local phosphoqtl higher?

par(mfrow=c(2,1),mar=c(2.5,2.5,3,0),cex=1)
boxplot(rowMeans(proteinLevelBatchCorrected[names(phrLog),],na.rm=T)~as.factor(phrLog))
boxplot(rowMeans(eBatchGeneLengthCorrected[names(phrLog),],na.rm=T)~as.factor(phrLog))

#transcript and protein levels for each set of local QTL
par(mfrow=c(2,1),mar=c(2.5,2.5,3,0),cex=1)
barplot(sapply(plotList,FUN=function(x){
  mean(rowMeans(eBatchGeneLengthCorrected[intersect(names(which(x)),rownames(eBatchGeneLengthCorrected)),],na.rm=T))
}),col=colTrait,main="RNA levels",names.arg = plotNames,space = 0)
abline(h=mean(rowMeans(eBatchGeneLengthCorrected)),lwd=2,lty=2)
barplot(sapply(plotList,FUN=function(x){
  mean(rowMeans(proteinLevelBatchCorrected[intersect(names(which(x)),rownames(proteinLevelBatchCorrected)),],na.rm=T))
}),col=colTrait,main="protein levels",names.arg = plotNames,space = 0)
abline(h=mean(rowMeans(proteinLevelBatchCorrected,na.rm=T)),lwd=2,lty=2)

#are genes with local eqtl enriched for having more snps upstream than downstream


sapply(list(eLog,ptLog),FUN=function(gs){
  upVsDown <- nSnpUp[names(gs)]>nSnpDown[names(gs)]
  fisher.test(upVsDown,gs,alternative="greater")
})

sapply(list(eLog,ptLog),FUN=function(gs){
  gs <- gs[utrLengths[[2]][names(gs)]>0&utrLengths[[1]][names(gs)]>0]
  upVsDown <- nSnpUp[names(gs)]>nSnpDown[names(gs)]
  fisher.test(upVsDown,gs,alternative="greater")
})

#gene length

par(mfrow=c(1,5))
invisible(sapply(1:5,FUN=function(i){
  gvec <- plotList[[i]]
  boxplot(geneLength[names(gvec)]~as.factor(gvec),outline=F,main=plotNames[i],ylab="CDS length",names=c("only distant","local"))
}))
par(mfrow=c(1,2))
boxplot(geneLength[phospho2prot[rownames(phosphoProtResidualsQTL$qv)[hasphosphoProtQTLTarget],2]]~as.factor(rownames(phosphoProtResidualsQTL$qv)[hasphosphoProtQTLTarget]%in%rownames(phosphoProtResidualsQTL$qv)[localphosphoProtQTLTarget]),outline=F)
wilcox.test(x=geneLength[phospho2prot[rownames(phosphoProtResidualsQTL$qv)[setdiff(hasphosphoProtQTLTarget,localphosphoProtQTLTarget)],2]],y=geneLength[phospho2prot[rownames(phosphoProtResidualsQTL$qv)[localphosphoProtQTLTarget],2]])
boxplot(geneLength[phospho2prot[rownames(phosphoLevelQTL$qv)[hasphosphoQTLTarget],2]]~as.factor(rownames(phosphoLevelQTL$qv)[hasphosphoQTLTarget]%in%rownames(phosphoLevelQTL$qv)[localphosphoQTLTarget]),outline=F)
wilcox.test(x=geneLength[phospho2prot[rownames(phosphoLevelQTL$qv)[setdiff(hasphosphoQTLTarget,localphosphoQTLTarget)],2]],y=geneLength[phospho2prot[rownames(phosphoLevelQTL$qv)[localphosphoQTLTarget],2]])
#UTR lengths
par(mfrow=c(1,5))
invisible(sapply(1:5,FUN=function(i){
  gvec <- plotList[[i]]
  boxplot(utrLengths[[2]][names(gvec)]~as.factor(gvec),outline=F)
}))
par(mfrow=c(1,5))
invisible(sapply(1:5,FUN=function(i){
  gvec <- plotList[[i]]
  boxplot(utrLengths[[1]][names(gvec)]~as.factor(gvec),outline=F)
}))
#
localEnoPT <- setdiff(names(which(eLog)),names(which(ptLog)))
distantEnoPT <- setdiff(names(which(!eLog)),names(which(ptLog)))
localPTnoE <- setdiff(names(which(ptLog)),names(which(eLog)))
distantPTnoE <- setdiff(names(which(!ptLog)),names(which(eLog)))

localEnoPT <- intersect(setdiff(names(which(eLog)),names(which(ptLog))),names(utrLengths[[2]])[utrLengths[[2]]>0])
distantEnoPT <- intersect(setdiff(names(which(!eLog)),names(which(ptLog))),names(utrLengths[[2]])[utrLengths[[2]]>0])
localPTnoE <- intersect(setdiff(names(which(ptLog)),names(which(eLog))),names(utrLengths[[2]])[utrLengths[[2]]>0])
distantPTnoE <- intersect(setdiff(names(which(!ptLog)),names(which(eLog))),names(utrLengths[[2]])[utrLengths[[2]]>0])

nBreaks <- 20
par(mfrow=c(2,2))
hist(unlist(upSnps[localEnoPT]),breaks=nBreaks,col=colTrait[1],xlab="distance to 5'UTR",main="SNPs upstream of local eQTL")
hist(unlist(upSnps[distantEnoPT]),breaks=nBreaks,col=colTrait[1],xlab="distance to 5'UTR",main="SNPs upstream of distant eQTL")
hist(unlist(downSnps[localEnoPT]),breaks=nBreaks,col=colTrait[1],xlab="distance to 3'UTR",main="SNPs downstream of local eQTL")
hist(unlist(downSnps[distantEnoPT]),breaks=nBreaks,col=colTrait[1],xlab="distance to 3'UTR",main="SNPs downstream of distant eQTL")
par(mfrow=c(2,2))
hist(unlist(upSnps[localPTnoE]),breaks=nBreaks,col=colTrait[2],xlab="distance to 5'UTR",main="SNPs upstream of local ptQTL")
hist(unlist(upSnps[distantPTnoE]),breaks=nBreaks,col=colTrait[2],xlab="distance to 5'UTR",main="SNPs upstream of distant ptQTL")
hist(unlist(downSnps[localPTnoE]),breaks=nBreaks,col=colTrait[2],xlab="distance to 3'UTR",main="SNPs downstream of local ptQTL")
hist(unlist(downSnps[distantPTnoE]),breaks=nBreaks,col=colTrait[2],xlab="distance to 3'UTR",main="SNPs downstream of distant ptQTL")
# local p non ptQTL
both <- intersect(rownames(pQTL$pv)[haspQTLTarget],rownames(ptQTL$pv)[hasptQTLTarget])

a<-intersect(setdiff(rownames(pQTL$pv)[haspQTLTarget],rownames(ptQTL$pv)[hasptQTLTarget]),
             rownames(ptQTL$pv))
length(a)/length(haspQTLTarget)

#Do genes with local pQTL have local eQTL more often than expected?
gs <- rownames(pQTL$qv)[haspQTLTarget]
localP <- rownames(pQTL$qv)[localpQTLTarget]
localE <- intersect(rownames(eQTL$qv)[localeQTLTarget],rownames(pQTL$qv))
contMatEP <- matrix(c(sum(!gs%in%localE&!gs%in%localP),
                      sum(gs%in%localE&!gs%in%localP),
                      sum(!gs%in%localE&gs%in%localP),
                      sum(gs%in%localE&gs%in%localP)),nrow=2,byrow=T)
fisher.test(contMatEP,alternative="greater")

#Do genes with local pQTL have local eQTL more often than expected?
gs <- intersect(rownames(pQTL$qv)[haspQTLTarget],rownames(ptQTL$qv))
localPT <- rownames(ptQTL$qv)[localptQTLTarget]
localE <- intersect(rownames(eQTL$qv)[localeQTLTarget],rownames(ptQTL$qv))
contMatEPT <- matrix(c(sum(!gs%in%localE&!gs%in%localPT),
                       sum(gs%in%localE&!gs%in%localPT),
                       sum(!gs%in%localE&gs%in%localPT),
                       sum(gs%in%localE&gs%in%localPT)),nrow=2,byrow=T)
fisher.test(contMatEPT,alternative="greater")

###compare genes with local pt QTL to the ones without###
##number of mutations

# ct <- 0
# inCoding <- apply(fullGenotype[,1:3],1,FUN=function(x){
#   ct <<- ct+1
#   subgff <- byGff[byGff==as.character(x[1]),,drop=F]
#   subgff <- subgff[subgff[,4]<=as.numeric(x[2]),,drop=F]
#   subgff <- subgff[subgff[,5]>=as.numeric(x[2]),,drop=F]
#   nrow(subgff)>0
# })
snpMatList <- lapply(rownames(ptQTL$qv),FUN=function(g){
  subgff <- byGff[byGff[,9]==g,,drop=F]
  chr <- subgff[1,1]
  cdsPos <- unlist(lapply(1:nrow(subgff),FUN=function(i){
    subgff[i,4]:subgff[i,5]
  }))
  if(subgff[1,7]=="-"){
    cdsPos <- cdsPos[length(cdsPos):1]
  }
  subSnps <- fullGenotype[fullGenotype[,1]==chr,1:4]
  subSnps <- subSnps[subSnps[,2]%in%cdsPos,]
  if(any(sapply(as.character(subSnps[,3]),nchar)>1)|any(sapply(as.character(subSnps[,4]),nchar)>1)|nchar(RMProts[[g]][[1]])!=nchar(BYProts[[g]][[1]])){
    return(NULL)
  }
  snpMat <- cbind(subSnps,rnaPos=sapply(subSnps[,2],FUN=function(pos){which(cdsPos==pos)}))
  rnaSeq <- byFasta[[chr]][cdsPos]
  if(nrow(snpMat)==0){
    return(NULL)
  }
  codons <- apply(snpMat,1,FUN=function(x){
    pos <- as.numeric(x[5])
    codonPos <- ((ceiling(pos/3)-1)*3+1):(ceiling(pos/3)*3)
    refCodon <- rnaSeq[codonPos]
    altCodon <- refCodon
    altCodon[which(codonPos==pos)] <- x[4]
    refCodon <- paste(refCodon,collapse="")
    altCodon <- paste(altCodon,collapse="")
    return(tolower(c(refCodon,altCodon)))
  })
  if(is.matrix(codons)){
    codons <- t(codons)
  }
  snpMat <- cbind(snpMat,codons)
  colnames(snpMat)[6:7] <- c("refCodon","altCodon")
  coding <- apply(snpMat[,6:7,drop=F],1,FUN=function(x){translate(strsplit(x[1],split="")[[1]])!=translate(strsplit(x[2],split="")[[1]])})
  snpMat <- cbind(snpMat,coding=as.numeric(coding))
  return(snpMat)
})
names(snpMatList) <- rownames(ptQTL$qv)
minMut <- t(sapply(snpMatList,FUN=function(mat){
  if(is.null(mat)){return(c(NA,NA))}
  if(any(mat[,"coding"]==0)){
    minSilent <- min(mat[mat[,"coding"]==0,"rnaPos"])
  }else{
    minSilent <- NA
  }
  if(any(mat[,"coding"]==1)){
    minLoud <- min(mat[mat[,"coding"]==1,"rnaPos"])
  }else{
    minLoud <- NA
  }
  return(c(minSilent,minLoud))
}))
colnames(minMut) <- c("silent","loud")
minMut <- cbind(minMut,combined=apply(minMut,1,min,na.rm=T))
minMut[is.na(minMut[,1])&is.na(minMut[,2]),3] <- NA

boxplot(minMut[rownames(minMut)[!is.na(minMut[,1])],1]~as.factor(rownames(minMut)[!is.na(minMut[,1])]%in%localPT),outline=F)
boxplot(minMut[rownames(minMut)[!is.na(minMut[,2])],2]~as.factor(rownames(minMut)[!is.na(minMut[,2])]%in%localPT),outline=F)
boxplot(minMut[rownames(minMut)[!is.na(minMut[,3])],3]~as.factor(rownames(minMut)[!is.na(minMut[,3])]%in%localPT),outline=F)

boxplot(log2(minMut[rownames(minMut)[!is.na(minMut[,1])],1])~as.factor(rownames(minMut)[!is.na(minMut[,1])]%in%localPT),outline=F)
boxplot(log2(minMut[rownames(minMut)[!is.na(minMut[,2])],2])~as.factor(rownames(minMut)[!is.na(minMut[,2])]%in%localPT),outline=F)
boxplot(log2(minMut[rownames(minMut)[!is.na(minMut[,3])],3])~as.factor(rownames(minMut)[!is.na(minMut[,3])]%in%localPT),outline=F)


snpCounts <- t(sapply(snpMatList,FUN=function(mat){c(sum(mat$coding==0),sum(mat$coding==1))}))
rownames(snpCounts) <- rownames(ptQTL$qv)
colnames(snpCounts) <- c("silent","missense")
boxplot(snpCounts[hasptQTLTarget,2]~as.factor(hasptQTLTarget%in%localptQTLTarget),outline=F)
boxplot(snpCounts[hasptQTLTarget,1]~as.factor(hasptQTLTarget%in%localptQTLTarget),outline=F)
par(mfrow=c(2,2))
hist(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),2],breaks = seq(0,42,1))
hist(snpCounts[localptQTLTarget,2],breaks = seq(0,42,1))
hist(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1],breaks = seq(0,42,1))
hist(snpCounts[localptQTLTarget,1],breaks = seq(0,42,1))

#local pt vs only distant pt
contSilent <- matrix(c(sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1]==0),sum(snpCounts[localptQTLTarget,1]==0),sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1]>0),sum(snpCounts[localptQTLTarget,1]>0)),byrow=T,ncol=2)
contMS <- matrix(c(sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),2]==0),sum(snpCounts[localptQTLTarget,2]==0),sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),2]>0),sum(snpCounts[localptQTLTarget,2]>0)),byrow=T,ncol=2)
contSNP <- matrix(c(sum(rowSums(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),])==0),sum(rowSums(snpCounts[localptQTLTarget,])==0),sum(rowSums(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),])>0),sum(rowSums(snpCounts[localptQTLTarget,])>0)),byrow=T,ncol=2)
fisher.test(contSNP,alternative = "greater")
#with qtl vs without
cont_hasqtl <- matrix(c(sum(snpCounts[setdiff(1:nrow(snpCounts),hasptQTLTarget),1]==0),sum(snpCounts[hasptQTLTarget,1]==0),sum(snpCounts[setdiff(1:nrow(snpCounts),hasptQTLTarget),1]>0),sum(snpCounts[hasptQTLTarget,1]>0)),byrow=T,ncol=2)
#with local qtl vs without any qtl
cont_localqtl <- matrix(c(sum(snpCounts[setdiff(1:nrow(snpCounts),hasptQTLTarget),1]==0),sum(snpCounts[localptQTLTarget,1]==0),sum(snpCounts[setdiff(1:nrow(snpCounts),hasptQTLTarget),1]>0),sum(snpCounts[localptQTLTarget,1]>0)),byrow=T,ncol=2)
#with only distant qtl vs without any qtl
cont_distantqtl <- matrix(c(sum(snpCounts[setdiff(1:nrow(snpCounts),hasptQTLTarget),1]==0),sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1]==0),sum(snpCounts[setdiff(1:nrow(snpCounts),hasptQTLTarget),1]>0),sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1]>0)),byrow=T,ncol=2)
#with local pt vs with only distant, both without missense
contSilentNoMS <- matrix(c(sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1]==0&snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),2]==0),
                           sum(snpCounts[localptQTLTarget,1]==0&snpCounts[localptQTLTarget,2]==0),
                           sum(snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),1]>0&snpCounts[setdiff(hasptQTLTarget,localptQTLTarget),2]==0),
                           sum(snpCounts[localptQTLTarget,1]>0&snpCounts[localptQTLTarget,2]==0)),byrow=T,ncol=2)

###check which codons are used
# ct <- 0
# snpMatList <- lapply(rownames(ptQTL$qv),FUN=function(g){
#   ct <<- ct+1
#   mat <- snpMatList[[g]]
#   if(is.null(mat)){
#     return(mat)
#   }
#   if(nrow(mat)==0){
#     return(mat)
#   }
#   subgff <- byGff[byGff[,9]==g,,drop=F]
#   chr <- subgff[1,1]
#   cdsPos <- unlist(lapply(1:nrow(subgff),FUN=function(i){
#     subgff[i,4]:subgff[i,5]
#   }))
#   if(subgff[1,7]=="-"){
#     cdsPos <- cdsPos[length(cdsPos):1]
#   }
#   rnaSeq <- byFasta[[chr]][cdsPos]
#   codons <- apply(mat,1,FUN=function(x){
#     pos <- as.numeric(x[5])
#     codonPos <- ((ceiling(pos/3)-1)*3+1):(ceiling(pos/3)*3)
#     refCodon <- rnaSeq[codonPos]
#     altCodon <- refCodon
#     altCodon[which(codonPos==pos)] <- x[4]
#     refCodon <- paste(refCodon,collapse="")
#     altCodon <- paste(altCodon,collapse="")
#     return(tolower(c(refCodon,altCodon)))
#   })
#   if(is.matrix(codons)){
#     codons <- t(codons)
#   }
#   mat <- cbind(mat,codons)
#   colnames(mat)[7:8] <- c("refCodon","altCodon")
#   return(mat)
# })
# names(snpMatList) <- rownames(ptQTL$qv)


# codons <- apply(expand.grid(c("a","c","g","t"),c("a","c","g","t"),c("a","c","g","t")),1,paste,collapse="")
# aas <- sapply(codons,FUN=function(s){translate(strsplit(x = s,split="")[[1]])})
# refCodonMat <- lapply(rownames(ptQTL$qv),FUN=function(g){
#   subgff <- byGff[byGff[,9]==g,,drop=F]
#   chr <- subgff[1,1]
#   cdsPos <- unlist(lapply(1:nrow(subgff),FUN=function(i){
#     subgff[i,4]:subgff[i,5]
#   }))
#   if(subgff[1,7]=="-"){
#     cdsPos <- cdsPos[length(cdsPos):1]
#   }
#   rnaSeq <- byFasta[[chr]][cdsPos]
#   if(length(rnaSeq)%%3!=0){return(NULL)}
#   rnaSeq <- matrix(rnaSeq,byrow = T,ncol=3)
#   rnaCodons <- apply(rnaSeq,1,paste,collapse="")
#   codonCount <- table(rnaCodons)[codons]
#   names(codonCount) <- codons
#   return(codonCount)
# })
# names(refCodonMat) <- rownames(ptQTL$qv)
# refCodonMat <- do.call("rbind",refCodonMat)
# refCodonMat[is.na(refCodonMat)] <- 0
# mutCodonMatSilent <- lapply(rownames(ptQTL$qv),FUN=function(g){
#   if(is.null(snpMatList[[g]])){
#     return(NULL)
#   }
#   refCodons <- snpMatList[[g]]$refCodon[snpMatList[[g]]$coding==0]
#   if(is.null(refCodons)){
#     out <- rep(0,length(codons))
#     names(out) <- codons
#     return(out)
#   }
#   out <- table(refCodons)[codons]
#   names(out) <- codons
#   return(out)
# })
# names(mutCodonMatSilent) <- rownames(ptQTL$qv)
# mutCodonMatSilent <- do.call("rbind",mutCodonMatSilent)
# mutCodonMatSilent[is.na(mutCodonMatSilent)] <- 0
# 
# mutCodonMatMS <- lapply(rownames(ptQTL$qv),FUN=function(g){
#   if(is.null(snpMatList[[g]])){
#     return(NULL)
#   }
#   refCodons <- snpMatList[[g]]$refCodon[snpMatList[[g]]$coding==1]
#   if(is.null(refCodons)){
#     out <- rep(0,length(codons))
#     names(out) <- codons
#     return(out)
#   }
#   out <- table(refCodons)[codons]
#   names(out) <- codons
#   return(out)
# })
# names(mutCodonMatMS) <- rownames(ptQTL$qv)
# mutCodonMatMS <- do.call("rbind",mutCodonMatMS)
# mutCodonMatMS[is.na(mutCodonMatMS)] <- 0
# 
# #test each ref codon for enrichment in the local pt qtl targets
# localPT <- intersect(rownames(ptQTL$qv)[localptQTLTarget],rownames(mutCodonMatSilent))
# distantPT <- intersect(rownames(ptQTL$qv)[setdiff(hasptQTLTarget,localptQTLTarget)],rownames(mutCodonMatSilent))
# codonContMat <- lapply(codons,FUN=function(codon){
#   nDistantNone <- sum(mutCodonMatSilent[distantPT,codon]==0)
#   nlocalNone <- sum(mutCodonMatSilent[localPT,codon]==0)
#   nDistantSome <- sum(mutCodonMatSilent[distantPT,codon]>0)
#   nlocalSome <- sum(mutCodonMatSilent[localPT,codon]>0)
#   out <- matrix(c(nDistantNone,nlocalNone,nDistantSome,nlocalSome),nrow=2,byrow=2)
#   return(out)
# })
# names(codonContMat) <- codons
# fisherPV <- sapply(codonContMat,FUN=function(mat){
#   fisher.test(mat,alternative = "greater")$p.value
# })
# aaContMat <- lapply(unique(aas),FUN=function(aa){
#   nDistantNone <- sum(rowSums(mutCodonMatSilent[distantPT,aas==aa,drop=F])==0)
#   nlocalNone <- sum(rowSums(mutCodonMatSilent[localPT,aas==aa,drop=F])==0)
#   nDistantSome <- sum(rowSums(mutCodonMatSilent[distantPT,aas==aa,drop=F])>0)
#   nlocalSome <- sum(rowSums(mutCodonMatSilent[localPT,aas==aa,drop=F])>0)
#   out <- matrix(c(nDistantNone,nlocalNone,nDistantSome,nlocalSome),nrow=2,byrow=2)
#   return(out)
# })
# names(aaContMat) <- unique(aas)
# fisherPVbyAA <- sapply(aaContMat,FUN=function(mat){
#   fisher.test(mat,alternative = "greater")$p.value
# })
# 
# library(randomForest)
# ptRF <- randomForest(x=cbind(mutCodonMatSilent,rowSums(mutCodonMatSilent))[c(distantPT,localPT),],y=as.factor(c(rep(0,length(distantPT)),rep(1,length(localPT)))),ntree=1000)
# 
# #collect pairwise codon switches per geneset
# localSilentPairs <- matrix(0,nrow = length(codons),ncol = length(codons))
# rownames(localSilentPairs) <- colnames(localSilentPairs) <- codons
# distantSilentPairs <- localSilentPairs
# for(g in localPT){
#   mat <- snpMatList[[g]]
#   if(is.null(mat)){next}
#   mat <- mat[mat[,"coding"]==0,,drop=F]
#   if(nrow(mat)==0){next}
#   for(i in 1:nrow(mat)){
#     ref <- as.character(mat[i,"refCodon"])
#     alt <- as.character(mat[i,"altCodon"])
#     localSilentPairs[ref,alt] <- localSilentPairs[ref,alt]+1
#     localSilentPairs[alt,ref] <- localSilentPairs[alt,ref]+1
#   }
# }
# for(g in distantPT){
#   mat <- snpMatList[[g]]
#   if(is.null(mat)){next}
#   mat <- mat[mat[,"coding"]==0,,drop=F]
#   if(nrow(mat)==0){next}
#   for(i in 1:nrow(mat)){
#     ref <- as.character(mat[i,"refCodon"])
#     alt <- as.character(mat[i,"altCodon"])
#     distantSilentPairs[ref,alt] <- distantSilentPairs[ref,alt]+1
#     distantSilentPairs[alt,ref] <- distantSilentPairs[alt,ref]+1
#   }
# }
# 
# silentMutsByGene <- lapply(snpMatList,FUN=function(mat){
#   if(is.null(mat)){return(NULL)}
#   mat <- mat[mat[,"coding"]==0&mat$rnaPos%%3==0,,drop=F]
#   if(nrow(mat)==0){return(rep(0,4))}
#   ref <- as.character(mat[,"REF"])
#   alt <- as.character(mat[,"ALT"])
#   nCT <- sum((ref=="C"&alt=="T")|(alt=="C"&ref=="T"))
#   nAG <- sum((ref=="A"&alt=="G")|(alt=="A"&ref=="G"))
#   nAT <- sum((ref=="A"&alt=="T")|(alt=="A"&ref=="T"))
#   nCG <- sum((ref=="C"&alt=="G")|(alt=="C"&ref=="G"))
#   return(c(nCT,nAG,nAT,nCG))
# })
# names(silentMutsByGene) <- names(snpMatList)
# silentMutsByGene <- do.call("rbind",silentMutsByGene)
# 
# #codon usage by AA
# codonUsage <- lapply(unique(aas),FUN=function(aa){
#   byCodon <- colSums(refCodonMat[,aas==aa,drop=F])
#   return(byCodon/sum(byCodon))
# })
# codonUsage <- unlist(codonUsage)
# 
# refCodonRF <- randomForest(x=refCodonMat[c(distantPT,localPT),],y=as.factor(c(rep(0,length(distantPT)),rep(1,length(localPT)))),ntree=1000)
# 
# codonBiasDiff <- sapply(snpMatList,FUN=function(mat){
#   if(is.null(mat)){return(NA)}
#   mat <- mat[mat$coding==0,,drop=F]
#   if(nrow(mat)==0){return(NA)}
#   refFreq <- codonUsage[as.character(mat$refCodon)]
#   altFreq <- codonUsage[as.character(mat$altCodon)]
#   freqMat <- cbind(refFreq,altFreq)
#   diffVec <- apply(freqMat,1,max)/apply(freqMat,1,min)
#   return(mean(diffVec))
# })
# 
# proteinChanges <- sapply(rownames(ptQTL$qv),FUN=function(g){
#   BYProts[[g]][[1]]!=RMProts[[g]]
# })

###Do local pQTL for kinases affect the phosphorylation of other proteins?
goTbl <- read.table("data/gene_association.sgd",comment.char = "!",sep="\t",quote="",as.is=T)
goTbl[,11] <- gsub("\\|{1}.*","",goTbl[,11],fixed = F)
library(topGO)
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

kinasesLocalPQTL <- gene[intersect(BPtermList[["GO:0006468"]],rownames(pQTL$qv)[localpQTLTarget]),]
kinasesLocalEQTL <- gene[intersect(BPtermList[["GO:0006468"]],rownames(eQTL$qv)[localeQTLTarget]),]

geneSub <- gene[which(gene[,5]%%1==0&gene[,6]%%1==0),]
phosphoHits <- lapply(1:nrow(geneSub),FUN=function(i){
  names(which(rowSums(phosphoProtResidualsQTL$qv[,as.numeric(geneSub[i,5]):as.numeric(geneSub[i,6]),drop=F]<=0.1)>0))
})
names(phosphoHits) <- rownames(geneSub)

###comparison of local and distant QTL###
geno <- t(genotype[,-(1:3)])
ePheno <- sapply(rownames(geno),FUN=function(s){
  rowMeans(eBatchGeneLengthCorrected[,meta[colnames(eBatchGeneLengthCorrected),"strain"]==s,drop=F],na.rm=T)
})
pPheno <- sapply(rownames(geno),FUN=function(s){
  rowMeans(proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]==s,drop=F],na.rm=T)
})
eEffect <- lapply(1:ncol(geno),FUN=function(m){
  rowMeans(ePheno[,geno[,m]==0,drop=F],na.rm=T)-rowMeans(ePheno[,geno[,m]==1,drop=F],na.rm=T)
})
eEffect <- do.call("cbind",eEffect)
rownames(eEffect) <- rownames(ePheno)
pEffect <- lapply(1:ncol(geno),FUN=function(m){
  rowMeans(pPheno[,geno[,m]==0,drop=F],na.rm=T)-rowMeans(pPheno[,geno[,m]==1,drop=F],na.rm=T)
})
pEffect <- do.call("cbind",pEffect)
rownames(pEffect) <- rownames(pPheno)

eQtlEffectMat <- lapply(1:length(eQTL$QtlList$FDR10),FUN=function(i){
  qtl <- eQTL$QtlList$FDR10[[i]]
  tgene <- rownames(eQTL$qv)[qtl$target]
  if(!tgene%in%rownames(ptQTL$qv)){return(NULL)}
  return(c(qtl$target,eEffect[tgene,qtl$mostSignificantPredictor],pEffect[tgene,qtl$mostSignificantPredictor],as.numeric(localeQTL[i])))
})
eQtlEffectMat <- do.call("rbind",eQtlEffectMat)
colnames(eQtlEffectMat) <- c("target","eEffect","pEffect","local")

#Are QTL enriched for local ones?
ct <- 0
nLocalByGene <- sapply(rownames(gene)[!is.infinite(gene[,"marker_1"])&!is.infinite(gene[,"marker_2"])],FUN=function(g){
  ct <<- ct+1
  corByMarker <- lapply(as.numeric(gene[g,"marker_1"]):as.numeric(gene[g,"marker_2"]),FUN=function(i){apply(geno,2,cor,y=geno[,i],use="pair")})
  corByMarker <- do.call("rbind",corByMarker)
  sum(colSums(abs(corByMarker)>=0.8)>0)
})

localProp <- c(mean(nLocalByGene)/ncol(geno),1-mean(nLocalByGene)/ncol(geno))
sapply(list(localeQTL,localptQTL,localpQTL,localphosphoProtQTL,localphosphoQTL),FUN=function(qtlvec){
  binom.test(x = sum(qtlvec),n = length(qtlvec),p = localProp[1],alternative = "greater")
})

centralMarker <- floor(apply(gene[,c("marker_1","marker_2")],1,median))
locMat <- matrix(F,nrow = length(centralMarker),ncol=3593)
for(i in 1:length(centralMarker)){
  locMat[i,centralMarker[i]] <- T
}
rownames(locMat) <- names(centralMarker)
fisher.test(as.vector(qv<=0.1),as.vector(locMat[rownames(qv),]))
fisher.test(as.vector(pQTL_results$pQTL$qv<=0.1),as.vector(locMat[rownames(pQTL_results$pQTL$qv),]))
fisher.test(as.vector(pQTL_results$ptQTL$qv<=0.1),as.vector(locMat[rownames(pQTL_results$ptQTL$qv),]))
fisher.test(as.vector(pQTL_results$phosphoLevelQTL$qv<=0.1),as.vector(locMat[phospho2prot[rownames(pQTL_results$phosphoLevelQTL$qv),2],]))
fisher.test(as.vector(pQTL_results$phosphoProtResidualsQTL$qv<=0.1),as.vector(locMat[phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv),2],]))


ct <- 0
localDistEffects <- t(sapply(rownames(ePheno),FUN=function(g){
  ct <<- ct+1
  corByMarker <- lapply(as.numeric(gene[g,"marker_1"]):as.numeric(gene[g,"marker_2"]),FUN=function(i){apply(geno,2,cor,y=geno[,i],use="pair")})
  corByMarker <- do.call("rbind",corByMarker)
  localMarkers <- which(colSums(abs(corByMarker)>=0.8)>0)
  localEffects <- mean(abs(eEffect[g,localMarkers]))
  distantEffects <- mean(abs(eEffect[g,-localMarkers]))
  return(c(localEffects,distantEffects))
}))

meanEffectDistant <- mean(2^abs(eqtlEffectSizes[eqtlEffectSizes[,2]<=0.1&eqtlEffectSizes[,"local"]==0,5])-1)
meanEffectLocal <- mean(2^abs(eqtlEffectSizes[eqtlEffectSizes[,2]<=0.1&eqtlEffectSizes[,"local"]==1,5])-1)
meanEffectLocal/meanEffectDistant

pQtlEffectMat <- lapply(1:length(pQTL$QtlList$FDR10),FUN=function(i){
  qtl <- pQTL$QtlList$FDR10[[i]]
  tgene <- rownames(pQTL$qv)[qtl$target]
  if(!tgene%in%rownames(ptQTL$qv)){return(NULL)}
  return(c(qtl$target,eEffect[tgene,qtl$mostSignificantPredictor],pEffect[tgene,qtl$mostSignificantPredictor],as.numeric(localpQTL[i])))
})
pQtlEffectMat <- do.call("rbind",pQtlEffectMat)
colnames(pQtlEffectMat) <- c("target","eEffect","pEffect","local")

ptQtlEffectMat <- lapply(1:length(ptQTL$QtlList$FDR10),FUN=function(i){
  qtl <- ptQTL$QtlList$FDR10[[i]]
  tgene <- rownames(ptQTL$qv)[qtl$target]
  if(!tgene%in%rownames(ptQTL$qv)){return(NULL)}
  return(c(qtl$target,eEffect[tgene,qtl$mostSignificantPredictor],pEffect[tgene,qtl$mostSignificantPredictor],as.numeric(localptQTL[i])))
})
ptQtlEffectMat <- do.call("rbind",ptQtlEffectMat)
colnames(ptQtlEffectMat) <- c("target","eEffect","pEffect","local")

phlQtlEffectMat <- lapply(1:length(phosphoLevelQTL$QtlList$FDR10),FUN=function(i){
  qtl <- phosphoLevelQTL$QtlList$FDR10[[i]]
  tgene <- phospho2prot[rownames(phosphoLevelQTL$qv)[qtl$target],2]
  if(!tgene%in%rownames(ptQTL$qv)){return(NULL)}
  return(c(qtl$target,eEffect[tgene,qtl$mostSignificantPredictor],pEffect[tgene,qtl$mostSignificantPredictor],as.numeric(localphosphoQTL[i])))
})
phlQtlEffectMat <- do.call("rbind",phlQtlEffectMat)
colnames(phlQtlEffectMat) <- c("target","eEffect","pEffect","local")



par(mfrow=c(2,2))
plot(density(eQtlEffectMat[eQtlEffectMat[,4]==0,2]),col="blue")
lines(density(eQtlEffectMat[eQtlEffectMat[,4]==1,2]),col="red")
plot(density(pQtlEffectMat[pQtlEffectMat[,4]==0,3]),col="blue")
lines(density(pQtlEffectMat[pQtlEffectMat[,4]==1,3]),col="red")
plot(density(abs(eQtlEffectMat[eQtlEffectMat[,4]==0,2])),col="blue")
lines(density(abs(eQtlEffectMat[eQtlEffectMat[,4]==1,2])),col="red")
plot(density(abs(pQtlEffectMat[pQtlEffectMat[,4]==0,3])),col="blue")
lines(density(abs(pQtlEffectMat[pQtlEffectMat[,4]==1,3])),col="red")

###Are local QTL easier to detect?###
#eQTL
localeQTLlists <- lapply(eQTL$QtlList,FUN=function(qtlList){
  localeQTL<-sapply(qtlList,function(qtl){
    target <- rownames(eLevelBatchCorrected)[qtl$target]
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M)))return(FALSE)
    COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
    COR>0.8
  })
})
eqtlHits <- cbind(local=sapply(localeQTLlists,sum),distant=sapply(localeQTLlists,length)-sapply(localeQTLlists,sum))

par(mfrow=c(1,2))
plot(y=eqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,max(eqtlHits)),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="number of eQTL")
points(y=eqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=eqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=eqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
plot(y=eqtlHits[,1]/max(eqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,1),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="proportion of eQTL at FDR<=25%")
points(y=eqtlHits[,1]/max(eqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=eqtlHits[,2]/max(eqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=eqtlHits[,2]/max(eqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
eqtlEffectSizes <- t(sapply(eQTL$QtlList$FDR25,FUN=function(qtl){
  target <- qtl$target
  marker <- qtl$mostSignificantPredictor
  pv <- eQTL$qv[target,marker]
  pvlevel <- which(c(0.01,0.05,0.1,0.15,0.2,0.25)>pv)[1]
  eEffect <- eEffect[target,marker]
  targetName <- rownames(ePheno)[target]
  if(targetName%in%rownames(ptQTL$qv)){
    pEffect <- pEffect[targetName,marker]
    inPT <- 1
  }else{
    pEffect <- NA
    inPT <- 0
  }
  return(c(target,pv,pvlevel,marker,eEffect,pEffect,inPT))
}))
colnames(eqtlEffectSizes) <- c("target","pv","pvlevel","marker","eEffect","pEffect","inPT")
eqtlEffectSizes <- cbind(eqtlEffectSizes,local=as.numeric(localeQTLlists$FDR25))



boxplot(abs(eqtlEffectSizes[,5])~as.factor(eqtlEffectSizes[,"local"]),outline=F)

corByLevelEQtl <- t(sapply(1:6,FUN=function(i){
  mat <- eqtlEffectSizes
  localI <- mat[mat[,"local"]==1&mat[,"inPT"]==1&mat[,"pvlevel"]==i,,drop=F]
  distantI <- mat[mat[,"local"]==0&mat[,"inPT"]==1&mat[,"pvlevel"]==i,,drop=F]
  nLocal <- nrow(localI)
  nDistant <- nrow(distantI)
  if(nLocal>=10){
    corLocal <- cor(localI[,"eEffect"],localI[,"pEffect"],method="spearman")
  }else{
    corLocal <- NA
  }
  if(nDistant>=10){
    corDistant <- cor(distantI[,"eEffect"],distantI[,"pEffect"],method="spearman")
  }else{
    corDistant <- NA
  }
  return(c(nLocal,nDistant,corLocal,corDistant))
}))

#pQTL
localpQTLlists <- lapply(pQTL$QtlList,FUN=function(qtlList){
  localpQTL<-sapply(qtlList,function(qtl){
    target <- rownames(pPheno)[qtl$target]
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M)))return(FALSE)
    COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
    COR>0.8
  })
})
pqtlHits <- cbind(local=sapply(localpQTLlists,sum),distant=sapply(localpQTLlists,length)-sapply(localpQTLlists,sum))

par(mfrow=c(1,2))
plot(y=pqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,max(pqtlHits)),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="number of pQTL")
points(y=pqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=pqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=pqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
plot(y=pqtlHits[,1]/max(pqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,1),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="proportion of pQTL at FDR<=25%")
points(y=pqtlHits[,1]/max(pqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=pqtlHits[,2]/max(pqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=pqtlHits[,2]/max(pqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)

pqtlEffectSizes <- t(sapply(pQTL$QtlList$FDR25,FUN=function(qtl){
  target <- qtl$target
  marker <- qtl$mostSignificantPredictor
  pv <- pQTL$qv[target,marker]
  pvlevel <- which(c(0.01,0.05,0.1,0.15,0.2,0.25)>pv)[1]
  pEffect <- pEffect[target,marker]
  targetName <- rownames(pPheno)[target]
  if(targetName%in%rownames(ptQTL$qv)){
    eEffect <- eEffect[targetName,marker]
    inPT <- 1
  }else{
    eEffect <- NA
    inPT <- 0
  }
  return(c(target,pv,pvlevel,marker,eEffect,pEffect,inPT))
}))
colnames(pqtlEffectSizes) <- c("target","pv","pvlevel","marker","eEffect","pEffect","inPT")
pqtlEffectSizes <- cbind(pqtlEffectSizes,local=as.numeric(localpQTLlists$FDR25))

boxplot(abs(pqtlEffectSizes[,6])~as.factor(pqtlEffectSizes[,"local"]),outline=F)

par(mfrow=c(1,2))
plot(y=pqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,max(pqtlHits)),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="number of pQTL")
points(y=pqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=pqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=pqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
plot(y=pqtlHits[,1]/max(pqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,1),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="proportion of pQTL at FDR<=25%")
points(y=pqtlHits[,1]/max(pqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=pqtlHits[,2]/max(pqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=pqtlHits[,2]/max(pqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)

corByLevelPQtl <- t(sapply(1:6,FUN=function(i){
  mat <- pqtlEffectSizes
  localI <- mat[mat[,"local"]==1&mat[,"inPT"]==1&mat[,"pvlevel"]==i,,drop=F]
  distantI <- mat[mat[,"local"]==0&mat[,"inPT"]==1&mat[,"pvlevel"]==i,,drop=F]
  nLocal <- nrow(localI)
  nDistant <- nrow(distantI)
  if(nLocal>=10){
    corLocal <- cor(localI[,"eEffect"],localI[,"pEffect"],method="spearman")
  }else{
    corLocal <- NA
  }
  if(nDistant>=10){
    corDistant <- cor(distantI[,"eEffect"],distantI[,"pEffect"],method="spearman")
  }else{
    corDistant <- NA
  }
  return(c(nLocal,nDistant,corLocal,corDistant))
}))

#pt 
localptQTLlists <- lapply(ptQTL$QtlList,FUN=function(qtlList){
  if(!is.list(qtlList[[1]])){return(NULL)}
  localpQTL<-sapply(qtlList,function(qtl){
    target <- rownames(ptQTL$qv)[qtl$target]
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M)))return(FALSE)
    COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
    COR>0.8
  })
})
ptqtlHits <- cbind(local=sapply(localptQTLlists,sum),distant=sapply(localptQTLlists,length)-sapply(localptQTLlists,sum))

par(mfrow=c(1,2))
plot(y=ptqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,max(ptqtlHits)),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="number of ptQTL")
points(y=ptqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=ptqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=ptqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
plot(y=ptqtlHits[,1]/max(ptqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,1),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="proportion of ptQTL at FDR<=25%")
points(y=ptqtlHits[,1]/max(ptqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=ptqtlHits[,2]/max(ptqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=ptqtlHits[,2]/max(ptqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)

#phl
localphlQTLlists <- lapply(phosphoLevelQTL$QtlList,FUN=function(qtlList){
  if(!is.list(qtlList[[1]])){return(NULL)}
  localpQTL<-sapply(qtlList,function(qtl){
    target <- phospho2prot[rownames(phosphoLevelQTL$qv)[qtl$target],2]
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M)))return(FALSE)
    COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
    COR>0.8
  })
})
phlqtlHits <- cbind(local=sapply(localphlQTLlists,sum),distant=sapply(localphlQTLlists,length)-sapply(localphlQTLlists,sum))

par(mfrow=c(1,2))
plot(y=phlqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,max(phlqtlHits)),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="number of phosphoLevelQTL")
points(y=phlqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=phlqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=phlqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
plot(y=phlqtlHits[,1]/max(phlqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,1),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="proportion of phosphoLevelQTL at FDR<=25%")
points(y=phlqtlHits[,1]/max(phlqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=phlqtlHits[,2]/max(phlqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=phlqtlHits[,2]/max(phlqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)

#phr
localphrQTLlists <- lapply(phosphoProtResidualsQTL$QtlList,FUN=function(qtlList){
  if(!is.list(qtlList[[1]])){return(NULL)}
  localpQTL<-sapply(qtlList,function(qtl){
    target <- phospho2prot[rownames(phosphoProtResidualsQTL$qv)[qtl$target],2]
    M<-c(gene[target,"marker_1"], gene[target,"marker_2"])
    if(any(is.infinite(M)))return(FALSE)
    COR<-max(cor.all[M[1]:M[2], unlist(lapply(1:nrow(qtl$predictors),function(i)qtl$predictors[i,1]:qtl$predictors[i,2]))] )
    COR>0.8
  })
})
phrqtlHits <- cbind(local=sapply(localphrQTLlists,sum),distant=sapply(localphrQTLlists,length)-sapply(localphrQTLlists,sum))

par(mfrow=c(1,2))
plot(y=phrqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,max(phrqtlHits)),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="number of phResQTL")
points(y=phrqtlHits[,1],x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=phrqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=phrqtlHits[,2],x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)
plot(y=phrqtlHits[,1]/max(phrqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),ylim=c(0,1),xlim=c(0,0.25),col="red",lty=2,type='l',xlab="FDR threshold",ylab="proportion of phResQTL at FDR<=25%")
points(y=phrqtlHits[,1]/max(phrqtlHits[,1]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),pch=20,col="red")
lines(y=phrqtlHits[,2]/max(phrqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",lty=2)
points(y=phrqtlHits[,2]/max(phrqtlHits[,2]),x=c(0.01,0.05,0.1,0.15,0.2,0.25),col="blue",pch=20)

#plots

par(mfrow=c(1,3))
boxplot(abs(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01,"eEffect"])~as.factor(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01,"local"]),outline=F,ylab="absolute RNA l2FC",names=c("distant","local"),main="FDR<=1%")
boxplot(abs(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.1,"eEffect"])~as.factor(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.1,"local"]),outline=F,ylab="absolute RNA l2FC",names=c("distant","local"),main="FDR<=10%")
boxplot(abs(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.25,"eEffect"])~as.factor(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.25,"local"]),outline=F,ylab="absolute RNA l2FC",names=c("distant","local"),main="FDR<=25%")

eLoc <- eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1&eqtlEffectSizes[,"local"]==1,c("eEffect","pEffect")]
eDis <- eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1&eqtlEffectSizes[,"local"]==0,c("eEffect","pEffect")]
pLoc <- pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1&pqtlEffectSizes[,"local"]==1,c("eEffect","pEffect")]
pDis <- pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1&pqtlEffectSizes[,"local"]==0,c("eEffect","pEffect")]


par(mfrow=c(2,2))
plot(eDis,pch=20,col=rgb(0,0,1,0.3),main=paste0("distant eQTL, rho = ",round(cor(eDis,method="spearman")[1,2],digits=2)),xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(eLoc,pch=20,col=rgb(1,0,0,0.3),main=paste0("local eQTL, rho = ",round(cor(eLoc,method="spearman")[1,2],digits=2)),xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(pDis,pch=20,col=rgb(0,0,1,0.3),main=paste0("distant pQTL, rho = ",round(cor(pDis,method="spearman")[1,2],digits=2)),xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(pLoc,pch=20,col=rgb(1,0,0,0.3),main=paste0("local pQTL, rho = ",round(cor(pLoc,method="spearman")[1,2],digits=2)),xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)

eLocM <- lm(eQtlEffectMat[eQtlEffectMat[,4]==1,3]~eQtlEffectMat[eQtlEffectMat[,4]==1,2]+0)
eDisM <- lm(eQtlEffectMat[eQtlEffectMat[,4]==0,3]~eQtlEffectMat[eQtlEffectMat[,4]==0,2]+0)

pLocM <- lm(pQtlEffectMat[pQtlEffectMat[,4]==1,3]~pQtlEffectMat[pQtlEffectMat[,4]==1,2]+0)
pDisM <- lm(pQtlEffectMat[pQtlEffectMat[,4]==0,3]~pQtlEffectMat[pQtlEffectMat[,4]==0,2]+0)

ptLocM <- lm(ptQtlEffectMat[ptQtlEffectMat[,4]==1,3]~ptQtlEffectMat[ptQtlEffectMat[,4]==1,2]+0)
ptDisM <- lm(ptQtlEffectMat[ptQtlEffectMat[,4]==0,3]~ptQtlEffectMat[ptQtlEffectMat[,4]==0,2]+0)


eLocM2 <- lm(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1&eqtlEffectSizes[,"local"]==1,"pEffect"]~eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1&eqtlEffectSizes[,"local"]==1,"eEffect"]+0)
eDisM2 <- lm(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1&eqtlEffectSizes[,"local"]==0,"pEffect"]~eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1&eqtlEffectSizes[,"local"]==0,"eEffect"]+0)

pLocM2 <- lm(pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1&pqtlEffectSizes[,"local"]==1,"pEffect"]~pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1&pqtlEffectSizes[,"local"]==1,"eEffect"]+0)
pDisM2 <- lm(pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1&pqtlEffectSizes[,"local"]==0,"pEffect"]~pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1&pqtlEffectSizes[,"local"]==0,"eEffect"]+0)


#local and distant together with individual slopes
par(mfrow=c(2,2))
plot(eQtlEffectMat[,2:3],pch=20,col=ifelse(eQtlEffectMat[,4]==0,rgb(0,0,1,0.3),rgb(1,0,0,0.3)),main="eQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(h=0,v=0,lty=2,col="black",lwd=1)
abline(eLocM,col="red")
abline(eDisM,col="blue")
legend(x = "topleft",legend = paste0("m = ",round(c(eLocM$coefficients[1],eDisM$coefficients[1]),digits=3)),lty=1,col=c("red","blue"))
plot(pQtlEffectMat[,2:3],pch=20,col=ifelse(pQtlEffectMat[,4]==0,rgb(0,0,1,0.3),rgb(1,0,0,0.3)),main="pQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(h=0,v=0,lty=2,col="black",lwd=1)
abline(pLocM,col="red")
abline(pDisM,col="blue")
legend(x = "topleft",legend = paste0("m = ",round(c(pLocM$coefficients[1],pDisM$coefficients[1]),digits=3)),lty=1,col=c("red","blue"))
plot(ptQtlEffectMat[,2:3],pch=20,col=ifelse(ptQtlEffectMat[,4]==0,rgb(0,0,1,0.3),rgb(1,0,0,0.3)),main="ptQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(h=0,v=0,lty=2,col="black",lwd=1)
abline(ptLocM,col="red")
abline(ptDisM,col="blue")
legend(x = "topleft",legend = paste0("m = ",round(c(ptLocM$coefficients[1],ptDisM$coefficients[1]),digits=3)),lty=1,col=c("red","blue"))

#only e at 10%FDR
par(mfrow=c(1,1))
plot(eQtlEffectMat[,2:3],pch=20,col=ifelse(eQtlEffectMat[,4]==0,rgb(0,0,1,0.3),rgb(1,0,0,0.3)),main="eQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(h=0,v=0,lty=2,col="black",lwd=1)
abline(eLocM,col="red")
abline(eDisM,col="blue")
legend(x = "topleft",legend = paste0("m = ",round(c(eLocM$coefficients[1],eDisM$coefficients[1]),digits=3)),lty=1,col=c("red","blue"))

#only e at 1%FDR
par(mfrow=c(1,1))
plot(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1,c("eEffect","pEffect")],pch=20,col=ifelse(eqtlEffectSizes[eqtlEffectSizes[,"pv"]<=0.01&eqtlEffectSizes[,"inPT"]==1,"local"]==0,rgb(0,0,1,0.3),rgb(1,0,0,0.3)),main="eQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(h=0,v=0,lty=2,col="black",lwd=1)
abline(eLocM2,col="red")
abline(eDisM2,col="blue")
legend(x = "topleft",legend = paste0("m = ",round(c(eLocM2$coefficients[1],eDisM2$coefficients[1]),digits=3)),lty=1,col=c("red","blue"))

#only p at 1%FDR
par(mfrow=c(1,1))
plot(pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1,c("eEffect","pEffect")],pch=20,col=ifelse(pqtlEffectSizes[pqtlEffectSizes[,"pv"]<=0.01&pqtlEffectSizes[,"inPT"]==1,"local"]==0,rgb(0,0,1,0.3),rgb(1,0,0,0.3)),main="pQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(h=0,v=0,lty=2,col="black",lwd=1)
abline(pLocM2,col="red")
abline(pDisM2,col="blue")
legend(x = "topleft",legend = paste0("m = ",round(c(pLocM2$coefficients[1],pDisM2$coefficients[1]),digits=3)),lty=1,col=c("red","blue"))


#all separately and with slope=1
par(mfrow=c(3,2))
plot(eQtlEffectMat[eQtlEffectMat[,4]==0,2:3],pch=20,col=rgb(0,0,1,0.3),main="distant eQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(eQtlEffectMat[eQtlEffectMat[,4]==1,2:3],pch=20,col=rgb(1,0,0,0.3),main="local eQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(pQtlEffectMat[pQtlEffectMat[,4]==0,2:3],pch=20,col=rgb(0,0,1,0.3),main="distant pQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(pQtlEffectMat[pQtlEffectMat[,4]==1,2:3],pch=20,col=rgb(1,0,0,0.3),main="local pQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(ptQtlEffectMat[ptQtlEffectMat[,4]==0,2:3],pch=20,col=rgb(0,0,1,0.3),main="distant ptQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)
plot(ptQtlEffectMat[ptQtlEffectMat[,4]==1,2:3],pch=20,col=rgb(1,0,0,0.3),main="local ptQTL",xlab="RNA l2FC",ylab="protein l2FC")
abline(a=0,b=1,h=0,v=0,lty=2,col="black",lwd=2)

###which proteins correlate well with their transcript
avGenes <- intersect(rownames(pPheno),rownames(ePheno))
corMat <- t(sapply(avGenes,FUN=function(g){
  gcor <- cor(ePheno[g,],pPheno[g,],use="pair",method="spearman")
  pv <- cor.test(ePheno[g,],pPheno[g,],alternative="greater",method="spearman")$p.value
  return(c(gcor,pv))
}))
colnames(corMat) <- c("cor","pv")
corMat <- cbind(corMat,fdr=p.adjust(corMat[,2],method="fdr"))

###compare the distance of phosphosites with local qtl to mutations to that for phosphosites with distant qtl
load(file = "data/phosphoProt.RData")
peptides<-gsub("\\(.+\\)","",metaPhosphoProtResiduals$peptide)
peptides<-gsub("\\[[^\\[]+\\]","",peptides)
peptides<-gsub("_[0-9]","",peptides)
metaPhosphoProtResiduals$phospho_pos<-sapply(1:nrow(metaPhosphoProtResiduals), function(i){
  regexpr(peptides[i],BYProts[[metaPhosphoProtResiduals$protein[i]]][1],fixed = T) + round(mean(gregexpr("S", peptides[i])[[1]]))
})

polymorphicProt<-!(unlist(RMProts)==unlist(BYProts))

protPolymorphismPos<-sapply(1:length(BYProts),function(p){
  if(!polymorphicProt[p]){return(NA)}
  which(do.call("!=",strsplit(c(BYProts[[p]], RMProts[[p]]), split = "")))
})
names(protPolymorphismPos)<-names(BYProts)

phrMeta <- cbind(metaPhosphoProtResiduals,minDist=apply(metaPhosphoProtResiduals,1,FUN=function(x){
  min(abs(as.numeric(x[3])-protPolymorphismPos[[x[2]]]))
}))
phrMeta <- phrMeta[phrMeta[,2]%in%unique(phospho2prot[rownames(pQTL_results$phosphoProtResidualsQTL$qv)[localphosphoProtQTLTarget],2]),]
boxplot(phrMeta[,"minDist"]~as.factor(phrMeta[,1]%in%rownames(pQTL_results$phosphoProtResidualsQTL$pv)[localphosphoProtQTLTarget]))
boxplot(phrMeta[,"minDist"]~as.factor(phrMeta[,1]%in%rownames(pQTL_results$phosphoProtResidualsQTL$pv)[localphosphoProtQTLTarget]),outline=F)
sigPeps <- rownames(pQTL_results$phosphoProtResidualsQTL$qv)[localphosphoProtQTLTarget]
intProts <- names(which(sapply(unique(phrMeta$protein),FUN=function(p){
  any(!phrMeta[phrMeta[,2]==p,1]%in%sigPeps)&any(phrMeta[phrMeta[,2]==p,1]%in%sigPeps)
})))
avDistsPerProt <- t(sapply(intProts,FUN=function(p){
  sigDist <- mean(phrMeta$minDist[phrMeta$protein==p&phrMeta$peptide%in%sigPeps])
  nonSigDist <- mean(phrMeta$minDist[phrMeta$protein==p&!phrMeta$peptide%in%sigPeps])
  return(c(sigDist,nonSigDist))
}))

pdfAndPng("graph/localVsDistantPhResSameProt",8,8,expr = expression({
  par(mfrow=c(1,1),cex=1.5)
  plot(avDistsPerProt,xlab="distance for significant linkages (aa)",ylab="distance for not significant linkages (aa)")
  abline(a=0,b=1)
  
}))


wilcox.test(x = avDistsPerProt[,1],y=avDistsPerProt[,2],paired=T,alternative="less")
phrMeta2 <- cbind(metaPhosphoProtResiduals,minDist=apply(metaPhosphoProtResiduals,1,FUN=function(x){
  min(abs(as.numeric(x[3])-protPolymorphismPos[[x[2]]]))
}))
avDistToPolys <- apply(phrMeta2[,c(2,3)],1,FUN=function(x){
  pos <- as.numeric(x[2])
  prot <- x[1]
  mean(abs(protPolymorphismPos[[prot]]-pos))
})
boxplot(avDistToPolys~as.factor(1:nrow(pQTL_results$phosphoProtResidualsQTL$qv)%in%localphosphoProtQTLTarget))
avDistToPolysCorr <- apply(phrMeta2[,c(2,3)],1,FUN=function(x){
  pos <- as.numeric(x[2])
  prot <- x[1]
  mean(abs(protPolymorphismPos[[prot]]-pos)/nchar(BYProts[[prot]]))
})

load("data/phosphoLevel.RData")
rownames(phospho2prot) <- phospho2prot[,1]
phMeta <- phospho2prot
peptides<-gsub("\\(.+\\)","",phMeta$peptide)
peptides<-gsub("\\[[^\\[]+\\]","",peptides)
peptides<-gsub("_[0-9]","",peptides)
phMeta$phospho_pos<-sapply(1:nrow(phMeta), function(i){
  regexpr(peptides[i],BYProts[[phMeta$protein[i]]][1],fixed = T) + round(mean(gregexpr("S", peptides[i])[[1]]))
})
phMeta <- phMeta[!phMeta[,2]%in%rownames(pQTL$qv),]
phMeta <- cbind(phMeta,minDist=apply(phMeta,1,FUN=function(x){
  min(abs(as.numeric(x[3])-protPolymorphismPos[[x[2]]]))
}))
sigPeps <- rownames(pQTL_results$phosphoLevelQTL$qv)[localphosphoQTLTarget]
intProts <- names(which(sapply(unique(phMeta$protein),FUN=function(p){
  any(!phMeta[phMeta[,2]==p,1]%in%sigPeps)&any(phMeta[phMeta[,2]==p,1]%in%sigPeps)&!p%in%rownames(eQTL$qv)[localeQTLTarget]
})))
avDistsPerProt2 <- t(sapply(intProts,FUN=function(p){
  sigDist <- mean(phMeta$minDist[phMeta$protein==p&phMeta$peptide%in%sigPeps])
  nonSigDist <- mean(phMeta$minDist[phMeta$protein==p&!phMeta$peptide%in%sigPeps])
  return(c(sigDist,nonSigDist))
}))

source("lib/general_function.R")
pdfAndPng("graph/phResQtlDistance",3,4,expr = expression({
  par(mfrow=c(1,1),cex=.9)
  boxplot(phrMeta[,"minDist"]~as.factor(!phrMeta[,1]%in%rownames(pQTL_results$phosphoProtResidualsQTL$pv)[localphosphoProtQTLTarget]),outline=F,names=c("local \nphResQTL","distant \nphResQTL"),ylab="closest polymorphism (aa)",xaxt='n',las=2,col=colTrait[4])
  axis(side = 1,tick = F,at = c(1,2),labels = c("local \nphResQTL","distant \nphResQTL"))
  
}))
