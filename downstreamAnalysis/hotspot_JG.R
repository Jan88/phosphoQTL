library(gtools) #JG
load(file = "data/phosphoProt.RData")
load("data/proteinLevel.RData")
load("data/pQTL_results170815.RData")
genotype <- read.table("data/genotype_for_mapping.tsv",header=T,check.names = F)
colnames(genotype)[1]<-"chr"
#meta <- read.table("metadata/metadata.csv", sep=",", header=T) 
meta <- read.table("metadata/metadata.tsv", sep="\t", header=T) #JG
rownames(meta) <- meta$culture
meta    <- meta[mixedsort(rownames(meta)), ]
#meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
meta <- meta[meta$matching_with_brem_genotype,] #JG

source("lib/genotyping_function.R")
source("lib/general_function.R")

genotypeChrLimit<-as.data.frame(t(sapply(unique(genotype[,1]),function(chr){range(which(genotype[,1]==chr))}))) 
#Consider that chr 9 is between 4 and 5. JG

colnames(genotypeChrLimit)<-c('start','end')


#################################
###### hotspot analysis #########
#################################

#make the bins
library(seqinr)
#genomeFasta <- read.fasta(file = "genomes/sacCer3/S228C_sorted.fa")
genomeFasta <- read.fasta(file = "Saccharomyces_cerevisiae/sacCer3/S228C_sorted.fa") #JG
chr_length <- sapply(genomeFasta,length)
chr_length <- chr_length[paste("chr",as.roman(1:16),sep="")]
rm(genomeFasta)

bins<-lapply(names(chr_length),function(chr){
  s<-seq(1,chr_length[chr],40000)
  s[length(s)]<-chr_length[chr]
  data.frame(chr=chr,start=s[1:(length(s)-1)],end=s[2:length(s)])
})
bins<- do.call("rbind",bins)
bins$abs <- mapply(getGenomicPosition,bins$chr,rowMeans(cbind(bins$start,bins$end)),MoreArgs=list(chromosome.size=chr_length))


# conts the phosphoProtResiduals falling in the bins -----
genotypeMeanPos<-rowMeans(genotype[,c("start","end")])

a<-lapply(pQTL_results, function(result){
  QTLList <- result$QtlList$FDR10
  QTL_genomic_pos <- data.frame( 
    pos=unlist(sapply(QTLList,function(qtl){
      genotypeMeanPos[qtl$mostSignificantPredictor]
    })),
    chr=unlist(sapply(QTLList,function(qtl){
      genotype$chr[qtl$mostSignificantPredictor]
    }))
  )
  #counts
  sapply(1:nrow(bins),function(i){
    sum(QTL_genomic_pos$chr==bins$chr[i] & QTL_genomic_pos$pos > bins$start[i] & QTL_genomic_pos$pos < bins$end[i])
  })
})

bins <- cbind(bins,a)

load("data/eQTL_results160831.RData")
QTLList <- QtlList$FDR10
QTL_genomic_pos <- data.frame( 
  pos=unlist(sapply(QTLList,function(qtl){
    genotypeMeanPos[qtl$mostSignificantPredictor]
  })),
  chr=unlist(sapply(QTLList,function(qtl){
    genotype$chr[qtl$mostSignificantPredictor]
  }))
)
#counts
bins$eQTL <- sapply(1:nrow(bins),function(i){
  sum(QTL_genomic_pos$chr==bins$chr[i] & QTL_genomic_pos$pos > bins$start[i] & QTL_genomic_pos$pos < bins$end[i])
})



sum(bins$phosphoProtResidualsQTL)/nrow(bins)
which(dpois(1:100,sum(bins$phosphoProtResidualsQTL)/nrow(bins))<1/nrow(bins))

M=0.05
max(which(ppois(1:50,sum(bins$eQTL)/nrow(bins),lower.tail = F)>M)) # 26
max(which(ppois(1:50,sum(bins$pQTL)/nrow(bins),lower.tail = F)>M)) # 11
max(which(ppois(1:50,sum(bins$ptQTL)/nrow(bins),lower.tail = F)>M)) # 7
max(which(ppois(1:50,sum(bins$phosphoLevelQTL)/nrow(bins),lower.tail = F)>M)) # 9
max(which(ppois(1:50,sum(bins$phosphoRnaResidualsQTL)/nrow(bins),lower.tail = F)>M)) # 7
max(which(ppois(1:50,sum(bins$phosphoProtResidualsQTL)/nrow(bins),lower.tail = F)>M)) # 3
M=0.01
max(which(ppois(1:50,sum(bins$eQTL)/nrow(bins),lower.tail = F)>M)) # 30
max(which(ppois(1:50,sum(bins$pQTL)/nrow(bins),lower.tail = F)>M)) # 13
max(which(ppois(1:50,sum(bins$ptQTL)/nrow(bins),lower.tail = F)>M)) # 9
max(which(ppois(1:50,sum(bins$phosphoLevelQTL)/nrow(bins),lower.tail = F)>M)) # 11
max(which(ppois(1:50,sum(bins$phosphoRnaResidualsQTL)/nrow(bins),lower.tail = F)>M)) # 9
max(which(ppois(1:50,sum(bins$phosphoProtResidualsQTL)/nrow(bins),lower.tail = F)>M)) # 4




binsChrLimit<-as.data.frame(t(sapply(unique(bins[,1]),function(chr){range(which(bins[,1]==chr))})))
colnames(binsChrLimit)<-c('start','end')

library(RColorBrewer)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]

names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

# pdfAndPng("graph/hotspotsJG",8,10,expression({
#   par(mfrow=c(5,1))
#   par(mar=c(2,4.2,1,1))
# 
#   plot(bins$eQTL, type="l", xaxt="n", xlab="", ylab=expression(italic(N)~~eQTL),las=2)
#   clip(1,nrow(bins),30,1000)
#   lines(bins$eQTL, type="l", col=colors[1],lwd=3.5)
#   do.call("clip", as.list(par("usr")))  # reset to plot region
#   abline(v=binsChrLimit$end[1:15],lty=2)
#   axis(1,at=rowMeans(binsChrLimit), as.roman(1:16), tick = F, line = -1)
# 
#   plot(bins$ptQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~ptQTL),las=2)
#   clip(1,nrow(bins),9,1000)
#   lines(bins$ptQTL, type="l", col=colors[2],lwd=3.5)
#   do.call("clip", as.list(par("usr")))  # reset to plot region
#   abline(v=binsChrLimit$end[1:15],lty=2)
#   axis(1,at=rowMeans(binsChrLimit), as.roman(1:16), tick = F, line = -1)
#   
#   plot(bins$pQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~pQTL),las=2)
#   clip(1,nrow(bins),13,1000)
#   lines(bins$pQTL, type="l", col=colors[3],lwd=3.5)
#   do.call("clip", as.list(par("usr")))  # reset to plot region
#   abline(v=binsChrLimit$end[1:15],lty=2)
#   axis(1,at=rowMeans(binsChrLimit), as.roman(1:16), tick = F, line = -1)
# 
#   plot(bins$phosphoProtResidualsQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~phResQTL),las=2)
#   clip(1,nrow(bins),4,1000)
#   lines(bins$phosphoProtResidualsQTL, type="l", col=colors[4],lwd=3.5)
#   do.call("clip", as.list(par("usr")))  # reset to plot region
#   abline(v=binsChrLimit$end[1:15],lty=2)
#   axis(1,at=rowMeans(binsChrLimit), as.roman(1:16), tick = F, line = -1)
# 
#   plot(bins$phosphoLevelQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~phQTL),las=2)
#   clip(1,nrow(bins),11,1000)
#   lines(bins$phosphoLevelQTL, type="l", col=colors[5],lwd=3.5)
#   do.call("clip", as.list(par("usr")))  # reset to plot region
#   abline(v=binsChrLimit$end[1:15],lty=2)
#   axis(1,at=rowMeans(binsChrLimit), as.roman(1:16), tick = F, line = -1)
# 
# 
#   
# }))

overLine <- function(binVec,limit,col,lwd,lty=1){
  showBins <- which(binVec>limit)
  plotMat <- matrix(c(showBins[1],showBins[1]),nrow=1)
  for(i in 2:length(showBins)){
    if(showBins[i]-showBins[i-1]==1){
      plotMat[nrow(plotMat),2] <- showBins[i]
    }else{
      plotMat <- rbind(plotMat,c(showBins[i],showBins[i]))
    }
  }
  clip(1,length(binVec),limit,1000)
  apply(plotMat,1,FUN=function(pl){
    lines(x=max(pl[1]-1,1):min(pl[2]+1,length(binVec)),y=binVec[(pl[1]-1):(pl[2]+1)],lwd=lwd,col=col,lty=lty)
  })
  do.call("clip", as.list(par("usr"))) 
}
labelFun <- function(values,roundingLevel){
  values <- min(0,min(values)):max(0,max(values))
  evenLevels <- values[values%%roundingLevel==0]
  evenLevels <- evenLevels[evenLevels!=0]
  return(c(range(evenLevels),length(evenLevels)-1))
}
pdfAndPng("graph/hotspotsJG",10,10,expression({
  par(mfrow=c(5,1))
  par(mar=c(1,4.2,1,1),cex=1.5,oma=c(2,0,0,0))
  lwd <- 3.5
  col <- "black"
  chrCex <- 1
  plot(bins$eQTL, type="l", xaxt="n", xlab="", ylab=expression(italic(N)~~eQTL),las=2,lwd=lwd,col=col)
  #clip(1,nrow(bins),30,1000)
  #plotMat <- cbind(bins$eQTL[bins$eQTL>30],which(bins$eQTL>30))
  #apply(plotMat,1,lines)
  #lines(y=bins$eQTL[bins$eQTL], type="l", col=colors[1],lwd=3.5)
  overLine(binVec = bins$eQTL,limit = 30,col = colors[1],lwd = 3.5)
  #do.call("clip", as.list(par("usr")))  # reset to plot region
  abline(v=binsChrLimit$end[1:15],lty=2,col="darkgrey")
  mtext(text = as.roman(1:16),side = 1,at = rowMeans(binsChrLimit),cex=chrCex)
  
  plot(bins$ptQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~ptQTL),las=2,lwd=lwd,col=col)
  # clip(1,nrow(bins),9,1000)
  # lines(bins$ptQTL, type="l", col=colors[2],lwd=3.5)
  # do.call("clip", as.list(par("usr")))  # reset to plot region
  overLine(binVec = bins$ptQTL,limit = 9,col = colors[2],lwd = 3.5)
  abline(v=binsChrLimit$end[1:15],lty=2,col="darkgrey")
  mtext(text = as.roman(1:16),side = 1,at = rowMeans(binsChrLimit),cex=chrCex)
  
  plot(bins$pQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~pQTL),las=2,lwd=lwd,col=col)
  # clip(1,nrow(bins),13,1000)
  # lines(bins$pQTL, type="l", col=colors[3],lwd=3.5)
  # do.call("clip", as.list(par("usr")))  # reset to plot region
  overLine(binVec = bins$pQTL,limit = 13,col = colors[3],lwd = 3.5)
  abline(v=binsChrLimit$end[1:15],lty=2,col="darkgrey")
  mtext(text = as.roman(1:16),side = 1,at = rowMeans(binsChrLimit),cex=chrCex)
  
  plot(bins$phosphoProtResidualsQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~phResQTL),las=2,lwd=lwd,col=col)
  # clip(1,nrow(bins),4,1000)
  # lines(bins$phosphoProtResidualsQTL, type="l", col=colors[4],lwd=3.5)
  # do.call("clip", as.list(par("usr")))  # reset to plot region
  overLine(binVec = bins$phosphoProtResidualsQTL,limit = 4,col = colors[4],lwd = 3.5)
  abline(v=binsChrLimit$end[1:15],lty=2,col="darkgrey")
  mtext(text = as.roman(1:16),side = 1,at = rowMeans(binsChrLimit),cex=chrCex)
  
  plot(bins$phosphoLevelQTL, type="l", xaxt="n", xlab="",ylab=expression(italic(N)~~phQTL),las=2,lwd=lwd,col=col)
  # clip(1,nrow(bins),11,1000)
  # lines(bins$phosphoLevelQTL, type="l", col=colors[5],lwd=3.5)
  # do.call("clip", as.list(par("usr")))  # reset to plot region
  overLine(binVec = bins$phosphoLevelQTL,limit = 11,col = colors[5],lwd = 3.5)
  abline(v=binsChrLimit$end[1:15],lty=2,col="darkgrey")
  mtext(text = as.roman(1:16),side = 1,at = rowMeans(binsChrLimit),cex=chrCex)
  
  mtext(text = "QTL-location",side = 1,line = 0.5,outer = T,at = 0.5,cex=2)
  
}))

#JG:
binMeanPos <- rowMeans(bins[,2:3])
binPerMarker <- sapply(1:length(genotypeMeanPos),FUN=function(i){
  which(bins[,1]==genotype[i,1]&genotypeMeanPos[i]>bins[,2]&genotypeMeanPos[i]<=bins[,3])
})
save(bins,binPerMarker, file = "data/binInfo.RData")

### hospot exploration --------------

bins[bins$phosphoProtResidualsQTL>4,]

RMprotsFasta   <- "genomes/RM/proteins_RM.fa"
BYProtsFasta   <- "genomes/sacCer3/REFINED/proteins_S288C_R6411.fa"
RMProts   <- read.fasta(RMprotsFasta, as.string=F, forceDNAtolower=FALSE)
BYProts   <- read.fasta(BYProtsFasta, as.string=F, forceDNAtolower=FALSE)

# get gene prosition (chr coordinates)
gff <- read.table(file = "genomes/RM/genes_RM_CDSonly_refined.gtf",stringsAsFactors = F)
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


genotypeRNA <- as.matrix(genotype[,as.character(meta[colnames(eLevelBatchCorrected),"strain"])])


load("data/expressionLevel.RData")
load(file = "data/phosphoProt.RData")
load("data/phosphoLevel.RData")
load("data/proteinLevel.RData")

proteinLevelBatchCorrected<-proteinLevelBatchCorrected[,meta[colnames(proteinLevelBatchCorrected),"strain"]%in%colnames(genotype)] # rmove sample without genotype
genotypeProtein <- as.matrix(genotype[,as.character(meta[colnames(proteinLevelBatchCorrected),"strain"])])

phosphoLevelBatchCorrected<-phosphoLevelBatchCorrected[,meta[colnames(phosphoLevelBatchCorrected),"strain"]%in%colnames(genotype)] # rmove sample without genotype
genotypePhospho <-  as.matrix(genotype[,as.character(meta[colnames(phosphoLevelBatchCorrected),"strain"])])



# hotspot chr II
load("data/eQTL_results160831.RData")
bins[bins$pQTL>14 & bins$chr=="chrII",]

sel<-sapply(pQTL_results$pQTL$QtlList$FDR10, function(qtl){
  "chrII" %in% genotype$chr[qtl$mostSignificantPredictor]
})

table(unlist(sapply(pQTL_results$pQTL$QtlList$FDR10[sel],function(qtl) sapply(nrow(qtl$predictors), function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))))
genotype[186:189,1:3]
genotype[250:264,1:3]

linkedII_1 <- rowSums(pQTL_results$pQTL$qv[,186:189]<0.1) >0
linkedII_2 <- rowSums(pQTL_results$pQTL$qv[,250:264]<0.1) >0

source("scripts/GO.R")
go<-topGO(names(which(linkedII_1)|which(linkedII_1)),gene2GO = gene2BP,algorithm = "weight01", background = rownames(proteinLevelBatchCorrected))
cat.table.redmine(go[1:15,])

go<-topGO(names(which(linkedII_2)),gene2GO = gene2BP,algorithm = "weight01", background = rownames(proteinLevelBatchCorrected))
cat.table.redmine(go[1:15,])

regulators<-rownames(gene)[gene$chr=="chrII" & gene$med >351335 & gene$med<361669]
qv[regulators[regulators%in%rownames(qv)],186:189]
sapply(regulators, function(prot){
  any(RMProts[[prot]][1:(length(RMProts[[prot]])-1)]==".")
})


regulators<-rownames(gene)[gene$chr=="chrII" & gene$med >513202 & gene$med<557246]
qv[regulators[regulators%in%rownames(qv)],250:264]
sapply(regulators, function(prot){
  any(RMProts[[prot]][1:(length(RMProts[[prot]])-1)]==".")
})




# hotspot chr XII

sel<-sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10, function(qtl){
  "chrXII" %in% genotype$chr[qtl$mostSignificantPredictor]
})
table(unlist(sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10[sel],function(qtl) sapply(nrow(qtl$predictors), function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))))
genotype[2385:2394,1:3]

linkedXII <- rowSums(pQTL_results$phosphoProtResidualsQTL$qv[,2385:2394]<0.1) >0

source("scripts/GO.R")
go<-topGO(unique(phospho2prot$protein[phospho2prot$peptide%in%(names(which(linkedXII)))]),gene2GO = gene2BP,algorithm = "weight01",background = unique(phospho2prot$protein[phospho2prot$peptide %in% rownames(phosphoProtResiduals)]))
cat.table.redmine(go[1:15,])


regulators<-rownames(gene)[gene$chr=="chrXII" & gene$med >646705 & gene$med<677005]
qv[regulators[regulators%in%rownames(qv)],2385:2394]
which(rowSums(qv[regulators[regulators%in%rownames(qv)],1589:1598]<.1)>0)



# hotspot chr XVI
sel<-sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10, function(qtl){
  "chrXVI" %in% genotype$chr[qtl$mostSignificantPredictor]
})
table(unlist(sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10[sel],function(qtl) qtl$predictors[,1]:qtl$predictors[,2])))
genotype[3469:347,1:3]

linkedXVI <- rowSums(pQTL_results$phosphoLevelQTL$qv[,3469:3477]<0.1) >0

source("scripts/GO.R")
go<-topGO(unique(phospho2prot$protein[linkedXVI]),gene2GO = gene2BP,algorithm = "weight01",background = unique(phospho2prot$protein))
cat.table.redmine(go[1:20,])



sel<-sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10, function(qtl){
  "chrXVI" %in% genotype$chr[qtl$mostSignificantPredictor]
})
table(unlist(sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10[sel],function(qtl) qtl$predictors[,1]:qtl$predictors[,2])))
genotype[3470:3478,1:3]

linkedXVI <- rowSums(pQTL_results$phosphoProtResidualsQTL$qv[,3470:3478]<0.1) >0

source("scripts/GO.R")
go<-topGO(unique(phospho2prot$protein[phospho2prot$peptide%in%(names(which(linkedXVI)))]),gene2GO = gene2BP,algorithm = "weight01",background = unique(phospho2prot$protein[phospho2prot$peptide %in% rownames(phosphoProtResiduals)]))
cat.table.redmine(go[1:15,])

unique(phospho2prot$protein[phospho2prot$peptide %in% rownames(phosphoProtResiduals)])

# hotspot chr V
load("data/eQTL_results160831.RData")
sel<-sapply(QtlList$FDR10, function(qtl){
  "chrV" %in% genotype$chr[qtl$mostSignificantPredictor]
})
table(unlist(sapply(QtlList$FDR10[sel],function(qtl) sapply(nrow(qtl$predictors), function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))))
genotype[1133:1141,1:3]

linkedV <- rowSums(qv[,1132:1142]<0.1) >0

source("scripts/GO.R")
go<-topGO(names(which(linkedV)),gene2GO = gene2BP,algorithm = "weight01")
cat.table.redmine(go[1:15,])

regulators<-rownames(gene)[gene$chr=="chrV" & gene$med >490549 & gene$med<524394]
qv[regulators[regulators%in%rownames(qv)],1133:1141]
sapply(regulators, function(prot){
  any(RMProts[[prot]][1:(length(RMProts[[prot]])-1)]==".")
})

sapply(regulators, function(prot){
  sel<-BYProts[[prot]]!=RMProts[[prot]]
  cbind(BYProts[[prot]][sel], RMProts[[prot]][sel])
})

boxplot(eLevelBatchCorrected["YER165W",]~(genotypeRNA[1137,]))







# hotspot chr VIII

bins[bins$phosphoProtResidualsQTL>4 & bins$chr=="chrVIII",]
sel<-sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10, function(qtl){
  "chrVIII" %in% genotype$chr[qtl$mostSignificantPredictor]
})
table(unlist(sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10[sel],function(qtl) sapply(nrow(qtl$predictors), function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))))
genotype[1589:1598,1:3]

linkedVIII <- rowSums(pQTL_results$phosphoLevelQTL$qv[,1589:1598]<0.1) >0
unique(phospho2prot$protein[linkedVIII])

source("scripts/GO.R")
go<-topGO(unique(phospho2prot$protein[linkedVIII]),gene2GO = gene2BP,algorithm = "weight01",background = unique(phospho2prot$protein))
cat.table.redmine(go[1:20,])

regulators<-rownames(gene)[gene$chr=="chrVIII" & gene$med >77502 & gene$med<118187]
qv[regulators[regulators%in%rownames(qv)],1589:1598]
which(rowSums(qv[regulators[regulators%in%rownames(qv)],1589:1598]<.1)>0)

sapply(regulators, function(prot){
  any(RMProts[[prot]][1:(length(RMProts[[prot]])-1)]==".")
})

sapply(regulators, function(prot){
  sel<-BYProts[[prot]]!=RMProts[[prot]]
  cbind(BYProts[[prot]][sel], RMProts[[prot]][sel])
})



sel<-sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10, function(qtl){
  "chrVIII" %in% genotype$chr[qtl$mostSignificantPredictor]
})
table(unlist(sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10[sel],function(qtl) sapply(nrow(qtl$predictors), function(i){qtl$predictors[i,1]:qtl$predictors[i,2]}))))

genotype[1592:1597,1:3]

linkedVIIIResidual <- rowSums(pQTL_results$phosphoProtResidualsQTL$qv[,1592:1597]<0.1) >0
linkedVIIIResidual_prot <- (unique(phospho2prot$protein[phospho2prot$peptide%in%(names(which(linkedVIIIResidual)))]))
linkedVIIIPhospho <- rowSums(pQTL_results$phosphoLevelQTL$qv[,1589:1598]<0.1) >0
linkedVIIIPhospho_prot <- unique(phospho2prot$protein[linkedVIIIPhospho])


unique(phospho2prot$protein[phospho2prot$peptide%in%(names(which(linkedVIII)))])

source("scripts/GO.R")
go<-topGO(unique(phospho2prot$protein[phospho2prot$peptide%in%(names(which(linkedVIII)))]),gene2GO = gene2BP,algorithm = "weight01",background = unique(phospho2prot$protein[phospho2prot$peptide %in% rownames(phosphoProtResiduals)]))
cat.table.redmine(go[1:20,])

# look for PPI interactor and for phosphorilatoin
biogrid <- read.table("data/interaction_data.tab.txt", sep="\t", quote="",comment.char = "",stringsAsFactors = F)
biogrid <- biogrid[,c(1,3,6)]
biogrid[,3]<- as.factor(biogrid[,3])

STE20_interactors <- unique(c(biogrid[biogrid[,2]=="YHL007C",1], biogrid[biogrid[,1]=="YHL007C",2]))
STE20_interactors_physical <- unique(c(biogrid[biogrid[,2]=="YHL007C" & biogrid[,3]=="physical interactions", 1], 
                                       biogrid[biogrid[,1]=="YHL007C" & biogrid[,3]=="physical interactions",2]))
STE20_interactors_genetic <- unique(c(biogrid[biogrid[,2]=="YHL007C" & biogrid[,3]=="genetic interactions", 1], 
                                       biogrid[biogrid[,1]=="YHL007C" & biogrid[,3]=="genetic interactions",2]))

fisher.test(unique(phospho2prot[,2])%in%linkedVIIIPhospho_prot, 
            unique(phospho2prot[,2])%in%STE20_interactors)

fisher.test(unique(phospho2prot[,2])%in%linkedVIIIPhospho_prot, 
            unique(phospho2prot[,2])%in%STE20_interactors_physical)

fisher.test(unique(phospho2prot[,2])%in%linkedVIIIPhospho_prot, 
            unique(phospho2prot[,2])%in%STE20_interactors_genetic)


# target of ste20
gene["YHL007C",]
STE20_RNA_target <- names(which(rowSums(qv[,1590:1592] < .1)>0))
STE20_protein_target <- names(which(rowSums(pQTL_results$pQTL$qv[,1590:1592] < .1)>0))
STE20_phospho_target <- rowSums(pQTL_results$phosphoLevelQTL$qv[,1590:1592]<0.1) >0
STE20_phospho_target <- unique(phospho2prot$protein[STE20_phospho_target])
STE20_phosphoProt_target <- rowSums(pQTL_results$phosphoLevelQTL$qv[,1590:1592]<0.1) >0
STE20_phosphoProt_target <- unique(phospho2prot$prot[phospho2prot$peptide%in%names(which(STE20_phosphoProt_target))])


load("data/expressionLevel.RData")
genotypeRNA <- as.matrix(genotype[,as.character(meta[colnames(eLevelBatchCorrected),"strain"])])

RNADirection <- sapply(STE20_RNA_target,function(g){
  mean(eLevelBatchCorrected[g,round(colMeans(genotypeRNA[1590:1592,]))==genotypeRNA[1591,"BY4716"]]) > 
    mean(eLevelBatchCorrected[g,round(colMeans(genotypeRNA[1590:1592,]))==genotypeRNA[1591,"RM11-1a"]])
})

protDirection <- sapply(STE20_protein_target,function(g){
  mean(proteinLevelBatchCorrected[g,round(colMeans(genotypeProtein[1590:1592,]))==genotypeProtein[1591,"BY4716"]],na.rm=T) > 
    mean(proteinLevelBatchCorrected[g,round(colMeans(genotypeProtein[1590:1592,]))==genotypeProtein[1591,"RM11-1a"]],na.rm=T)
})


phosphoDirection<- sapply(STE20_phospho_target,function(g){
  if(length(phospho2prot$peptide[phospho2prot$protein==g])>1) {
  rowMeans(phosphoLevelBatchCorrected[phospho2prot$peptide[phospho2prot$protein==g], 
           round(colMeans(genotypePhospho[1590:1592,]))==genotypePhospho[1591,"BY4716"]],na.rm=T) >
    rowMeans(phosphoLevelBatchCorrected[phospho2prot$peptide[phospho2prot$protein==g], 
             round(colMeans(genotypePhospho[1590:1592,]))==genotypePhospho[1591,"RM11-1a"]],na.rm=T) 
  } else {
    mean(phosphoLevelBatchCorrected[phospho2prot$peptide[phospho2prot$protein==g],round(colMeans(genotypePhospho[1590:1592,]))==genotypePhospho[1591,"BY4716"]],na.rm=T) > 
      mean(phosphoLevelBatchCorrected[phospho2prot$peptide[phospho2prot$protein==g],round(colMeans(genotypePhospho[1590:1592,]))==genotypePhospho[1591,"RM11-1a"]],na.rm=T)
    }
})


pQTL_results$phosphoLevelQTL$qv[phospho2prot$peptide[phospho2prot$protein=="YKL062W"],1590:1592]


linkedVIIIResidual <- rowSums(pQTL_results$phosphoProtResidualsQTL$qv[,1592:1597]<0.1) >0
linkedVIIIResidual_prot <- (unique(phospho2prot$protein[phospho2prot$peptide%in%(names(which(linkedVIIIResidual)))]))
linkedVIIIPhospho <- rowSums(pQTL_results$phosphoLevelQTL$qv[,1589:1598]<0.1) >0
linkedVIIIPhospho_prot <- unique(phospho2prot$protein[linkedVIIIPhospho])



qv["YLR362W",1589:1598]

boden<-read.table("data/table_bodenmiller.txt", sep="\t", header=T, stringsAsFactors = F,comment.char = "",quote = "")[,c(1,2,4)]
sgd <- read.table("data/SGD_kinase_phosphatase.txt", sep="\t", header=T, stringsAsFactors = F,comment.char = "",quote = "")
colnames(sgd) <- colnames(boden)
cat(sgd[!sgd[,1]%in%boden[,1],1],sep="\n")
sgd[,1]%in%boden[,1]
boden[,1]%in%sgd[,1]
kinasePhosphatase <- rbind(boden, sgd[!sgd[,1]%in%boden[,1],][c(56,57,58,65,138),]) 
write.table(kinasePhosphatase,"data/kinase_phosphatase_yeast.tsv", sep="\t",row.names = F,quote = F)
kinasePhosphatase <- read.table("data/kinase_phosphatase_yeast.tsv", sep="\t", header=T, stringsAsFactors = F,comment.char = "",quote = "")


unique(phospho2prot$protein[linkedVIII])










# load data -------
meta <- read.table("metadata/metadata.csv", sep=",", header=T)
rownames(meta) <- meta$culture
meta    <- meta[mixedsort(rownames(meta)), ]
meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
load("data/proteinLevel.RData")
load("data/phosphoLevel.RData")
load("data/expressionLevel.RData")
sharedSamples <- mixedsort(intersect(intersect(intersect(colnames(proteinLevelBatchCorrected), colnames(phosphoLevelBatchCorrected)),rownames(meta)),colnames(eLevelBatchCorrected)))
meta <- meta[mixedsort(sharedSamples),]

prot <- as.matrix(proteinLevelBatchCorrected[,sharedSamples])
phospho <- as.matrix((phosphoLevelBatchCorrected[,sharedSamples]))
rna <- as.matrix(eLevelBatchCorrected[,sharedSamples])
sel <- intersect(rownames(prot), rownames(rna))
prot <- prot[sel,]
rna <- rna[sel,]


# Scale the data for propar comparison of RNA and prot level
rnaScale <- protScale <- matrix(NA,nrow(prot),150)
for(i in 1:nrow(prot)){ 
  p <- prot[i,]
  r <- rna[i,]
  sel <- !is.na(p) & !is.na(r)
  p <- (p - mean(p[sel])) / sd(p[sel])
  r <- (r - mean(r[sel])) / sd(r[sel])
  rnaScale[i,] <- r
  protScale[i,] <- p
}

library(MASS)
mmProtRna<-sapply(1:nrow(protScale),function(i){
  rlm( protScale[i,]~rnaScale[i,],method = "MM")$coefficients
})
mmProtRna<-t(mmProtRna)
colnames(mmProtRna) <- c("Intersecpt", "Slope")
rownames(mmProtRna) <- sel


linked<-mmProtRna[intersect(names(which(linkedV)),sel),2]
unlinked<-mmProtRna[intersect(names(which(!linkedV)),sel),2]

pdfAndPng("graph/corelationProtRnaV", 10,5,expression({
  par(mfrow=c(1,2))
  boxplot(linked,unlinked, col="grey", names=c("hotspot chrV\ntargets","non hotspot chrV \ntargets"),varwidth=T,
          ylab="Prot to RNA pearson's r")
  hist(unlinked, col="grey",breaks=20, main="", xlab="Protein-RNA resgression not hotspot targets")
  abline(v=0, col="red",lwd=2)
}))



#counts
sapply(1:nrow(bins),function(i){
  sum(QTL_genomic_pos$chr==bins$chr[i] & QTL_genomic_pos$pos > bins$start[i] & QTL_genomic_pos$pos < bins$end[i])
})

















