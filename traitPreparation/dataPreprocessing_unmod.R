## Load required libraries
library(data.table)
library(gtools)
library(reshape2)
library(xlsx)
library(sva)
library(RColorBrewer)
library(seqinr)
library(dplyr)
library(tidyr)
source("lib/PlotHclustWithColClasses.R")
source("lib/general_function.R")

## Inputs
input      <- "data/unmod/lgillet_feature_aligned_Nonenriched_targetfdr001_lowess_141231.tsv"
RTsub_iRT  <- "data/unmod/lookup_RTsubgroup_iRT_unmod.tsv"
metadat    <- "metadata/metadata.xlsx"

RMprotsFasta   <- "genomes/RM/proteins_RM.fa"
BYProtsFasta   <- "genomes/sacCer3/REFINED/proteins_S288C_R6411.fa"

## Outputs
trs        <- "data/unmod/transitions.csv"
trsmod     <- "data/unmod/transition_mod.csv"
output_trs <- "data/unmod/SWATH-QTL_unmod_trslevel"
output_mapDIA <- "proteinInference/mapDIA_transition/mapDIA_SWATH-QTL_unmod_trslevel.tsv"


## Read in data file
raw <- fread(input, sep="\t", header=TRUE)
raw <- as.data.frame(raw)


## Read in metadata file and order by culture number
meta    <- read.xlsx(metadat, sheetName="Main", row.names=1, header=TRUE,
                     as.data.frame=TRUE, stringsAsFactors=TRUE)
meta    <- meta[mixedsort(rownames(meta)), ]



### Refinement of raw data -----------------------------------------------------

## Remove rows containing decoys
rawRef <- raw[-grep("DECOY", raw$ProteinName), ]


## Remove rows containing features identified after 7000 sec
rawRef <- rawRef[rawRef$RT < 7000, ]


## Refine transition_group_id and look up iRTs and protein names including 
## subgroup information from RT-based splitting in assay library
rawRef$transition_group_id <- gsub("_run0", "", rawRef$transition_group_id)
RTiRTlookup                <- read.table(RTsub_iRT, header=TRUE)
rawRef                     <- merge(rawRef, RTiRTlookup, by="transition_group_id")
rawRef$RTsubgroup          <- gsub("Subgroup_", "Subgroup.", rawRef$RTsubgroup)
rawRef$RTsubgroup          <- gsub("_.*", "", rawRef$RTsubgroup)
rawRef$RTsubgroup          <- gsub("Subgroup\\.", "", rawRef$RTsubgroup)


## Refine ProteinName: if peptide belongs to a normal and a reverse protein, 
## remove reverse protein name from ProteinName and fix counter
for(i in 1:5){
  rawRef$ProteinName <- gsub("2(/.*)(/reverse_.*)", "1\\1", rawRef$ProteinName, perl=TRUE)
  rawRef$ProteinName <- gsub("3(/.*)(/.*)(/reverse_.*)", "2\\1\\2", rawRef$ProteinName, perl=TRUE)
  rawRef$ProteinName <- gsub("4(/.*)(/.*)(/*)(/reverse_.*)", "3\\1\\2\\3", rawRef$ProteinName, perl=TRUE)
  rawRef$ProteinName <- gsub("5(/.*)(/.*)(/*)(/*)(/reverse_.*)", "4\\1\\2\\3\\4", rawRef$ProteinName, perl=TRUE)
  rawRef$ProteinName <- gsub("6(/.*)(/.*)(/*)(/*)(/*)(/reverse_.*)", "5\\1\\2\\3\\4\\5", rawRef$ProteinName, perl=TRUE)
}
# unique(rawRef$ProteinName[grep("reverse", rawRef$ProteinName)])
# rawRef$ProteinName <- gsub("2(/.*)(/reverse_.*)", "1\\1", rawRef$ProteinName, perl=TRUE)
# rawRef$ProteinName <- gsub("3(/.*)(/.*)(/reverse_.*)", "2\\1\\2", rawRef$ProteinName, perl=TRUE)
# rawRef$ProteinName <- gsub("4(/.*)(/.*)(/reverse_.*)(/reverse_.*)", "2\\1\\2", rawRef$ProteinName, perl=TRUE)
# 


## Refine ProteinName: remove 1/ for unique proteins
rawRef$ProteinName <- gsub("^1/", "", rawRef$ProteinName)



## Remove non-unique peptides (contain a "/" in ProteinName)
## CAUTION 1: SpectraST annotates also Ile/Leu variants of a peptide as non-unique!
## CAUTION 2: If non-unique peptides are to be kept, need to refine ProteinName 
## such that protein IDs occur in alphabetical order (same two IDs can be swapped!)
rawRef <- rawRef[-grep("/", rawRef$ProteinName), ]


## Remove features with >=90 missing values (count entries in list - no NAs!)
notMissing <- tapply(rawRef$transition_group_id, rawRef$transition_group_id, length)
rawRefRed  <- rawRef[notMissing[rawRef$transition_group_id] >= 90, ]

# remove oxidesed M and all petides which sequence is the same as a peptides with oxidized M
oxiSeq <- unique(rawRefRed$FullPeptideName[grep("M\\(UniMod:35\\)", rawRefRed$FullPeptideName)])
oxiSeq <- gsub( " *\\(.*?\\) *", "", oxiSeq)
peptideSeq<- gsub( " *\\(.*?\\) *", "", rawRefRed$FullPeptideName)
rawRefRed<-rawRefRed[!(peptideSeq%in%oxiSeq),]
  
print(paste(length(unique(rawRef$ProteinName)), "proteins before filtering", sep=" "))
print(paste(length(unique(rawRefRed$ProteinName)), "proteins after filtering", sep=" "))
# "2937 proteins before filtering"
# "2773 proteins after filtering"

## Refine file names and lookup corresponding culture IDs
rawRefRed$filename <- gsub(".*/", "", rawRefRed$filename) 
rawRefRed$filename <- gsub("\\..*", "", rawRefRed$filename) 
cultureLookup      <- cbind(filename=as.character(meta$MSrun_unmod), culture=rownames(meta))
rawRefRed          <- merge(rawRefRed, cultureLookup, by="filename")


write.table(rawRefRed[, c("ProteinName", "FullPeptideName", "Charge", 
                           "RTsubgroup", "iRT", "culture", "Intensity", 
                           "aggr_Fragment_Annotation", "aggr_Peak_Area")], 
            file=trs, sep=",", row.names=FALSE, quote=FALSE)
bashcommand   <- paste("lib/precursor2trs.bash",trs,trsmod,sep=" ")
system(bashcommand)
rawRefTrs <- fread(trsmod, sep=",", header=TRUE)
rawRefTrs <- as.data.frame(rawRefTrs)

#convert OpenSWATH format to wide format 
dat <- dcast(rawRefTrs, 
             ProteinName + FullPeptideName + Charge + RTsubgroup + iRT + aggr_Fragment_Annotation ~ culture, 
             value.var="aggr_Peak_Area")

colnames(dat)[1:6] <- c("protein", "peptide", "charge", "RTsubgroup", "iRT", "fragment")

datsort       <- dat[, mixedsort(colnames(dat[, 7:ncol(dat)]))]
dat           <- cbind(dat[1:6], datsort)

##### Remove peptides that vary in BY / RM --------

RMProts   <- read.fasta(RMprotsFasta, as.string=TRUE, forceDNAtolower=FALSE)
BYProts   <- read.fasta(BYProtsFasta, as.string=TRUE, forceDNAtolower=FALSE)

peptideSeq <- dat$peptide
peptideSeq<- gsub( " *\\(.*?\\) *", "", peptideSeq) # remove text in parenthesis like [9991] "VIANVFSYFKPNMDDYC(UniMod:4)NR"
uniquePeptideSeq <- unique(peptideSeq)


inRM<-sapply(uniquePeptideSeq,function(s){
  length(grep(s,RMProts,fixed = TRUE))
})
inBY<-sapply(uniquePeptideSeq,function(s){
  length(grep(s,BYProts,fixed = TRUE))
})

sel<- peptideSeq%in%(uniquePeptideSeq[inBY==1 & inRM ==1])

dat<-dat[sel,]

write.table(dat, 
            file=paste(output_trs, ".tsv", sep=""), 
            sep="\t", row.names=FALSE, col.names=TRUE)


### Write mapDIA input files (slight reformatting of data, nothing more) -------
dat <- read.table(paste(output_trs, ".tsv", sep=""), sep="\t", header = T)
dat$fragment   <- paste(dat$fragment, dat$RTsubgroup,  sep="_")
dat$charge     <- NULL
dat$RTsubgroup <- NULL
dat <- dat[, c(1,2,4,5:ncol(dat),3)]
write.table(dat, file=output_mapDIA,sep="\t", row.names=FALSE, col.names=TRUE,quote=F)

# Run mapDIA ------

bashcommand   <- paste("cd proteinInference/mapDIA_transition/ && mapDIA input.params")
bashcommand   <- paste("cd proteinInference/mapDIA_transition2/ && mapDIA input.params")

system(bashcommand)


# protein level for --


# prtein data batch removal --------------



pLevel <- read.table("proteinInference/mapDIA_transition/protein_level.txt",sep="\t", header=T)
p <- as.matrix(pLevel[,rownames(meta)])
p[p==0] <- NA
rownames(p)<-pLevel$Protein
p<-log2(p)



batch <- meta[colnames(p),"MSbatch_unmod"]
a<-apply(p,1,function(x){
  sapply(split(x,batch),sd,na.rm=T)
})
sel1 <- colSums(is.na(a))==0
batch <- meta[colnames(p),"cultureBatch"]
a<-apply(p,1,function(x){
  sapply(split(x,batch),sd,na.rm=T)
})
sel2 <- colSums(is.na(a))==0


p<-p[sel1&sel2,]
mod <- model.matrix(~1, data=meta)

batch <- meta$MSbatch_unmod
pComBat = ComBat(dat=p, batch=batch, mod=mod, par.prior=F, prior.plots=F)
batch2 <- meta$cultureBatch
pComBatComBat = ComBat(dat=pComBat, batch=batch2, mod=mod, par.prior=F, prior.plots=F)

## verifying the influence of culture batch
a<-apply(pComBat,1,function(x){
  summary(aov(x~batch2))[[1]][1,5]
})

b<-apply(pComBat,1,function(x){
  summary(aov(sample(x)~batch2))[[1]][1,5]
})

pdfAndPng("graph/test_culture_batch",8,8, expression({
  hist(a, col=rgb(1,0,0,0.5),xlim=c(0,1), xlab="anova p.val", main="")
  hist(b, col=rgb(0,0,1,0.5), add=T)
  legend("topright",fill = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),legend = c("real","permuted"),bty="n",cex=2,xpd = T )
}))



pal <- c("#00B395", "#F035A9", "#00C35B", "#505FE7", "#90C111", "#E4AFFF", "#CB8900", "#35518B", "#F7BB78", "#007449", "#FF7A67", "#8D3252", "#784704", "#9A2911")

colBatchMS      <- pal[as.numeric(meta[colnames(p),]$MSbatch_unmod)]
colBatchCul     <- pal[as.numeric(meta[colnames(p),]$cultureBatch)]
#colCol     <- rainbow(length(unique(meta[colnames(p),]$LCcolumn_phospho)))[as.numeric(meta[colnames(p),]$LCcolumn_phospho)]


pdfAndPng(file = "graph/unmod_batch",10,10,expression({
  par(mfrow=c(3,1),mar=c(4,2,1,1))
  pTemp<-p
  colnames(pTemp)<-meta$strain
  PlotHclustWithColClasses(pTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
  title(ylab = "raw",line = -1.5)
  pComBatTemp<-pComBat
  colnames(pComBatTemp)<-meta$strain
  PlotHclustWithColClasses(pComBatTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
  title(ylab = "MS batch corrected",line = -1.5)
  pComBatTemp<-pComBatComBat
  colnames(pComBatTemp)<-paste(meta[colnames(p),]$strain, colnames(p), sep=" - ")
  PlotHclustWithColClasses(pComBatTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
  title(ylab = "MS batch + culture batch corrected",line = -1.5)
}))


proteinLevelBatchCorrected <- pComBatComBat
proteinLevelNonBatch <- p



meta    <- read.xlsx("metadata/metadata.xlsx", sheetName="Main", row.names=1, header=TRUE,as.data.frame=TRUE, stringsAsFactors=F,na.strings = "NA")
meta$id <- rownames(meta) 
meta[meta=="NA"]<-NA

meta <- meta[!is.na(meta$MSrun_unmod) & meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]

# remove problemaitic samples
proteinLevelNonBatch <- proteinLevelNonBatch[,meta$id]
proteinLevelBatchCorrected <- proteinLevelBatchCorrected[,meta$id]
save(proteinLevelNonBatch,proteinLevelBatchCorrected, file="data/proteinLevel.RData")
write.table(x=proteinLevelBatchCorrected, file="data/protein_level.tsv", sep="\t", quote = F,row.names = T)

write.table(x=phosphoLevelBatchCorrected, file="data/phospho_level.tsv", sep="\t", quote = F,row.names = T)

