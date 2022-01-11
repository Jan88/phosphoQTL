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
input      <- "data/phospho/final_out_300616.tsv"

metadat    <- "metadata/metadata.xlsx"
RTsub_iRT  <- "data/phospho/lookup_RTsubgroup_iRT_phospho.tsv"

RMprotsFasta   <- "genomes/RM/proteins_RM.fa"
BYProtsFasta   <- "genomes/sacCer3/REFINED/proteins_S288C_R6411.fa"

## Outputs
trs        <- "data/phospho/transitions.csv"
trsmod     <- "data/phospho/transition_mod.csv"
output_trs <- "data/phospho/SWATH-QTL_unmod_trslevel"
output_mapDIA <- "phosphoInference/mapDIA_SWATH-QTL_phospho_trslevel.tsv"


## Read in data file
raw <- fread(input, sep="\t", header=TRUE)
raw <- as.data.frame(raw)


## Read in metadata file and order by culture number
meta    <- read.xlsx(metadat, sheetName="Main", row.names=1, header=TRUE,
                     as.data.frame=TRUE, stringsAsFactors=TRUE)
meta    <- meta[mixedsort(rownames(meta)), ]

meta <- read.table("metadata/metadata.csv", sep=",", header=T)
rownames(meta) <- meta$culture
meta    <- meta[mixedsort(rownames(meta)), ]



### Refinement of raw data -----------------------------------------------------
#old <- fread("data/phospho/lgillet_feature_aligned_Phospho_targetfdr001_lowess_141227.tsv", sep="\t", header=TRUE)

## Remove rows containing features identified after 6000 sec
rawRef <- raw[raw$RT < 6000, ]

## Refine file names and remove runs which are of bad quality
rawRef$filename <- gsub(".*/", "", rawRef$filename) 
rawRef$filename <- gsub("\\..*", "", rawRef$filename) 
rawRef          <- rawRef[rawRef$filename %in% as.character(meta$MSrun_phospho), ]


## Refine transition_group_id and look up iRTs and protein names including 
## subgroup information from RT-based splitting in assay library
rawRef$collapsed_peptide_cluster <- gsub("_run0", "", rawRef$collapsed_peptide_cluster)
rawRef$transition_group_id <- gsub("_run0", "", rawRef$transition_group_id)
RTiRTlookup                <- read.table(RTsub_iRT, header=TRUE,stringsAsFactors = F)
rawRef                     <- merge(rawRef, RTiRTlookup, by="transition_group_id")
rawRef$RTsubgroup          <- gsub("Subgroup_", "Subgroup.", rawRef$RTsubgroup)
rawRef$RTsubgroup          <- gsub("_.*", "", rawRef$RTsubgroup)
rawRef$RTsubgroup          <- gsub("Subgroup\\.", "", rawRef$RTsubgroup)

## Refine ProteinName: remove subgroup information and 1/ for unique proteins
rawRef$ProteinName <- gsub(".*_", "", rawRef$ProteinName)
rawRef$ProteinName <- gsub("^1/", "", rawRef$ProteinName)

## Remove non-unique peptides (contain a "/" in ProteinName)
## CAUTION 1: SpectraST annotates also Ile/Leu variants of a peptide as non-unique!
## CAUTION 2: If non-unique peptides are to be kept, need to refine ProteinName 
## such that protein IDs occur in alphabetical order (same two IDs can be swapped!)
rawRef <- rawRef[-grep("/", rawRef$ProteinName), ]


## Remove features with >=90 missing values (need to count entries in list - no NAs!)
notMissing <- tapply(rawRef$collapsed_peptide_cluster, rawRef$collapsed_peptide_cluster, length)
rawRefRed  <- rawRef[notMissing[rawRef$collapsed_peptide_cluster] >= 100, ]
rawRefRed <- rawRefRed[-grep("(0P)",rawRefRed$collapsed_peptide_cluster,fixed=T),]
rawRefRed <- rawRefRed[-grep("UniMod:35",rawRefRed$collapsed_peptide_cluster,fixed=T),]

print(paste(length(unique(rawRef$FullPeptideName)), "peptides before filtering", sep=" "))
print(paste(length(unique(rawRefRed$FullPeptideName)), "peptides after filtering", sep=" "))
# "5020 peptides before filtering"
# "4039 peptides after filtering"



## Refine file names and lookup corresponding culture IDs
rawRefRed$filename <- gsub(".*/", "", rawRefRed$filename) 
rawRefRed$filename <- gsub("\\..*", "", rawRefRed$filename) 
cultureLookup      <- cbind(filename=as.character(meta$MSrun_phospho), culture=rownames(meta))
rawRefRed          <- merge(rawRefRed, cultureLookup, by="filename")

temp <- rawRefRed[, c("ProteinName", "collapsed_peptide_cluster", "Charge", 
                      "RTsubgroup", "iRT", "culture", "Intensity", 
                      "aggr_Fragment_Annotation", "aggr_Peak_Area")]
colnames(temp)[2] <- "FullPeptideName"

write.table(temp, 
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
peptideSeq<- gsub( " *\\[.*?\\] *", "", peptideSeq) # remove text in braket like [9991] "VIANVFSYFKPNMDDYC(UniMod:4)NR"
peptideSeq<- gsub( "_[0-9]", "", peptideSeq) 
uniquePeptideSeq <- unique(peptideSeq)

inRM<-sapply(uniquePeptideSeq,function(s){
  length(grep(s,RMProts,fixed = TRUE))
})
inBY<-sapply(uniquePeptideSeq,function(s){
  length(grep(s,BYProts,fixed = TRUE))
})

sel<- peptideSeq%in%(uniquePeptideSeq[inBY==1 & inRM ==1])
sum(sel) #[1] 26592
length(sel) #[1] 27252

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
dat$protein <- paste(dat$protein, dat$peptide, sep="|")

write.table(dat, file=output_mapDIA,sep="\t", row.names=FALSE, col.names=TRUE,quote=F)

# Run mapDIA ------
bashcommand   <- paste("cd phosphoInference/ && mapDIA input.params")
system(bashcommand)


pLevel <- read.table("phosphoInference/peptide_level.txt",sep="\t", header=T,stringsAsFactors = F)
p <- as.matrix(pLevel[,intersect(rownames(meta[!is.na(meta$MSbatch_phospho),]),colnames(pLevel))])
p[p==0] <- NA
rownames(p)<-pLevel$Peptide
p<-log2(p)

# removing the problematic ppetide (no variance in on of the subgroup)
batch <- meta[colnames(p),"MSbatch_phospho"]
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

mod <- model.matrix(~1, data=meta[colnames(p),])

batch <- meta[colnames(p),"MSbatch_phospho"]
pComBat = ComBat(dat=p, batch=batch, mod=mod, par.prior=F, prior.plots=F)
batch2 <- meta[colnames(p),"cultureBatch"]
pComBatComBat = ComBat(dat=pComBat, batch=batch2, mod=mod, par.prior=F, prior.plots=F)

pal <- c("#00B395", "#F035A9", "#00C35B", "#505FE7", "#90C111", "#E4AFFF", "#CB8900", "#35518B", "#F7BB78", "#007449", "#FF7A67", "#8D3252", "#784704", "#9A2911")

colBatchMS      <- pal[as.numeric(meta[colnames(p),]$MSbatch_phospho)]
colBatchCul     <- pal[as.numeric(meta[colnames(p),]$cultureBatch)]
colCol     <- rainbow(length(unique(meta[colnames(p),]$LCcolumn_phospho)))[as.numeric(meta[colnames(p),]$LCcolumn_phospho)]


pdfAndPng(file = "graph/phospho_batch",10,10,expression({
  par(mfrow=c(3,1),mar=c(4,2,1,1))
  pTemp<-p
  colnames(pTemp)<-meta[colnames(p),]$strain
  PlotHclustWithColClasses(pTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
  title(ylab = "raw",line = -1.5)
  pComBatTemp<-pComBat
  colnames(pComBatTemp)<-meta[colnames(p),]$strain
  PlotHclustWithColClasses(pComBatTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
  title(ylab = "MS batch corrected",line = -1.5)
  pComBatTemp<-pComBatComBat
  colnames(pComBatTemp)<-paste(meta[colnames(p),]$strain, colnames(p), sep=" - ")
  PlotHclustWithColClasses(pComBatTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
  title(ylab = "MS batch + culture batch corrected",line = -1.5)
}))


# test the batch effects

batchMS <- as.character(meta[colnames(p),]$MSbatch_phospho)
batchCul  <- as.character(meta[colnames(p),]$cultureBatch)

ms<-apply(as.matrix(p),1,function(x){
  summary(aov(x~batchMS))[[1]][1,5]
})

batchCul <- sample(batchCul)
cul<-apply(as.matrix(p),1,function(x){
  summary(aov(x~batchCul))[[1]][1,5]
})


phosphoLevelBatchCorrected <- pComBatComBat
phosphoLevelNonBatch <- p


meta <- read.table("metadata/metadata.csv", sep=",", header=T)
rownames(meta) <- meta$culture
meta    <- meta[mixedsort(rownames(meta)), ]


meta <- meta[!is.na(meta$MSrun_phospho) & meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
#meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]

# remove problemaitic samples
phosphoLevelNonBatch <- phosphoLevelNonBatch[,rownames(meta)]
phosphoLevelBatchCorrected <- phosphoLevelBatchCorrected[,rownames(meta)]

prot<-sapply(rownames(phosphoLevelBatchCorrected), function(x){
  pLevel$Protein[which(pLevel[,2]==x)[1]]
})
prot <- sub("\\|.*","",prot)
phospho2prot <- data.frame(peptide=rownames(phosphoLevelBatchCorrected), protein= prot,stringsAsFactors=F)
rownames(phospho2prot) <-NULL
save(phosphoLevelNonBatch,phosphoLevelBatchCorrected,phospho2prot, file="data/phosphoLevel.RData")

# 
# load("data/phosphoLevel_old.RData")
# 
# non_phospho <- grepl("(0P)",phospho2prot$peptide,fixed = T)
# 
# oxi<-grep("UniMod:35",phospho2prot[,1])
# peptideSeq<- phospho2prot$peptide
# #peptideSeq <- gsub( " *\\(.*?\\) *", "", phospho2prot$peptide)
# peptideSeq<- gsub( " *\\[.*?\\] *", "", peptideSeq) # remove text in braket like [9991] "VIANVFSYFKPNMDDYC(UniMod:4)NR"
# #peptideSeq<- gsub( "_[0-9]", "", peptideSeq) 
# oxi<-sapply(peptideSeq[oxi],function(n) which(peptideSeq==n))
# oxi <- unlist(oxi)
# 
# sel <- !non_phospho
# sel[oxi] <- F
# 
# phosphoLevelNonBatch <- phosphoLevelNonBatch[sel,]
# phosphoLevelBatchCorrected <- phosphoLevelBatchCorrected[sel,]
# phospho2prot <- phospho2prot[sel,]
# 
# save(phosphoLevelNonBatch,phosphoLevelBatchCorrected,phospho2prot, file="data/phosphoLevel.RData")
# 
# 
# dim(phosphoLevelBatchCorrected)
# 
# 
# 
