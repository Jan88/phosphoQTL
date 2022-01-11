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
input      <- "ste20Validation/proteinData/E2001080903_feature_alignment.tsv"
RTsub_iRT  <- "data/phospho/lookup_RTsubgroup_iRT_phospho.tsv"
metadat    <- "ste20Validation/proteinData/SampleAnnotation.txt"

RMprotsFasta   <- "Saccharomyces_cerevisiae/RM/proteins_RM.fa"
BYProtsFasta   <- "Saccharomyces_cerevisiae/sacCer3/REFINED/proteins_S288C_R6411.fa"

## Outputs
trs        <- "ste20Validation/proteinData/phosphoNoCor/transitions.csv"
trsmod     <- "ste20Validation/proteinData/phosphoNoCor/transition_mod.csv"
output_trs <- "ste20Validation/proteinData/phosphoNoCor/trslevel"
output_mapDIA <- "ste20Validation/proteinData/phosphoNoCor/mapDIA/mapDIA_trslevel.tsv"


## Read in data file
raw <- fread(input, sep="\t", header=TRUE)
raw <- as.data.frame(raw)


## Read in metadata file and order by culture number
meta    <- read.table(metadat,header=T,as.is=T)
rownames(meta) <- meta$CondRep
meta    <- meta[mixedsort(rownames(meta)), ]



### Refinement of raw data -----------------------------------------------------
#old <- fread("data/phospho/lgillet_feature_aligned_Phospho_targetfdr001_lowess_141227.tsv", sep="\t", header=TRUE)

## Remove rows containing features identified after 6000 sec
rawRef <- raw[raw$RT < 6000, ]

## Refine file names and remove runs which are of bad quality
rawRef$filename <- gsub(".*/", "", rawRef$filename) 
rawRef$filename <- gsub("\\..*", "", rawRef$filename) 
#rawRef          <- rawRef[rawRef$filename %in% as.character(meta$MSrun_phospho), ]


## Refine transition_group_id and look up iRTs and protein names including 
## subgroup information from RT-based splitting in assay library
#rawRef$collapsed_peptide_cluster <- gsub("_run0", "", rawRef$collapsed_peptide_cluster)
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
rawRef <- rawRef[!grepl("/", rawRef$ProteinName), ]


## Remove features with >=6 missing values (need to count entries in list - no NAs!)
notMissing <- tapply(rawRef$transition_group_id, rawRef$transition_group_id, length) #same as table(rawRef$transition_group_id)
rawRefRed  <- rawRef[notMissing[rawRef$transition_group_id] >= 6, ]
#rawRefRed <- rawRefRed[!grepl("(0P)",rawRefRed$transition_group_id,fixed=T),]
rawRefRed <- rawRefRed[!grepl("UniMod:35",rawRefRed$transition_group_id,fixed=T),]
rawRefRed <- rawRefRed[grepl("UniMod:21",rawRefRed$transition_group_id,fixed=T),]
print(paste(length(unique(rawRef$FullPeptideName)), "peptides before filtering", sep=" "))
print(paste(length(unique(rawRefRed$FullPeptideName)), "peptides after filtering", sep=" "))
# "8865 peptides before filtering"
# "6873 peptides after filtering"



## Refine file names and lookup corresponding culture IDs
rawRefRed$filename <- gsub(".*/", "", rawRefRed$filename) 
rawRefRed$filename <- gsub("\\..*", "", rawRefRed$filename) 
cultureLookup      <- cbind(filename=as.character(meta$FileName.Phospho), culture=rownames(meta))
rawRefRed          <- merge(rawRefRed, cultureLookup, by="filename")

temp <- rawRefRed[, c("ProteinName", "FullPeptideName", "Charge", 
                      "RTsubgroup", "iRT", "culture", "Intensity", 
                      "aggr_Fragment_Annotation", "aggr_Peak_Area")]
#colnames(temp)[2] <- "FullPeptideName"

write.table(temp,file=trs, sep=",", row.names=FALSE, quote=FALSE)

bashcommand   <- paste("lib/precursor2trs_no_extra_columns2.bash",trs,trsmod,sep=" ")
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
peptideSeq<- gsub( ".*_", "", peptideSeq)
uniquePeptideSeq <- unique(peptideSeq)

#sel<- peptideSeq%in%(uniquePeptideSeq[inBY==1 & inRM ==1])

#dat<-dat[sel,]

pepsTotest <- unique(peptideSeq[dat$protein%in%c("YHL007C","YHR005C")])

inRM<-sapply(pepsTotest,function(s){
  length(grep(s,RMProts,fixed = TRUE))
})
inBY<-sapply(pepsTotest,function(s){
  length(grep(s,BYProts,fixed = TRUE))
})
all(inRM==1)&all(inBY==1) #one polymorph peptide
pepsRem <- pepsTotest[which(inRM!=1|inBY!=1)]

dat <- dat[!peptideSeq%in%pepsRem,]

# inRM<-sapply(uniquePeptideSeq,function(s){
#   length(grep(s,RMProts,fixed = TRUE))
# })
# inBY<-sapply(uniquePeptideSeq,function(s){
#   length(grep(s,BYProts,fixed = TRUE))
# })
# 
# sel<- peptideSeq%in%(uniquePeptideSeq[inBY==1 & inRM ==1])
# sum(sel) #[1] 26592
# length(sel) #[1] 27252
# 
# dat<-dat[sel,]

write.table(dat, file=paste(output_trs, ".tsv", sep=""), 
            sep="\t", row.names=FALSE, col.names=TRUE)


### Write mapDIA input files (slight reformatting of data, nothing more) -------
dat <- read.table(paste(output_trs, ".tsv", sep=""), sep="\t", header = T,as.is=T)

# dat$fragment   <- paste(dat$fragment, dat$RTsubgroup,  sep="_")
# dat$charge     <- NULL
# dat$RTsubgroup <- NULL
# dat <- dat[, c(1,2,4,5:ncol(dat),3)]
# dat$protein <- paste(dat$protein, dat$peptide, sep="|")
#count phospho modifications per peptide
phMod <- "(UniMod:21)"
nPhospho <- sapply(dat$peptide,FUN=function(pep){
  origLength <- nchar(pep)
  newLength <- nchar(gsub(phMod,"",pep,fixed = T))
  out <- (origLength-newLength)/nchar(phMod)
  return(out)
})
strippedPep <- gsub(phMod,"",dat$peptide,fixed = T)
newPepName <- paste0(strippedPep,"(",nPhospho,"P)")
dat$peptide <- newPepName
dat$charge     <- NULL
dat$RTsubgroup <- NULL
dat <- dat[, c(1,2,4,5:ncol(dat),3)]
write.table(dat, file=output_mapDIA,sep="\t", row.names=FALSE, col.names=TRUE,quote=F)

# Run mapDIA ------
bashcommand   <- paste("cd ste20Validation/proteinData/phosphoNoCor/mapDIA/ && ./mapDIA input.params")
system(bashcommand)


pLevel <- read.table("ste20Validation/proteinData/phosphoNoCor/mapDIA/peptide_level.txt",sep="\t", header=T,stringsAsFactors = F)
p <- as.matrix(pLevel[,grepl("X",colnames(pLevel))])
p[p==0] <- NA
rownames(p)<-pLevel$Peptide
p<-log2(p)

# removing the problematic peptide (no variance in on of the subgroup)
batch <- gsub("^.*_","",colnames(p))
a<-apply(p,1,function(x){
  sapply(split(x,batch),sd,na.rm=T)
})
sel1 <- colSums(is.na(a))==0 #variance in each batch
# batch <- meta[colnames(p),"cultureBatch"]
# a<-apply(p,1,function(x){
#   sapply(split(x,batch),sd,na.rm=T)
# })
# sel2 <- colSums(is.na(a))==0

p<-p[sel1,]

mod <- model.matrix(~1, data=meta[colnames(p),])

#batch <- meta[colnames(p),"MSbatch_phospho"]
pComBat = ComBat(dat=p, batch=batch, mod=mod, par.prior=F, prior.plots=F)
# batch2 <- meta[colnames(p),"cultureBatch"]
# pComBatComBat = ComBat(dat=pComBat, batch=batch2, mod=mod, par.prior=F, prior.plots=F)

# pal <- c("#00B395", "#F035A9", "#00C35B", "#505FE7", "#90C111", "#E4AFFF", "#CB8900", "#35518B", "#F7BB78", "#007449", "#FF7A67", "#8D3252", "#784704", "#9A2911")
# 
# colBatchMS      <- pal[as.numeric(meta[colnames(p),]$MSbatch_phospho)]
# colBatchCul     <- pal[as.numeric(meta[colnames(p),]$cultureBatch)]
# colCol     <- rainbow(length(unique(meta[colnames(p),]$LCcolumn_phospho)))[as.numeric(meta[colnames(p),]$LCcolumn_phospho)]
# 
# 
# pdfAndPng(file = "graph/phospho_batch",10,10,expression({
#   par(mfrow=c(3,1),mar=c(4,2,1,1))
#   pTemp<-p
#   colnames(pTemp)<-meta[colnames(p),]$strain
#   PlotHclustWithColClasses(pTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
#   title(ylab = "raw",line = -1.5)
#   pComBatTemp<-pComBat
#   colnames(pComBatTemp)<-meta[colnames(p),]$strain
#   PlotHclustWithColClasses(pComBatTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
#   title(ylab = "MS batch corrected",line = -1.5)
#   pComBatTemp<-pComBatComBat
#   colnames(pComBatTemp)<-paste(meta[colnames(p),]$strain, colnames(p), sep=" - ")
#   PlotHclustWithColClasses(pComBatTemp,rbind(colBatchMS,colBatchCul),classLabel = c("MS","Culture"))
#   title(ylab = "MS batch + culture batch corrected",line = -1.5)
# }))
# 
# 
# # test the batch effects
# 
# batchMS <- as.character(meta[colnames(p),]$MSbatch_phospho)
# batchCul  <- as.character(meta[colnames(p),]$cultureBatch)
# 
# ms<-apply(as.matrix(p),1,function(x){
#   summary(aov(x~batchMS))[[1]][1,5]
# })
# 
# batchCul <- sample(batchCul)
# cul<-apply(as.matrix(p),1,function(x){
#   summary(aov(x~batchCul))[[1]][1,5]
# })
# 
# 
# phosphoLevelBatchCorrected <- pComBatComBat
# phosphoLevelNonBatch <- p
# 
# 
# meta <- read.table("metadata/metadata.csv", sep=",", header=T)
# rownames(meta) <- meta$culture
# meta    <- meta[mixedsort(rownames(meta)), ]
# 
# 
# meta <- meta[!is.na(meta$MSrun_phospho) & meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
# #meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
# 
# # remove problemaitic samples
# phosphoLevelNonBatch <- phosphoLevelNonBatch[,rownames(meta)]
# phosphoLevelBatchCorrected <- phosphoLevelBatchCorrected[,rownames(meta)]

phosphoLevel <- pComBat
# phospho2ProtRep <- t(sapply(rownames(phosphoLevel),FUN=function(pep){
#   #pep <- strsplit(pep,"_")[[1]][2]
#   pos <- grep(pep,raw$FullPeptideName,fixed=T)[1]
#   prot <- gsub(".*/","",raw$ProteinName[pos])
#   return(c(pep,prot))
# }))
# rownames(phospho2ProtRep) <- rownames(phosphoLevel) <- phospho2ProtRep[,1]
phospho2ProtRep <- pLevel[,2:1]
rownames(phospho2ProtRep) <- phospho2ProtRep[,1]
phospho2ProtRep <- phospho2ProtRep[rownames(phosphoLevel),]

save(phosphoLevel,phospho2ProtRep,file="data/phosphoLevelReplacementStrainsNoCor.RData")
# 
# prot<-sapply(rownames(phosphoLevelBatchCorrected), function(x){
#   pLevel$Protein[which(pLevel[,2]==x)[1]]
# })
# prot <- sub("\\|.*","",prot)
# phospho2prot <- data.frame(peptide=rownames(phosphoLevelBatchCorrected), protein= prot,stringsAsFactors=F)
# rownames(phospho2prot) <-NULL
# save(phosphoLevelNonBatch,phosphoLevelBatchCorrected,phospho2prot, file="data/phosphoLevelReplacementStrains.RData")

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



###test phospho between strains###
comparisons <- rbind(c("X1473","X1475"),c("X1473","X1476"),c("X1473","X1477"),c("X1475","X1477"),c("X1476","X1477"))
phosphoLevel <- phosphoLevel[,colnames(phosphoLevel)!="X1473_r3"]
diffPhospho <- lapply(1:nrow(comparisons),FUN=function(i){
  s1 <- comparisons[i,1]
  s2 <- comparisons[i,2]
  fc <- rowMeans(phosphoLevel[,grepl(s1,colnames(phosphoLevel))],na.rm=T)-rowMeans(phosphoLevel[,grepl(s2,colnames(phosphoLevel))],na.rm=T)
  pv <- sapply(rownames(phosphoLevel),FUN=function(p){
    s1Samples <- grepl(s1,colnames(phosphoLevel))
    s2Samples <- grepl(s2,colnames(phosphoLevel))
    if(sum(is.na(phosphoLevel[p,s1Samples]))>1|sum(is.na(phosphoLevel[p,s2Samples]))>1){return(NA)}
    t.test(phosphoLevel[p,s1Samples],phosphoLevel[p,s2Samples])$p.value
  })
  fdr <- p.adjust(pv,method="BH")
  out <- cbind(l2FC=fc,p=pv,fdr=fdr)
  rownames(out) <- rownames(phosphoLevel)
  return(out)
})
names(diffPhospho) <- apply(comparisons,1,paste,collapse="_")
save(diffPhospho,comparisons,file="data/diffPhosphoReplacementStrainsNoCor.RData")