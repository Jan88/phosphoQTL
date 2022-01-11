library(RFQTL)
library(xlsx)
source("lib/general_function.R")

# load data #####
load("data/Nogami2007/Nogami_mappingData.RData")
info = read.xlsx("data/Nogami2007/journal.pgen.0030031.st001_curated.XLS",
                 sheetIndex=1, startRow=5, header=T,stringsAsFactors=F)
info$ID = gsub(x=info$ID, pattern="-",replacement=".")
load("data/Nogami2007/Nogami_genotypes.RData")
markerPos = read.table("data/genotype_for_mapping.tsv", header = T, as.is=T)[,1:3]
chrVec <- markerPos[,1]
corMat <- cor(genotype, use = "pair")^2
corMatSq <- corMat*corMat
group2genotype = mappingData$group2genotype
corThreshold = 0.8
#####


# get QTL on all traits #####
scores = t(sapply(1:nrow(mappingData$phenotype), function(i){
   load(paste0("data/Nogami2007/realScores/trait_", i, ".RData"))
   return(realScores)
}))
pValues <- pEst(path="data/Nogami2007/permScores/",
                scores=scores,
                markersPerIteration = 350,
                printProg = T,
                pCorrection = "fdr")
save(pValues, file = "data/Nogami2007/pValues.RData")
pValuesX <- pValues[,mappingData$genotype2group]
QTL_list <- QTLgrouper(pmat = pValuesX,
                       sigThreshold = 0.1,
                       corThreshold = corThreshold,
                       distThreshold = 9,
                       genotype = genotype,
                       chrVec = chrVec)
save(QTL_list,
     file="data/Nogami2007/QTL_list.RData")
writeQTL(QTLlist = QTL_list,
         traitNames = "",
         markerPositions = markerPositions,
         path = "data/Nogami2007/QTL_list.qtl")
# make human-readable version
QTL_list_HR = do.call(rbind,lapply(QTL_list, function(x){
  ph = rownames(mappingData$phenotype)[x$target]
  return(data.frame(targetID = x$target,
                    info[info$ID == ph,],
                    markerRegionID = x$predictors,
                    peakID = x$mostSignificantPredictor,
                    markerPos[x$predictors[,"start"],1:2],
                    end = markerPos[x$predictors[,"end"],3],
                    peak = markerPos[x$mostSignificantPredictor,2:3], 
                    minP = x$minP))
}))
save(QTL_list_HR,
     file="data/Nogami2007/QTL_list_HR.RData")
#####

#  compare with their results #####
Nogami2007 = read.xlsx("data/otherQTLstudies/Nogami2007_suplTable4.xls",
                       sheetIndex=1, header = T ,startRow=2, stringsAsFactors = F)
Nogami2007$Marker_chr = paste0("chr", as.roman(Nogami2007$Marker_chr))
Nogami2007$Trait = gsub(pattern="-", replacement=".", 
                        x=Nogami2007$Trait, fixed=T)
oursInTheirs =  apply(QTL_list_HR,1, function(x){
  sum(Nogami2007$Marker_chr == x["CHR"] & 
        Nogami2007$Trait == x["ID"] &
        Nogami2007$Marker_pos >= as.numeric(x["start"]) & 
        Nogami2007$Marker_pos <= as.numeric(x["end"]))
})
theirsInOurs = apply(Nogami2007, 1, function(x){
  sum(QTL_list_HR$CHR == x["Marker_chr"] & 
    QTL_list_HR$ID == x["Trait"] & 
    QTL_list_HR$start <= as.numeric(x["Marker_pos"]) & 
    QTL_list_HR$end >= as.numeric(x["Marker_pos"]))
})
table(theirsInOurs)
table(oursInTheirs)
######


# get the nearby genes for each gQTL ######
gff <- read.table("/cellnet/phosphoQTL/Saccharomyces_cerevisiae/sacCer3/S288C_R6411_CDSonly.gff",
                  sep="\t",comment.char = "#",quote = "",as.is=T)
colnames(gff) = c("chr", "source", "type", "start", "end", 
                  "something", "strand", "score", "name")
anno = read.table("Saccharomyces_cerevisiae/SGD_features.tab",
                  header = F,sep = "\t",quote="",as.is=T)
anno = anno[anno$V4 != "",]
rownames(anno) = anno$V4
colnames(anno) = c("SGD_ID", "type", "status", "name_systematic", 
                   "name_standard", "description", "chromosome", 
                   "ID_alternative", "chromosome_nr", "start", "end",
                   "strand", "score", "date", "old_date", "function")
inPeak = do.call(rbind,apply(QTL_list_HR,1,function(x){
  s = as.integer(x["peak.start"])
  e = as.integer(x["peak.end"])
  genes = gff[gff$chr == x["CHR"] & 
                ((gff$start >= s & gff$start <= e) |
                   (gff$end >= s & gff$end <= e) | 
                   (gff$start <= s & gff$end >= e)),"name"]
  anno[genes,c(4,5,6,16)]
}))
save(inPeak, file = "data/Nogami2007/GenesInPeak.RData")
######

