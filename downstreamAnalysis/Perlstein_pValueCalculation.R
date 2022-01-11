library(RFQTL)
source("lib/general_function.R")

# load some general data #####
load("data/PerlsteinData/Perlstein_mappingData.RData")
conditions = rownames(mappingData$phenotype)
temp = sapply(conditions, function(x){
  x = strsplit(x,split=" ")[[1]]
  if(length(x) == 5){
    return(c(x[1], paste(x[2:5], collapse = "_")))
  } else if(length(x) == 6){
    return(c(paste(x[1:2], collapse = " "), paste(x[3:6], collapse = "_")))
  } else
    return(x)
})
compounds = temp[1,]
exposures = temp[2,]
load("data/PerlsteinData/Perlstein_genotypes.RData")
geno = read.table("data/genotype_for_mapping.tsv", header = T, as.is=T)
chrVec <- geno[,1]
corMat <- cor(genotype, use = "pair")^2
corMatSq <- corMat*corMat
markerPositions <- geno[,1:3]
group2genotype = mappingData$group2genotype
corThreshold = 0.8
#####


# get p-Values #####
scores = t(sapply(1:307, function(i){
   load(paste0("data/PerlsteinData/realScores/trait_", i, ".RData"))
   return(realScores)
}))
pValues <- pEst(path="data/PerlsteinData/permScores/",
                scores=scores,
                markersPerIteration = 350,
                printProg = T,
                pCorrection = "fdr")
save(pValues, file = "data/PerlsteinData/pValues.RData")
#####
load("data/PerlsteinData/pValues.RData")

# merge QTLs #####
# merge QTLS that are correlated or close together
pValuesX <- pValues[,mappingData$genotype2group]
QTL_list <- QTLgrouper(pmat = pValuesX,
                       sigThreshold = 0.1,
                       corThreshold = corThreshold,
                       distThreshold = 9,
                       genotype = genotype,
                       chrVec = chrVec)
# merge QTLs with the same compound (but in different concentrations)
QTLcompounds = compounds[sapply(QTL_list, function(x){x$target})]
corThresholdSq = corThreshold^2
QTL_list_merged = lapply(unique(QTLcompounds), function(compound){
   qtls = QTL_list[QTLcompounds == compound]
   loci = do.call(rbind, lapply(qtls, function(x){
      cbind(x$predictors, x$target, x$mostSignificantPredictor, x$minP)
   }))
   colnames(loci) = c("start", "end", "targets", 
                      "mostSignificantPredictor", "minP")
   if(nrow(loci) == 1){
      return(loci)
   }
   loci = loci[order(loci[,1], loci[,2]),]
   col = loci[1,1:2]
   tar = loci[1,3]
   Peak = loci[1,4]
   P = loci[1,5]
   outMat = NULL
   for (i in 2:nrow(loci)) {
      sofar = col[1]:col[2]
      new = loci[i,1]:loci[i,2]
      sofarChr = unique(geno[col[1]:col[2],1])
      newChr = unique(geno[loci[i,1]:loci[i,2],1])
      sofarEnd = geno[col[2],3]
      newStart = geno[loci[i,1],2]
      if(sofarChr == newChr & 
         (newStart - sofarEnd < 50000 | any(corMatSq[sofar,new] > corThresholdSq))){
         col = c(min(c(sofar,new)), 
                 max(c(sofar,new)))
         if(loci[i,5]<P){
           Peak = loci[i,4]
           P = loci[i,5]
         }
         tar = sort(unique(c(tar,loci[i,3])))
      } else{
         outMat <- rbind(outMat, 
                         c(col[1:2],paste(tar, collapse=","),Peak,P))
         tar = loci[i,3]
         col = loci[i,1:2]
         Peak = loci[i,4]
         P = loci[i,5]
      }
   }
   outMat <- rbind(outMat, 
                   c(col[1:2],paste(tar, collapse=","),Peak,P))
   colnames(outMat) = c("start", "end", "targets", 
                           "mostSignificantPredictor", "minP")
   return(outMat)
})
names(QTL_list_merged)= unique(QTLcompounds)
save(QTL_list_merged, 
     file="data/PerlsteinData/Perlstein_QTL_list_merged.RData")

# Make human readable version, i.e., translate marker and 
# phenotype IDs into positions and phenotype descriptions  
QTL_list_merged_HR = lapply(names(QTL_list_merged), function(compound){
   x = QTL_list_merged[[compound]]
   pos = cbind(geno[as.numeric(x[,1]),1:2], 
               end = geno[as.numeric(x[,2]),3])
   peakPos = cbind(geno[as.numeric(x[,4]),2:3])
   t = strsplit(as.character(x[,3]), split=",")
   pheno = sapply(t, function(y){
      paste(unname(exposures[as.numeric(y)]), collapse = ",")
   })
   res = cbind(pos, peakPos, x[,5], pheno, compound)
   colnames(res) = c("chr", "start", "end", "peakStart", "peakEnd",
                     "minP", "exposure", "compound")
   return(res)
})
names(QTL_list_merged_HR) = names(QTL_list_merged)
QTL_list_merged_HR_table = do.call(rbind,QTL_list_merged_HR)
QTL_list_merged_HR_table$exposure = as.character(QTL_list_merged_HR_table$exposure)
QTL_list_merged_HR_table$compound = as.character(QTL_list_merged_HR_table$compound)
QTL_list_merged_HR_table$minP = as.numeric(as.character(QTL_list_merged_HR_table$minP))
save(QTL_list_merged_HR, QTL_list_merged_HR_table,
     file="data/PerlsteinData/Perlstein_QTL_list_merged_HR.RData")
#####

# compare with the results from Perlstein #####
load("data/PerlsteinData/Perlstein_QTL_list_merged_HR.RData")
theirQTLs = read.table("data/PerlsteinData/Perlstein_SupFile5_collapsedCI_2.csv", 
                     header = T, sep = ",", as.is = T)
theirCompounds = sapply(theirQTLs$FullCompoundName, function(x){
  x = strsplit(x,split=" ")[[1]]
  if(length(x) <= 5){
    return(x[1])
  } else if(length(x) == 6){
    return(paste(x[1:2], collapse = " "))
  } else if(length(x) == 7){
    return(paste(x[1:3], collapse = " "))
  }
})
theirCompounds[theirCompounds == "DFI"] = "diphenyleneiodonium"
theirQTLs$compound = theirCompounds
theirQTLs$chrom = paste0("chr", as.roman(theirQTLs$chrom))
c1 = sort(c("benzethonium chloride", "clomitrazole", "fendiline", 
            "resveratrol",  "SK&F 96366",  "SK&F 96367", 
            "trifluoperazine", "trimeprazine", "sodium chloride",
            "thiram", unique(theirCompounds)))
c2 = sort(c("chlorpromazine", "SK&F 96365", "tamoxifen",
            unique(QTL_list_merged_HR$compound)))
compoundLegend = cbind(theirs = c1, ours = c2)
oursInPerlstein = sapply(1:nrow(QTL_list_merged_HR), function(i){
  x = unlist(QTL_list_merged_HR[i,])
  s = as.integer(x["start"])
  e = as.integer(x["end"])
  theircomp = compoundLegend[compoundLegend[,"ours"] == x["compound"],"theirs"]
  theirs = theirQTLs[theirQTLs$compound == theircomp,]
  dists = abs(c(s-theirs$ci_start,s-theirs$ci_end., 
                e-theirs$ci_start,e-theirs$ci_end.) )
  sum(theirs$chrom == x["chr"] & 
        ((theirs$ci_start > s & theirs$ci_start < e) |
           (theirs$ci_end. > s & theirs$ci_end. < e) | 
           (theirs$ci_start < s & theirs$ci_end. > e) | 
           any(dists < 10000)))
})
PerlsteinInOurs = apply(theirQTLs, 1, function(x){
  ourcomp = compoundLegend[compoundLegend[,"theirs"] == x["compound"],"ours"]
  ours = QTL_list_merged_HR[QTL_list_merged_HR$compound == ourcomp,]
  s = as.numeric(x["ci_start"])
  e = as.numeric(x["ci_end."])
  dists = abs(c(s-ours$start,s-ours$end, 
                e-ours$start,e-ours$end) )
  sum(ours$chr == x["chrom"] & 
        ((ours$start > s & ours$start < e) |
           (ours$end > s & ours$end < e) | 
           (ours$start < s & ours$end > e) | 
           any(dists < 10000)))
})
#####



# get the nearby genes for each gQTL ####
load("data/PerlsteinData/Perlstein_QTL_list_merged_HR.RData") #QTL_list_merged_HR
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
inPeak = lapply(QTL_list_merged_HR,function(x){
  res = unique(do.call(rbind,apply(x,1,function(y){
    s = as.integer(y["peakStart"])
    e = as.integer(y["peakEnd"])
    genes = gff[gff$chr == y["chr"] & 
          ((gff$start > s & gff$start < e) |
             (gff$end > s & gff$end < e) | 
             (gff$start < s & gff$end > e)),"name"]
    anno[genes,c(4,5,6,16)]
  })))
  # compound = c(unique(as.character(x$compound)), rep("",nrow(res)-1))
  # cbind(compound, res)
})
inPeakTable = do.call(rbind,lapply(names(inPeak), function(compound){
  res = inPeak[[compound]]
  cbind(compound = c(compound, rep("",nrow(res)-1)), res)
}))
save(inPeak, file = "data/PerlsteinData/GenesInPeak.RData")
write.table(inPeakTable, file="data/PerlsteinData/GenesInPeak.tsv",
            sep = "\t", row.names=F, col.names=T)
#####
