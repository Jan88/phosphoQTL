
# response to small molecule phentype data from Perlstein et a.
molecules = read.table("data/PerlsteinData/Perlstein2007_resistanceToChemicals.csv",
                       sep = ",", skip = 2, header = T, as.is = T)
compound = sapply(molecules$Compound.response,function(x){
   strsplit(x,split=" ")[[1]][1]
}, USE.NAMES=F)
molecules = molecules[compound != "DMSO",]
# take mean of compounds that were measured more than once
# in exactly the same concentration and time
double = molecules$Compound.response[duplicated(molecules$Compound.response)]
for(i in double){
   sub = as.data.frame(t(molecules[molecules$Compound.response == i,-(1:4)]))
   plot(sub, main = i)
   new = rowMeans(sub)
   rows = which(molecules$Compound.response == i)
   molecules[rows[1],5:ncol(molecules)] = new
   molecules = molecules[-rows[2],]
}
rownames(molecules) = molecules[,1]
molecules = molecules[,-1]
molecules = as.matrix(molecules)
# get only strains that are present both in our and their data
# load("data/eQtlMappingData160817.RData")
geno <- read.table("data/genotype_for_mapping.tsv",header=T)
genotype = geno[,-(1:3)]
colnames(molecules) = gsub(x=colnames(molecules), 
                           pattern="_", replacement=".")
both = intersect(colnames(molecules), colnames(genotype))
molecules = molecules[,both]
genotype = t(genotype[,both])
save(genotype, file = "data/PerlsteinData/Perlstein_genotypes.RData")

# mapping
library(RFQTL)
sampleInfo = 1:nrow(genotype)
mappingData = preMap(genotype = genotype, phenotype=molecules, 
                     sampleInfo = sampleInfo, propVar=0.6)
save(mappingData, file = "data/PerlsteinData/Perlstein_mappingData.RData")

