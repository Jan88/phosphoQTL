args <- commandArgs(trailingOnly = TRUE)
nTrait <- as.numeric(args)
# library(RandomForestExtended, lib.loc="/home/cschmal1/R/3.1.1")
library(randomForest) #, lib.loc="/home/cschmal1/R/3.1.1"
library(snow)
library(RFQTL)

###specify parameters for rf
ntree <- 200
nforest <- 100
nPermutations <- 400
nCl <- 8


load("/home/cschmal1/phQTL/data/Perlstein_mappingData.RData")
print("real scores")
realScores = rfMapper(mappingData=mappingData, nTrait=nTrait,
                      permute=F, nforest=nforest, ntree=ntree)
print("saving real scores")
save(realScores, file = paste0("/home/cschmal1/phQTL/data/realScores/trait_",nTrait,".RData"))
print("perm scores")
permutedScores <- rfMapper(mappingData = mappingData, nTrait=nTrait,
                           permute = T, nforest = nforest,
                           ntree = ntree, nPermutations=nPermutations,
                           file=paste0("/home/cschmal1/phQTL/data/permutations/trait_",nTrait,".RData"),
                           nCl=nCl,clType="SOCK")