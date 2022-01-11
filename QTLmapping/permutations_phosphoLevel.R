###get arguments
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
nCl <- 8
permsToMap <- args[1]*10
permsToMap <- (permsToMap-9):permsToMap

###load job data and libraries
library(RandomForestExtended,lib.loc = "/home/jgrossb1/R/x86_64-redhat-linux-gnu-library/3.1")
library(snow)
source("/home/jgrossb1/phosphoQTL/mappingFunctions.R")
load("/home/jgrossb1/phosphoQTL/phosphoLevelQTL/phosphoLevelQtlMappingData170831.RData")


###specify parameters for rf
ntree <- 200
nforest <- 100
exclude <- 3594:3600
nPermutations <- 600
phenotype <- phenotype4perm

###find the missing genotypic values
NAlist <- extractNAs(genotype)

###create the cluster, share information
cl <- makeCluster(nCl,type="SOCK")
clusterExport(cl=cl,list = c("genotype",
                       "phenotype",
                       "rfMapperPar",
                       "replaceGenoNAs",
                       "cScore",
                       "ntree",
                       "nforest",
                       "NAlist",
                       "exclude",
                       "randomForest",
                       "combine"))

###loop through the trait-job vector
for(i in permsToMap){
  a <- proc.time()
  strains <- rownames(genotype)
  pMat_new <- lapply(1:nPermutations,FUN=function(x){
    pScheme <- sample(1:nrow(genotype))
    return(pScheme)
  })
  pMat_new <- do.call("rbind",pMat_new)
  pMat <- pMat_new
  clusterExport(cl=cl,list=c("pMat"))
  #work the trait
  out <- rfMapperPar(phenotype=phenotype[i,], 
                     genotype=genotype, 
                     NAlist=NAlist, 
                     exclude=exclude, 
                     nPermutations=nPermutations, 
                     cl=cl, 
                     nforest=nforest, 
                     ntree=ntree,
		                 pMat=pMat)
  save(out,file=paste("/scratch/jgrossb1/phosphoQTL/phosphoLevelQTL/perms/perm",i,".RData",sep=""))
  b <- proc.time()
  print((b-a)[3]/60)
}
stopCluster(cl)
mpi.exit()
