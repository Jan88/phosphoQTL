###get arguments
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
nCl <- 8

###parameters
ntree <- 200
nforest <- 100
exclude <- 3594:3600
load("/home/jgrossb1/phosphoQTL/ptQTL/ptQtlMappingData170831.RData")


###load job data and libraries
library(RandomForestExtended,lib.loc = "/home/jgrossb1/R/x86_64-redhat-linux-gnu-library/3.1")
library(snow)
source("/home/jgrossb1/phosphoQTL/mappingFunctions.R")

cl <- makeCluster(nCl,type="SOCK")
clusterExport(cl=cl,list = c("rfMapper",
                             "replaceGenoNAs",
                             "cScore",
                             "ntree",
                             "nforest",
                             "exclude",
                             "randomForest",
                             "combine",
                             "phenotype",
                             "genotype",
                             "extractNAs"))


AllScores <- parSapply(cl=cl,1:nrow(phenotype),FUN=function(x){
  phenotype <- phenotype[x,]
  if(any(is.na(phenotype))){
    NApos <- which(is.na(phenotype))
    phenotype <- phenotype[-NApos]
    genotype <- genotype[-NApos,]
  }
  
  NAlist <- extractNAs(genotype)
  
  scores <- rfMapper(phenotype = phenotype,
                     genotype = genotype,
                     NAlist = NAlist,
                     exclude = exclude,
                     permute = F,
                     nforest = nforest,
                     ntree = ntree)
  save(scores,file=paste("/scratch/jgrossb1/phosphoQTL/ptQTL/real/trait",x,".RData",sep=""))
  return(scores)
})
stopCluster(cl)
mpi.exit()
