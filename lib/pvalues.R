###pEst###########################################
#function that uses finished permutations to estimate an empirical p-value
#Arguments:
#           path: folder where the .RData files with the matrices of empirical values can be found. They have to
#           conform to predictors (rows) X permutations (columns).
#           markersPerIteration: total number of null distributions that are loaded at the same time. The possible
#           number is dependent on available working memory and size of the null distributions.
#           scores: real predictor-scores to be compared to the empirical null distributions. If scores is a 
#           matrix, it has to conform to traits (rows) X predictors (columns). Scores for the same predictor in
#           different traits (rows) are compared to the same null distribution.
#           printProg: logical, if TRUE the last finished predictor is printed
#Value:     
#           A matrix with empirical p-values is returned.

pEst <- function(path,markersPerIteration,scores,printProg=T){
  #if scores is a vetor (i.e. only one trait), it is transformed into a matrix
  if(is.vector(scores)){
    scores <- matrix(scores,nrow=1)
  }
  
  #identify data objects
  files <- list.files(path)
  files <- files[grep(files,pattern = ".RData")] 
  
  #separate the predictors into packages according to markersPerIteration
  blocks <- seq(1,ncol(scores),markersPerIteration)
  
  #get null distributions per block and compare it to the actual scores
  pMatBlocks <- lapply(1:length(blocks),FUN=function(blockN){
    if(blockN<length(blocks)){
      preds <- blocks[blockN]:(blocks[blockN+1]-1)
    }else{
      preds <- blocks[blockN]:ncol(scores)
    }
    #load the null distributions and bind them into a single matrix
    nullList <- lapply(files,FUN=function(file){
      load(paste0(path,file))
      out <- out[preds,]
      return(out)
    })
    nullMat <- do.call(cbind, nullList)
    rm(nullList)
    gc()
    
    #compute the size of the null distribution to calculate the smallest possible p-value
    nullSize <- ncol(nullMat)
    
    #compare scores to the empirical null distributions
    pMat <- lapply(1:length(preds),FUN=function(predN){
      null <- ecdf((-1)*nullMat[predN,])
      pPred <- null((-1)*scores[,preds[predN]])
      pPred <- pmax(pPred,1/nullSize) #p-values cannot be zero
      rm(null)
      gc()
      return(pPred)
    })
    pMat <- do.call(cbind,pMat)
    print(preds[length(preds)])
    return(pMat)
  })
  #bind all matrices to a single one
  pMatBlock <- do.call(cbind,pMatBlocks)
  rm(pMatBlocks)
  gc()
  
  #return results
  return(pMatBlock)
}

###pEstByFile###########################################
#function that uses finished permutations to estimate an empirical p-value
#Arguments:
#           path: folder where the .RData files with the matrices of empirical values can be found. They have to
#           conform to predictors (rows) X permutations (columns).
#           fileBatch: Number of files with permutations that are loaded at the same time. The possible
#           number is dependent on available working memory.
#           scores: real predictor-scores to be compared to the empirical null distributions. If scores is a 
#           matrix, it has to conform to traits (rows) X predictors (columns). Scores for the same predictor in
#           different traits (rows) are compared to the same null distribution.
#Value:     
#           A matrix with empirical p-values is returned.

pEstByFile <- function(scores,path,fileBatch=1){
  files <- list.files(path)
  files <- files[grep(files,pattern = ".RData")] 
  nFiles <- length(files)
  steps <- seq(1,nFiles,fileBatch)
  if(is.vector(scores)){
    scores <- matrix(scores,nrow=1)
  }
  nPerm <- 0
  scores <- t(scores)
  scoreRanks <- t(apply(scores,1,rank,ties.method="min"))
  if(nrow(scoreRanks)==1){
    scoreRanks <- t(scoreRanks)
  }
  biggerEqual <- matrix(0,nrow = nrow(scores),ncol = ncol(scores))
  for(i in steps){
    allMat <- lapply(i:min(i+fileBatch-1,nFiles),FUN=function(fI){
      load(paste0(path,files[fI]))
      return(out)
    })
    allMat <- do.call("cbind",allMat)
    nPerm <- nPerm+ncol(allMat)
    allMat <- cbind(scores,allMat)
    allRanks <- t(apply(allMat,1,rank,ties.method="min"))
    allRanks <- allRanks[,1:ncol(scores),drop=F]
    allRanks <- allRanks-scoreRanks
    biggerEqual <- biggerEqual+(ncol(allMat)-ncol(scores)-allRanks)
    print(i)
  }
  biggerEqual[biggerEqual==0] <- 1
  pv <- biggerEqual/nPerm
  pv <- t(pv)
  return(pv)
}
