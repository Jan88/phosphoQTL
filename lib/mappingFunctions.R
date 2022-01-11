### extractNAs##############
###         function that calculates the allele frequencies (1 or 0) for each marker (column).
##input:    genotype-matrix with one row per individual, one column per marker and binary marker information (1,0,NA)
##output:   list with one vector for each marker. The first element is the frequency of allele 1. If there
#           are missing values the indices of the respective samples follow.
extractNAs <- function(genotype){
  apply(genotype,2,FUN=function(predictor){
    out <- length(which(predictor==1))/length(which(!is.na(predictor)))
    if(any(is.na(predictor))){
      out <- c(out,which(is.na(predictor)))
    }
    return(out)
  })
}

### rfMapper########################
###           wrapper function for qtl mapping with rf
##input:      phenotype: trait values as a numerical vector. NAs are removed in this function.
#             genotype: genotype as a matrix with 1,0,NA. Predictors that are indicated in exclude may be quantitative
#             permute: logical, should permutations be performed (if false the real trait-values are used)?
#             nPermutations: number of permutations to be performed. If permute is false, nPermutations is irrelevant.
#             nforest: number of randomizations for missing genotypes
#             pMat: optional, predetermined permutation scheme in form of a matrix. one row should contain 
#             the new indices for all individuals
##output:     for unpermuted traits: a numerical vector with importance scores for all predictors aside from those in exclude
#             for permuted traits: a matrix of scores with one column per permutation
rfMapper <- function(phenotype, genotype, NAlist, exclude, permute, nPermutations=1000, nforest, ntree, pMat=NULL){
  if(any(is.na(phenotype))){ #remove individuals with missing trait-values
    stop("missing outcomes")
    #missing <- which(is.na(phenotype))
    #phenotype <- phenotype[-missing]
    #genotype <- genotype[-missing,]
  }
  if(permute){ #do permutations
    if(is.null(pMat)){
      out <- sapply(1:nPermutations,FUN=function(nPerm){
        pScheme <- sample(1:length(phenotype)) #generate a permutation scheme
        permPhenotype <- phenotype[pScheme] #permute the trait vector
        permGenotype <- genotype
        if(!is.null(exclude)){
          permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
        }
        rfs <- lapply(1:nforest,FUN=function(nF){
          permGenotype <- replaceGenoNAs(genotype=permGenotype,NAlist=NAlist)
          rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE,keep.forest=F)
          return(rf)
        })
        bigRf <- combine(rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
        scores <- cScore(bigRf)
        if(!is.null(exclude)){
          scores <- scores[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
        }
        scores
      })
    }else{
      if(class(pMat)!="matrix"){
        stop("permutation-matrix has to be a matrix")
      }
      if(nrow(pMat)<nPermutations|ncol(pMat)!=length(phenotype)){
        stop("permutation-matrix has wrong dimensions")
      }
      out <- sapply(1:nPermutations,FUN=function(nPerm){
        pScheme <- pMat[nPerm,] #generate a permutation scheme
        permPhenotype <- phenotype[pScheme] #permute the trait vector
        permGenotype <- genotype 
        if(!is.null(exclude)){
          permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
        }
        rfs <- lapply(1:nforest,FUN=function(nF){
          permGenotype <- replaceGenoNAs(genotype=permGenotype,NAlist=NAlist)
          rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE,keep.forest=F)
          return(rf)
        })
        bigRf <- combine(rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
        scores <- cScore(bigRf)
        if(!is.null(exclude)){
          scores <- scores[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
        }
        scores
        })
    }
  }else{ #map real traits
    rfs <- lapply(1:nforest,FUN=function(nF){
      genotype <- replaceGenoNAs(genotype=genotype,NAlist=NAlist)
      rf <- randomForest(x=genotype, y=phenotype, ntree=ntree, importance=TRUE)
      return(rf)
    })
    bigRf <- combine(rfs)
    out <- cScore(bigRf)
    if(!is.null(exclude)){
      out <- out[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
    }
    
  }
  
  return(out) #return the scores as a vector for unpermuted traits and as a matrix for permuted traits
}

### rfMapperPar########################
###           wrapper function for qtl mapping with rf while parallelizing permutations
##input:      phenotype: trait values as a numerical vector. NAs are removed in this function.
#             genotype: genotype as a matrix with 1,0,NA. Predictors that are indicated in exclude may be quantitative
#             permute: logical, should permutations be performed (if false the real trait-values are used)?
#             nPermutations: number of permutations to be performed. If permute is false, nPermutations is irrelevant.
#             nforest: number of randomizations for missing genotypes
#             pMat: optional, predetermined permutation scheme in form of a matrix. one row should contain 
#             the new indices for all individuals
##output:     for unpermuted traits: a numerical vector with importance scores for all predictors aside from those in exclude
#             for permuted traits: a matrix of scores with one column per permutation
rfMapperPar <- function(phenotype, genotype, NAlist, exclude, nPermutations=1000, cl, nforest, ntree, pMat=NULL){
  if(any(is.na(phenotype))){ #remove individuals with missing trait-values
    stop("missing outcomes")
    #missing <- which(is.na(phenotype))
    #phenotype <- phenotype[-missing]
    #genotype <- genotype[-missing,]
  }
  if(is.null(pMat)){  
    out <- parSapply(cl=cl,1:nPermutations,FUN=function(nPerm){
        pScheme <- sample(1:length(phenotype)) #generate a permutation scheme
        permPhenotype <- phenotype[pScheme] #permute the trait vector
        permGenotype <- genotype 
        permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
        rfs <- lapply(1:nforest,FUN=function(nF){
          permGenotype <- replaceGenoNAs(genotype=permGenotype,NAlist=NAlist)
          rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE, keep.forest=F)
          return(rf)
        })
        bigRf <- combine(rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
        scores <- cScore(bigRf)[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
        return(scores)
      })
  }else{
    if(class(pMat)!="matrix"){
      stop("permutation-matrix has to be a matrix")
    }
    if(nrow(pMat)<nPermutations|ncol(pMat)!=length(phenotype)){
      stop(paste(dim(pMat),"permutation-matrix has wrong dimensions"))
    }
    out <- parSapply(cl=cl,1:nPermutations,FUN=function(nPerm){
      pScheme <- pMat[nPerm,] #generate a permutation scheme
      permPhenotype <- phenotype[pScheme] #permute the trait vector
      permGenotype <- genotype 
      permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
      rfs <- lapply(1:nforest,FUN=function(nF){
        permGenotype <- replaceGenoNAs(genotype=permGenotype,NAlist=NAlist)
        rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE, keep.forest=F)
        return(rf)
      })
      bigRf <- combine(rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
      scores <- cScore(bigRf)[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
      return(scores)
    })
  }
  
  return(out) #return the scores as a vector for unpermuted traits and as a matrix for permuted traits
}

### replaceGenoNAs#########################################
###function that replaces missing genotypes according to the local allele frequencies
###input:
#       genotype:   genotype as a matrix with 1,0,NA. 
#                   Predictors that are indicated in exclude may be quantitative
#       NAlist:     list object with one vector per predictor. The vector has to contain 
#                   the frequency of the 1-allele first. Indices of individuals with missing values
#                   follow.
###output:
#       Matrix with the same dimensions as the input-genotype. NAs are replaced with alleles (0,1)
#       according to the allele-frequencies in the input (NAlist)
replaceGenoNAs <- function(genotype,NAlist){
  sapply(1:ncol(genotype),FUN=function(x){ #per predictor
    lNA <- length(NAlist[[x]])-1
    if(lNA==0){
      return(genotype[,x]) #if no values are missing, return the original vector
    }
    genotype[NAlist[[x]][2:(lNA+1)],x] <- sample(c(1,0),size=lNA,prob = c(NAlist[[x]][1],1-NAlist[[x]][1]),replace=T) #sample NAs according to NAlist
    return(genotype[,x]) #return the extended genotype
  })
}

### cScore#################       
###extracts a combined importance score from a randomForest object.
### input:
#  rf:                      A randomForest-object that contains both permutation importance (PI) and the increase in node purity (RSS).
#                           Random Forests that were generated with importance=TRUE are suited for this.
#  predictorsOfInterest:    Predictors for which the combined score should be computed. As a default all predictors are selected.
#  normalize:               Should the scores be normalized with the trait variance? This makes the importance scores comparable across different traits.
### output:
#                           The output consists of a numeric vector containing a score for each predictor of interest.
#                           The score is the product of RSS and PI. Negative RSS- or PI-scores are set to zero first.
#                           These scores are influenced by correlation between predictors and therefore contain a marker-specific bias.
cScore <- function(rf=NULL,predictorsOfInterest=NULL,normalize=TRUE){
  if(class(rf)!="randomForest"){
    stop("A randomForest object is required as input!")
  }
  if(is.null(predictorsOfInterest)){
    predictorsOfInterest <- 1:nrow(rf$importance)
  }else{
    if((!is.numeric(predictorsOfInterest))|(!is.vector(predictorsOfInterest))){
      stop("Predictors of interest have to be entered as a numeric vector!")
    }
  }
  if(ncol(rf$importance)<2){
    stop("A randomForest object with both permutation importance and increase in node purity (importance=TRUE) is required as input!")
  }
  if(normalize){ #divide scores by the trait variance. default
    if(is.null(rf$y)|rf$type!="regression"){
      stop("Normalization can only be performed if the randomForest object was created by supervised regression!")
    }
    tVar <- var(rf$y)
    PI <- rf$importance[predictorsOfInterest,1]/tVar
    PI[which(PI<0)] <- 0
    RSS <- rf$importance[predictorsOfInterest,2]/tVar
    RSS[which(RSS<0)] <- 0
  }else{ #don't normalize
    PI <- rf$importance[predictorsOfInterest,1]
    PI[which(PI<0)] <- 0
    RSS <- rf$importance[predictorsOfInterest,2]
    RSS[which(RSS<0)] <- 0
  }
  combinedScore <- PI*RSS
  names(combinedScore) <- rownames(rf$importance)[predictorsOfInterest]
  return(combinedScore)
}

###condID####################
condID <- function(y ,x1 ,x2){
  m1 <- lm(y~x1*x2)
  m2 <- lm(y~x1+x2)
  return(anova(m1,m2)[2,6])
}

###condIDall############
condIDall <- function(phenotype, genotype, predictorsOfInterest, condInd, nCl, clType="SOCK"){
  require(snow)
  cl <- makeCluster(nCl,clType)
  clusterExport(cl,
                list = c("phenotype",
                         "genotype",
                         "predictorsOfInterest",
                         "condInd",
                         "condID"),
                envir = environment())
  out <- parApply(cl,phenotype[1:100,],1,FUN=function(trait){
    apply(genotype[,predictorsOfInterest],2,FUN=condID,y=trait,x2=genotype[,condInd])
  })
  stopCluster(cl)
  return(t(out))
}

speedTest <- function(){
  a <- proc.time()
  work  <- 0
  for(i in 1:10000000){
    work <- 1*5
  }
  b <- proc.time()
  return(as.numeric(b[3]-a[3]))
}
