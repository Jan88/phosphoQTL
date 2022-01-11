###QTLgrouper##################
#function that identifies significant marker-trait associations and joins close markers to a single QTL
#if they are consecutive and/or highly correlated and joins distant QTL to QTL-groups if they 
#contain highly correlated markers
#Arguments:
#           pmat: matrix with p-values, traits (rows) X loci (columns)
#           sigThreshold: significance threshold to applied on the p-values
#           corThreshold: threshold for absolute correlation values above which signifcant loci are treated as 
#           linked to each other. If they are close they are joined to a single region, otherwise they are grouped
#           together without including inbetween loci
#           distThreshold: distance threshold up until which loci containing highly correlated markers are linked
#           genotype: allele-information for all samples (rows) and markers (columns)
#           chrVec: character-vector containing the chromosome on which a given marker is located
#Value:     
#           list with QTLgroups. Each entry is a list of the target and a matrix with the QTL-regions involved.
QTLgrouper <- function(pmat, sigThreshold, corThreshold, distThreshold, genotype, chrVec){
  
  QTLmat <- which(pmat<=sigThreshold,arr.ind=T)
  
  #check if any significant QTL are present
  if(nrow(QTLmat)==0){
    return("no significant loci")
  }
  colnames(QTLmat) <- c("target","predictor")
  
  #square threshold and correlation matrix to account for extreme negative correlation
  corThreshold <- corThreshold^2
  corMat <- cor(genotype,use="pair")^2
  
  #join consecutive markers regulating the same trait to a single QTL
  QTLlist <- joinConsecutive(QTLmat, chrVec)
  
  #join markers to QTL if they are close and correlated
  QTLlist <- joinNear(QTLlist,corMat,corThreshold,distThreshold,chrVec)
  
  #link distant QTL that are correlated
  QTLlist <- joinCorrelated(QTLlist,corMat,corThreshold)
  
  #find the predictor with the smallest p-value in each qtl
  QTLlist <- minPV(QTLlist,pmat)
  return(QTLlist)
}

###joinConsecutive######################
#function that joins consecutive significantly linked markers to QTL
#Arguments:
#           QTLmat: matrix consisting of two columns (target, predictor) X one row per linked predictor
#           chrVec: character-vector containing the chromosome on which a given marker is located
#Value:
#           list with one entry for each regulated trait. Each entry contains the target and a matrix with the start
#           and end of the QTL-regions
joinConsecutive <- function(QTLmat, chrVec){
  #identify traits with at least one qtl
  targets <- sort(unique(QTLmat[,1]))
  
  #join consecutive significant predictors per trait
  QTLperTarget <- lapply(targets,FUN=function(target){
    targetMat <- QTLmat[which(QTLmat[,"target"]==target),] #subset for the trait
    if(is.vector(targetMat)){ #if there is only a single significant predictor, there is no joining to do
      outMat <-matrix(c(targetMat[2],targetMat[2]),nrow=1)
      colnames(outMat) <- c("start", "end")
      return(list(target=target,predictors=outMat))
    }
    predictors <- sort(targetMat[,2])
    col <- predictors[1]
    outMat <- NULL
    for(i in 2:length(predictors)){
      #join loci if they are on the same chromosome and are direct neighbors
      if(chrVec[col[1]]==chrVec[predictors[i]]&(predictors[i]-col[length(col)])==1){
        col <- c(col,predictors[i])
      }else{
        outMat <- rbind(outMat,c(min(col),max(col)))
        col <- predictors[i] #start with a new region
      }
    }
    outMat <- rbind(outMat,c(min(col),max(col)))
    colnames(outMat) <- c("start","end")
    #return one matrix per trait containing start and end points of all qtl
    return(list(target=target,predictors=outMat)) 
  })
  #return a list containg these matrices for all regulated traits
  return(QTLperTarget)
}

###joinNear#################################
#function that joins different QTL to a single one if they contain highly correlated markers 
#and are close to each other
#Arguments:
#           QTLlist: output from joinConsecutive
#           corMat: squared matrix of predictor-correlation coefficients
#           corThreshold: threshold indicating how strong QTL-regions have to be correlated to be joined
#           distThreshold: threshold indicating how close QTL-regions have to be to be joined
#           chrVec: character vector containing the information on which chromosome a marker is located
#Value:
#           list with an entry for each regulated trait. Each entry contains the target and a matrix with the start
#           and end of the QTL-regions after combining close and correlated QTL
joinNear <- function(QTLlist, corMat, corThreshold, distThreshold, chrVec){
  out <- lapply(QTLlist,FUN=function(targetList){
    qtl <- targetList[[2]] #matrix with QTL-regions for this trait
    if(nrow(qtl)==1){ #if there is only one region, no regions can be joined
      return(list(target=targetList[[1]],predictors=qtl))
    }
    prog <- T #progress-variable
    while(prog){ #as long as there is progress
      nextProg <- F #loop stops if there is no change
      #each qtl is tested to its right neighbor
      for(i in 1:(nrow(qtl)-1)){
          #are the QTL correlated and close enough?
          if(max(corMat[qtl[i,1]:qtl[i,2],qtl[i+1,1]:qtl[i+1,2]])>=corThreshold&(qtl[i+1,1]-qtl[i,2])<=distThreshold){
            qtl[i,] <- c(qtl[i,1],qtl[i+1,2]) #end of the combined QTL is the end of the second QTL
            qtl <- qtl[-(i+1),] #the second QTL is deleted as it is part of the combined QTL now
            qtl <- matrix(as.vector(qtl),ncol=2,byrow=F) #if there is one QTL remaining it still has the same format
            nextProg <- T
            break #break for loop because indices and matrix size are changed
          }
      }
      prog <- nextProg #prepare for the next iteration
      if(nrow(qtl)==1){ #if only one QTL remains, no further joining can be done
        prog <- F 
      }
    }
    colnames(qtl) <- c("start","end")
    #return the target and the QTL-matrix with start and end points
    return(list(target=targetList[[1]],predictors=qtl)) 
  })
  return(out)
}

###joinCorrelated###############################################
#join QTL-regions that are distant from each other but contain highly correlated markers to QTL-groups using
#hierarchical clustering
#Arguments:
#           QTLlist: output of joinNear
#           corMat: squared matrix of predictor-correlation coefficients
#           corThreshold: threshold indicating how strong QTL-regions have to be correlated to be joined
#Values:    
#           list with an entry for each QTL. Each entry contains the target and a matrix with the start
#           and end of the QTL-regions after grouping highly correlated distant QTL.
joinCorrelated <- function(QTLlist, corMat, corThreshold){
  out <- lapply(QTLlist,FUN=function(targetList){
    qtl <- targetList[[2]] #matrix with QTL-regions for this trait
    if(nrow(qtl)==1){ #if there is only one region, no regions can be joined
      return(list(list(target=targetList[[1]],predictors=qtl)))
    }
    mat <- sapply(1:nrow(qtl), function(i) { #correlation matrix
      sapply(1:nrow(qtl), function(j) {
        return(max(corMat[qtl[i,1]:qtl[i,2],qtl[j,1]:qtl[j,1]], na.rm=T))
      })
    })
    colnames(mat) <- rownames(mat) <- 1:nrow(qtl)
    diag(mat) <- 0 # remove the diagonal
    mat[is.na(mat)] <- 0
    #here we convert the correlation-matrix into distances for clustering 
    dcorm <- as.dist((1 - mat)/2) 
    dcor_cut <- (1-corThreshold)/2 #convert the threshold accordingly
    dcorm_hclust <- hclust(dcorm)
    group <- cutree(hclust(dcorm), h=dcor_cut) #apply the threshold to the tree
    group <- split(as.numeric(names(group)), group)
    out <- lapply(group,FUN=function(gr){ #join regions to QTL-groups
      outMat <- matrix(qtl[gr,],nrow=length(gr),byrow = F)
      colnames(outMat) <- c("start","end")
      list(target=targetList[[1]],predictors=outMat)
    })
    return(out)
  })
  out <- do.call(c,out) #parse list so that one QTL-group is one entry regardless of the target
  names(out) <- NULL
  return(out)
}

###minPV######################
#identify the predictor with the smallest p-value in each QTL
#Arguments:
#           QTLlist: list-object with an entry for each QTL as returned by joinCorrelated
#           pmat: matrix containing the p-values for each predictor(columns)/trait(row)
#Value:
#           object similar to QTLlist with added information regarding the most significant 
#           predictor in each QTL

minPV <- function(QTLlist,pmat){
  out <- lapply(QTLlist,FUN=function(qtl){
    target <- qtl$target
    predictors <- apply(qtl$predictors,1,FUN=function(preds){return(preds[1]:preds[2])})
    if(is.list(predictors)){
      predictors <- do.call("c",predictors)
    }
    sigPred <- predictors[which.min(pmat[target,predictors])]
    qtl <- c(qtl,sigPred)
    names(qtl)[length(qtl)] <- "mostSignificantPredictor"
    return(qtl)
  })
  return(out)
}

###testConditional#################################
#test QTL-groups as returned by QTLgrouper for conditional effects
#Arguments:
#           QTLgroups: list-object as returned by QTLgrouper, with each element being a list of the target
#           and a matrix with predictor regions
#           phenotype: matrix with trait-values for different targets (rows) and individuals (columns). If
#           genotype2group is not NULL, colnames have to be assigned.
#           genotype: genotype containing predictor status for different individuals (rows) 
#           and different markers (columns)
#           covariates: covariates for all individuals (rows)
#           condition: vector with condition-information for each sample
#           condThresholdLocus: significance threshold for the calling of conditional QTL
#           condThresholdCondition: significance threshold for the calling of the condition the QTL is active in
#           genotype2group: optional. If the genotype contains non-unique columns a vector can be used to translate 
#           a genotype with redundant columns into a smaller one with unique columns. This can decrease the burden
#           of multiple testing correction. Predictors are transformed as well. If the default is used, all 
#           predictors and genotypes are used as they are.
#Value:
#           list object like QTLgroups. In addition each entry contains the QTL-status.
testConditional <- function(QTLgroups, phenotype, genotype, covariates, condition, condThresholdLocus, condThresholdCondition){
  #if not all predictors were used for the mapping individually, they can be reformatted to avoid additional testing
  
  if(is.null(colnames(phenotype))){return("phenotypes have to be matched to strains via colnames")}
  
  pList <- lapply(1:length(QTLgroups),FUN=function(nGroup){
    QTLgroup <- QTLgroups[[nGroup]] 
    target <- QTLgroup$target
    genotype <- genotype[colnames(phenotype),]
    minPred <- QTLgroup$mostSignificantPredictor
    
    #test the most significant predictor
    pVals <- testConditionalSub(expres=phenotype[target,],
                       genotype=genotype[,minPred],
                       covariates=covariates,
                       condition=as.factor(condition))
    
    #pVals <- lapply(predictors,FUN=function(predictor){
    #  testConditionalSub(expres=phenotype[target,],
    #                     genotype=genotype[,predictor],
    #                    covariates=covariates,
    #                     condition=as.factor(condition))
    #})
    #pVals <- do.call(rbind, pVals)
    pVals <- c(nGroup, pVals)
    return(pVals)
  })
  #rbind all QTL-groups to facilitate FDR-estimation
  pList <- do.call(rbind,pList)
  
  #estimate the FDR only for the interaction test between condition and predictor, not the tests aiming to find the 
  #affected condition-level
  pList[,2] <- p.adjust(pList[,2],method="fdr") 
  ts1 <<- pList
  
  #separate the QTL-groups again and assign a status regarding the conditional regulation according to the thresholds
  pList <- lapply(unique(pList[,1]),FUN=function(nOut){
    subMat <- subset(pList,pList[,1]==nOut)
    qtlStatus <- assignStatus(subMat[,2:4],
                              condThreshold1=condThresholdLocus,
                              condThreshold2=condThresholdCondition)
    return(list(target=QTLgroups[[nOut]]$target,
                predictors=QTLgroups[[nOut]]$predictors,
                pValues=subMat[,2:4],
                status=qtlStatus))
  })
}

###testConditionalSub###############################
#tests for an interaction between a marker and a condition and tests in which condition the predictor is more 
#meaningful
#Arguments:
#           expres: trait-vector
#           genotype: genotype containing predictor (columns) status for different individuals (rows)
#           covariates: population covariates (row/sample) which are included as variables in the models
#           condition: vector indicating the condition for each sample
#Value:     
#           vector with a probability for the observed interaction between predictor and condition and probabilities
#           for the effect of the predictor in condition-levels assuming a linear relationship
testConditionalSub <- function(expres, genotype, covariates, condition){
  full<-lm(expres~genotype+condition+covariates+genotype:condition)
  cond <- anova(full)[4,5] #p-value for the interaction
  PH<-sapply(levels(condition),function(i){ #p-values for the predictor for a specific condition-level
    anova(lm(expres[condition==i]~covariates[condition==i,]+genotype[condition==i]))[2,5] 
  })
  return(c(cond,p.adjust(PH,method="bonferroni"))) #only the p-values for the levels are corrected with bonferroni
}

###assignStatus##################################
#assigns a status regarding the interaction between condition and predictor according to the thresholds based
#on the predictor with the smallest p-value for an interaction with the condition
#Arguments:
#           mat: matrix with the p-values for each tested predictor
#           condThreshold1: significance threshold for the calling of conditional QTL
#           condThreshold2: significance threshold for the calling of the condition the QTL is active in
#Value:
#           vector containing the status and optionally the condition the QTL is active in
assignStatus <- function(mat,condThreshold1,condThreshold2){
  if(is.null(nrow(mat))){
    mat <- t(mat)
  }
  rowOI <- mat[which.min(mat[,1]),] #use the predictor with the smallest p-value
  
  #if the interaction is significant but the qtl is not significant in any condition, return static
  if(rowOI[1]<=condThreshold1&sum(rowOI[2:length(rowOI)]<=condThreshold2)==0){
    return("static")
  }
  
  #if the interaction is not significant, return static
  if(rowOI[1]>condThreshold1){
    return("static")
  }
  
  #if the interaction is significant and at least one affected condition can be identified, return conditional and
  #the conditions
  if(rowOI[1]<=condThreshold1&sum(rowOI[2:length(rowOI)]<=condThreshold2)>0){
    conds <- colnames(mat)[which(rowOI[2:length(rowOI)]<=condThreshold2)+1]
    return(c("conditional",conds))
  }
}

###QTLplotter#########################################################################################
#plot a QTLmatrix generated with QTLgrouper
#Arguments:
#           QTLlist: output of QTLgrouper
#           targetLocs: location of targets in the genome in absolute cummulative positions
#           predictorLocs: location of predictors in the genome in absolute cummulative positions
#           chrLen: vector containing the length of each chromosome with names
#           xlab: desired label for the x-axis
#           ylab: desired label for the y-axis
#           main: desired main-title
#           qtl.lwd: width of the lines for the qtl
QTLplotter <- function(QTLlist, targetLocs, predictorLocs, chrLen, xlab, ylab, main, col="black",labcex=0.6,qtl.lwd=3,chrLty=2,chrCol="red"){
  plot(NULL,
       ylim=c(0,max(targetLocs,na.rm=T)),
       xlim=c(0,max(predictorLocs,na.rm=T)),
       xlab=xlab,
       ylab=ylab,
       main=main,
       xaxt="n",
       yaxt="n")
  sapply(QTLlist,FUN=function(qtl){
    if(!is.na(targetLocs[qtl$target])){
      apply(qtl$predictors,1,FUN=function(row){
        lines(y=rep(targetLocs[qtl$target],2),x=predictorLocs[row],col=col,lwd=qtl.lwd)
      })
    }
  })
  chrLenCum <- sapply(1:length(chrLen),FUN=function(x){
    sum(chrLen[1:x])
  })
  for(i in 1:(length(chrLenCum)-1)){
    abline(h=chrLenCum[i],lty=chrLty,col=chrCol)
    abline(v=chrLenCum[i],lty=chrLty,col=chrCol)
  }
  borders <- c(0,chrLenCum)
  for(i in 1:length(chrLen)){
    mtext(names(chrLen[i]),side=1,at=mean(borders[i:(i+1)]),cex=par()$cex*labcex)
    mtext(names(chrLen[i]),side=2,at=mean(borders[i:(i+1)]),cex=par()$cex*labcex)
  }
}

###plotSplit#######################################
#function that plots individuals according to their trait values in different conditions. Points are colored
#according to the allele the individual carries.
#Arguments: 
#           traitC1: trait matrix for condition 1. Each trait is a row.
#           traitC2: trait matrix for condition 2. Each trait is a row.
#           genotype: Allele information (0 or 1) as a matrix in the form of individual (row) X markers (columns)
#           xlab: desired label for the x-axis
#           ylab: desired label for the y-axis
#           main: desired main-title
#           nT: number of the trait of interest
#           nG: number of the marker of interest
plotSplit <- function(traitC1,traitC2,genotype,main="",xlab="",ylab="",nT,nG){
  traitC1 <- traitC1[nT,]
  traitC2 <- traitC2[nT,]
  genotype <- genotype[,nG]
  missing <- which(is.na(traitC1)|is.na(traitC2))
  if(length(missing)>=1){
    traitC1 <- traitC1[-missing]
    traitC2 <- traitC2[-missing]
    genotype <- genotype[-missing]
  }
  plot(NULL,xlim=range(traitC1),ylim=range(traitC2),xlab=xlab,ylab=ylab,main=main)
  abline(a=0,b=1,col="grey",lty=2)
  cols <- c("red","blue")
  invisible(sapply(c(0,1),FUN=function(allele){
    toPlot <- which(genotype==allele)
    points(x=traitC1[toPlot],y=traitC2[toPlot],col=cols[allele+1],pch=20)
  }))
}

###plotSplit2#######################################
#function that plots individuals according to their trait values in different conditions. Points are colored
#according to the allele the individual carries.
#Arguments: 
#           traitC1: trait matrix for condition 1. Each trait is a row.
#           traitC2: trait matrix for condition 2. Each trait is a row.
#           genotype: Allele information (0 or 1) as a matrix in the form of individual (row) X markers (columns)
#           arg.names: desired label for the x-axis
#           ylab: desired label for the y-axis
#           main: desired main-title
#           nT: number of the trait of interest
#           nG: number of the marker of interest
plotSplit2 <- function(traitC1,traitC2,genotype,arg.names="",main=c("",""),nT,nG,ex_par=F){
  require(beeswarm)
  traitC1 <- traitC1[nT,]
  traitC2 <- traitC2[nT,]
  genotype <- genotype[,nG]
  if(!ex_par){
    par(mfrow=c(1,2))
  }
  cols <- c("red","blue")
  toPlot <- list(traitC1[which(genotype==0)],traitC1[which(genotype==1)])
  names(toPlot) <- arg.names
  beeswarm(toPlot,col=cols,pch=20)
  title(main[1])
  toPlot <- list(traitC2[which(genotype==0)],traitC2[which(genotype==1)])
  names(toPlot) <- arg.names
  beeswarm(toPlot,col=cols,pch=20)
  title(main[2])
}
