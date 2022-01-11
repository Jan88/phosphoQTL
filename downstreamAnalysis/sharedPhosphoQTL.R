source("lib/general_function.R")
library(RColorBrewer)
library(topGO)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
###load qtl results###
load("data/pQTL_results170815.RData")
load("data/phosphoNew/phosphoLevel.RData")
rownames(phospho2prot) <- phospho2prot[,1]

###prepare matrix of shared QTL###
phResPeps <- rownames(pQTL_results$phosphoProtResidualsQTL$qv)
phResProts <- unique(phospho2prot[phResPeps,2])

sharedProtein <- matrix(0,ncol=length(phResProts),nrow=length(phResProts))
rownames(sharedProtein) <- colnames(sharedProtein) <- phResProts
diag(sharedProtein) <- 1

###load complex data###
uniProt2geneSymbol <- read.table(file = "Saccharomyces_cerevisiae/uniprot2ens.tab",header = T,sep = "\t",quote="",as.is=T)
rownames(uniProt2geneSymbol) <- uniProt2geneSymbol[,3]
complexAnnotation <- read.table("Saccharomyces_cerevisiae/saccharomyces_cerevisiae_complexes.tsv",sep="\t",as.is=T,quote="",comment.char = "",header = T)
complexMembers <- lapply(complexAnnotation[,5],FUN=function(idString){
  ids <- strsplit(x = idString,split = "|",fixed=T)[[1]]
  ids <- gsub(pattern = "\\(.*",replacement = "",x = ids)
  geneSymbols <- uniProt2geneSymbol[ids,1]
  geneSymbols <- geneSymbols[!is.na(geneSymbols)]
  return(geneSymbols)
})
names(complexMembers) <- complexAnnotation[,1]
complexMembers <- complexMembers[sapply(complexMembers,length)>1]
sharedComplex <- matrix(0,ncol=length(phResProts),nrow=length(phResProts))
rownames(sharedComplex) <- colnames(sharedComplex) <- phResProts
for(i in 1:length(complexMembers)){
  gvec <- intersect(complexMembers[[i]],phResProts)
  if(length(gvec)>0){
    sharedComplex[gvec,gvec] <- 1
  }
}
###prepare compartment info###
uniGenes <- phResProts
goTbl <- read.table("data/gene_association.sgd",comment.char = "!",sep="\t",quote="",as.is=T)
CC_annotation <- lapply(uniGenes,FUN=function(x){
  sub_go <- goTbl[which(goTbl[,11]==x&goTbl[,9]=="C"),]
  return(as.character(sub_go[,5]))
})
names(CC_annotation) <- uniGenes

allGenes <- rep(0,length(uniGenes))
names(allGenes) <- uniGenes
allGenes[1:(length(allGenes)-1)] <- 1
allGenes <- as.factor(allGenes)

GOdataCC <- new("topGOdata",
                ontology = "CC",
                allGenes = allGenes,
                annot = annFUN.gene2GO,
                gene2GO = CC_annotation,
                nodeSize=1)
resultFisherCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")

CCterms <- c("GO:0005634", #Nucleus
             "GO:0005829", #Cytosol
             "GO:0005886", #Membrane
             "GO:0005739", #Mitochondrion
             "GO:0005794", #Golgi
             "GO:0005783", #ER
             "GO:0005618") #Cell wall
CCtermList <- lapply(CCterms,FUN=function(term){
  unlist(genesInTerm(object = GOdataCC,whichGO = term))
})
names(CCtermList) <- CCterms

sharedCompartment <- matrix(0,ncol=length(phResProts),nrow=length(phResProts))
rownames(sharedCompartment) <- colnames(sharedCompartment) <- phResProts
for(i in CCterms){
  sharedCompartment[CCtermList[[i]],CCtermList[[i]]] <- 1
}

###investigate which phospho sites have overlapping qtl
#phr
sharedPhrQTL <- matrix(0,nrow = nrow(pQTL_results$phosphoProtResidualsQTL$qv),ncol = nrow(pQTL_results$phosphoProtResidualsQTL$qv))
phrTargets <- sapply(pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10,FUN=function(x){x$target})
for(i in (1:(length(phResPeps)-1))){
  for(j in (i+1):length(phResPeps)){
    
    if(all(c(i,j)%in%phrTargets)){
      iTargets <- which(phrTargets==i)
      jTargets <- which(phrTargets==j)
      iQTL <- lapply(iTargets,FUN=function(n){pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10[[n]]$predictors})
      jQTL <- lapply(jTargets,FUN=function(n){pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10[[n]]$predictors})
      ovMat <- matrix(NA,nrow=length(iQTL),ncol = length(jQTL))
      for(nI in 1:length(iQTL)){
        for(nJ in 1:length(jQTL)){
          ovMat[nI,nJ] <- length(intersect(as.vector(unlist(apply(iQTL[[nI]],1,FUN=function(x){x[1]:x[2]}))),as.vector(unlist(apply(jQTL[[nJ]],1,FUN=function(x){x[1]:x[2]})))))>0
        }
      }
      overlap <- min(sum(rowSums(ovMat)),sum(colSums(ovMat)),1)
    }else{
      overlap <- 0
    }
    sharedPhrQTL[i,j] <- sharedPhrQTL[j,i] <- overlap
  }
  print(i)
}
rownames(sharedPhrQTL) <- colnames(sharedPhrQTL) <- phResPeps
diag(sharedPhrQTL) <- NA
#ph
sharedPhQTL <- matrix(0,nrow = nrow(pQTL_results$phosphoProtResidualsQTL$qv),ncol = nrow(pQTL_results$phosphoProtResidualsQTL$qv))
phTargets <- sapply(pQTL_results$phosphoLevelQTL$QtlList$FDR10,FUN=function(x){x$target})
for(i in (1:(length(phResPeps)-1))){
  for(j in (i+1):length(phResPeps)){
    phI <- which(rownames(pQTL_results$phosphoLevelQTL$qv)==phResPeps[i])
    phJ <- which(rownames(pQTL_results$phosphoLevelQTL$qv)==phResPeps[j])
    if(all(c(phI,phJ)%in%phTargets)){
      iTargets <- which(phTargets==phI)
      jTargets <- which(phTargets==phJ)
      iQTL <- lapply(iTargets,FUN=function(n){pQTL_results$phosphoLevelQTL$QtlList$FDR10[[n]]$predictors})
      jQTL <- lapply(jTargets,FUN=function(n){pQTL_results$phosphoLevelQTL$QtlList$FDR10[[n]]$predictors})
      ovMat <- matrix(NA,nrow=length(iQTL),ncol = length(jQTL))
      for(nI in 1:length(iQTL)){
        for(nJ in 1:length(jQTL)){
          ovMat[nI,nJ] <- length(intersect(as.vector(unlist(apply(iQTL[[nI]],1,FUN=function(x){x[1]:x[2]}))),as.vector(unlist(apply(jQTL[[nJ]],1,FUN=function(x){x[1]:x[2]})))))>0
        }
      }
      overlap <- min(sum(rowSums(ovMat)),sum(colSums(ovMat)),1)
    }else{
      overlap <- 0
    }
    sharedPhQTL[i,j] <- sharedPhQTL[j,i] <- overlap
  }
  print(i)
}
rownames(sharedPhQTL) <- colnames(sharedPhQTL) <- phResPeps
diag(sharedPhQTL) <- NA

###compute non exclusive genesets and proportions of shared qtl
#No relation
noRel <- sharedCompartment[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0&sharedComplex[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0&sharedProtein[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0

#same compartment
sameCC <- sharedCompartment[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==1

#same complex
sameComplex <- sharedComplex[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==1

#same protein
sameProtein <- sharedProtein[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==1
rownames(sameProtein) <- rownames(sameComplex) <- rownames(sameCC) <- rownames(noRel) <- phResPeps
colnames(sameProtein) <- colnames(sameComplex) <- colnames(sameCC) <- colnames(noRel) <- phResPeps
nonExList <- list(noRel=noRel,sameCC=sameCC,sameComplex=sameComplex,sameProtein=sameProtein)

nonExRes <- lapply(nonExList,FUN=function(logMat){
  diag(logMat) <- NA
  keepPh <- colSums(logMat,na.rm=T)>0
  logMat <- logMat[keepPh,keepPh]
  resMat <- cbind(as.vector(sharedPhQTL[colnames(logMat),colnames(logMat)]),as.vector(logMat))
  resPh <- c(sum(!resMat[,1]&!resMat[,2],na.rm=T),sum(resMat[,1]&!resMat[,2],na.rm=T),sum(!resMat[,1]&resMat[,2],na.rm=T),sum(resMat[,1]&resMat[,2],na.rm=T))
  
  resMat <- cbind(as.vector(sharedPhrQTL[colnames(logMat),colnames(logMat)]),as.vector(logMat))
  resPhr <- c(sum(!resMat[,1]&!resMat[,2],na.rm=T),sum(resMat[,1]&!resMat[,2],na.rm=T),sum(!resMat[,1]&resMat[,2],na.rm=T),sum(resMat[,1]&resMat[,2],na.rm=T))
  
  resMat <- rbind(resPh,resPhr)
  rownames(resMat) <- c("ph","phr")
  colnames(resMat) <- c("NN","PN","NP","PP") #shared qtl, shared annotation
  return(resMat/2)
})
###compute exclusive genesets and proportions of shared qtl
#No relation
noRel <- sharedCompartment[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0&sharedComplex[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0&sharedProtein[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0

#same compartment but not same complex or protein
sameCC <- sharedCompartment[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==1&sharedComplex[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0&sharedProtein[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0

#same complex but not same protein
sameComplex <- sharedComplex[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==1&sharedProtein[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==0

#same protein
sameProtein <- sharedProtein[phospho2prot[phResPeps,2],phospho2prot[phResPeps,2]]==1
rownames(sameProtein) <- rownames(sameComplex) <- rownames(sameCC) <- rownames(noRel) <- phResPeps
colnames(sameProtein) <- colnames(sameComplex) <- colnames(sameCC) <- colnames(noRel) <- phResPeps

exList <- list(noRel=noRel,sameCC=sameCC,sameComplex=sameComplex,sameProtein=sameProtein)
exRes <- lapply(exList,FUN=function(logMat){
  diag(logMat) <- NA
  keepPh <- colSums(logMat,na.rm=T)>0
  logMat <- logMat[keepPh,keepPh]
  resMat <- cbind(as.vector(sharedPhQTL[colnames(logMat),colnames(logMat)]),as.vector(logMat))
  resPh <- c(sum(!resMat[,1]&!resMat[,2],na.rm=T),sum(resMat[,1]&!resMat[,2],na.rm=T),sum(!resMat[,1]&resMat[,2],na.rm=T),sum(resMat[,1]&resMat[,2],na.rm=T))
  
  resMat <- cbind(as.vector(sharedPhrQTL[colnames(logMat),colnames(logMat)]),as.vector(logMat))
  resPhr <- c(sum(!resMat[,1]&!resMat[,2],na.rm=T),sum(resMat[,1]&!resMat[,2],na.rm=T),sum(!resMat[,1]&resMat[,2],na.rm=T),sum(resMat[,1]&resMat[,2],na.rm=T))
  
  resMat <- rbind(resPh,resPhr)
  rownames(resMat) <- c("ph","phr")
  colnames(resMat) <- c("NN","PN","NP","PP") #shared qtl, shared annotation
  return(resMat/2)
})
###compute proportions of pairs of phosphosites that share a qtl
sharedQtlEx <- t(sapply(exList,FUN=function(mat){
  ph <- sum(sharedPhQTL[mat],na.rm=T)/sum(!is.na(sharedPhQTL[mat]))
  phr <- sum(sharedPhrQTL[mat],na.rm=T)/sum(!is.na(sharedPhrQTL[mat]))
  return(c(ph,phr))
}))
colnames(sharedQtlEx) <- c("ph","phr")

sharedQtlNonEx <- t(sapply(nonExList,FUN=function(mat){
  ph <- sum(sharedPhQTL[mat],na.rm=T)/sum(!is.na(sharedPhQTL[mat]))
  phr <- sum(sharedPhrQTL[mat],na.rm=T)/sum(!is.na(sharedPhrQTL[mat]))
  return(c(ph,phr))
}))
colnames(sharedQtlNonEx) <- c("ph","phr")

###panel pdf for all data as an ovierview
pdfAndPng("graph/sharedPhosphoQTL", width = 10, height = 10, expr = expression({
  par(mfrow=c(2,2))
  barplot(sharedQtlNonEx[,1],col=colors[5],names.arg = c("no \nrelation","same \ncompartment","same \ncomplex","same \nprotein"),ylab="% of pairs with a shared QTL",main="phQTL, non exclusive conditions")
  barplot(sharedQtlNonEx[,2],col=colors[4],names.arg = c("no \nrelation","same \ncompartment","same \ncomplex","same \nprotein"),ylab="% of pairs with a shared QTL",main="phResQTL, non exclusive conditions")
  barplot(sharedQtlEx[,1],col=colors[5],names.arg = c("no \nrelation","same \ncompartment","same \ncomplex","same \nprotein"),ylab="% of pairs with a shared QTL",main="phQTL, exclusive conditions")
  barplot(sharedQtlEx[,2],col=colors[4],names.arg = c("no \nrelation","same \ncompartment","same \ncomplex","same \nprotein"),ylab="% of pairs with a shared QTL",main="phResQTL, exclusive conditions")
}))

###shared QTL for phosphosites for the same and different proteins
pdfAndPng("graph/sharedPhosphoQTLSameProtein", width = 7.5, height = 10, expr = expression({
  par(mfrow=c(1,1),cex=2)
  barplot(c(exRes[[4]][2,2]/sum(exRes[[4]][2,1:2]),exRes[[4]][2,4]/sum(exRes[[4]][2,3:4]))*100,ylab="affected by the same QTL (%)",names.arg = c("in different\nproteins","in the\nsame protein"),xlab="phosphosites",col=colors[4])
}))
