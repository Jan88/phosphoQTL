load("data/phosphoLevelReplacementStrainsNoCor.RData")
load("data/protLevelsReplacementStrainsNoCor.RData")
colnames(prot) <- gsub("X","s",colnames(prot))
colnames(prot) <- gsub("r","",colnames(prot))
colnames(phosphoLevel) <- gsub("X","s",colnames(phosphoLevel))
colnames(phosphoLevel) <- gsub("r","",colnames(phosphoLevel))
sharedSamples <- intersect(colnames(prot),colnames(phosphoLevel))
sharedSamples <- setdiff(sharedSamples,"s1473_3")
strainLabels <- gsub("_.*","",sharedSamples)

usePep <- phospho2ProtRep[,2]%in%rownames(prot) 
regPeps <- phospho2ProtRep[usePep,1]

corByPep <- sapply(regPeps,FUN=function(pep){
  cor(phosphoLevel[pep,sharedSamples],prot[phospho2ProtRep[pep,2],sharedSamples],use="pair")
})

prelimLMs <- t(sapply(regPeps,FUN=function(pep){
  m <- lm(phosphoLevel[pep,sharedSamples]~prot[phospho2ProtRep[pep,2],sharedSamples])
  coeffs <- m$coefficients
  p <- summary(m)$coefficients["prot[phospho2ProtRep[pep, 2], sharedSamples]","Pr(>|t|)"]
  return(c(coeffs,p))
}))

plot(prelimLMs[,2],-log10(prelimLMs[,3]))
abline(v=c(-1,1),col="red")

scaledP <- t(scale(t(prot[,sharedSamples]),scale=T,center=T))
scaledPh <- t(scale(t(phosphoLevel[,sharedSamples]),scale=T,center=T))


#estimate coefficients
library(MASS)
mmPhosphoProt<-sapply(regPeps,function(pep){
  if(sum(is.na(scaledPh[pep,]+scaledP[phospho2ProtRep[pep,2],]))>6){return(c(NA,NA))}
  rlm(scaledPh[pep,]~scaledP[phospho2ProtRep[pep,2],],method = "MM")$coefficients
})
mmPhosphoProt<-t(mmPhosphoProt)
colnames(mmPhosphoProt) <- c("Intercept", "Slope")

phosphoProtResiduals <- sapply(regPeps,function(pep){
  scaledPh[pep,] - (mmPhosphoProt[pep,1] + mmPhosphoProt[pep,2]*scaledP[phospho2ProtRep[pep,2],])
})
phosphoProtResiduals<-t(phosphoProtResiduals)
rownames(phosphoProtResiduals) <- regPeps
save(phosphoProtResiduals,file="data/phosphoProtResidualsReplacementStrains.RData")

#quick look at heritability
varByStrain <- sapply(unique(strainLabels),FUN=function(s){
  apply(phosphoProtResiduals[,strainLabels==s],1,var,na.rm=T)
})

plot(apply(varByStrain,1,mean,na.rm=T),apply(phosphoProtResiduals,1,var,na.rm=T))
abline(a=0,b=1,col="red")
varByStrainPerm <- matrix(0,ncol=4,nrow=nrow(varByStrain))
for(i in 1:200){
  set.seed(i)
  strainLabelsPerm <- strainLabels[sample(1:length(strainLabels))]
  varByStrainPerm <- varByStrainPerm+sapply(unique(strainLabels),FUN=function(s){
    apply(phosphoProtResiduals[,strainLabelsPerm==s],1,var,na.rm=T)
  })
}
varByStrainPerm <- varByStrainPerm/i
varByStrainPerm[varByStrainPerm==0] <- NA
plot(apply(varByStrainPerm,1,mean,na.rm=T),apply(phosphoProtResiduals,1,var,na.rm=T))
abline(a=0,b=1,col="red")

#compute phospho prot ratios
useSamples <- setdiff(intersect(colnames(prot),colnames(phosphoLevel)),"s1473_3")
phosphoProtRatio <- phosphoLevel[regPeps,useSamples]-prot[phospho2ProtRep[regPeps,2],useSamples]

#quick look at heritability
varByStrain <- sapply(unique(strainLabels),FUN=function(s){
  apply(phosphoProtRatio[,strainLabels==s],1,var,na.rm=T)
})

plot(apply(varByStrain,1,mean,na.rm=T),apply(phosphoProtRatio,1,var,na.rm=T))
abline(a=0,b=1,col="red")
varByStrainPerm <- matrix(0,ncol=4,nrow=nrow(varByStrain))
for(i in 1:200){
  set.seed(i)
  strainLabelsPerm <- strainLabels[sample(1:length(strainLabels))]
  varByStrainPerm <- varByStrainPerm+sapply(unique(strainLabels),FUN=function(s){
    apply(phosphoProtRatio[,strainLabelsPerm==s],1,var,na.rm=T)
  })
}
varByStrainPerm <- varByStrainPerm/i
varByStrainPerm[varByStrainPerm==0] <- NA
plot(apply(varByStrainPerm,1,mean,na.rm=T),apply(phosphoProtRatio,1,var,na.rm=T))
abline(a=0,b=1,col="red")

#test comparisons
comparisons <- rbind(c("X1473","X1475"),c("X1473","X1476"),c("X1473","X1477"),c("X1475","X1477"),c("X1476","X1477"))
comparisons <- gsub("X","s",comparisons)

diffPhosphoProt <- lapply(1:nrow(comparisons),FUN=function(i){
  samples1 <- which(grepl(comparisons[i,1],colnames(phosphoProtResiduals)))
  samples2 <- which(grepl(comparisons[i,2],colnames(phosphoProtResiduals)))
  fcs <- rowMeans(phosphoProtResiduals[,samples1],na.rm=T)-rowMeans(phosphoProtResiduals[,samples2],na.rm=T)
  p <- apply(phosphoProtResiduals,1,FUN=function(x){
    if(sum(!is.na(x[samples1]))<3|sum(!is.na(x[samples2]))<3){return(NA)}
    t.test(x[samples1],x[samples2])$p.value
  })
  fdr <- p.adjust(p,method="BH")
  out <- cbind(p=p,fdr=fdr,fc=fcs)
  rownames(out) <- regPeps
  return(out)
})
names(diffPhosphoProt) <- apply(comparisons,1,paste,collapse="_")

diffPhosphoProtRatio <- lapply(1:nrow(comparisons),FUN=function(i){
  samples1 <- which(grepl(comparisons[i,1],colnames(phosphoProtRatio)))
  samples2 <- which(grepl(comparisons[i,2],colnames(phosphoProtRatio)))
  fcs <- rowMeans(phosphoProtRatio[,samples1],na.rm=T)-rowMeans(phosphoProtRatio[,samples2],na.rm=T)
  p <- apply(phosphoProtRatio,1,FUN=function(x){
    if(sum(!is.na(x[samples1]))<3|sum(!is.na(x[samples2]))<3){return(NA)}
    t.test(x[samples1],x[samples2])$p.value
  })
  fdr <- p.adjust(p,method="BH")
  out <- cbind(p=p,fdr=fdr,fc=fcs)
  rownames(out) <- regPeps
  return(out)
})
names(diffPhosphoProtRatio) <- apply(comparisons,1,paste,collapse="_")
save(diffPhosphoProt,diffPhosphoProtRatio,file="data/diffPhosphoProtReplacementStrainsNoCor.RData")
