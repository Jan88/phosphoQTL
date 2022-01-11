library(gtools)
library(beeswarm)
source("lib/general_function.R")

# load data -------
meta <- read.table("metadata/metadata.csv", sep=",", header=T)
rownames(meta) <- meta$culture
meta    <- meta[mixedsort(rownames(meta)), ]
meta <- meta[meta$matching_with_brem_genotype & !meta$suspect_multiple_strain,]
load("data/phosphoLevel.RData")
load("data/proteinLevel.RData")
load("data/expressionLevel.RData")

#write(phospho2prot[,1], file="data/PhosphoPeptide.txt")


# indetify usable shared samples -----------
sharedSamples <- mixedsort(intersect(intersect(intersect(colnames(proteinLevelBatchCorrected), colnames(phosphoLevelBatchCorrected)),rownames(meta)),colnames(eLevelBatchCorrected)))
meta <- meta[mixedsort(sharedSamples),]

prot <- as.matrix(proteinLevelBatchCorrected[,sharedSamples])
phospho <- as.matrix((phosphoLevelBatchCorrected[,sharedSamples]))
rna <- as.matrix(eLevelBatchCorrected[,sharedSamples])


###  explore correlation between phosopho prot and RNA ----------
sum(unique(phospho2prot[,2]) %in% rownames(prot))  # 402
sum(unique(phospho2prot[,2])%in% rownames(prot))/length(unique(phospho2prot[,2])) # 0.4068826
sum((phospho2prot[,2]) %in% rownames(prot))  #  879
sum((phospho2prot[,2]) %in% rownames(prot))/nrow(phospho2prot) #  0.4154064
sum((phospho2prot[,2]) %in% rownames(rna))/nrow(phospho2prot)
sum((rownames(prot) %in% rownames(rna))/nrow(prot))

phosphoProtCor<- sapply(1:nrow(phospho),function(i){
  if(phospho2prot[i,2] %in% rownames(prot)) {
    cor(phospho[i,],prot[phospho2prot[i,2],],method = "pear",use = "pair" )
  } else (NA)
})

phosphoPhosphoCor<- sapply(unique(phospho2prot[,2]), function(p) {
  sel <- phospho2prot[,2]==p
  if(sum(sel)==1) return(NA) 
  out <- cor(t(phospho[sel,]), method = "pear",use = "pair")
  out <- out[upper.tri(out,diag = F)]
})
phosphoPhosphoCor<-unlist(phosphoPhosphoCor)

phosphoRnaCor<- sapply(1:nrow(phospho),function(i){
  if(phospho2prot[i,2] %in% rownames(rna)) {
    cor(phospho[i,],rna[phospho2prot[i,2],],method = "pear",use = "pair" )
  } else (NA)
})

protRnaCor <- sapply(rownames(prot),function(i){
  if(i %in% rownames(rna)) {
    cor(prot[i,],rna[i,],method = "pear",use = "pair" )
  } else (NA)
})


phosphoProtCorRandom<- sapply(1:nrow(phospho),function(i){
  if(phospho2prot[i,2] %in% rownames(prot)) {
  cor(phospho[sample(1:nrow(phospho))[1],],prot[phospho2prot[i,2],],method = "pear",use = "pair" )
  } else (NA)
})

phosphoPhosphoCorRandom<- sapply(unique(phospho2prot[,2]), function(p) {
  sel <- sample(phospho2prot[,2]==p)
  if(sum(sel)==1) return(NA) 
  out <- cor(t(phospho[sel,]), method = "pear",use = "pair")
  out <- out[upper.tri(out,diag = F)]
})
phosphoPhosphoCorRandom<-unlist(phosphoPhosphoCorRandom)

phosphoRnaCorRandom<- sapply(1:nrow(phospho),function(i){
  if(phospho2prot[i,2] %in% rownames(rna)) {
    cor(phospho[sample(1:nrow(phospho))[1],],rna[phospho2prot[i,2],],method = "pear",use = "pair" )
  } else (NA)
})

protRnaCorRandom <- sapply(rownames(prot),function(i){
  if(i %in% rownames(rna)) {
    cor(prot[sample(1:nrow(prot))[1],],rna[i,],method = "pear",use = "pair" )
  } else (NA)
})

pdfAndPng("graph/phospho_prot_correlation_2", 8,8, expression({      
  beeswarm(list(phosphoProtCor, phosphoPhosphoCor), pch=16, col=rgb(.3,.3,.3,.8), cex = .4, ylim=c(-1,1),
           labels = c('phospho-petide\nto matching protein', "phospho-petide\nto matching phospho-petide"), ylab=expression(paste("pearson's corretaltion ",rho))) 
  boxplot(phosphoProtCorRandom, phosphoPhosphoCorRandom,   add=T,col="#0000ff44",boxwex=.10,
          xaxt="n", yaxt="n",notch = F,outline = F)
  abline(h=0,lty=3)
}))

pdfAndPng("graph/phospho_prot_correlation_4", 10,8, expression({      
  beeswarm(list(phosphoPhosphoCor,phosphoProtCor,phosphoRnaCor,protRnaCor), pch=16, col=rgb(.3,.3,.3,.8), cex = .4, ylim=c(-1,1),
           labels = c("phospho-petide\nto matching phospho-petide",'phospho-petide\nto matching protein',"phospho-petide\nto matching RNA","protein\nto matching RNA" ),cex.axis=.9, ylab=expression(paste("pearson's corretaltion ",rho))) 
  boxplot(phosphoPhosphoCorRandom,phosphoProtCorRandom,  phosphoRnaCorRandom,protRnaCorRandom,  add=T,col="#0000ff44",boxwex=.10,
          xaxt="n", yaxt="n",notch = F,outline = F)
  abline(h=0,lty=3)
}))

pdfAndPng("graph/phospho_prot_correlation_paper", 10,8, expression({      
beeswarm(list(phosphoPhosphoCor,phosphoProtCor,protRnaCor), pch=16, col=rgb(0,0,0,.8), cex = .4, ylim=c(-1,1),
         labels = c("phospho-petide\nto matching phospho-petide",'phospho-petide\nto matching protein',"protein\nto matching RNA" ),cex.axis=.9, ylab=expression(paste("pearson's corretaltion ",rho))) 
abline(h=0,lty=3)
}))

wilcox.test(phosphoProtCor,phosphoRnaCor, paired = T)


require(VennDiagram)
pdfAndPng("graph/phospho_prot_overlap_2", 8,8, expression({      
grid.newpage()
temp<-venn.diagram(list(phospho = unique(phospho2prot[,2]), rownames(prot)),fill = c("red", "green"),
             alpha = c(0.5, 0.5), cex = 2,lty =2, lwd=0, filename = NULL,category.names = c("proteins with\nphosphoPeptides","measured proteins"), cat.cex=1.5, margin=.1, cat.fontface=2,cat.fontfamily = rep("sans", 2), ,fontfamily=rep("sans",3))
grid.draw(temp)
}))





###### Generate phospho to prot ratio AND phospho-prot regression--------
# indetify usable shared samples 
sharedSamples <- mixedsort(intersect(intersect(intersect(colnames(proteinLevelBatchCorrected), colnames(phosphoLevelBatchCorrected)),rownames(meta)),colnames(eLevelBatchCorrected)))
meta <- meta[mixedsort(sharedSamples),]

sel <- phospho2prot[,2] %in% rownames(prot)
metaPhosphoProtRatio <- phospho2prot[sel,]

phospho <- as.matrix((phosphoLevelBatchCorrected[,sharedSamples]))[sel,]
prot <- as.matrix(proteinLevelBatchCorrected[,sharedSamples])[metaPhosphoProtRatio[,2],]




# phospho to prot ratio
phosphoProtRatio<- sapply(1:nrow(phospho),function(i){
    # Normalization by Z-transformation of the phospho and protein data.
    # The Z transformation parameter (ie mean and sd) are calculalted by 
    # only taking into account the values that not missing in both p and phospho 
    ph <- phospho[i,]
    p <- prot[i,]
    sel <- !is.na(ph) & !is.na(p)
    ph <- (ph - mean(ph[sel])) / sd(ph[sel])
    p <- (p - mean(p[sel])) / sd(p[sel])
    return(ph -p)
})
phosphoProtRatio <- t(phosphoProtRatio)



###### regression 

# Scale the data for propar comparison of RNA and prot level
protScale <- phosphoScale <- matrix(NA,nrow(phospho),150)
for(i in 1:nrow(phospho)){ 
  ph <- phospho[i,]
  p <- prot[i,]
  sel <- !is.na(ph) & !is.na(p)
  ph <- (ph - mean(ph[sel])) / sd(ph[sel])
  p <- (p - mean(p[sel])) / sd(p[sel])
  protScale[i,] <- p
  phosphoScale[i,] <- ph
}

# corPhosphoProt<-sapply(1:nrow(phosphoScale),function(i){
#   unlist(cor.test(phosphoScale[i,], protScale[i,],use = "pair")[c("estimate","p.value")])
# })
# corPhosphoProt <- t(corPhosphoProt)
# # select the phospho and prot significaltly correlated at at least 15% FRD
# sel <- p.adjust(corPhosphoProt[,2],method = "BH")<0.15

library(MASS)
mmPhoshoProt<-sapply(1:nrow(phosphoScale),function(i){
  rlm( phosphoScale[i,]~protScale[i,],method = "MM")$coefficients
})
mmPhoshoProt<-t(mmPhoshoProt)
colnames(mmPhoshoProt) <- c("Intersecpt", "Slope")

pdf("graph//phospho_prot_regression_mm.pdf",20,10)
par(mfrow=c(4,8),mar=c(4,4,3,1),cex.main = 0.9)
sapply(1:nrow(phosphoScale),function(i){
  plot(x=protScale[i,], y=phosphoScale[i,],  xlab="prot", ylab="phospho",
       main=paste(metaPhosphoProtRatio[i,1],metaPhosphoProtRatio[i,2],sep="-"), 
       pch=21, cex=1, bg="grey", xlim=range(c(protScale[i,], y=phosphoScale[i,]),na.rm=T))
  abline(mmPhoshoProt[i,], col="red", lwd=3)
  
})
dev.off()

# Compute the residuals
phosphoProtResiduals <- sapply(1:nrow(phosphoScale),function(i){
  phosphoScale[i,] - (mmPhoshoProt[i,1] + mmPhoshoProt[i,2]*protScale[i,])
})
phosphoProtResiduals<-t(phosphoProtResiduals)

rownames(phosphoProtRatio) <- rownames(phosphoProtResiduals) <- metaPhosphoProtRatio$peptide
colnames(phosphoProtRatio) <- colnames(phosphoProtResiduals) <- sharedSamples

metaPhosphoProtResiduals<-metaPhosphoProtRatio


save(phosphoProtRatio,phosphoProtResiduals, metaPhosphoProtRatio,metaPhosphoProtResiduals, file="data/phosphoProt.RData")
write(metaPhosphoProtResiduals[,1], file="data/phosphoProtResisualPeptide.txt")

###### Generate phospho to RNA regression--------
sharedSamples <- mixedsort(intersect(intersect(intersect(colnames(proteinLevelBatchCorrected), colnames(phosphoLevelBatchCorrected)),rownames(meta)),colnames(eLevelBatchCorrected)))
meta <- meta[mixedsort(sharedSamples),]

sel <- phospho2prot[,2] %in% rownames(eLevelBatchCorrected)
metaPhosphoRnaResiduals <- phospho2prot[sel,]

phospho <- as.matrix((phosphoLevelBatchCorrected[,sharedSamples]))[sel,]
rna <- as.matrix(eLevelBatchCorrected[,sharedSamples])[metaPhosphoRnaResiduals[,2],]


# Scale the data for propar comparison of RNA and prot level
rnaScale <- phosphoScale <- matrix(NA,nrow(phospho),150)
for(i in 1:nrow(phospho)){ 
  ph <- phospho[i,]
  p <- rna[i,]
  sel <- !is.na(ph) & !is.na(p)
  ph <- (ph - mean(ph[sel])) / sd(ph[sel])
  p <- (p - mean(p[sel])) / sd(p[sel])
  rnaScale[i,] <- p
  phosphoScale[i,] <- ph
}

library(MASS)
mmPhoshoRna<-sapply(1:nrow(phosphoScale),function(i){
  rlm( phosphoScale[i,]~rnaScale[i,],method = "MM")$coefficients
})
mmPhoshoRna<-t(mmPhoshoRna)
colnames(mmPhoshoRna) <- c("Intersecpt", "Slope")

phosphoRnaResiduals <- sapply(1:nrow(phosphoScale),function(i){
  phosphoScale[i,] - (mmPhoshoRna[i,1] + mmPhoshoRna[i,2]*rnaScale[i,])
})
phosphoRnaResiduals<-t(phosphoRnaResiduals)

rownames(phosphoRnaResiduals) <- metaPhosphoRnaResiduals$peptide
colnames(phosphoRnaResiduals) <- sharedSamples



save(phosphoRnaResiduals, metaPhosphoRnaResiduals, file="data/phosphoRna.RData")

write(metaPhosphoRnaResiduals[,1], file="data/PhosphoRnaResidualPeptide.txt")

###### Generate prot to RNA regression--------

sharedSamples <- mixedsort(intersect(intersect(intersect(colnames(proteinLevelBatchCorrected), colnames(phosphoLevelBatchCorrected)),rownames(meta)),colnames(eLevelBatchCorrected)))
meta <- meta[mixedsort(sharedSamples),]

prot <- as.matrix(proteinLevelBatchCorrected[,sharedSamples])
phospho <- as.matrix((phosphoLevelBatchCorrected[,sharedSamples]))
rna <- as.matrix(eLevelBatchCorrected[,sharedSamples])


sel <- intersect(rownames(prot), rownames(rna))

prot <- prot[sel,]
rna <- rna[sel,]


# Scale the data for propar comparison of RNA and prot level
rnaScale <- protScale <- matrix(NA,nrow(prot),150)
for(i in 1:nrow(prot)){ 
  p <- prot[i,]
  r <- rna[i,]
  sel <- !is.na(p) & !is.na(r)
  p <- (p - mean(p[sel])) / sd(p[sel])
  r <- (r - mean(r[sel])) / sd(r[sel])
  rnaScale[i,] <- r
  protScale[i,] <- p
}



library(MASS)
mmProtRna<-sapply(1:nrow(protScale),function(i){
  rlm( protScale[i,]~rnaScale[i,],method = "MM")$coefficients
})
mmProtRna<-t(mmProtRna)
colnames(mmProtRna) <- c("Intersecpt", "Slope")

plot(rnaScale[18,], protScale[18,])
abline(mmProtRna[18,])


pdfAndPng("graph/prot_rna_regression",6,6,expression({
hist(mmProtRna[,2], col="grey",breaks=20, main="", xlab="Protein-RNA resgression slope")
abline(v=0, col="red",lwd=2)
}))

protRnaResiduals <- sapply(1:nrow(protScale),function(i){
  protScale[i,] - (mmProtRna[i,1] + mmProtRna[i,2]*rnaScale[i,])
})
protRnaResiduals<-t(protRnaResiduals)

rownames(protRnaResiduals) <- rownames(prot)
colnames(protRnaResiduals) <- sharedSamples



save(protRnaResiduals, file="data/protRna.RData")



