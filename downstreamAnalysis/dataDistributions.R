library(RColorBrewer)


# load mapping data #####
load("data/PerlsteinData/Perlstein_mappingData.RData")
strains = rownames(mappingData$genotype)
load("data/eQtlMappingData160817.RData")
ePheno = phenotype[,strains]
load("data/ptQtlMappingData170801.RData")
ptPheno = phenotype[,strains]
load("data/pQtlMappingData161117.RData")
pPheno = phenotype[,strains]
load("data/phosphoProtResidualsQtlMappingData170727.RData")
phResPheno = phenotype[,strains]
load("data/phosphoLevelQtlMappingData170727.RData")
phPheno = phenotype[,strains]
rm(phenotype4perm, phenotype, genotype)
#####



# get set of common genes #####
transcripts = rownames(ePheno)
ptTraits = rownames(ptPheno)
proteins = rownames(pPheno)
load("data/phosphoLevel.RData")
phospho2protHash = phospho2prot[,2]
names(phospho2protHash) = phospho2prot[,1]
prot2phosphoHash = split(phospho2prot[,1], phospho2prot[,2])
rm(phosphoLevelBatchCorrected, phosphoLevelNonBatch)
phosphoTargets = rownames(phPheno)
phosphoTargets = phospho2protHash[phosphoTargets]
phosphoTargets = unique(phosphoTargets)
phosphoResTargets = rownames(phResPheno)
phosphoResTargets = phospho2protHash[phosphoResTargets]
phosphoResTargets = unique(phosphoResTargets)
genesAvail = Reduce(intersect, list(transcripts,
                                    ptTraits,
                                    proteins,
                                    phosphoTargets,
                                    phosphoResTargets))
ePhenoSub = ePheno[genesAvail,]
ptPhenoSub = ptPheno[genesAvail,]
pPhenoSub = pPheno[genesAvail,]
pepAvail = unlist(prot2phosphoHash[genesAvail])
phResPhenoAverage = t(sapply(genesAvail, function(x){
   pepAvail = unlist(prot2phosphoHash[x])
   sub = phResPheno[rownames(phResPheno) %in% pepAvail,, drop = F]
   return(colMeans(sub))
}))
phPhenoAverage = t(sapply(genesAvail, function(x){
   pepAvail = unlist(prot2phosphoHash[x])
   sub = phPheno[rownames(phPheno) %in% pepAvail,, drop = F]
   return(colMeans(sub))
}))
phResPhenoSub = phResPheno[rownames(phResPheno) %in% pepAvail,]
phPhenoSub = phPheno[rownames(phPheno) %in% pepAvail,]
#####



# visualize distributions of molecular layers ######
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
png("graph/compareTraitDistributions.png", height=1500, width=2000, pointsize=25)
par(mfrow = c(8,8), mar = c(1.5,2,0.9,0.1), oma = c(2,2,0,0))
for(i in sample(1:nrow(ePhenoSub),63)){
   d1 = density(ePhenoSub[i,], na.rm=T)
   d2 = density(ptPhenoSub[i,], na.rm=T)
   d3 = density(pPhenoSub[i,], na.rm=T)
   d4 = density(phResPhenoAverage[i,], na.rm=T)
   d5 = density(phPhenoAverage[i,], na.rm=T)
   maxval = max(c(d1$y, d2$y, d3$y, d4$y, d5$y))
   plot(d1, main = rownames(ePhenoSub)[i], mgp = c(1.4,0.5,0), 
        xlab = "", ylab = "", ylim = c(0,maxval+0.01),
        col = colors["eQTL"], lwd = 3, las = 1)
   lines(d2, col = colors["ptQTL"], lwd = 3)
   lines(d3, col = colors["pQTL"], lwd = 3)
   lines(d4, col = colors["phResQTL"], lwd = 3)
   lines(d5, col = colors["phQTL"], lwd = 3)
}
plot.new()
legend("topleft", legend=c("transcript", "pt", "protein", "phRes", "phospho"), 
       col = colors[1:5], lwd = 3, xpd = NA)
mtext(text="trait value", side=1, outer = T, cex = 1.3, line= 0.5)
mtext(text="density", side=2, outer = T, cex = 1.3, line= 0.2)
dev.off()
######

# visualize (and test) normal distribution of growth traits. #####
load("data/PerlsteinData/Perlstein_mappingData.RData")
pheno = mappingData$phenotype
png("graph/PerlsteinDistributions.png", 
    height=1500, width=2000, pointsize=25)
par(mfrow = c(10,10), mar = c(1,1,0.9,0.1), oma = c(2,2,0,0))
for(i in sample(1:nrow(pheno),100)){
   plot(density(pheno[i,], na.rm=T), main = rownames(pheno)[i], 
        cex.main = 0.7,cex.axis = 0.5,tcl = -0.2,
        mgp = c(1.4,0.2,0),xlab = "", ylab = "",  lwd = 3, las = 1)
}
mtext(text="trait value", side=1, outer = T, cex = 1.3, line= 0.5)
mtext(text="density", side=2, outer = T, cex = 1.3, line= 0.2)
dev.off()
#####
