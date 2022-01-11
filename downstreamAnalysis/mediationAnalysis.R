# library(data.table)
library(RColorBrewer)
# library(randomForest)
# library(ranger)
# source("lib/general_function.R")
library(mediation)
library(xlsx)
library(ggplot2)

# functions for mediation #####
# calcMediation = function(marker, mediator, outcomes){
#   res = lapply(1:nrow(outcomes), function(i){
#     outcome = outcomes[i,]
#     dat = data.frame(outcome, marker , mediator)
#     dat = dat[!is.na(dat$outcome),]
#     fitM <- lm(mediator ~ marker, data=dat)
#     fitComb <- lm(outcome ~ marker + mediator, data=dat) 
#     fitMedBoot <- mediate(fitM, fitComb, boot=T, sims=999, 
#                           treat="marker", mediator="mediator")
#     temp = fitMedBoot[c("d0", "d0.ci", "d0.p",
#                         "z0", "z0.ci", "z0.p",
#                         "tau.coef", "tau.ci",  "tau.p",
#                         "n0", "n0.ci", "n0.p")]
#     matrix(unlist(temp), 
#            ncol = 4, byrow=T, 
#            dimnames=list(c("ACME", "ADE", "Total Effect", "Prop.Mediated"),
#                          c("Estimate", "2.5% CI", "97.5% CI", "p-value")))
#   })
#   names(res) = rownames(outcomes)
#   return(res)
# }
# # a version without bootstraps (i.e. Sobel test) and with parallelization
# calcMediationFast = function(marker, mediator, outcomes, nthread){
  require(parallel)
  cl = makeCluster(nthread, type="FORK")
  res = parLapply(cl, 1:nrow(outcomes), function(i){
    outcome = outcomes[i,]
    dat = data.frame(outcome, marker , mediator)
    dat = dat[!is.na(dat$outcome),]
    fitM <- lm(mediator ~ marker, data=dat)
    fitComb <- lm(outcome ~ marker + mediator, data=dat) 
    fitMedBoot <- mediate(fitM, fitComb, boot=F, 
                          treat="marker", mediator="mediator")
    temp = fitMedBoot[c("d0", "d0.ci", "d0.p",
                        "z0", "z0.ci", "z0.p",
                        "tau.coef", "tau.ci",  "tau.p",
                        "n0", "n0.ci", "n0.p")]
    matrix(unlist(temp), 
           ncol = 4, byrow=T, 
           dimnames=list(c("ACME", "ADE", "Total Effect", "Prop.Mediated"),
                         c("Estimate", "2.5% CI", "97.5% CI", "p-value")))
  })
  names(res) = rownames(outcomes)
  stopCluster(cl)
  return(res)
#}
#####

# prepare data #####
# mapping input data/phenotypes 
fullGeno = read.table("data/genotype_for_mapping.tsv", header = T, as.is=T)
load("data/PerlsteinData/Perlstein_genotypes.RData")
strainGeno = genotype
load("data/PerlsteinData/Perlstein_mappingData.RData")
gPheno = mappingData$phenotype
strains = rownames(mappingData$genotype)
gGeno = mappingData$genotype
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
rm(phenotype4perm, phenotype, genotype, mappingData)

# QTLs 
load("data/eQTL_results160831.RData")
eQTL <- QtlList$FDR10
load("data/pQTL_results170815.RData")
ptQTL <- pQTL_results$ptQTL$QtlList$FDR10
pQTL <- pQTL_results$pQTL$QtlList$FDR10
phResQTL <- pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10
phQTL <- pQTL_results$phosphoLevelQTL$QtlList$FDR10
rm(pQTL_results, pv, qv,QtlList)
load("data/PerlsteinData/Perlstein_QTL_list_merged.RData")
gQTL = lapply(QTL_list_merged, function(x){
   x = data.frame(x,stringsAsFactors=F)
   x$start = as.integer(x$start)
   x$end = as.integer(x$end)
   x$mostSignificantPredictor = as.integer(x$mostSignificantPredictor)
   x$minP = as.numeric(x$minP)
   return(x)
})


# hotspots
# load("data/otherQTLstudies/Albert2018Candidates.RData")
load("data/hotspotLocations.RData") # "leader" and "markerLocations"
load("data/hsInfo.RData") # "hsOvMat", "targetsByHS", "GOanalysis" 
hsOvMat=data.frame(hsOvMat[,1:2],
                   apply(hsOvMat[,3:11],2,as.integer),
                   apply(hsOvMat[,12:16],2,as.logical), 
                   stringsAsFactors=F)
markerLocations$CHR = as.character(markerLocations$CHR)
rm(GOanalysis)                   

# table to look up gene positions
anno = read.table("Saccharomyces_cerevisiae/SGD_features.tab",
                  header = F,sep = "\t",quote="",as.is=T)
colnames(anno) = c("ID", "type", "status", "name_systematic", 
                       "name_standard", "alias", "parent", "altID", 
                       "chr", "start", "end", "strand", "genetic_pos",
                       "coord_version", "seq_version", "description")
id2name = setNames(anno$name_standard[anno$name_systematic != ""], 
                   anno$name_systematic[anno$name_systematic != ""])
id2name[id2name == ""] = names(id2name[id2name == ""])
name2id = setNames(anno$name_systematic[anno$name_systematic != ""], 
                   anno$name_standard[anno$name_systematic != ""])
name2id[name2id == ""] = names(name2id[name2id == ""])
geneAnno = anno[anno$type == "ORF" & anno$chr != "2-micron",]
geneAnno$chr = paste0("chr", as.roman(geneAnno$chr))

# create translation between peptide and protein names
load("data/phosphoLevel.RData")
phospho2protHash = setNames(phospho2prot[,2], phospho2prot[,1]) # setNames can be used to create a named vector
prot2phosphoHash = split(phospho2prot[,1], phospho2prot[,2])
rm(phosphoLevelBatchCorrected, phosphoLevelNonBatch, phospho2prot)
# for phospho traits: average over peptides belonging to the same protein
phResPhenoAverage = t(sapply(names(prot2phosphoHash), function(x){
   sub = phResPheno[rownames(phResPheno) %in% unlist(prot2phosphoHash[x]),, drop = F]
   return(colMeans(sub, na.rm= T))
}))
phPhenoAverage = t(sapply(names(prot2phosphoHash), function(x){
   sub = phPheno[rownames(phPheno) %in% unlist(prot2phosphoHash[x]),, drop = F]
   return(colMeans(sub, na.rm=T))
}))

# get set of common genes/peptides
# genesAvail = Reduce(intersect, list(rownames(ePheno),
#                                     rownames(ptPheno),
#                                     rownames(pPheno),
#                                     unique(phospho2protHash[rownames(phPheno)]),
#                                     unique(phospho2protHash[rownames(phResPheno)])))
# pepAvail = unlist(prot2phosphoHash[genesAvail])

# subset phenotype tables
# ePhenoSub = ePheno[genesAvail,]
# ptPhenoSub = ptPheno[genesAvail,]
# pPhenoSub = pPheno[genesAvail,]
# phResPhenoSub = phResPheno[rownames(phResPheno) %in% pepAvail,]
# phPhenoSub = phPheno[rownames(phPheno) %in% pepAvail,]
# phResPhenoAverageSub =  phResPhenoAverage[genesAvail,]
# phPhenoAverageSub =  phPhenoAverage[genesAvail,]


#####

# mediation analysis for Ste20 hs (8chrVIII:2) through gpa1 and ste20 #######
hs = "chrVIII:1"
hsInfo = unlist(hsOvMat[hsOvMat$name == hs,])
leadMarker_e = leader[which(hsOvMat$name == hs),"e"]
marker_e = strainGeno[,leadMarker_e]
targets = targetsByHS[[which(hsOvMat$name == hs)]]
lay2Pheno = c("e" = "ePheno", "pt" = "ptPheno", "p" = "pPheno", 
              "ph" = "phPhenoAverage", "phr" = "phResPhenoAverage")
outcomes = do.call(rbind,sapply(names(targets), function(x){
  temp = get(lay2Pheno[x])[targets[[x]],,drop = F]
  rownames(temp) = paste(rownames(temp),x, sep="_")
  return(temp)
}))

# GPA1
GPA1 = "YHR005C" #GPA1
mediations_GPA1_e = calcMediation(marker=marker_e,
                                  mediator=ePheno[GPA1,], 
                                  outcomes=outcomes)
mediations_GPA1_p = calcMediation(marker=marker_e,
                                  mediator=pPheno[GPA1,], 
                                  outcomes=outcomes)
mediations_GPA1_ph = calcMediation(marker=marker_e,
                                   mediator=phPhenoAverage[GPA1,], 
                                   outcomes=outcomes)
save(mediations_GPA1_e, mediations_GPA1_p, file = "data/mediations_GPA1.RData")
# STE20
STE20 = "YHL007C"
mediations_STE20_e = calcMediation(marker=marker_e,
                                   mediator=ePheno[STE20,], 
                                   outcomes=outcomes)
mediations_STE20_p = calcMediation(marker=marker_e,
                                   mediator=pPheno[STE20,], 
                                   outcomes=outcomes)
mediations_STE20_ph = calcMediation(marker=marker_e,
                                    mediator=phPhenoAverage[STE20,], 
                                    outcomes=outcomes)
save(mediations_STE20_e, mediations_STE20_p, mediations_STE20_ph, 
     file = "data/mediations_STE20.RData")
load("data/mediations_GPA1.RData")
load("data/mediations_STE20.RData")
mediation_pvals = data.frame("GPA1_transcript" = sapply(mediations_GPA1_e, function(x){x[1,4]}),
                             "GPA1_protein" = sapply(mediations_GPA1_p, function(x){x[1,4]}),
                             "STE20_transcript" = sapply(mediations_STE20_e, function(x){x[1,4]}),
                             "STE20_protein" = sapply(mediations_STE20_p, function(x){x[1,4]}),
                             "STE20_phospho" = sapply(mediations_STE20_ph, function(x){x[1,4]}))
mediation_pvals[mediation_pvals == 0] = 2*10^(-16)
png("graph/mediation_Gpa1hotspot.png", height=800, width=1400, pointsize=28)
par(mfrow = c(2,3), mar = c(4,4,0,0), oma = c(2,2,3,0.1))
for(i in 1:2){
  for(j in 3:5){
    plot(-log10(mediation_pvals[,j]), 
         -log10(mediation_pvals[,i]),
         xlab = "-log10(p-value)",
         ylab = "-log10(p-value)", 
         cex.lab = 1, las = 1, mgp = c(2.2,1,0),
         xlim = c(0,max(-log10(mediation_pvals[,3:5]))))
    abline(0,1)
    if(i == 2){
      mtext(1,line = 4.5, cex = 0.8,font = 2,
            text = paste0("mediation through\n",
                          colnames(mediation_pvals)[j]))}
    if(j == 3){
      mtext(2, line = 3.5, cex = 0.8, font = 2,
            text = paste0("mediation through\n",
                          colnames(mediation_pvals)[i]))}
  }
}
mtext(3,outer=T, text = "mediation of chrVIII:1 hotspot targets", 
      font=2, line = 1)
dev.off()
#######

# Ste20/Gpa1 hotspot: mediation of intermediate ph/p/e traits on 
# downstream targets (transcripts) #####
hs = "chrVIII:1"
hsInfo = unlist(hsOvMat[hsOvMat$name == hs,])
leadMarker_e = leader[which(hsOvMat$name == hs),"e"]
marker_e = strainGeno[,leadMarker_e]
# first, targets of the mating response
matTargets = c("Apa2" = "YDR530C", "Aga1" = "YNR044W", 
               "Fus1" = "YCL027W", "Fus2" = "YMR232W", 
               "Hym1" = "YKL189W","Kar4" = "YCL055W", 
               "Mfa1" = "YDR461W", "Prm1" = "YNL279W",
               "Prm4" = "YPL156C", "Prm5" = "YIL117C", 
               "Prp39" = "YML046W", "Hal1" = "YPR005C",
               "Pst1" = "YDR055W")
matCandidates = c("Dig1" = "YPL049C", "Dig2" = "YDR480W",
                  "Fus3" = "YBL016W", "Mnr2" = "YKL064W",
                  "Far1" = "YJL157C", "Tec1" = "YBR083W",
                  "Ste20" = "YHL007C", "Gpa1" = "YHR005C",
                  "Ste12"="YHR084W", "Ste7"="YDL159W", 
                  "Ste11"="YLR362W", "Ste5"="YDR103W")
colLayer <-brewer.pal(6,name="Paired")[-1]
names(colLayer) <-c("ePheno","ptPheno","pPheno","phResPhenoAverage", "phPhenoAverage")
comb = expand.grid(target = names(matTargets),
                   mediator = names(matCandidates), 
                   layer = names(colLayer),
                   stringsAsFactors=F)
comb$targetID = matTargets[comb$target]
comb$mediatorID = matCandidates[comb$mediator]
medRes = t(apply(comb,1, function(x){
  layer = x["layer"]
  if(!x["mediatorID"] %in% rownames(get(layer))){
    return(c(NA, NA))
  }
  dat = data.frame(outcome = ePheno[x["targetID"],], 
                   marker = marker_e, 
                   mediator = get(layer)[x["mediatorID"],])
  dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
  if(nrow(dat) < 3){
    return(c(NA,NA))
  }
  fitM <- lm(mediator ~ marker, data=dat)
  fitComb <- lm(outcome ~ marker + mediator, data=dat) 
  fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                         treat="marker", mediator="mediator"))
  res = unlist(fitMedBoot[c("d0", "d0.p")])
  return(res)
  # fitMedBoot["d0.p"]
}))
medRes = cbind(comb, medRes)
colnames(medRes)[6:7] = c("ACME", "ACME_pVal")
medRes$ACME_pVal[medRes$ACME_pVal == 0] = 2e-16
medRes$layer = factor(medRes$layer,levels=names(colLayer))
medRes$neg_log10_pVal = -log10(medRes$ACME_pVal)
medRes$significant = medRes$ACME_pVal < 0.05
save(medRes, file="data/mediations_matingresponse.RData")


# switch mediator and outcome variables:
medResRev = t(apply(comb,1, function(x){
  layer = x["layer"]
  if(!x["mediatorID"] %in% rownames(get(layer))){
    return(c(NA, NA))
  }
  dat = data.frame(mediator = ePheno[x["targetID"],], 
                   marker = marker_e, 
                   outcome  = get(layer)[x["mediatorID"],])
  dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
  if(nrow(dat) < 3){
    return(c(NA,NA))
  }
  fitM <- lm(mediator ~ marker, data=dat)
  fitComb <- lm(outcome ~ marker + mediator, data=dat) 
  fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                         treat="marker", mediator="mediator"))
  res = unlist(fitMedBoot[c("d0", "d0.p")])
  return(res)
  # fitMedBoot["d0.p"]
}))
medResRev = cbind(comb, medResRev)
colnames(medResRev)[6:7] = c("ACME", "ACME_pVal")
medResRev$ACME_pVal[medResRev$ACME_pVal == 0] = 2e-16
medResRev$layer = factor(medResRev$layer,levels=names(colLayer))
medResRev$neg_log10_pVal = -log10(medResRev$ACME_pVal)
medResRev$significant = medResRev$ACME_pVal < 0.05
save(medResRev, file="data/mediations_matingresponse_reversed.RData")

# now plot
load("data/mediations_matingresponse.RData")
load("data/mediations_matingresponse_reversed.RData")
colLayer <-brewer.pal(6,name="Paired")[-1]
names(colLayer) <-c("ePheno","ptPheno","pPheno","phResPhenoAverage", "phPhenoAverage")
g1 = ggplot(medRes, aes(x=mediator, y=neg_log10_pVal, fill=layer)) +
  geom_boxplot() +
  scale_fill_manual(values=colLayer)+
  geom_jitter()+
  ggtitle("true mediations")
g2 = ggplot(medResRev, aes(x=mediator, y=neg_log10_pVal, fill=layer)) +
  geom_boxplot() +
  scale_fill_manual(values=colLayer)+
  geom_jitter() +
  ggtitle("switched mediations")
medRes$mediator = factor(medRes$mediator, levels=c("Gpa1","Ste20",
                                        "Ste5", "Ste11","Ste7",
                                        "Fus3","Mnr2",
                                        "Far1","Dig1",  "Dig2",  
                                        "Ste12", "Tec1"))
g3 = ggplot(medRes, aes(x=layer, fill = significant))+
  geom_bar(stat="count") +
  facet_wrap(~ mediator, )+
  scale_x_discrete(labels=c("e", "pt", "p", "phRes", "ph")) + 
  theme_classic() +
  theme(#axis.text.x = element_text(angle = 45, hjust=1),
    strip.text=element_text(face = "bold"),
        strip.background = element_rect(size=1.5, linetype="blank")); g3
g4 = ggplot(medResRev, aes(x=layer, fill = significant))+
  geom_bar(stat="count") +
  facet_wrap(~ mediator)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_discrete(labels=c("e", "pt", "p", "phRes", "ph"))
d = data.frame(true_ACME = medRes$ACME, rev_ACME = medResRev$ACME, 
               true_pVal = medRes$ACME_pVal, rev_pVal = medResRev$ACME_pVal)
g5 = ggplot(d,aes(x=true_ACME, y=rev_ACME))+
  geom_point()
g6 = ggplot(d,aes(x=true_pVal, y=rev_pVal))+
  geom_point()
gcomb = gridExtra::grid.arrange(g1, g2, g3, g4,g5,g6)
ggsave(gcomb, filename= "graph/medRes.png", width = 8, height=8)
ggsave(g3, filename= "graph/medRes_barplot.png", width = 8, height=8)

# targets of the osmotic response:
osmTargets = c("Crh1"="YGR189C", "Vhs3" = "YOR054C",
               "Kar5" = "YMR065W", "Cdc20" = "YGL116W",
               "Rgc1" = "YPR115W")
osmCandidates = c("Msn2" = "YMR037C", "Msn4" = "YKL062W", 
                  "Rsc9" = "YML127W", "Ste20"= "YHL007C", 
                  "Gpa1" = "YHR005C", "Ssk1" = "YLR006C")
combOsm = expand.grid(target = names(osmTargets),
                      mediator = names(osmCandidates), 
                      layer = names(colLayer),
                      stringsAsFactors=F)
combOsm$targetID = osmTargets[combOsm$target]
combOsm$mediatorID = osmCandidates[combOsm$mediator]
osmRes = t(apply(combOsm,1, function(x){
  layer = x["layer"]
  if(!x["mediatorID"] %in% rownames(get(layer))){
    return(c(NA, NA))
  }
  dat = data.frame(outcome = ePheno[x["targetID"],], 
                   marker = marker_e, 
                   mediator = get(layer)[x["mediatorID"],])
  dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
  if(nrow(dat) < 3){
    return(c(NA,NA))
  }
  fitM <- lm(mediator ~ marker, data=dat)
  fitComb <- lm(outcome ~ marker + mediator, data=dat) 
  suppressMessages(fitMedBoot <- mediate(fitM, fitComb, boot=T, sims=999, 
                        treat="marker", mediator="mediator"))
  res = unlist(fitMedBoot[c("d0", "d0.p")])
  return(res)
}))
osmRes = cbind(combOsm, osmRes)
colnames(osmRes)[6:7] = c("ACME", "ACME_pVal")
osmRes$ACME_pVal[osmRes$ACME_pVal == 0] = 2e-16
osmRes$layer = factor(osmRes$layer,
                      levels=names(colLayer))
osmRes$neg_log10_pVal = -log10(osmRes$ACME_pVal)
osmRes$significant = osmRes$ACME_pVal < 0.05
save(osmRes, file="data/mediations_osmoticresponse.RData")
load("data/mediations_osmoticresponse.RData")
osmResRev = t(apply(combOsm,1, function(x){
  layer = x["layer"]
  if(!x["mediatorID"] %in% rownames(get(layer))){
    return(c(NA, NA))
  }
  dat = data.frame(mediator = ePheno[x["targetID"],], 
                   marker = marker_e, 
                   outcome  = get(layer)[x["mediatorID"],])
  dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
  if(nrow(dat) < 3){
    return(c(NA,NA))
  }
  fitM <- lm(mediator ~ marker, data=dat)
  fitComb <- lm(outcome ~ marker + mediator, data=dat) 
  fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                         treat="marker", mediator="mediator"))
  res = unlist(fitMedBoot[c("d0", "d0.p")])
  return(res)
  # fitMedBoot["d0.p"]
}))
osmResRev = cbind(combOsm, osmResRev)
colnames(osmResRev)[6:7] = c("ACME", "ACME_pVal")
osmResRev$ACME_pVal[osmResRev$ACME_pVal == 0] = 2e-16
osmResRev$layer = factor(osmResRev$layer,levels=names(colLayer))
osmResRev$neg_log10_pVal = -log10(osmResRev$ACME_pVal)
osmResRev$significant = osmResRev$ACME_pVal < 0.05
save(osmResRev, file="data/mediations_osmoticresponse_reversed.RData")
load("data/mediations_osmoticresponse_reversed.RData")
g1Osm = ggplot(osmRes, aes(x=mediator, y=neg_log10_pVal, fill=layer)) +
  geom_boxplot() +
  scale_fill_manual(values=colLayer)+
  geom_jitter()+
  ggtitle("true mediations")
g2Osm = ggplot(osmResRev, aes(x=mediator, y=neg_log10_pVal, fill=layer)) +
  geom_boxplot() +
  scale_fill_manual(values=colLayer)+
  geom_jitter() +
  ggtitle("switched mediations")
g3Osm = ggplot(osmRes, aes(x=layer, fill = significant))+
  geom_bar(stat="count") +
  facet_wrap(~ mediator)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_discrete(labels=c("e", "pt", "p", "phRes", "ph"))
g4Osm = ggplot(osmResRev, aes(x=layer, fill = significant))+
  geom_bar(stat="count") +
  facet_wrap(~ mediator)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_discrete(labels=c("e", "pt", "p", "phRes", "ph"))
dOsm = data.frame(true_ACME = osmRes$ACME, rev_ACME = osmResRev$ACME, 
               true_pVal = osmRes$ACME_pVal, rev_pVal = osmResRev$ACME_pVal)
g5Osm = ggplot(dOsm,aes(x=true_ACME, y=rev_ACME))+
  geom_point()
g6Osm = ggplot(dOsm,aes(x=true_pVal, y=rev_pVal))+
  geom_point()
gcombOsm = gridExtra::grid.arrange(g1Osm, g2Osm, g3Osm, 
                                   g4Osm,g5Osm,g6Osm)
ggsave(gcombOsm, filename= "graph/osmRes.png", width = 8, height=8)
#####

# Ste20 hotspot: mediation only for phResQTL #####
hs = "chrVIII:1"
hsInfo = unlist(hsOvMat[hsOvMat$name == hs,])
leadMarker_e = leader[which(hsOvMat$name == hs),"e"]
marker_e = strainGeno[,leadMarker_e]
# first, targets of the mating response
outcomes = setNames(targetsByHS[[which(hsOvMat$name==hs)]]$e,
                    id2name[targetsByHS[[which(hsOvMat$name==hs)]]$e])
candidates = c("GIP4", "FUS1", "KAR4", "LRG1", "RGA2", "MFA1", "APA2", 
               "RGD2", "CRH1", "GPA1", "PRM5", "EXO70", "FAR1",
               "STE18", "HYM1", "SST2", "PRP39", "KAR5", "FUS2", 
               "PRM1", "AGA1", "VHS3", "STE13", "DIG1", "PRM4",
               "RGC1", "TEC1", "PST1", "AFR1", "CDC20", "SPO11",
               "RIM101", "HAL1", "FUS3", "ERG11", "SNL1", "SKT5",
               "AKL1", "REG1", "GIC2", "DIG2", "BOI2", "AVT1", 
               "MSN4", "MNR2", "SRP40", "SSK1", "RCK2", "NUP2",
               "BDF1", "RSC9", "MLF3", "NBA1", "RTS1", "MSC3",
               "STE4", "STE18", "CDC42", "CDC24", "STE50", "STE5",
               "STE7", "STE11", "STE12", "STE20")
candidates = setNames(name2id[candidates], candidates)
mediators = candidates[candidates %in% rownames(phResPhenoAverage)]
mediators = mediators[!apply(phResPhenoAverage[mediators,],1,function(x){all(is.na(x))})]
outcomes = outcomes[!outcomes %in% mediators]
comb = expand.grid(outcome = names(outcomes),
                   mediator = names(mediators), 
                   stringsAsFactors=F)
comb$outcomeID = outcomes[comb$outcome]
comb$mediatorID = mediators[comb$mediator]
medphRes = t(apply(comb,1, function(x){
  dat = data.frame(outcome = ePheno[x["outcomeID"],], 
                   marker = marker_e, 
                   mediator = phResPhenoAverage[x["mediatorID"],])
  dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
  fitM <- lm(mediator ~ marker, data=dat)
  fitComb <- lm(outcome ~ marker + mediator, data=dat) 
  fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                         treat="marker", mediator="mediator"))
  res = unlist(fitMedBoot[c("d0", "d0.p")])
  # switch mediator and outcome:
  dat = data.frame(mediator = ePheno[x["outcomeID"],], 
                   marker = marker_e, 
                   outcome = phResPhenoAverage[x["mediatorID"],])
  dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
  fitM <- lm(mediator ~ marker, data=dat)
  fitComb <- lm(outcome ~ marker + mediator, data=dat) 
  fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                         treat="marker", mediator="mediator"))
  rev = unlist(fitMedBoot[c("d0", "d0.p")])
  return(c(res,rev))
}))
medphRes = cbind(comb, medphRes)
colnames(medphRes)[5:8] = c("ACME", "ACME_pVal", "ACME_rev", "ACME_pVal_rev")
medphRes$ACME_pVal[medphRes$ACME_pVal == 0] = 2e-16
medphRes$neg_log10_pVal = -log10(medphRes$ACME_pVal)
medphRes$significant = medphRes$ACME_pVal < 0.05
save(medphRes, file="data/mediations_Ste20_phRes.RData")
ggplot(medphRes, aes(x=mediator, fill = significant))+
  geom_bar(stat="count") +
  theme_classic() 
ggsave(filename="graph/mediations_Ste20_phRes.png")
plot(medphRes$ACME, medphRes$ACME_rev, 
     col = (medphRes$ACME_pVal<0.05)+1, 
     pch = (medphRes$ACME_pVal_rev < 0.05)+1)
png("graph/mediations_Ste20_phRes_truevsReversed.png")
plot(medphRes$ACME_pVal, medphRes$ACME_pVal_rev, xlab = "mediation p-value",
     ylab = "reversed mediation p-value")
dev.off()
######


# mediation for Psk2-hotspot #####
hs = "chrXII:2"
hsInfo = unlist(hsOvMat[hsOvMat$name == hs,])
leadMarker = leader[which(hsOvMat$name == hs),"ph"]
marker = strainGeno[,leadMarker]

phRestargets = targetsByHS[[which(hsOvMat$name == hs)]]$phr
phResoutcomes = phResPhenoAverage[phRestargets,]
Psk2 = "YOL045W"
lay2Pheno = c("e" = "ePheno", "pt" = "ptPheno", "p" = "pPheno",
              "ph" = "phPhenoAverage", "phr"= "phResPhenoAverage")

medPsk2 = lapply(names(lay2Pheno), function(layer){
  mediator = get(lay2Pheno[layer])[Psk2,]
  if(all(is.na(mediator))){
    return(NA)
  }
  med = t(apply(phResoutcomes,1, function(x){
    dat = data.frame(outcome = x, 
                     marker = marker, 
                     mediator = mediator)
    dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
    if(nrow(dat) < 3){
      return(c(NA,NA))
    }
    fitM <- lm(mediator ~ marker, data=dat)
    fitComb <- lm(outcome ~ marker + mediator, data=dat) 
    fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                           treat="marker", mediator="mediator"))
    res = unlist(fitMedBoot[c("d0", "d0.p")])
    return(res)
  }))
  rev = t(apply(phResoutcomes,1, function(x){
    dat = data.frame(outcome = mediator, 
                     marker = marker, 
                     mediator = x)
    dat = dat[apply(dat,1,function(y){!any(is.na(y))}),]
    if(nrow(dat) < 3){
      return(c(NA,NA))
    }
    fitM <- lm(mediator ~ marker, data=dat)
    fitComb <- lm(outcome ~ marker + mediator, data=dat) 
    fitMedBoot <- suppressMessages(mediate(fitM, fitComb, boot=T, sims=99, 
                                           treat="marker", mediator="mediator"))
    res = unlist(fitMedBoot[c("d0", "d0.p")])
    return(res)
  }))
  colnames(med) = c("ACME", "ACME_pval")
  colnames(rev) = c("ACME_rev", "ACME_rev_pval")
  return(cbind(med, rev))
  
  # fitMedBoot["d0.p"]
})
names(medPsk2) = names(lay2Pheno)
plot(ACME_pval ~ ACME_rev_pval, medPsk2[[1]])
barplot(sapply(medPsk2, function(x){
  table(factor(x[,"ACME_pval"] < 0.05,levels=c(TRUE, FALSE)))
}), ylab = "count", legend.text=T)

######


# mediation analyses for the other hotspots ######
lay2Pheno = c("e" = "ePheno", "pt" = "ptPheno", "p" = "pPheno", 
              "ph" = "phPhenoAverage", "phr" = "phResPhenoAverage")
hotspotCandidates = read.xlsx("data/hotspot_table.xlsx", 
                              sheetIndex=2, stringsAsFactors = F)
load("data/otherQTLstudies/Albert2018Candidates.RData")
candidates = sapply(hsOvMat$name, function(hs){
  candCurated = hotspotCandidates[hotspotCandidates$name == hs,2]
  candCurated = setNames(name2id[candCurated], candCurated)
  candAlbert = names(Albert2018Candidates[[hs]][1:5])
  candAlbert = setNames(name2id[candAlbert], candAlbert)
  return(list(candCurated = candCurated, candAlbert = candAlbert))
}, simplify=F)

curatedMediations = sapply(names(candidates), function(hs){
  print(hs)
  candCurated = candidates[[hs]]$candCurated
  if(is.na(candCurated)){
    return(NULL)
  }
  targets = targetsByHS[[which(hsOvMat$name == hs)]]
  outcomes = do.call(rbind,sapply(names(targets), function(x){
    temp = get(lay2Pheno[x])[targets[[x]],,drop = F]
    if(nrow(temp)>0){rownames(temp) = paste(rownames(temp),x, sep="_")}
    return(temp)
  }))
  leadMarker = round(mean(leader[which(hsOvMat$name == hs),], na.rm = T))
  marker = strainGeno[,leadMarker]
  CuratedRes = sapply(c("ePheno", "pPheno", "phPhenoAverage"), function(lay){
    print(lay)
    if(!candCurated %in% rownames(get(lay))){
      return(NULL)
    }
    mediator = get(lay)[candCurated,]
    med = calcMediationFast(marker=marker,
                            mediator=mediator, 
                            outcomes=outcomes, 
                            nthread=8)
    return(med)
  }, simplify=F)
  return(CuratedRes)
}, simplify = F)
save(curatedMediations, file = "data/mediations_hotspots_curatedCandidates.RData")
png("graph/mediations_hotspots.png", width=1200, height=1000, pointsize=25)
par(mfrow = c(3,4))
sapply(names(curatedMediations), function(hs){
  print(hs)
  x = curatedMediations[[hs]]
  if(is.null(x) | all(sapply(x,is.null))){
    return(NULL)
  }
  pVals = lapply(x, function(y){
    sapply(y, function(z)z[1,4])
  })
  pVals = pVals[sapply(pVals, length)>1]
  if(length(pVals)<2){
    hist(unlist(pVals), 
         main = paste0(hs, "\nmediation through ", 
                       names(pVals), "\nof ",
                       id2name[candidates[[hs]]$candCurated]), 
         xlab = "p-Values")
  } else{
    pValsLog = lapply(pVals, function(y) {
      if(length(y) < 1) {return(NULL)}
      y[y ==0] = 2.220446e-16
      -log10(y)})
    boxplot(pValsLog, ylab = "-log10(p-Values)", las = 1,
            main = paste0(hs, "\nmediation through ",
                          id2name[candidates[[hs]]$candCurated]))
  }
})
dev.off()
# png("graph/mediations_hotspots.png", 
#     width= 1000, height=600, pointsize = 25)
# dev.off()
#####


# Psk1/Hap1 locus #####


#####