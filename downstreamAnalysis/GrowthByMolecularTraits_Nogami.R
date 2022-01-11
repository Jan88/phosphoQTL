# general preparation ######
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(randomForest)
library(ranger)
library(xlsx)
library(sinaplot)
library(mediation)
source("lib/general_function.R")
geno = read.table("data/genotype_for_mapping.tsv", header = T, as.is=T)
anno = read.table("Saccharomyces_cerevisiae/SGD_features.tab",
                  header = F,sep = "\t",quote="",as.is=T)
id2name = anno$V5[anno$V4 != ""]
names(id2name) = anno$V4[anno$V4 != ""]
id2name[id2name == ""] = names(id2name[id2name == ""])
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
calcMediation = function(marker, mediator, outcomes){
  res = lapply(1:nrow(outcomes), function(i){
    outcome = outcomes[i,]
    dat = data.frame(outcome, marker , mediator)
    dat = dat[!is.na(dat$outcome),]
    fitM <- lm(mediator ~ marker, data=dat)
    fitComb <- lm(outcome ~ marker + mediator, data=dat) 
    fitMedBoot <- mediate(fitM, fitComb, boot=T, sims=999, 
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
  return(res)
}
######

# load the mapping data #####
load("data/Nogami2007/Nogami_mappingData.RData")
gPheno = mappingData$phenotype
info = read.xlsx("data/Nogami2007/journal.pgen.0030031.st001_curated.XLS",
                 sheetIndex=1, startRow=5, header=T,stringsAsFactors=F)
info$ID = gsub(info$ID, pattern="-",replacement=".", fixed=T)
ID2Description = setNames(info$Description, info$ID)
strains = rownames(mappingData$genotype)
strains[strains == "RM11.1a"] = "RM11-1a"
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


# load QTL data #####
# load the relevant data
load("data/Nogami2007/QTL_list.RData")
gQTL = QTL_list
load("data/Nogami2007/Nogami_genotypes.RData")
load("data/eQTL_results160831.RData")
eQTL <- QtlList$FDR10
load("data/pQTL_results170815.RData")
ptQTL <- pQTL_results$ptQTL$QtlList$FDR10
pQTL <- pQTL_results$pQTL$QtlList$FDR10
phResQTL <- pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10
phQTL <- pQTL_results$phosphoLevelQTL$QtlList$FDR10
rm(pQTL_results, pv, qv,QtlList, QTL_list)
allMarkers = seq(ncol(genotype))
#####






# correlate growth traits with each molecular trait #####
# prepare data
layersSub = c("transcripts" = "ePhenoSub",
              "pt" = "ptPhenoSub",
              "protein" = "pPhenoSub", 
              "phosphoRes" = "phResPhenoAverage",
              "phospho" = "phPhenoAverage")

# compute all correlations
correlationsSub = sapply(names(layersSub), function(type){
   predictors = t(get(layersSub[type]))
   temp = abs(cor(predictors, t(gPheno), use = "pair"))
   return(temp)
}, simplify=F)
# get only maximum correlations for each layer
maxCorrelationsSub = sapply(names(layersSub), function(type){
   predictors = t(get(layersSub[type]))
   temp = abs(cor(predictors, t(gPheno), use = "pair"))
   maxs = apply(temp,2,max)
   return(maxs)
})
# get the corresponding gene names
maxCorrelationsGeneSub = sapply(names(layersSub), function(type){
   predictors = t(get(layersSub[type]))
   temp = abs(cor(predictors, t(gPheno), use = "pair"))
   maxs = apply(temp,2,function(x){
      id2name[names(x)[which.max(x)]]
   })
   return(maxs)
})
maxCorrelationsGeneSubTable = cbind("compound" = rownames(maxCorrelationsGeneSub), 
                                    maxCorrelationsGeneSub, deparse.level=0)
write.table(maxCorrelationsGeneSubTable, 
            file="data/Nogami2007/maxCorrelationsGeneSub.tsv",
            sep = "\t", quote=F, col.names=T, row.names=F)

png("graph/Nogami2007/MaximumCorrelateGrowthWithMolecular_CustomplotNogami.png",
    width = 1380, height = 1200,pointsize=32)
temp = maxCorrelationsSub
colnames(temp) = c("transcripts", "pt", "protein", "phRes", "phospho")
par(mfrow = c(4,4), mar = c(0,0,0,0), oma = c(5,10,0.1,0))
for(i in 1:4){
  for(j in 1:4){
    if(i < j){
      plot.new()
    }else{
      plot(NA,xlim = c(min(temp),max(temp)), 
           ylim = c(min(temp),max(temp)),ylab = "", xlab = "",
           mgp = c(2,0.9,0), las = 1,
           xaxt = ifelse(i==4,"s","n"), yaxt = ifelse(j==1,"s","n"))
      abline(0,1, col = "grey", lwd = 1.5)
      points(temp[,j],temp[,i+1], 
             pch = 19, cex = 0.8, col = rgb(0,0,0,0.3),)
      res = wilcox.test(temp[,j],temp[,i+1],paired=T,exact=T)
      txt = paste0("p=", format(res$p.value, digits=3))
      legend("bottomright", legend=txt,col=NA, bty="n", cex = 1.2)
    }
    if(i==4){
      mtext(side=1,text=colnames(temp)[j], line=3.7,
            col = colors[j], font=2)
    }
    if(j==1){
      mtext(side=2, text=colnames(temp)[i+1], line = 4.5,las = 1,
            col = colors[i+1], font=2) #,
    }
  }
}
mtext(side=2,text="maximum absolute R", outer=T, line = 2.5)
mtext(side=1,text="maximum absolute R", outer=T, line = 2.5)
par( fig=c(0.5,1,0.5,1), new=TRUE, mar=c(2,3.5,0,0) )
sinaplot(as.data.frame(maxCorrelationsSub), las = 1, col=colors[1:5], 
         labels="",
         ylab = "maximum absolute R", mgp = c(2.1,0.8,0), 
         cex.lab = 1.5)
axis(side=1,at=1:5,cex.axis = 1.4,mgp = c(2.1,0.8,0), 
     labels=c("e", "pt", "p", "phRes", "ph"))
boxplot(maxCorrelationsSub, add = T, col = NA, boxwex=0.2, 
        axes = F, outline = F,lwd = 1.5)
dev.off()




# get the outliers, i.e. the best growth trait-molecular trait combos
outliers = sapply(names(layersSub), function(type){
   ord = order(maxCorrelationsSub[,type], decreasing=T)[1:15]
   temp = maxCorrelationsGeneSub[ord, type]
   cbind(condition = names(maxCorrelationsGeneSub[ord, type]), 
         gene = unname(maxCorrelationsGeneSub[ord, type]), 
         absR = format(unname(maxCorrelationsSub[ord, type]), digits=3))
}, simplify=F)
for(x in names(outliers)){print(x); cat.table.redmine(outliers[[x]])}
# get the genes that are top correlators with most growth traits
topGenes = apply(maxCorrelationsGeneSub,2, function(x){
   temp = table(x)
   temp = temp[order(temp, decreasing = T)]
   data.frame(gene = names(temp[1:10]), times = as.integer(temp[1:10]))
})
cat.table.redmine(do.call(cbind, topGenes))


# here comes all the plotting #
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
#custom panel function for matrixplots, with diagonal line
panelfunction= function(x,y){
   points(x,y, pch = 19, cex = 0.8, col = rgb(0,0,0,0.3))
   abline(0,1, col = "red")
   abline(h=0,v=0)
}
# matrix plot of MAX correlations
png("graph/Nogami2007/MaximumCorrelateGrowthWithMolecular_Matrixplot_Nogami.png",
    width = 800, height = 800,pointsize=20)
par(oma = c(0,0,0,0), mar = c(1,1,0,0))
pairs(maxCorrelationsSub, panel = panelfunction, cex.labels=1.6,  gap=0.5,
      labels=c("transcripts", "pt", "protein", "phRes", "phospho"))
mtext("maximum absolute R", side = 1, line=-0.5)
mtext("maximum absolute R", side = 2, line=-0.5)
dev.off()
# boxplot of MAX correlations
png("graph/Nogami2007/MaximumCorrelateGrowthWithMolecular_Nogami.png")
boxplot(maxCorrelationsSub, ylab = "maximum abs. Pearson's R",
        col = colors, border = colors, las = 1)
boxplot(maxCorrelationsSub, col = rgb(0,0,0,0), add = T,
        outline = F, names = NA,  axes = F)
dev.off()
# sinaplot of MAX correlations
png("graph/Nogami2007/MaximumCorrelateGrowthWithMolecular_sinaplot_Nogami.png",
    width = 300,height=300, pointsize=20)
par(mar = c(3,4,0.1,0.1))
sinaplot(as.data.frame(maxCorrelationsSub), ylab = "maximum absolute R",
        col = colors, las = 1. mgp=c(2,1,0))
boxplot(maxCorrelationsSub, col = rgb(0,0,0,0), add = T,
        outline = F, names = NA,  axes = F, boxwex =0.2)
dev.off()
# violin plot of MAX correlations
temp = melt(maxCorrelationsSub)
ggplot(temp, aes(x = Var2, y = value, fill=Var2)) + 
   guides(fill=FALSE)+
   labs(y="maximum\nabsolute Pearson's R", x = "")+
   theme(axis.title = element_text(size=20),
         axis.text.x = element_text(size=20, face = "plain",  
                                    colour = "black", 
                                    angle=45, hjust = 1, vjust = 1),
         axis.text.y = element_text(size=15, colour = "black",
                                    margin=ggplot2::margin(l=10)),
         axis.ticks =  element_blank(),
         axis.ticks.y = element_line(size = 0.5),
         axis.ticks.length = unit(x = 4,units = "mm"),
         axis.line.y = element_line(size=0.5))+
   geom_violin(adjust = 1)+ 
   stat_summary(fun.y = "median", geom = "point")+
   geom_boxplot(width = 0.2)+
   theme(panel.background = element_rect(fill = NA),
         panel.grid.major.x = element_blank())+
   scale_fill_manual(values=unname(colors))
ggsave(filename = "graph/Nogami2007/MaximumCorrelateGrowthWithMolecular_violin_Nogami.png", 
       width=6, height=6)

# and a test:
temp = melt(maxCorrelationsSub)
pairwise.wilcox.test(x=temp$value, g=temp$Var2, paired=T)
# matrix plot of all correlations
png("graph/Nogami2007/CorrelateGrowthWithMolecular_Matrixplot_Nogami.png")
correlationsSubTemp = sapply(correlationsSub, as.vector)
pairs(correlationsSubTemp, panel = panelfunction) 
mtext("abs. Pearson R", side = 1, line=4)
mtext("abs. Pearson R", side = 2, line=2.5)
dev.off()
temp = melt(correlationsSub)
temp$L1 = factor(temp$L1, levels=unique(temp$L1))
ggplot(temp, aes(x = L1, y = value, fill=L1)) + 
   guides(fill=FALSE)+
   labs(y="absolute Pearson's R", x = "")+
   theme(axis.title = element_text(size=20),
         axis.text.x = element_text(size=20, face = "plain",  
                                    colour = "black", 
                                    angle=45, hjust = 1, vjust = 1),
         axis.text.y = element_text(size=15, colour = "black",
                                    margin=ggplot2::margin(l=10)),
         axis.ticks =  element_blank(),
         axis.ticks.y = element_line(size = 0.5),
         axis.ticks.length = unit(x = 4,units = "mm"),
         axis.line.y = element_line(size=0.5))+
   geom_violin(adjust = 1)+ 
   stat_summary(fun.y = "median", geom = "point")+
   geom_boxplot(width = 0.2)+
   theme(panel.background = element_rect(fill = NA),
         panel.grid.major.x = element_blank())+
   scale_fill_manual(values=unname(colors))
ggsave(filename = "graph/Nogami2007/CorrelateGrowthWithMolecular_Nogami.png")
# test difference between correlations
pairwise.wilcox.test(x=temp$value, g=temp$Var2, paired=T)
colSums(correlationsSub)
######


# extract which layer provides highest correlation for each gene-growthTrait pair #####
# compute all correlations
layersSub = c("transcripts" = "ePhenoSub",
              "pt" = "ptPhenoSub",
              "protein" = "pPhenoSub", 
              "phosphoRes" = "phResPhenoAverage",
              "phospho" = "phPhenoAverage")
correlationsSub = sapply(names(layersSub), function(type){
  predictors = t(get(layersSub[type]))
  temp = abs(cor(predictors, t(gPheno), use = "pair"))
  return(temp)
}, simplify=F)
topCorPerGrowthTraitExtended = sapply(rownames(gPheno), function(x){
  temp2 = sapply(correlationsSub, function(y){
    y[,x]
  })
  index = arrayInd(ind = which.max(temp2),.dim=dim(temp2))
  setNames(colnames(temp2)[index[,2]],rownames(temp2)[index[,1]])
}, USE.NAMES=F)
topCorCounts = table(topCorPerGrowthTraitExtended)
png("graph/Nogami2007/topCorPerGrowthTrait_Nogami.png")
barplot(topCorCounts[5:1], col = colors, las = 1, ylab = "n growth Traits")
dev.off()

# which genes are top correlators?
topGeneCounts = table(id2name[names(topCorPerGrowthTraitExtended)])
topGeneCounts[order(topGeneCounts,decreasing=T)[1:10]]
# which genes are top correlators for each layer?
topGeneCountsSplit = sapply(split(id2name[names(topCorPerGrowthTraitExtended)],
                                  topCorPerGrowthTraitExtended),
                            function(x){
                              temp = table(x)
                              temp[order(temp, decreasing=T)]
                            })

pdfAndPng("graph/Nogami2007/topCorPerGrowthTrait_topGenes_over4_Nogami", 
          width=6, height=6, expr=({
            par(mar = c(2,4,0.1,0))
            for(i in 1:5){
              temp = do.call(cbind,sapply(5:1, function(j){
                if(j==i){
                  return(topGeneCountsSplit[[i]])
                }else{return(NA)}
              }))
              n = if(i==1){c("transcripts", "pt", "protein", "phRes", "phospho")}else{rep(NA,5)}
              pos = barplot(temp, add = (i!=1), las = 1,yaxt = "n",
                            col = colors[5-i+1],xpd = NA,names.arg=n,
                            mgp = c(2,0.5,0), ylim = c(0,max(table(topCorPerGrowthTraitExtended))))
              toMark = topGeneCountsSplit[[i]][topGeneCountsSplit[[i]]>4]
              if(length(toMark)>0){
                bold = (names(toMark) %in% c("AVT1", "SPE4", "MKT1")) +1
                xPos = pos[6-i]
                yPos = c(0,cumsum(toMark))+c(toMark/2,-sum(toMark))
                yPos = yPos[-length(yPos)]
                text(x = xPos, y = yPos, names(toMark), font = bold)
              }
            }
            axis(2, las = 1)
            mtext(side=2, text="n traits", line = 2.5)https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4017318/
          }))
# png("graph/Nogami2007/topCorPerGrowthTrait_topGenes_over4_Nogami.png", 
#     height=800, width=800, pointsize=27)
# par(mar = c(2,4,0.1,0))
# for(i in 1:5){
#   temp = do.call(cbind,sapply(5:1, function(j){
#     if(j==i){
#       return(topGeneCountsSplit[[i]])
#     }else{return(NA)}
#   }))
#   n = if(i==1){c("transcripts", "pt", "protein", "phRes", "phospho")}else{rep(NA,5)}
#   pos = barplot(temp, add = (i!=1), las = 1,yaxt = "n",
#                 col = colors[5-i+1],xpd = NA,names.arg=n,
#                 mgp = c(2,0.5,0), ylim = c(0,max(table(topCorPerGrowthTraitExtended))))
#   toMark = topGeneCountsSplit[[i]][topGeneCountsSplit[[i]]>4]
#   if(length(toMark)>0){
#     bold = (names(toMark) %in% c("AVT1", "SPE4", "MKT1")) +1
#     xPos = pos[6-i]
#     yPos = c(0,cumsum(toMark))+c(toMark/2,-sum(toMark))
#     yPos = yPos[-length(yPos)]
#     text(x = xPos, y = yPos, names(toMark), font = bold)
#   }
# }
# axis(2, las = 1)
# mtext(side=2, text="n traits", line = 2.5)
# dev.off()

topCorPerGrowthTraitGene = sapply(rownames(gPheno), function(x){
  temp2 = sapply(correlationsSub, function(y){
    y[,x]
  })
  index = arrayInd(ind = which.max(temp2),.dim=dim(temp2))
  unname(id2name[rownames(temp2)[index[,1]]])
})
# names(topCorPerGrowthTraitGene) = ID2Description[names(topCorPerGrowthTraitGene)]

candidateVals = cbind("AVT1" = phResPhenoAverage[names(id2name[id2name == "AVT1"]),],
                      "SPE4" = phResPhenoAverage[names(id2name[id2name == "SPE4"]),],
                      "PSK1" = phPhenoAverage[names(id2name[id2name == "PSK1"]),], 
                      "YLR257W" = phPhenoAverage["YLR257W",])
candidateCorsAll = sapply(colnames(candidateVals), function(x){
  targets = names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == x]
  pheno = candidateVals[,x]
  cors = cor(pheno, t(gPheno), use = "pair")[1,]
})
candidateCorsTargets = sapply(colnames(candidateVals), function(x){
  targets = names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == x]
  pheno = candidateVals[,x]
  cors = cor(pheno, t(gPheno[targets,]), use = "pair")[1,]
})

# boxplot(candidateCorsAll, las = 1, ylab = "Pearson's R")
temp = do.call(rbind,lapply(colnames(candidateCorsAll),function(x){
  vals = candidateCorsAll[,x]
  trts = rownames(candidateCorsAll)
  targts = names(candidateCorsTargets[[x]])
  data.frame(gene = x,cor = vals, target = as.integer(trts %in% targts))
}))
png("graph/Nogami2007/candidatesCorrelationNogami.png")
temp2 = sinaplot(cor ~ gene, data = temp, las = 1,
                 ylab = "Pearson's R", xlab = "", pch = NA, 
                 labels = paste0(unique(temp$gene), "\n",
                                 c("phosphoRes", "phosphoRes", 
                                   "phospho", "phospho")), 
                 mgp= c(3,2,0), yaxt = "n")
axis(2, las = 1)
points(temp2$scaled,temp2$y,col = temp$target+1)
points(temp2$scaled[temp$target == 1],temp2$y[temp$target==1],
       col = "red", ticks = F)
abline(h =0)
dev.off()

# AVT1: Vacuolar transporter; imports large neutral amino acids into the vacuole 
unname(ID2Description[names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == "AVT1"]])
# SPE4: Spermine synthase,  also involved in biosynthesis of pantothenic acid
unname(ID2Description[names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == "SPE4"]])
# PSK1: no comment needed. 
unname(ID2Description[names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == "PSK1"]])
# YLR257W: unknown function; protein abundance increases in response to DNA replication stress
unname(ID2Description[names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == "YLR257W"]])

# closer look at YLR257W
YLR257Wpeptide = unlist(prot2phosphoHash["YLR257W"])
phQTLtargets = sapply(phQTL,function(x){x$target})
phQTLtargets = rownames(phPheno)[phQTLtargets]
YLR257Wqtl = lapply(phQTL[phQTLtargets %in% YLR257Wpeptide], function(x){x$predictors[1,]})
YLR257Wpeaks = sapply(phQTL[phQTLtargets %in% YLR257Wpeptide], function(x){x$mostSignificantPredictor})
YLR257WpeakPos = markerPos[YLR257Wpeaks,] 
YLR257WpeakPos[order(YLR257WpeakPos[,1], YLR257WpeakPos[,2]),]
# there are multiple QTL, including multiple loci in the AMN1-hostpot on chrII, 
# the HAP1-locus on chrXII, the MKT1-locus and the IRA2 locus

# is level of YLR257 correlated with HAP1, IRA2, MKT1, or AMN1?
otherGenes = sapply(c("HAP1", "IRA2", "MKT1", "AMN1"), function(x){
  id = names(id2name[id2name == x])
  sapply(c("ePheno", "ptPheno", "pPheno", 
           "phResPhenoAverage", "phPhenoAverage"), function(y){
             if(id %in% rownames(get(y))){
               return(get(y)[id,])
             } else{ return(rep(NA, ncol(get(y))))}
           })
}, simplify=F)
YLR257WphosphoCorOtherGenes = sapply(names(otherGenes), function(x){
  cor(otherGenes[[x]], phPhenoAverage["YLR257W",], use = "pair")[,1]
})
YLR257WtranscriptCorOtherGenes = sapply(names(otherGenes), function(x){
  cor(otherGenes[[x]], ePheno["YLR257W",], use = "pair")[,1]
}
png("graph/Nogami2007/corr_YLR257W_otherGene.png")
par(mfrow = c(2,1), mar = c(4,5,0.1,1))
barplot(YLR257WphosphoCorOtherGenes, beside = T, col = colors[1:5], 
        las= 1, ylab = "correlation with\nYLR257W phopho", mgp=c(2.5,1,0))
abline(h = 0)
barplot(YLR257WtranscriptCorOtherGenes, beside = T, col = colors[1:5], 
        las= 1, ylab = "correlation with\nYLR257W transcript", mgp=c(2.5,1,0))
abline(h = 0)
dev.off()


# closer look at AVT1:
# AVT1 lies at 	chrX:436802-438610 
AVT1 = "YJR001W"
AVT1peptide = unlist(prot2phosphoHash[AVT1])
# peptides are at protein position 180-191 and 2-12, both known
phResQTLtargets = sapply(phResQTL,function(x){x$target})
phResQTLtargets = rownames(phResPheno)[phResQTLtargets]
AVT1qtl = lapply(phResQTL[phResQTLtargets %in% AVT1peptide], function(x){x$predictors[1,]})
AVT1peaks = sapply(phResQTL[phResQTLtargets %in% AVT1peptide], function(x){x$mostSignificantPredictor})
markerPos[AVT1peaks,] # chrVIII 97536-98643, so within GPA1 hotspot
# two synonymous SNPs in AVT1 in our cross:
rawGeno = read.table("data/strain_genotype.tsv", header = T)
rawGeno[(rawGeno$CHROM == "chrX" & 
        rawGeno$POS >= 436802 & rawGeno$POS <= 438610),1:4]
# shape traits affected by the same QTL?
AVT1gTraits = names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == "AVT1"]
affectedByAVT1qtl = t(sapply(gQTL, function(x){
  sameT = rownames(gPheno)[x$target] %in% AVT1gTraits
  sameQ = any((min(unlist(AVT1qtl))-10):(max(unlist(AVT1qtl))+10) %in%
                (min(x$predictors)-10):(max(x$predictors)+10))
  return(c(sameT=sameT, sameQ = sameQ))
}))
affectedByAVT1qtl[affectedByAVT1qtl[,1],] # all relevant growth traits are also affected by QTL
AVT1mediation = calcMediation(marker = genotype[,AVT1peaks],
                              mediator=phResPhenoAverage[AVT1,],
                              outcomes=gPheno[AVT1gTraits,])
save(AVT1mediation, file = "data/Nogami2007/AVT1mediation.RData")
AVT1mediationPvals = sapply(AVT1mediation, function(x){x[1,4]})
AVT1mediationPvals[AVT1mediationPvals == 0] = 0.000000001



# closer look at Spe4:
# Spe4 is located at chrXII:433725-432823
SPE4 = "YLR146C"
SPE4peptide = unlist(prot2phosphoHash[SPE4])
# peptide at position 2-11 (known)
SPE4qtl = lapply(phResQTL[phResQTLtargets %in% SPE4peptide], function(x){x$predictors[1,]})
SPE4peaks = sapply(phResQTL[phResQTLtargets %in% SPE4peptide], function(x){x$mostSignificantPredictor})
markerPos[SPE4peaks,] # chrVIII 380673-388664, not within any hotspot
rawGeno[(rawGeno$CHROM == "chrXII" & 
           rawGeno$POS >= 433725-500 & rawGeno$POS <= 432823+500),1:4]
# no known polymorphisms
# mediation:
SPE4gTraits = names(topCorPerGrowthTraitGene)[topCorPerGrowthTraitGene == "SPE4"]
affectedBySPE4qtl = t(sapply(gQTL, function(x){
  sameT = rownames(gPheno)[x$target] %in% SPE4gTraits
  sameQ = any((min(unlist(SPE4qtl))-20):(max(unlist(SPE4qtl))+20) %in%
                (min(x$predictors)-20):(max(x$predictors)+20))
  return(c(sameT=sameT, sameQ = sameQ))
}))
affectedBySPE4qtl[affectedBySPE4qtl[,1],] # all relevant growth traits are also affected by QTL
affectedBySPE4qtlTest = apply(gPheno[SPE4gTraits,],1,function(x){
  mod = summary(lm(x ~ genotype[,SPE4peaks]))$coefficients[2,4]
})
png("graph/Nogami2007/affectedBySPE4qtlTest.png")
hist(affectedBySPE4qtlTest, xlab = "p-value", xlim = c(0,1), 
     main="Effect of SPE4-phResQTL on actin traits")
abline(v = 0.05, lty = 2)
dev.off()


SPE4mediation = calcMediation(marker = genotype[,SPE4peaks],
                              mediator=phResPheno[SPE4peptide,],
                              outcomes=gPheno[SPE4gTraits,]                              )
save(SPE4mediation, file = "data/Nogami2007/SPE4mediation.RData")
SPE4mediationPvals = sapply(SPE4mediation, function(x){x[1,4]})
SPE4mediationPvals[SPE4mediationPvals == 0] = 0.000000001
SPE4mediationEstimates = sapply(SPE4mediation, function(x){x[,1]})

png("graph/Nogami2007/mediation_AVT1_SPE4_pvals.png",
    width=800, height=800, pointsize=25)
par(mfrow = c(2,1), mar = c(4,3,1,0))
hist(p.adjust(AVT1mediationPvals, method="fdr"), xlim = c(0,1), 
     las = 1, xlab = "p-value", main = "AVT1", mgp = c(2,1,0))
abline(v = 0.05, lty = 2)
hist(p.adjust(SPE4mediationPvals, method="fdr"), xlim = c(0,1), 
     las = 1, xlab = "p-value", main = "SPE4", mgp = c(2,1,0))
abline(v = 0.05, lty = 2)
dev.off()

png("graph/Nogami2007/topCorPerGrowthTrait_Nogami.png")
barplot(topCorCounts[5:1], col = colors, las = 1, ylab = "n growth Traits")
dev.off()

# compare the actual correlation values
topCorPerGrowthTraitValues = sapply(rownames(gPheno), function(x){
  temp2 = sapply(correlationsSub, function(y){
    y[,x]
  })
  index = arrayInd(ind = which.max(temp2),.dim=dim(temp2))
  setNames(temp2[index],colnames(temp2)[index[,2]])
}, simplify = T, USE.NAMES=F)

res = pairwise.wilcox.test(topCorPerGrowthTraitValues,
                           names(topCorPerGrowthTraitValues), p.adjust.method="none")
cat.table.redmine(cbind(rownames(res$p.value), 
                        format(res$p.value, digits = 2)))


library(sinaplot)
png("graph/Perlstein/topCorPerGrowthTrait_corValues.png",
    height=800, width=800, pointsize=27)
par(mar = c(2,4,0.1,0.1))
boxplot(split(topCorPerGrowthTraitValues, 
              names(topCorPerGrowthTraitValues))[5:1], 
        col = colors, las = 1, outline=F, 
        ylim = c(min(topCorPerGrowthTraitValues),
                 max(topCorPerGrowthTraitValues)),
        ylab = "maximum\nabsolute Pearson's R")
temp = sinaplot(topCorPerGrowthTraitValues, xaxt = "n",yaxt = "n",
                groups = factor(names(topCorPerGrowthTraitValues), 
                                levels=c("transcripts", "pt", "protein", 
                                         "phosphoRes", "phospho")), 
                pch = 21, col = colors[1:5], plot = F)
points(temp[,"scaled"], temp[,"y"], bg = temp[,"col"],  pch = 21)
dev.off()




# Vise versa: for each gene, on which layer 
# does it correlate most with growth traits?
correlationsPerGene = sapply(genesAvail, function(x){
   sapply(correlationsSub, function(y){y[x,]})
}, simplify=F)
correlationsPerGeneMax = sapply(genesAvail, function(x){
   names(which.max(apply(abs(correlationsPerGene[[x]]), 2, median)))
})
png("graph/Nogami2007/topCorPerGene_Nogami.png")
barplot(table(correlationsPerGeneMax)[5:1], col = colors, las = 1, ylab = "n growth Traits")
dev.off()
#####



# compare with the QTL for each layer #####
gQTL_peaks = unlist(sapply(gQTL, function(x){
   x$mostSignificantPredictor}))
is_gQTL = allMarkers %in% gQTL_peaks
layers = c("eQTL", "ptQTL", "pQTL",  "phResQTL", "phQTL")
pVals = t(sapply(layers, function(type){
   qtl = get(type)
   peaks = sapply(qtl, function(x){x$mostSignificantPredictor})
   is_QTL = allMarkers %in% peaks
   tab = matrix(c(sum( is_QTL &  is_gQTL), 
                  sum( is_QTL & !is_gQTL), 
                  sum(!is_QTL &  is_gQTL),
                  sum(!is_QTL & !is_gQTL)), nrow = 2)
   test = fisher.test(tab, alternative = "greater")
   return(c(pval = test$p.value, test$estimate))
}))
cat.table.redmine(cbind("layer" = rownames(pVals),
                        "p-value" = format(pVals[,1], digits = 3), 
                        "Odds ratio" = format(pVals[,2], digits = 3)))
library(eulerr)
PeakMarkersPerLayer = sapply(c("eQTL", "ptQTL", "pQTL",  
                               "phResQTL", "phQTL", "gQTL"),function(x){
                                  unique(unlist(sapply(get(x),function(y){
                                     y$mostSignificantPredictor})))
                               })
sapply(c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), function(x){
   png(paste0("graph/Nogami2007/Venn_gQTL_",x,"_Nogami.png"))
   fit = euler(PeakMarkersPerLayer[c(x, "gQTL")])
   temp = plot(fit, quantities =T,
               fills = list(fill = c(colors[x], "grey"), alpha = 0.9),
               edges = F, labels=list(cex=2),
               main=paste("p-value =", format(pVals[x],digits=2)))
   plot(temp)
   dev.off()
})

#####


# map growth traits based on molecular QTL #####
layers = c("gQTL", "eQTL", "ptQTL", "pQTL",  "phResQTL", "phQTL")
predictions = sapply(layers, function(type){
   print(type)
   qtl = get(type)
   peaks = sort(unique(unlist(sapply(qtl, function(x){
      x$mostSignificantPredictor}))))
   predictors = genotype[,peaks]
   predictors = predictors[,!apply(predictors,2,function(x){any(is.na(x))})]
   res = apply(gPheno,1,function(trait){
      dat = data.frame(trait, predictors)
      rf = ranger(trait ~ ., data=dat)
      return(c(rf$r.squared))
   })
})

temp = melt(predictions)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
ggplot(temp, aes(x = factor(Var2), y = value, fill=Var2)) + 
   guides(fill=FALSE)+
   labs(y="RF R-squared", x = "")+
   theme(axis.title = element_text(size=20),
         axis.text.x = element_text(size=20,face = "plain",  colour = "black", 
                                    angle=45, hjust = 1, vjust = 1),
         axis.text.y = element_text(size=15, colour = "black",
                                    margin=ggplot2::margin(l=10)),
         axis.ticks =  element_blank(),
         axis.ticks.y = element_line(size = 0.5),
         axis.ticks.length = unit(x = 4,units = "mm"),
         axis.line.y = element_line(size=0.5))+
   geom_violin(adjust = 1)+ 
   stat_summary(fun.y = "median", geom = "point")+
   geom_boxplot(width = 0.2)+
   theme(panel.background = element_rect(fill = NA),
         panel.grid.major.x = element_blank())+
   scale_fill_manual(values=unname(c("grey", colors[])))
ggsave(filename = "graph/Nogami2007/GrowthByMolecularQTL_Nogami.png", height=6,width=6)


png("graph/Nogami2007/GrowthByMolecularQTL_Matrixplot_Nogami.png",
    width = 800, height = 800,pointsize=20)
panelfunction= function(x,y){
   points(x,y, pch = 19, cex = 0.8, col = rgb(0,0,0,0.3))
   abline(0,1, col = "red")
   abline(h=0,v=0)
}
par(oma = c(0,0,0,0), mar = c(1,1,0,0))
pairs(predictions, panel = panelfunction, cex.labels=1.6,  gap=0.5)
mtext("RF R-squared", side = 1, line=-0.5)
mtext("RF R-squared", side = 2, line=-0.5)
dev.off()




t.test(predictions[,"eQTL"], predictions[,"pQTL"], paired=T)$p.value
#####

# ira2 ########
# chrXV 160001-200001
ira2Start = 160001
ira2End = 200001
gQTLpos = t(sapply(gQTL, function(x){
  m = x$predictors
  mChr = as.character(unique(markerPos[m,1]))
  mStart = min(markerPos[m,2])
  mEnd = max(markerPos[m,3])
  c(mChr, mStart, mEnd)
}))
  
table(gQTLpos[,1] == "chrXV" & 
  ((as.integer(gQTLpos[,2]) >= ira2Start & as.integer(gQTLpos[,2]) <= ira2End) |
  (as.integer(gQTLpos[,3]) >= ira2Start & as.integer(gQTLpos[,3]) <= ira2End) |
  (as.integer(gQTLpos[,2]) <= ira2Start & as.integer(gQTLpos[,3]) >= ira2End)))
#######


# hotspots #####
load("/cellnet/phosphoQTL/data/hotspotLocations.RData") #"markerLocations" "leader"   
load("/cellnet/phosphoQTL/data/hsInfo.RData") #"hsOvMat" "targetsByHS" "GOanalysis" 
gQTLpos = do.call(rbind,lapply(gQTL, function(x){
  m = x$predictors
  mChr = as.character(unique(markerPos[m,1]))
  mStart = min(markerPos[m,2])
  mEnd = max(markerPos[m,3])
  data.frame(mChr, mStart, mEnd, stringsAsFactors=F)
}))
hs_gQTL = apply(hsOvMat, 1, function(x){
  sum(apply(gQTLpos, 1,function(y){
    y["mChr"] == x["chr"] & 
      as.numeric(y["mStart"]) >= as.numeric(x["startPos"]) &
      as.numeric(y["mEnd"]) <= as.numeric(x["endPos"])
  }))
})
names(hs_gQTL) = hsOvMat[,"name"]
t(t(hs_gQTL))
######
