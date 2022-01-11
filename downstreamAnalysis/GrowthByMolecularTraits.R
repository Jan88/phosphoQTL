# prepare data  ######
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(randomForest)
library(ranger)
source("lib/general_function.R")
geno = read.table("data/genotype_for_mapping.tsv", header = T, as.is=T)
anno = read.table("Saccharomyces_cerevisiae/SGD_features.tab",
                  header = F,sep = "\t",quote="",as.is=T)
id2name = anno$V5[anno$V4 != ""]
names(id2name) = anno$V4[anno$V4 != ""]
id2name[id2name == ""] = names(id2name[id2name == ""])
library(mediation)
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
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
######


# load the Perlstein data #####
load("data/PerlsteinData/Perlstein_mappingData.RData")
conditions = rownames(mappingData$phenotype)
gPheno = mappingData$phenotype
temp = sapply(conditions, function(x){
  x = strsplit(x,split=" ")[[1]]
  if(length(x) == 5){
    return(c(x[1], paste(x[2:5], collapse = "_")))
  } else if(length(x) == 6){
    return(c(paste(x[1:2], collapse = " "), paste(x[3:6], collapse = "_")))
  } else
    return(x)
})
compounds = temp[1,]
exposures = temp[2,]
rm(temp)
load("data/PerlsteinData/Perlstein_QTL_list_merged.RData")
load("data/PerlsteinData/Perlstein_QTL_list_merged_HR.RData")
strains = rownames(mappingData$genotype)
#####

# load mapping data #####
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
load("data/PerlsteinData/Perlstein_QTL_list_merged.RData")
load("data/PerlsteinData/Perlstein_genotypes.RData")
load("data/eQTL_results160831.RData")
eQTL <- QtlList$FDR10
load("data/pQTL_results170815.RData")
ptQTL <- pQTL_results$ptQTL$QtlList$FDR10
pQTL <- pQTL_results$pQTL$QtlList$FDR10
phResQTL <- pQTL_results$phosphoProtResidualsQTL$QtlList$FDR10
phQTL <- pQTL_results$phosphoLevelQTL$QtlList$FDR10
rm(pQTL_results, pv, qv,QtlList)
allMarkers = seq(ncol(genotype))
gQTL = lapply(QTL_list_merged, function(x){
  x = data.frame(x,stringsAsFactors=F)
  x$start = as.integer(x$start)
  x$end = as.integer(x$end)
  x$mostSignificantPredictor = as.integer(x$mostSignificantPredictor)
  x$minP = as.numeric(x$minP)
  return(x)
})
#####





# map growth traits based on molecular traits #####
layers = c("transcripts" = "ePheno",
           "pt" = "ptPheno",
           "protein" = "pPheno", 
           "phosphoRes" = "phResPheno",
           "phospho" = "phPheno")
load("data/PerlsteinData/Perlstein_mappingData.RData")
pheno = mappingData$phenotype

predictions = sapply(names(layers), function(type){
  print(type)
  predictors = t(get(layers[type]))
  predictors = apply(predictors,2,function(x){
    missing = which(is.na(x))
    if(length(missing)>0){
      x[missing] = mean(x,na.rm=T)
    }
    return(x)
  })
  res = apply(pheno[1:100,],1,function(trait){
    samp = sample(length(trait), size=10)
    trainPredictors = predictors[-samp,]
    testPredictors = predictors[samp,]
    trainTrait = trait[-samp]
    testTrait = trait[samp]
    rf = randomForest(x=trainPredictors, y=trainTrait)
    pred = predict(object=rf, newdata=testPredictors)
    return(cor(pred,testTrait))
  })
})
save(predictions, 
     file ="data/PerlsteinData/predictions_GrowthByMolecularTraits.RData")
layersSub = c("transcripts" = "ePhenoSub",
              "pt" = "ptPhenoSub",
              "protein" = "pPhenoSub", 
              "phosphoRes" = "phResPhenoAverage",
              "phospho" = "phPhenoAverage")
predictionsSub = sapply(names(layersSub), function(type){
  print(type)
  predictors = t(get(layersSub[type]))
  predictors = apply(predictors,2,function(x){
    missing = which(is.na(x))
    if(length(missing)>0){
      x[missing] = mean(x,na.rm=T)
    }
    return(x)
  })
  res = apply(pheno[1:100,],1,function(trait){
    samp = sample(length(trait), size=10)
    trainPredictors = predictors[-samp,]
    testPredictors = predictors[samp,]
    trainTrait = trait[-samp]
    testTrait = trait[samp]
    rf = randomForest(x=trainPredictors, y=trainTrait)
    pred = predict(object=rf, newdata=testPredictors)
    return(cor(pred,testTrait))
  })
}) 
save(predictionsSub, 
     file ="data/PerlsteinData/predictionsSub_GrowthByMolecularTraits.RData")
apply(predictionsSub, 2, function(x){
  wilcox.test(x,mu=0)$p.value
})


# colors = c("gQTL" = "grey", colors)
temp = melt(predictions)
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
  scale_fill_manual(values=unname(colors))
ggsave(filename = "graph/Perlstein/GrowthByMolecularTraits.png")

temp = melt(predictionsSub)
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
  scale_fill_manual(values=unname(colors))
ggsave(filename = "graph/Perlstein/GrowthByMolecularTraits_subset.png")
#####


# correlate growth traits with each molecular trait #####
# prepare data
layersSub = c("transcripts" = "ePhenoSub",
              "pt" = "ptPhenoSub",
              "protein" = "pPhenoSub", 
              "phosphoRes" = "phResPhenoAverage",
              "phospho" = "phPhenoAverage")
#custom panel function for matrixplots, with diagonal line
panelfunction= function(x,y){
  points(x,y, pch = 19, cex = 0.8, col = rgb(0,0,0,0.3))
  abline(0,1, col = "red")
  abline(h=0,v=0)
}

# compute all correlations
correlationsSub = sapply(names(layersSub), function(type){
  predictors = t(get(layersSub[type]))
  temp = abs(cor(predictors, t(gPheno), use = "pair"))
  return(temp)
}, simplify = F)
# plots of all correlations
# matrix plot of all correlations
{png("graph/Perlstein/CorrelateGrowthWithMolecular_Matrixplot.png")
correlationsSubTemp = sapply(correlationsSub, as.vector)
pairs(correlationsSubTemp, panel = panelfunction) 
mtext("abs. Pearson R", side = 1, line=4)
mtext("abs. Pearson R", side = 2, line=2.5)
dev.off()
# violin plot of all correlations
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
ggsave(filename = "graph/Perlstein/CorrelateGrowthWithMolecular.png")}
# test difference between correlations
pairwise.wilcox.test(x=temp$value, g=temp$Var2, paired=T)
colSums(correlationsSub)

# extract which layer provides highest correlation for each gene-growthTrait pair
topCorPerGrowthTraitExtended = sapply(conditions, function(x){
  temp2 = sapply(correlationsSub, function(y){
    y[,x]
  })
  index = arrayInd(ind = which.max(temp2),.dim=dim(temp2))
  setNames(colnames(temp2)[index[,2]],rownames(temp2)[index[,1]])
}, USE.NAMES=F)
topCorCounts = table(topCorPerGrowthTraitExtended)
png("graph/Perlstein/topCorPerGrowthTrait.png")
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

# png("graph/Perlstein/topCorPerGrowthTrait_topGenes_over3.png", 
#     height=800, width=800, pointsize=27)
pdfAndPng("graph/Perlstein/topCorPerGrowthTrait_topGenes_over3", 
          width=6, height=6, expr=({
par(mar = c(2,4,0.2,0))
for(i in 1:5){
  temp = do.call(cbind,sapply(5:1, function(j){
    if(j==i){
      return(topGeneCountsSplit[[i]])
    }else{return(NA)}
  }))
  n = if(i==1){c("transcripts", "pt", "protein", "phRes", "phospho")}else{rep(NA,5)}
  pos = barplot(temp, add = (i!=1), las = 1,yaxt = "n",
                col = colors[5-i+1],xpd = NA,names.arg=n,
                mgp = c(2,0.5,0))
  toMark = topGeneCountsSplit[[i]][topGeneCountsSplit[[i]]>3]
  if(length(toMark)>0){
    bold = names(toMark) %in% c("MKT1", "RPS31")+1
    xPos = pos[6-i]
    yPos = c(0,cumsum(toMark))+c(toMark/2,-sum(toMark))
    yPos = yPos[-length(yPos)]
    text(x = xPos, y = yPos, names(toMark), font=bold)
  }
}
axis(2, las = 1)
mtext(side=2, text="n traits", line = 2.5)
}))
# dev.off()
png("graph/Perlstein/topCorPerGrowthTrait_topGenes_over2.png",
    height=800, width=800, pointsize=27)
par(mar = c(2,4,0.2,0))
for(i in 1:5){
  temp = do.call(cbind,sapply(5:1, function(j){
    if(j==i){
      return(topGeneCountsSplit[[i]])
    }else{return(NA)}
  }))
  n = if(i==1){c("transcripts", "pt", "protein", "phRes", "phospho")}else{rep(NA,5)}
  pos = barplot(temp, add = (i!=1), las = 1,yaxt = "n",
                col = colors[5-i+1],
                names.arg=n, mgp=c(2,0.5,0))
  toMark = topGeneCountsSplit[[i]][topGeneCountsSplit[[i]]>2]
  if(length(toMark)>0){
    xPos = pos[6-i]
    yPos = c(0,cumsum(toMark))+c(toMark/2,-sum(toMark))
    yPos = yPos[-length(yPos)]
    text(x = xPos, y = yPos, names(toMark))
  }
}
axis(2, las = 1)
mtext(side=2, text="n growth traits", line = 2.5)
dev.off()


topCorPerGrowthTraitValues = sapply(conditions, function(x){
  temp2 = sapply(correlationsSub, function(y){
    y[,x]
  })
  index = arrayInd(ind = which.max(temp2),.dim=dim(temp2))
  setNames(temp2[index],colnames(temp2)[index[,2]])
}, simplify = T, USE.NAMES=F)

res = pairwise.wilcox.test(topCorPerGrowthTraitValues,
                     names(topCorPerGrowthTraitValues), p.adjust.method="none")
cat.table.redmine(cbind(rownames(res$p.value), format(res$p.value, digits = 2, scientific = T)))


library(sinaplot)
pdfAndPng("graph/Perlstein/topCorPerGrowthTrait_corValues", 
          width=6, height=6, expr=({
# png("graph/Perlstein/topCorPerGrowthTrait_corValues.png",
#     height=800, width=800, pointsize=27)
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
}))
# dev.off()

# temp = data.frame(R = topCorPerGrowthTraitValues, 
#                   layer = names(topCorPerGrowthTraitValues))
# ggplot(data=temp, aes(x = layer, y = R, fill = layer))+
#   geom_violin(adjust = 1)+ 
#   scale_fill_manual(values=unname(colors[5:1]))+
#   scale_x_discrete(limits=c("transcripts", "pt","protein", "phosphoRes", "phospho")) +
#   geom_boxplot(width=0.1) +
#   labs(y="maximum\nabsolute Pearson's R", x = "")+
#   guides(fill=FALSE)+
#   theme(legend.position="none",
#         axis.title = element_text(size=20),
#         axis.text.x = element_text(size=20, face = "plain",  
#                                    colour = "black", 
#                                    angle=45, hjust = 1, vjust = 1),
#         axis.text.y = element_text(size=15, colour = "black",
#                                    margin=ggplot2::margin(l=10)),
#         axis.ticks =  element_blank(),
#         axis.ticks.y = element_line(size = 0.5),
#         axis.ticks.length = unit(x = 4,units = "mm"),
#         axis.line.y = element_line(size=0.5),
#         panel.background = element_rect(fill = NA),
#         panel.grid.major.x = element_blank())



# Vise versa: for each gene, on which layer 
# does it correlate most with growth traits?
correlationsPerGene = sapply(genesAvail, function(x){
  sapply(correlationsSub, function(y){y[x,]})
}, simplify=F)
correlationsPerGeneMax = sapply(genesAvail, function(x){
  names(which.max(apply(abs(correlationsPerGene[[x]]), 2, median)))
})
png("graph/Perlstein/topCorPerGene.png")
barplot(table(correlationsPerGeneMax)[5:1], col = colors, las = 1, ylab = "n growth Traits")
dev.off()
correlationsPerGeneMax2 = sapply(genesAvail, function(x){
  temp = apply(abs(correlationsPerGene[[x]]), 1, which.max)
  table(colnames(correlationsPerGene[[x]])[temp])
})
png("graph/Perlstein/topCorPerGeneDetail.png", 
    width=1200, height=800, pointsize=20)
par(mfrow =c(1,2), mar = c(3,4,0,0))
barplot(correlationsPerGeneMax2[5:1,], col = colors, las = 1, 
        space=0,border=NA, cex.names=0.2, horiz=T)
mtext(1,text="number of growth traits", line=2)
mtext(2,text="genes", line=2)
barplot(rowSums(correlationsPerGeneMax2)[5:1], 
        col = colors, las = 1, 
        names.arg=c("transcripts", "pt", "protein", "phRes", "phospho"))
dev.off()
# get only maximum correlations for each layer
maxCorrelationsSub = sapply(names(layersSub), function(type){
  predictors = t(get(layersSub[type]))
  temp = abs(cor(predictors, t(gPheno), use = "pair"))
  maxs = apply(temp,2,max)
  return(maxs)
})
# plots 
{png("graph/Perlstein/MaximumCorrelateGrowthWithMolecular_Matrixplot.png",
    width = 800, height = 800,pointsize=20)
par(oma = c(0,0,0,0), mar = c(1,1,0,0))
pairs(maxCorrelationsSub, 
      panel = panelfunction,
      cex.labels=1.6,  gap=0.5,
      labels=c("transcripts", "pt", "protein", "phRes", "phospho"))
mtext("maximum absolute R", side = 1, line=-0.5)
mtext("maximum absolute R", side = 2, line=-0.5)
dev.off()

png("graph/Perlstein/MaximumCorrelateGrowthWithMolecular_Customplot.png",
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


# boxplot of MAX correlations
png("graph/Perlstein/MaximumCorrelateGrowthWithMolecular.png")
boxplot(maxCorrelationsSub, ylab = "maximum abs. Pearson's R",
        col = colors, border = colors, las = 1)
boxplot(maxCorrelationsSub, col = rgb(0,0,0,0), add = T,
        outline = F, names = NA,  axes = F)
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
ggsave(filename = "graph/Perlstein/MaximumCorrelateGrowthWithMolecular_violin.png", 
       width=6, height=6)
}
# and a test:
temp = melt(maxCorrelationsSub)
pairwise.wilcox.test(x=temp$value, g=temp$Var2, paired=T)


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
            file="data/PerlsteinData/maxCorrelationsGeneSub.tsv",
            sep = "\t", quote=F, col.names=T, row.names=F)
# get the outliers, i.e. the best growth trait-molecular trait combos
outliers = lapply(names(layersSub), function(type){
  ord = order(maxCorrelationsSub[,type], decreasing=T)[1:15]
  temp = maxCorrelationsGeneSub[ord, type]
  cbind(condition = names(maxCorrelationsGeneSub[ord, type]), 
        gene = unname(maxCorrelationsGeneSub[ord, type]), 
        absR = format(unname(maxCorrelationsSub[ord, type]), digits=3))
})
names(outliers) = names(layersSub)
for(x in names(outliers)){print(x); cat.table.redmine(outliers[[x]])}
# get the genes that are top correlators with most growth traits
topGenes = apply(maxCorrelationsGeneSub,2, function(x){
  temp = split(x,compounds)
  temp = sapply(temp, unique)
  temp = table(unlist(temp))
  temp = temp[order(temp, decreasing = T)]
  data.frame(gene = names(temp[1:5]), times = as.integer(temp[1:5]))
})
cat.table.redmine(do.call(cbind, topGenes))
save(maxCorrelationsGeneSub, topGenes,
     file = "data/PerlsteinData/maxCorrelationsGeneSub.RData")
# other way to compute it for the phospho-traits: compute correlation
# for each peptide, then average for each protein
phlayersSub = c("phosphoRes" = "phResPhenoSub",
                "phospho" = "phPhenoSub")
correlationsPhSubAv = lapply(names(phlayersSub), function(type){
  predictors = t(get(phlayersSub[type]))
  temp = abs(cor(predictors, t(gPheno), use = "pair"))
  prots = split(colnames(predictors),
                phospho2protHash[colnames(predictors)])
  temp2 = t(sapply(prots, function(p){
    sub = abs(temp[p,,drop = F])
    apply(sub,2,mean)
  }))
  return(temp2)
})

correlationsPhSubMax = lapply(names(phlayersSub), function(type){
  predictors = t(get(phlayersSub[type]))
  temp = abs(cor(predictors, t(gPheno), use = "pair"))
  prots = split(colnames(predictors),
                phospho2protHash[colnames(predictors)])
  temp2 = t(sapply(prots, function(p){
    sub = abs(temp[p,,drop = F])
    apply(sub,2,max)
  }))
  return(temp2)
})
MKT1 = "YNL085W"
RPS31 = "YLR167W"

rowInds = row(correlationsPhSubAv[[2]])
rowN = rownames(correlationsPhSubAv[[2]])
cols = as.vector((rowInds == which(rowN == MKT1)) +
                   (rowInds == which(rowN == RPS31))*2 + 1)
png(paste0("graph/Perlstein/CorrelateGrowthWithPhospho",
           "_compareMethodsAverage_corrected.png"),
    width = 1000, height = 500)
par(mfrow = c(1,2))
plot(correlationsPhSubAv[[1]], correlationsSub[[4]], 
     xlab = "average individual abs. corr.", 
     ylab = "abs. corr. of averaged peptides", 
     main = "phosphoRes", col = cols)
plot(correlationsPhSubAv[[2]], correlationsSub[[5]], 
     xlab = "average individual abs. corr.", 
     ylab = "abs. corr. of averaged peptides", 
     main = "phospho", col = cols)
legend("topleft", col = unique(cols), pch = 1,
       legend=c("all genes", "RPS31", "MKT1"))
dev.off()
png(paste0("graph/Perlstein/CorrelateGrowthWithPhospho",
           "_compareMethodsMax_corrected.png"),
    width = 1000, height = 500)
par(mfrow = c(1,2))
plot(correlationsPhSubMax[[1]], correlationsSub[[4]], 
     xlab = "max abs. corr.", 
     ylab = "abs. corr. of averaged peptides", 
     main = "phosphoRes")
plot(correlationsPhSubMax[[2]], correlationsSub[[5]], 
     xlab = "max abs. corr.", 
     ylab = "abs. corr. of averaged peptides", main = "phospho")
dev.off()
wilcox.test(correlationsPhSubMax[[2]], correlationsSub[[5]], paired = T)
######


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
cat.table.redmine(cbind(layer = rownames(pVals),
                        "p-value" = format(pVals[,1], digits = 3), 
                        "Odds ratio" = format(pVals[,2], digits = 3)))
# make Venn diagrams
library(eulerr)
PeakMarkersPerLayer = sapply(c("eQTL", "ptQTL", "pQTL",  
                   "phResQTL", "phQTL", "gQTL"),function(x){
                     unique(unlist(sapply(get(x),function(y){
                       y$mostSignificantPredictor})))
                   })
sapply(c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"), function(x){
  png(paste0("graph/Perlstein/Venn_gQTL_",x,".png"))
  fit = euler(PeakMarkersPerLayer[c(x, "gQTL")])
  temp = plot(fit, quantities =T,
       fills = list(fill = c(colors[x], "grey"), alpha = 0.9),
       edges = F, labels=list(cex=2),
       main=paste("p-value =", format(pVals[x,"pval"],digits=2)))
  plot(temp)
  dev.off()
})
# for all markers, correlate number of growth traits affected with 
# number of each molecular type affected
nTraitsAffected = cbind(sapply(c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL"),
                               function(x){
                                 peaks = unlist(sapply(get(x),function(y){
                                   y$mostSignificantPredictor
                                 }))
                                 table(factor(peaks,levels=allMarkers))
                               }),
                        gQTL = sapply(allMarkers, function(i){
                          sum(sapply(gQTL, function(x){
                            i %in% x$mostSignificantPredictor}))
                        }))
nTraitsAffected2 = nTraitsAffected[rowSums(nTraitsAffected)>0,]
png("graph/Perlstein/nTraitsPerMarker_layerComparison.png", 
    width = 1000, height=200, pointsize=20)
par(mfrow = c(1,5), mar = c(3,1,1,0), oma = c(0,3,0,0))
sapply(colnames(nTraitsAffected2)[1:5], function(x){
  plot(nTraitsAffected2[,x],nTraitsAffected2[,"gQTL"],
       las = 1, cex.lab = 1.3,
       xlab = paste0("n Traits ",x), yaxt = "n",
       ylab = "",mgp =c(2,0.8,0))
  m = cor(nTraitsAffected2[,"gQTL"],nTraitsAffected2[,x], 
          method="spearman")
  legend("topright", bty ="n",
         legend=paste0("Spearman corr.: ",format(m, digits = 2)))
})
axis(2,outer=T,line=-1, las = 1)
mtext(side=2,text="n Traits gQTL",outer=T, line = 1)
dev.off()
#####


# compare with the local QTL for each layer #####
load("data/localQTL.RData")
eQTL_local = eQTL[localeQTL]
ptQTL_local = ptQTL[localptQTL]
pQTL_local = pQTL[localpQTL]
phResQTL_local = phResQTL[localphosphoProtQTL]
phQTL_local = phQTL[localphosphoQTL]

layers = c("eQTL_local", "ptQTL_local", "pQTL_local", 
           "phResQTL_local", "phQTL_local")
pVals_local = sapply(layers, function(type){
  qtl = get(type)
  peaks = sapply(qtl, function(x){x$mostSignificantPredictor})
  is_QTL = allMarkers %in% peaks
  tab = matrix(c(sum( is_QTL &  is_gQTL), 
                 sum( is_QTL & !is_gQTL), 
                 sum(!is_QTL &  is_gQTL),
                 sum(!is_QTL & !is_gQTL)), nrow = 2)
  test = fisher.test(tab, alternative = "greater")
  return(test$p.value)
})
t(t(pVals_local))
#####


# look at MKT1 #####
MKT1 = "YNL085W"
phQTLtargets = sapply(phQTL,function(x){x$target})
phQTLtargets = rownames(phPheno)[phQTLtargets]
MKT1peptide = unlist(prot2phosphoHash[MKT1])
MKT1_QTL = phQTL[phQTLtargets == MKT1peptide][[1]]$predictors[1,]
# there is only one MKT1 QTL. Where is it?
geno[MKT1_QTL,1:3]
# Other growth trait affected by the same QTL hotspot
affectedByMKT1QTL = sapply(gQTL, function(x){
  any(apply(x, 1,function(y){
    any(y["start"]:y["end"] %in% MKT1_QTL["start"]:MKT1_QTL["end"])
  }))
})
# mediation analysis
MKT1_peakQTL = phQTL[phQTLtargets == MKT1peptide][[1]]$mostSignificantPredictor
affectedByMKT1QTL_detail = unique(unlist(sapply(names(gQTL), function(x){
  x =gQTL[[x]]
  overlap = apply(x, 1,function(y){
    any(y["start"]:y["end"] %in% MKT1_QTL["start"]:MKT1_QTL["end"])
  })
  if(!any(overlap)){
    return(NULL)
  }else{
    return(unlist(strsplit(as.character(x[overlap,"targets"]),
                           split=",")))
  }
})))
mediationRes = calcMediation(marker=genotype[,MKT1_peakQTL], 
                             mediator=phPheno[MKT1peptide,],
                             outcomes=gPheno[as.integer(affectedByMKT1QTL_detail),])
MKT1layers = cbind("transript" = ePheno[MKT1,], 
                   "protein" = pPheno[MKT1,],
                   "phospho" = phPheno[MKT1peptide,])
mediationResPerLayer = apply(MKT1layers, 2,function(mediator){
  calcMediation(marker=genotype[,MKT1_peakQTL], 
                mediator=mediator,
                outcomes=gPheno[as.integer(affectedByMKT1QTL_detail),])
})
pVals = sapply(mediationRes, function(x){x[1,4]})
save(pVals, file = "data/PerlsteinData/mediation_MKT1_pVals.RData")
png("graph/Perlstein/mediationPvals_MKT1.png", pointsize=20, width=800, height=500)
par(mfrow = c(1,2))
plot(sapply(mediationRes, function(x){x[4,1]}), 
     -log10(pVals), las = 1,
     ylab= "-log10(p-value)", xlab = "% mediated")
abline(h = -log10(0.05), lty = 2)
plot(sapply(mediationRes, function(x){x[4,1]}), 
     -log10(p.adjust(pVals, method="fdr")), las = 1,
     ylab= "-log10(FDR)", xlab = "% mediated")
abline(h = -log10(0.1), lty = 2)
dev.off()
png("graph/Perlstein/mediationPvals_MKT1_finalplot.png", pointsize=20, width=800, height=500)
par(mfrow = c(1,2))
plot(sapply(mediationRes, function(x){x[4,1]}), 
     -log10(pVals), las = 1,
     ylab= "-log10(p-value)", xlab = "% mediated")
abline(h = -log10(0.05), lty = 2)
plot(sapply(mediationRes, function(x){x[4,1]}), 
     -log10(p.adjust(pVals, method="fdr")), las = 1,
     ylab= "-log10(FDR)", xlab = "% mediated")
abline(h = -log10(0.1), lty = 2)
dev.off()

# how many traits are affected at the other mol. layers.
layers = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL")
affectingOtherLayers = sapply(layers, function(layer){
  qtl = get(layer)
  affectedByMKT1QTL = sapply(qtl, function(x){
    x=x$predictors
    any(apply(x, 1,function(y){
      any(y["start"]:y["end"] %in% MKT1_QTL["start"]:MKT1_QTL["end"])
    }))
  })
  table(affectedByMKT1QTL)
})
# is MKT1 also regulated on the transcript,
# protein abundance, and phosphoRes level by this hotspot?
eID = which(rownames(ePheno) == MKT1)
eQTLtargets = sapply(eQTL,function(x){x$target})
eID %in% eQTLtargets
pID = which(rownames(pPheno) == MKT1)
pQTLtargets = sapply(pQTL,function(x){x$target})
pID %in% pQTLtargets
phResQTLtargets = sapply(phResQTL,function(x){x$target})
phResQTLtargets = rownames(phResPheno)[phResQTLtargets]
phResID = unlist(prot2phosphoHash[MKT1])
phResID %in% phResQTLtargets

# look at PBP1:
PBP1 = "YGR178C"
PBP1peptide = unlist(prot2phosphoHash[PBP1])
layers = c("transcripts" = "ePheno",  "pt" = "ptPheno", "protein" = "pPheno",
           "phRes" = "phResPhenoAverage", "phospho" = "phPhenoAverage")
png("graph/Perlstein/MKT1_PBP1_correlation.png",width=800, height=600, pointsize=15)
par(mfrow = c(2,3), mar = c(3,3,2,1))
for(i in names(layers)){
  phe = get(layers[i])
  plot(t(phe[c(PBP1, MKT1),]), pch = 19, 
       col = colors[i], main = i, las=1)
  abline(v = 0, h = 0,col = "grey")
  abline(0,1,col = "grey")
  abline(lm(phe[MKT1,] ~ phe[PBP1,]))
  mtext(side=1,text="PBP1", line = 2.1, cex = 0.8)
  mtext(side=2,text="MKT1", line = 1.7, cex = 0.8)
  cr = cor(phe[PBP1,],phe[MKT1,], use="pair")
  legend("topleft", paste0("R^2 = ", format(cr^2,digits=2)))
}
dev.off()

load("data/PerlsteinData/Perlstein_mappingData.RData")
pheno = mappingData$phenotype
PBP1pheno = sapply(names(layers), function(i){
  get(layers[i])[PBP1,]
})
colnames(PBP1pheno) = paste0("PBP1_", colnames(PBP1pheno))
MKT1pheno = sapply(names(layers), function(i){
  get(layers[i])[MKT1,]
})
colnames(MKT1pheno) = paste0("MKT1_", colnames(MKT1pheno))

PBP1_MKT1_cors = cor(PBP1pheno, MKT1pheno, use = "pair")
library(corrplot)
png("graph/Perlstein/correlate_MKT1_PBP1_corrplot.png")
corrplot(cor(cbind(PBP1pheno, MKT1pheno), use = "pair"))
dev.off()


PBP1cors = cor(t(pheno), PBP1pheno, use = "pair")
colnames(PBP1cors) = names(layers)
MKT1cors = cor(t(pheno), MKT1pheno, use = "pair")
colnames(MKT1cors) = names(layers)

png("graph/Perlstein/correlate_PBP1_histograms.png", 
    height=800, width= 650, pointsize = 20)
par(mfrow = c(5,1), mar = c(4,4,2,0), oma = c(0,0,3,0))
lim = c(-(max(abs(PBP1cors))+0.05),max(abs(PBP1cors))+0.05)
for (i in 1:ncol(PBP1cors)){
  hist(PBP1cors[,i], xlim = lim, 
       xlab = "",ylab="",las = 1, 
       col = colors[i], main = colnames(PBP1cors)[i])
  mtext(1,text="Pearson's R", line = 2.5, cex = 0.8)
  mtext(2,text="Frequency", line = 2.5, cex = 0.8)
  abline(v=0,  lwd = 3)
}
mtext(3,text="PBP1", outer=T, font = 2, adj=0.55)
dev.off()
png("graph/Perlstein/correlate_MKT1_histograms.png", 
    height=800, width= 650, pointsize = 20)
par(mfrow = c(5,1), mar = c(4,4,2,0), oma = c(0,0,3,0))
lim = c(-(max(abs(MKT1cors))+0.05),max(abs(MKT1cors))+0.05)
for (i in 1:ncol(MKT1cors)){
  hist(MKT1cors[,i], xlim = lim, 
       xlab = "",ylab="",las = 1, 
       col = colors[i], main = colnames(MKT1cors)[i])
  mtext(1,text="Pearson's R", line = 2.5, cex = 0.8)
  mtext(2,text="Frequency", line = 2.5, cex = 0.8)
  abline(v=0,  lwd = 3)
}
mtext(3,text="MKT1", outer=T, font = 2, adj=0.55)
dev.off()

PBP1cors = abs(cor(t(pheno), PBP1pheno, use = "pair"))
MKT1cors = abs(cor(t(pheno), MKT1pheno, use = "pair"))
png("graph/Perlstein/correlate_MKT1_PBP1_growth_boxplot.png",
    width=850, height=450, pointsize = 14)
par(mfrow = c(1,2), mar = c(3,4,3,1))
boxplot(MKT1cors, main = "MKT1", 
        ylab = "abs. R with growth traits", 
        line = 0.5, las = 1, col = colors)
boxplot(PBP1cors, main = "PBP1", 
        ylab = "abs. R with growth traits", 
        line = 0.5, las = 1, col = colors)
dev.off()
panelfunction= function(x,y){
  points(x,y)
  abline(0,1, col = "red")
  abline(h=0,v=0)
}
png("graph/Perlstein/correlate_PBP1_growth_matrixplot.png",
    width=800, height=800, pointsize=15)
pairs(PBP1cors, panel=panelfunction, main = "PBP1")
dev.off()
png("graph/Perlstein/correlate_MKT1_growth_matrixplot.png",
    width=800, height=800, pointsize=15)
pairs(MKT1cors, panel=panelfunction, main = "MKT1")
dev.off()

#####


# look at RPS31 #####
RPS31 = "YLR167W"
phQTLtargets = sapply(phQTL,function(x){x$target})
phQTLtargets = rownames(phPheno)[phQTLtargets]
RPS31peptide = unlist(prot2phosphoHash[RPS31])
RPS31_QTL = phQTL[phQTLtargets == RPS31peptide][[1]]$predictors[1,]
RPS31_peaks = phQTL[phQTLtargets == RPS31peptide][[1]]$mostSignificantPredictor

# there is only one QTL. Where is it?
RPS31_QTL_region = geno[RPS31_QTL,1:3]  # chrXIII:28519-29394 and chrXIII:42127-42196
RPS31_peak_region = (geno[RPS31_peaks,1:3]) # chrXIII:42127-42196
# since it is a transQTL: which genes lie in the QTL?
anno[anno$V2 == "ORF" &
       !is.na(anno$V11) &
       anno$V9 == as.integer(as.roman(substr(RPS31_peak_region[[1]],4,8))) & 
       ((anno$V10 > RPS31_peak_region$start & anno$V10 < RPS31_peak_region$end) | 
          (anno$V11 > RPS31_peak_region$start & anno$V11 < RPS31_peak_region$end) |
          (anno$V10 < RPS31_peak_region$start & anno$V11 > RPS31_peak_region$end)),]
withinRPS31QTL = (anno[anno$V2 == "ORF" &
                         !is.na(anno$V11) &
                         anno$V9 == as.integer(as.roman(substr(RPS31_QTL_region[[1]][1],4,8))) & 
                         ((anno$V10 > RPS31_QTL_region$start[1] & anno$V10 < RPS31_QTL_region$end[2]) | 
                            (anno$V11 > RPS31_QTL_region$start[1] & anno$V11 < RPS31_QTL_region$end[2]) |
                            (anno$V10 < RPS31_QTL_region$start[1] & anno$V11 > RPS31_QTL_region$end[2])),])
RPS31QTLgenes = withinRPS31QTL$V5
RPS31QTLgenes[RPS31QTLgenes == ""] = withinRPS31QTL$V4[RPS31QTLgenes == ""]
# Other growth trait affected by the same QTL 
affectedByRPS31QTL = sapply(gQTL, function(x){
  any(apply(x, 1,function(y){
    any(y["start"]:y["end"] %in% RPS31_QTL["start"]:RPS31_QTL["end"])
  }))
})
names(affectedByRPS31QTL[affectedByRPS31QTL])

# mediation analysis
RPS31_peakQTL = phQTL[phQTLtargets == RPS31peptide][[1]]$mostSignificantPredictor
affectedByMKT1QTL_detail = unique(unlist(sapply(names(gQTL), function(x){
  x =gQTL[[x]]
  overlap = apply(x, 1,function(y){
    any(y["start"]:y["end"] %in% RPS31_QTL["start"]:RPS31_QTL["end"])
  })
  if(!any(overlap)){
    return(NULL)
  }else{
    return(unlist(strsplit(as.character(x[overlap,"targets"]),
                           split=",")))
  }
})))
mediationRes = calcMediation(marker=genotype[,RPS31_peakQTL], 
                             mediator=phPheno[RPS31peptide,],
                             outcomes=gPheno[as.integer(affectedByMKT1QTL_detail),])
save(mediationRes, file = "data/PerlsteinData/mediationRes_PRS31.RData")
pVals = sapply(mediationRes, function(x){x[1,4]})
png("graph/Perlstein/mediationPvals_RPS31.png", 
    pointsize=20, width=800, height=500)
par(mfrow = c(1,2))
plot(sapply(mediationRes, function(x){x[4,1]}), 
     -log10(pVals), las = 1,
     ylab= "-log10(p-value)", xlab = "% mediated")
abline(h = -log10(0.05), lty = 2)
plot(sapply(mediationRes, function(x){x[4,1]}), 
     -log10(p.adjust(pVals, method="fdr")), las = 1,
     ylab= "-log10(FDR)", xlab = "% mediated")
abline(h = -log10(0.1), lty = 2)
dev.off()



# how many traits are affected at the other mol. layers.
layers = c("eQTL", "ptQTL", "pQTL", "phResQTL", "phQTL")
affectingOtherLayers = sapply(layers, function(layer){
  qtl = get(layer)
  affectedByRPS31QTL = sapply(qtl, function(x){
    x=x$predictors
    any(apply(x, 1,function(y){
      any(y["start"]:y["end"] %in% RPS31_QTL["start"]:RPS31_QTL["end"])
    }))
  })
  table(affectedByRPS31QTL)
})
cat.table.redmine(cbind(affectingLayer = c("FALSE", "TRUE"), affectingOtherLayers))
# which genes are regulated by this hotspot, and are they the candidates?
layers = c("eQTL" = "ePheno", "ptQTL" = "ptPheno", "pQTL" = "pPheno",
           "phResQTL" = "phResPheno", "phQTL" = "phPheno")
genesAffectedByRPS31qtl = sapply(names(layers), function(layer){
  qtl = get(layer)
  affectedByRPS31QTL = sapply(qtl, function(x){
    x=x$predictors
    any(apply(x, 1,function(y){
      any(y["start"]:y["end"] %in% RPS31_QTL["start"]:RPS31_QTL["end"])
    }))
  })
  targets = sapply(qtl[affectedByRPS31QTL], function(x){
    x$target
  })
  ph = rownames(get(layers[layer]))[targets]
  if(layer %in% c("phResQTL", "phQTL")){
    ph = phospho2protHash[ph]
  }
  id2name[ph]
})
RPS31QTLgenes %in% unlist(genesAffectedByRPS31qtl)
# is RPS31 also regulated on the transcript,
# protein abundance, and phosphoRes level by this hotspot?
eID = which(rownames(ePheno) == RPS31)
eQTLtargets = sapply(eQTL,function(x){x$target})
eID %in% eQTLtargets
pID = which(rownames(pPheno) == RPS31)
pQTLtargets = sapply(pQTL,function(x){x$target})
pID %in% pQTLtargets
phResQTLtargets = sapply(phResQTL,function(x){x$target})
phResQTLtargets = rownames(phResPheno)[phResQTLtargets]
phResID = unlist(prot2phosphoHash[RPS31])
phResID %in% phResQTLtargets


RPS31pos = as.integer(anno[grep(pattern=RPS31, x = anno$V4),9:11])
RPS31var = geno[geno$CHR == paste0("chr", as.roman(RPS31pos[1])) & 
                  ((geno$start > RPS31pos[2] & geno$start < RPS31pos[3]) |
                     (geno$end > RPS31pos[2] & geno$end < RPS31pos[3]) | 
                     (geno$start < RPS31pos[2] & geno$end > RPS31pos[3])),1:3]

straingeno = read.table("data/strain_genotype.tsv",
                        sep="\t", as.is = T, header = T)
RPS31snps = straingeno[straingeno$CHR == paste0("chr", as.roman(RPS31pos[1])) & 
                         straingeno$POS > (RPS31pos[2]-50) & 
                         straingeno$POS < (RPS31pos[3]+50),1:4]

load("data/PerlsteinData/Perlstein_mappingData.RData")
pheno = mappingData$phenotype
layers = c("transcripts" = "ePheno",  "pt" = "ptPheno", "protein" = "pPheno",
           "phRes" = "phResPhenoAverage", "phospho" = "phPhenoAverage")
RPS31pheno = sapply(names(layers), function(i){
  get(layers[i])[RPS31,]
})
RPS31cors = cor(RPS31pheno, t(pheno), use = "pair")

# png("graph/Perlstein/correlate_RPS31_histograms.png", 
#     height=500, width= 400, pointsize = 15)
pdfAndPng("graph/Perlstein/correlate_RPS31_histograms", 
          width=6, height=6, expr=({
par(mfrow = c(5,1), mar = c(0,0,0,0), oma = c(4,4,1,0.1), lab = c(4,1,7))
lim = c(-(max(abs(RPS31cors))+0.05),max(abs(RPS31cors))-.05)
for (i in 1:nrow(RPS31cors)){
  res = hist(RPS31cors[i,], xlim = lim, 
       xlab = "",ylab="",las = 1, main = "",
       col = colors[i], xaxt = "n") #, yaxp = c(0,max(res$counts),2)
  p = wilcox.test(RPS31cors[i,],mu=0)$p.value
  abline(h = 0)
  text(x=0.28, y = max(res$counts)/2, 
       labels=paste0(rownames(RPS31cors)[i], "\np=", format(p,digits=2)), 
       cex = 1.2, adj=0)
  abline(v=0,  lwd = 2, lty = 2)
}
axis(1,pos=0)
# mtext(3,text="RPS31", outer=T, font = 2, adj=0.55)
mtext(1,text="Pearson's R", outer=T, cex = 1, line = 2)
mtext(2,text="Frequency", outer=T, cex = 1, line = 2.5)
}))
# dev.off()

library(sinaplot)
pdfAndPng("graph/Perlstein/correlate_RPS31_sinaplot", 
          width=6, height=6, expr=({
# png("graph/Perlstein/correlate_RPS31_sinaplot.png",
#     height=500, width=500, pointsize=15)
par(mar = c(3,4,0.1,0.1))
temp = as.data.frame(t(RPS31cors))
pvals = apply(temp,2,function(x){
  wilcox.test(x,mu=0)$p.value
})
sinaplot(temp, col = colors,las = 1, xlab = "", tck = 0, mgp = c(3,0.5,0),
         yaxt = "n",ylim = c(-0.48, 0.48))
axis(2, las = 1)
title(ylab = "Correlation with growth traits",mgp = c(2.5,1,0))
abline(h=0)
boxplot(temp, col = colors, bty = "n",
        add = T,outline=F, boxwex = 0.2, ann = F, xaxt ="n", yaxt = "n")
text(x=1:5, y=0.45, labels=paste0("p=",format(pvals,digits=2)), cex = 0.8)
}))
# dev.off()
#####


# map growth traits based on molecular QTL #####
load("data/PerlsteinData/Perlstein_mappingData.RData")
pheno = mappingData$phenotype
layers = c("gQTL", "eQTL", "ptQTL", "pQTL",  "phResQTL", "phQTL")
predictions = sapply(layers, function(type){
  print(type)
  qtl = get(type)
  peaks = sort(unique(unlist(sapply(qtl, function(x){
    x$mostSignificantPredictor}))))
  predictors = genotype[,peaks]
  predictors = predictors[,!apply(predictors,2,function(x){any(is.na(x))})]
  res = apply(pheno,1,function(trait){
    dat = data.frame(trait, predictors)
    rf = ranger(trait ~ ., data=dat)
    return(c(rf$r.squared))
  })
})

temp = melt(predictions)
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
ggsave(filename = "graph/Perlstein/GrowthByMolecularQTL.png", height=6,width=6)


png("graph/Perlstein/GrowthByMolecularQTL_Matrixplot.png",
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

# hotspots #####
load("/cellnet/phosphoQTL/data/hotspotLocations.RData") #"markerLocations" "leader"   
load("/cellnet/phosphoQTL/data/hsInfo.RData") #"hsOvMat" "targetsByHS" "GOanalysis" 
hs_gQTL = apply(hsOvMat, 1, function(x){
  sum(unlist(sapply(QTL_list_merged_HR, function(y){
    y$chr == x["chr"] & 
      y$peakStart >= as.numeric(x["startPos"]) &
      y$peakEnd <= as.numeric(x["endPos"])
  })))
})
names(hs_gQTL) = hsOvMat[,"name"]
t(t(hs_gQTL))
png("graph/Perlstein/HS_nTraits_layerComparison.png", 
    width = 1000, height=200, pointsize=20)
par(mfrow = c(1,5), mar = c(3,1,1,0), oma = c(0,3,0,0))
sapply(colnames(hsOvMat[,7:11]),function(x){
  plot(as.integer(hsOvMat[,x]),hs_gQTL, las = 1, yaxt = "n",
       xlab = x, mgp = c(2,1,0))
  m = cor(as.integer(hsOvMat[,x]),hs_gQTL, method = "pearson")
  legend("topright", bty ="n",
         legend=paste0("Spearman corr.: ",format(m, digits = 2)))
})
axis(2,outer=T,line=-1, las = 1)
mtext(side=2,text="gQTL targets",outer=T, line = 1, cex = 0.7)
hs_cors = apply(hsOvMat[,7:11],2, function(x){
  cor(as.integer(x),hs_gQTL, method = "pearson")
})
dev.off()
#####


# recreate old figures about hotspots ####
load("/cellnet/phosphoQTL/data/hotspotLocations.RData") 
load("/cellnet/phosphoQTL/data/hsInfo.RData")
load("data/PerlsteinData/Perlstein_QTL_list_merged_HR.RData")
load("data/PerlsteinData/Perlstein_mappingData.RData")
pheno = mappingData$phenotype
hotspots=data.frame(hsOvMat, 
                    eStart=markerLocations[leader[,1],2],
                    eEnd=markerLocations[leader[,1],3], 
                    ptStart=markerLocations[leader[,2],2], 
                    ptEnd=markerLocations[leader[,2],3],
                    pStart=markerLocations[leader[,3],2], 
                    pEnd=markerLocations[leader[,3],3],
                    phResStart=markerLocations[leader[,4],2],
                    phResEnd=markerLocations[leader[,4],3],
                    phStart=markerLocations[leader[,5],2],
                    phEnd=markerLocations[leader[,5],3], 
                    stringsAsFactors=F)
row.names(hotspots)= hotspots[,1]
hotspots[,3:11] = apply(hotspots[,3:11],2,as.integer)
hotspots[,12:16] = apply(hotspots[,12:16],2,as.logical)
hotspots = cbind(hotspots, 
                 size=as.numeric(hotspots[,4])-as.numeric(hotspots[,3]),
                 peakStart = apply(hotspots[,17:26], 1,
                                   function(x) min(as.numeric(x),na.rm = T)),
                 peakEnd = apply(hotspots[,17:26], 1,
                                 function(x) max(as.numeric(x), na.rm = T)))
chrLengths = c(230218,813184,316620,1531933,576874,270161,
               1090940,562643,439888,745751,666816,1078177,
               924431,784333,1091291,948066)
names(chrLengths) = c("chrI","chrII","chrIII","chrIV","chrV","chrVI",
                      "chrVII","chrVIII","chrIX","chrX","chrXI",
                      "chrXII","chrXIII","chrXIV","chrXV","chrXVI")
chrPos = cumsum(c(1,chrLengths[-length(chrLengths)]))
names(chrPos) = names(chrLengths)
qtlStarts = QTL_list_merged_HR_table$start
qtlEnds = QTL_list_merged_HR_table$end

# 20 kb bins:
bins20kb = data.frame(do.call(rbind,sapply(names(chrLengths), function(chr){
  x = chrLengths[chr]
  start = seq(1,x,by=20000)
  end = c(start[-1] - 1, x)
  cbind(chr,start, end, size = end-start)
})), row.names=NULL, stringsAsFactors=F)
bins20kb$start = as.integer(bins20kb$start)
bins20kb$end = as.integer(bins20kb$end)
bins20kb$size = as.integer(bins20kb$size)
growthPer20kbbin = apply(bins20kb, 1, function(x){
  s = as.integer(x["start"])
  e = as.integer(x["end"])
  overlap = which(QTL_list_merged_HR_table$chr == x["chr"] & 
                    ((qtlStarts > s & qtlStarts < e) |
                       (qtlEnds > s & qtlEnds < e) |
                       (qtlStarts < s & qtlEnds > e)))
  rownames(pheno)[overlap]
})
nTraitsPerBin = sapply(growthPer20kbbin, length)/ bins20kb$size * 20000
positions = chrPos[bins20kb[,"chr"]] + 
  rowMeans(bins20kb[,c("start","end")])
png("graph/Perlstein/growthHotspots.png", 
    width=1000, height=400, pointsize=20)
par(mar = c(5,5,0.1,0.1))
plot(positions, nTraitsPerBin,  
     ylab = "n growth traits", xaxt = "n", las = 1, 
     xlim = c(0,sum(chrLengths)),xaxs = "i",type = "l",
     ylim=c(0,max(nTraitsPerBin)), lwd = 1.5)
abline(v = c(chrPos[-1]), lty = 1, col = "grey")
text(x=(chrPos + cumsum(chrLengths))/2, y = 0,
     pos=3,offset = -2.3,xpd = T,
     labels=names(chrLengths), srt = 45)
dev.off()
# identify which growth traits are affected by each hotspot

hotspotGrowthPeaks=apply(hotspots, 1, function(x){
  s = as.integer(x["peakStart"])
  e = as.integer(x["peakEnd"])
  QTL_list_merged_HR_table[QTL_list_merged_HR_table$chr == x["chr"] & 
                             ((qtlStarts > s & qtlStarts < e) |
                                (qtlEnds > s & qtlEnds < e) |
                                (qtlStarts < s & qtlEnds > e)),"compound"]
  # rownames(growth)[overlap]
  # return(overlap)
})
hotspotGrowthPeaks_table = do.call(rbind,sapply(names(hotspotGrowthPeaks), function(x){
  a = hotspotGrowthPeaks[[x]]
  if(length(a) < 1){
    return(c(hotspot=x, compounds = "none"))
  }
  b = c(x,rep("", length(a)-1))
  cbind(hotspot = b, compounds = a)
}))
cat.table.redmine(hotspotGrowthPeaks_table)
region = 500
hotspotGrowthPeaks_region=apply(hotspots, 1, function(x){
  s = as.integer(x["peakStart"])
  e = as.integer(x["peakEnd"])
  overlap = which(QTL_list_merged_HR_table$chr == x["chr"] & 
                    ((qtlStarts > s-region & qtlStarts < e+region) |
                       (qtlEnds > s-region & qtlEnds < e+region) |
                       (qtlStarts < s & qtlEnds > e)))
  rownames(pheno)[overlap]
  # return(overlap)
})
# how many per hotspot?
hotspots$nGrowthTraits = sapply(hotspotGrowthPeaks, length)
barplot(hotspots$nGrowthTraits, las = 2, names.arg= hotspots$name)
cor(hotspots[,7:11], hotspots$nGrowthTraits)
# density of growth traits in each hotspot
hotspotDensity = hotspots$nGrowthTraits / (hotspots$endPos-hotspots$startPos)
names(hotspotDensity) = hotspots$name
# density of growth traits outside of hotspots
# nonhotspots = sapply(QTL_list_merged_HR, function(x){})
nonhotspots = seq(nrow(QTL_list_merged_HR))[-unlist(hotspotGrowthPeaks)]
nonhotspotGrowth = QTL_list_merged_HR[nonhotspots,]
nonhotspotDensity = sapply(names(chrLengths), function(chr){
  Nqtl = sum(nonhotspotGrowth$chr == chr)
  hsregion = sum(hotspots$endPos[hotspots[,"chr"] == chr] - 
                   hotspots$startPos[hotspots[,"chr"] == chr])
  region=chrLengths[chr] - hsregion
  Nqtl/region
})
nonhotspotDensity = nonhotspotDensity*1000
hotspotDensity = hotspotDensity*1000
pdfAndPng(file = "graph/Perlstein/HS_QTLdensity", 
          width = 6,height = 6, expr = expression({
            boxplot(nonhotspotDensity,hotspotDensity, NA,
                    xlab = "number of growth QTL per kb",
                    las = 1, xaxt='n',yaxt = "n", 
                    frame=FALSE, outline = F)
            axis(side = 1,at = 0:5,  lwd.ticks = FALSE, labels=NA)
            axis(side = 1,at = 1:2,
                 labels=c("outside\nhotspots","within\nhotspots"),
                 lwd.ticks = FALSE, line = 0.5, lty = 0, las = 1)
            axis(side = 2, las = 1)
            mtext(2,text="number of gQTL per kb", line=3)
            axis(side = 2,lwd.ticks = FALSE, labels=NA, at = -1:10)
            points(rep(2, length(hotspotDensity)),
                   hotspotDensity, col = "red")
            points(rep(1, length(nonhotspotDensity)),
                   nonhotspotDensity, col = "blue")
            bottom5=names(hotspotDensity[order(hotspotDensity)[1:5]])
            legend("bottomright",legend = bottom5, bty = "n", 
                   title = expression(bold("bottom 5")),
                   x.intersp = 0, title.adj = 0.5, inset = 0.02)
            top5=names(sort(hotspotDensity, decreasing = T)[1:5])
            legend("topright",legend = top5, bty = "n", 
                   title = expression(bold("top 5   ")),
                   title.adj = 0.5, x.intersp = 0, inset = 0.02)
            legend("topleft", legend = c("hotspots", "chromosomes"), 
                   col = c("red", "blue"), pch = 1, inset = 0.05)
          }))
wilcox.test(hotspotDensity, nonhotspotDensity)
pdfAndPng(file = "graph/Perlstein/HS_nGrowthQTL_vs_nHsTargets", 
          width = 8,height = 8, Cairo = T,expr = expression({
            xpoints = apply(hotspots[,7:11],1,
                            function(x)sum(as.numeric(x)))
            plot(xpoints,hotspots$nGrowthTraits,
                 xlab = "number of hotspot targets", 
                 ylab = "number of growth QTL", las = 1, 
                 pch = 21, bg = ggplot2::alpha("black", 0.3), 
                 lwd = 1.5)
            abline(lm(hotspotDensity ~ xpoints, lwd = 1.5))
            colorsTemp = c("black", "dark blue", "light green",
                       "dark green", "pink", "red")
            alphas = c(0.3,0.3,0.5,0.3,0.6,0.3)
            for(i in 7:11){
              points(as.numeric(hotspots[,i]), hotspots$nGrowthTrait,
                     col = colorsTemp[(i-5)],
                     pch = 21, lwd = 1.5,
                     bg = ggplot2::alpha(colorsTemp[(i-5)], alphas[(i-5)]))
              abline(lm(hotspots$nGrowthTrait ~ as.numeric(hotspots[,i])),
                     col = colorsTemp[(i-5)], lwd = 1.5)
            }
            cors = format(c(cor(xpoints,hotspots$nGrowthTrait),
                            sapply(7:11,function(i){ 
                              cor(as.numeric(hotspots[,i]),
                                  hotspots$nGrowthTrait)})),
                          digits = 2)
            text = c("all QTL", "eQTL", "ptQTL", 
                     "pQTL", "phResQTL", "phQTL")
            text2 = paste0("(R=",cors,")")
            legend(x=1250,y=5,lty=1,lwd = 1.5,pch = 21,
                   pt.bg = ggplot2::alpha(colors, alphas),
                   col = c("black", "dark blue", "light green", 
                           "dark green", "pink", "red"),
                   legend =text, bty = 'n', inset = 0.02)
            legend(x=1620,y=5,
                   legend =text2, bty = 'n', inset = 0.02)
          }))
# are the correlations significant?
corp=c(cor.test(apply(hotspots[,7:11],1,
                      function(x)sum(as.numeric(x))),
                hotspots$nGrowthTrait)$p.value,
       sapply(7:11,function(i){ 
         cor.test(as.numeric(hotspots[,i]),
                  hotspots$nGrowthTrait)$p.value
       }))
names(corp) = c("all QTL", "eQTL", "ptQTL", 
                "pQTL", "phResQTL", "phQTL")
corp = p.adjust(corp)

#####



