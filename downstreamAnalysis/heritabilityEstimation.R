asRGB <- function(vec){
  rgb(vec[1],vec[2],vec[3])
}

# load data: #####
source("lib/general_function.R")
meta <- read.table("metadata/metadata.tsv",header=T,sep="\t", as.is = T)
dontuse <- c("RM11-1-1", "RM11-1b", "BY4724") 
strainsX <- sapply(meta$strain,FUN=function(s){
  s1 <- substr(s,start = 1,stop=1)
  ifelse(is.na(as.numeric(s1)),s,paste0("X",s))
})
meta$strainX = strainsX
genotype4mapping <- read.table("data/genotype_for_mapping.tsv",header=T)
genotype4mapping = genotype4mapping[,-(1:3)]
library(RColorBrewer)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
#####


# set of proteins where we have expression, pQTL and phosphorylation #####
load("data/proteinLevel.RData")# use proteinLevelBatchCorrected
load("data/expressionLevel.RData")#use eBatchGeneLengthCorrected
load("data/phosphoLevel.RData") # use phosphoLevelBatchCorrected
rm(eLevelNonBatch, eLevelBatchCorrected, phosphoLevelNonBatch, proteinLevelNonBatch)
proteins = rownames(proteinLevelBatchCorrected)
transcripts = rownames(eBatchGeneLengthCorrected)
phosphoTargets = rownames(phosphoLevelBatchCorrected)
phospho2protHash = phospho2prot[,2]
names(phospho2protHash) = phospho2prot[,1]
phosphoTargets = phospho2protHash[phosphoTargets]
phosphoTargets = unique(phosphoTargets)
genesAvail = intersect(intersect(transcripts, proteins),phosphoTargets)
#####


# Venn diagram of the Sets of genes  #####
library(eulerr)
library(RColorBrewer)
library(grid)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
colorMat <- col2rgb(colors)/255
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")
fit = euler(list("proteins" = proteins,
                  "transcripts" = transcripts, 
                  "phospho-peptides" = phosphoTargets))
#colvec <- c(colors[c(3,1,5)],rep("grey",20))
colvec <- c(asRGB(rowMeans(colorMat[,c(1,3)])),
            asRGB(rowMeans(colorMat[,c(5,1)])),
            asRGB(rowMeans(colorMat[,c(5,1,3)])),
            asRGB(colorMat[,c(3)]),
            asRGB(colorMat[,c(1)]),
            asRGB(colorMat[,c(5)]))

eg = plot(fit, fills = colors[c(3,1,5)], 
     alpha = 0.5, edges = F,
     quantities = list(col = colvec, font = 2, cex = 2.5),
     labels = list(col = colors[c(3,1,5)], font = 2, cex = 2.5))
#eg$gp <- gpar(mar=c(10,0,0,0))
eg$vp <- viewport(w = .7, h = .7, gp = gpar(col="blue"))
eg
# eg$quantities$gp$col = eg$fills$gp$fill[c(4,6,7,1:3)]
# eg$data$labels$x = eg$data$labels$x + c(-10, 0, -12.3)
# eg$data$labels$y = eg$data$labels$y + c(-0.7, -0.2, -1)
# eg$data$quantities$x = eg$data$quantities$x + c(0,0,0,-10,0,-12.3)
# eg$data$quantities$y = eg$data$quantities$y + c(0,0,0,-0.7, -0.2, -1)
# eg$data$xlim = c(-60, 44)
pdfAndPng(file = "graph/Venn_availableMeasurements",8,8, expression({
  plot(eg)
}))
rm(eg, fit)
#####


# functions/packages needed for heritability estimation #####
library(lme4)
library(bootstrap)
theta <- function(x, xdata){ # function for the estimation of heritabiliy, to be used within jackknife
  model = lmer(y~0+(1|strain), data=xdata[x,])
  res = as.data.frame(VarCorr(model))
  if(!all(dim(res) == c(2,5)))
    return(NA)
  g = res$vcov[1]
  e = res$vcov[2]
  return(g/(g+e))
}
#####


# protein levels heritability #####
proteinLevelBatchCorrected = proteinLevelBatchCorrected[genesAvail,]
subMeta <- subset(meta,culture %in% colnames(proteinLevelBatchCorrected))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))

H2protein = t(sapply(1:nrow(proteinLevelBatchCorrected), function(i){
  cat(i,' ')
  y = proteinLevelBatchCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  
  # get H2
  model = lmer(values~1+(1|strain), data=df)
  res = as.data.frame(VarCorr(model))
  if(!all(dim(res) == c(2,5)))
    return(NA)
  g = res$vcov[1]
  e = res$vcov[2]
  h = g/(g+e)
  
  # get SE
  jkn <- jackknife(1:32, theta, xdata=df)
 
  return(c(H = mean(results$jack.values), se = jkn$jack.se))
}))
save(H2protein, file = "data/heritabilityProtein.RData")


# compare real and permuted estimates
H2proteinPermute = t(sapply(1:nrow(proteinLevelBatchCorrected), function(i){
  cat(i,' ')
  y = proteinLevelBatchCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = sample(strainNames))
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
# H2proteinPermute2 = t(sapply(1:nrow(proteinLevelBatchCorrected), function(i){
#   cat(i,' ')
#   sapply(1:100, function(j){
#     y = proteinLevelBatchCorrected[i, unlist(CultureByReplicates)]
#     y = scale(y)
#     df = data.frame(y = y, strain = sample(strainNames))
#     model = lmer(y~0+(1|strain), data=df)
#     res = as.data.frame(VarCorr(model))
#     if(!all(dim(res) == c(2,5)))
#       return(NA)
#     g = res$vcov[1]
#     e = res$vcov[2]
#     return(g/(g+e))
#   })
# }))
# boxplot(t(H2proteinPermute2[1:100,]))
# points(H2protein[1:100,1], col = "red", pch = 19)


H2proteinTest = wilcox.test(H2protein[,1], H2proteinPermute[,1], alternative = "two.sided")
pdfAndPng("graph/heritability/heritabilityHistogramsProtein", 8,8,expression({
  par(mfrow = c(2,1))
  hist(H2protein[,1], col = colors["pQTL"], xlim = c(0,1),breaks = seq(0,1,0.05),
       xlab = "broad-sense heritability", ylab = "number of proteins", main = "")
  hist(H2proteinPermute[,1], col = colors["pQTL"], xlim = c(0,1),breaks = seq(0,1,0.05),
       xlab = "broad-sense heritability of Permutations", ylab = "number of proteins", main = "")
}))


# compare how it looks if only parent replicates are used:
CultureByParents = CultureByReplicates[1:2]
parentNames = unlist(sapply(names(CultureByParents), function(x) 
  rep(x, length(CultureByParents[[x]]))))
H2proteinParents = t(sapply(1:nrow(proteinLevelBatchCorrected), function(i){
  cat(i,' ')
  y = proteinLevelBatchCorrected[i, unlist(CultureByParents)]
  y = scale(y)
  df = data.frame(y = y, strain = parentNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
pdfAndPng("graph/heritability/heritabilityProteinParentsVsAllstrains", 8,8,expression({
  x = H2protein[,1]
  SEx = H2protein[,2]/2
  y = H2proteinParents[,1]
  SEy = H2proteinParents[,2]/2
  plot (x,y, xlim = c(0,1), ylim = c(0,1), pch = 19, cex = 0.5,
        xlab = "H2 expression", ylab = "H2 protein")
  arrows(x,y-SEy,x,y+SEy, code=3, length=0, angle = 90, col = "grey")
  arrows(x-SEx,y,x+SEx,y, code=3, length=0, angle = 90, col = "grey")
  plot (x,y, xlim = c(0,1), ylim = c(0,1), pch = 19, cex = 0.5,
        xlab = "H2 all strains", ylab = "H2 parent strains")
  arrows(x,y-SEy,x,y+SEy, code=3, length=0, angle = 90, col = "grey")
  arrows(x-SEx,y,x+SEx,y, code=3, length=0, angle = 90, col = "grey")
  points (x,y, xlim = c(0,1), ylim = c(0,1), pch = 19, cex = 0.5)
  abline(a = 0, b = 1)
}))# no big difference

#####

# protein Levels Foss method #####
proteinLevelBatchCorrected = proteinLevelBatchCorrected[genesAvail,]
subMeta <- subset(meta,culture %in% colnames(proteinLevelBatchCorrected))
strains = subMeta$strainX
VossHdisregardRepl = apply(proteinLevelBatchCorrected,1,function(x){
  substrain = strains[!is.na(x)]
  x = x[!is.na(x)]
  gm = mean(x)
  indVar = sum((x-gm)^2)
  xbystrains = split(x,substrain)
  Ns = sapply(xbystrains, length)
  Xs = sapply(xbystrains, mean)
  groupVar = sum((Xs-gm)^2*Ns)
  h = groupVar/indVar
})

replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
withRepl = subMeta$culture[subMeta$strainX %in% replicates]
SubproteinLevelBatchCorrected = proteinLevelBatchCorrected[genesAvail,withRepl]
subMeta <- subset(meta,culture %in% colnames(SubproteinLevelBatchCorrected))
strains = subMeta$strainX
VossH = apply(SubproteinLevelBatchCorrected,1,function(x){
  substrain = strains[!is.na(x)]
  x = x[!is.na(x)]
  x = scale(x)[,1]
  gm = mean(x)
  indVar = sum((x-gm)^2)
  xbystrains = split(x,substrain)
  Ns = sapply(xbystrains, length)
  Xs = sapply(xbystrains, mean)
  groupVar = sum((Xs-gm)^2*Ns)
  h = groupVar/indVar
})
# plot(H2protein[,1],VossH, col = apply(SubproteinLevelBatchCorrected,1,function(x) any(is.na(x)))+1)
# plot(H2protein[,1],VossHdisregardRepl, col = apply(SubproteinLevelBatchCorrected,1,function(x) any(is.na(x)))+1)
# abline(0,1)
save(VossHdisregardRepl, VossH, file = "data/heritabilityProtein_Fossmethod.RData")
#####


# eqtl #####
eBatchGeneLengthCorrected = eBatchGeneLengthCorrected[genesAvail,]
subMeta <- subset(meta,culture %in% colnames(eBatchGeneLengthCorrected))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))


H2expression = t(sapply(1:nrow(eBatchGeneLengthCorrected), function(i){
  cat(i,' ')
  y = eBatchGeneLengthCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2expression, file = "data/heritabilityExpression.RData")

# test difference between real and permuted
H2expressionPermute = t(sapply(1:nrow(eBatchGeneLengthCorrected), function(i){
  cat(i,' ')
  y = eBatchGeneLengthCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = sample(strainNames))
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
H2expressionTest = wilcox.test(H2expression[,1], H2expressionPermute[,1], alternative = "two.sided")
pdfAndPng("graph/heritability/heritabilityHistogramsExpression", 8,8,expression({
  par(mfrow = c(2,1))
  hist(H2expression[,1], col = colors["eQTL"], xlim = c(0,1),breaks = seq(0,1,0.05),
       xlab = "broad-sense heritability", ylab = "number of transcripts", main = "")
  hist(H2expressionPermute[,1], col = colors["eQTL"], xlim = c(0,1),breaks = seq(0,1,0.05),
       xlab = "broad-sense heritability of Permutations", ylab = "number of transcripts", main = "")
}))
#####


# phospholevel #####
peptides2Use = phospho2prot[phospho2prot[,2] %in% genesAvail,1]
phosphoLevelBatchCorrected = phosphoLevelBatchCorrected[peptides2Use,]
subMeta <- subset(meta,culture %in% colnames(phosphoLevelBatchCorrected))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))


H2phospho = t(sapply(1:nrow(phosphoLevelBatchCorrected), function(i){
  cat(i,' ')
  y = phosphoLevelBatchCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2phospho, file = "data/heritabilityPhospholevel.RData")

# test difference between real and permuted
H2phosphoPermute = t(sapply(1:nrow(phosphoLevelBatchCorrected), function(i){
  cat(i,' ')
  y = phosphoLevelBatchCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = sample(strainNames))
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
H2phosphoTest = wilcox.test(H2phospho[,1], H2phosphoPermute[,1], alternative = "two.sided")
pdfAndPng("graph/heritability/heritabilityHistogramsPhospho", 8,8,expression({
  par(mfrow = c(2,1))
  hist(H2phospho[,1], col = colors["phQTL"], xlim = c(0,1),breaks = seq(0,1,0.04),
       xlab = "broad-sense heritability", ylab = "number of phospho-peptides", main = "")
  hist(H2phosphoPermute[,1], col = colors["phQTL"], xlim = c(0,1),breaks = seq(0,1,0.04),
       xlab = "broad-sense heritability of Permutations", ylab = "number of phospho-peptides", main = "")
}))
#####

# phospho Foss method ######
peptides2Use = phospho2prot[phospho2prot[,2] %in% genesAvail,1]
phosphoLevelBatchCorrected = phosphoLevelBatchCorrected[peptides2Use,]
subMeta <- subset(meta,culture %in% colnames(phosphoLevelBatchCorrected))
strains = subMeta$strainX
VossH_phospho_disregardRepl = apply(phosphoLevelBatchCorrected,1,function(x){
  substrain = strains[!is.na(x)]
  x = x[!is.na(x)]
  gm = mean(x)
  indVar = sum((x-gm)^2)
  xbystrains = split(x,substrain)
  Ns = sapply(xbystrains, length)
  Xs = sapply(xbystrains, mean)
  groupVar = sum((Xs-gm)^2*Ns)
  h = groupVar/indVar
})

replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
withRepl = subMeta$culture[subMeta$strainX %in% replicates]
SubphosphoLevelBatchCorrected = phosphoLevelBatchCorrected[peptides2Use,withRepl]
subMeta <- subset(meta,culture %in% colnames(SubphosphoLevelBatchCorrected))
strains = subMeta$strainX
VossH_phospho = apply(SubphosphoLevelBatchCorrected,1,function(x){
  substrain = strains[!is.na(x)]
  x = x[!is.na(x)]
  x = scale(x)[,1]
  gm = mean(x)
  indVar = sum((x-gm)^2)
  xbystrains = split(x,substrain)
  Ns = sapply(xbystrains, length)
  Xs = sapply(xbystrains, mean)
  groupVar = sum((Xs-gm)^2*Ns)
  h = groupVar/indVar
})
load("data/heritabilityPhospholevel.RData")
plot(H2phospho[,1],VossH_phospho, 
     col = apply(SubphosphoLevelBatchCorrected,1,function(x) any(is.na(x)))+1)
plot(H2phospho[,1],VossH_phospho_disregardRepl, 
     col = apply(SubphosphoLevelBatchCorrected,1,function(x) any(is.na(x)))+1)
abline(0,1)
save(VossH_phospho_disregardRepl, VossH_phospho,
     file = "data/heritabilityPhospho_Fossmethod.RData")

#######

# phosphoProtResiduals #####
load("data/phosphoNew/phosphoProt.RData") # phosphoProtResiduals
peptides2Use = phospho2prot[phospho2prot[,2] %in% genesAvail,1]
phosphoProtResiduals = phosphoProtResiduals[peptides2Use,]
subMeta <- subset(meta,culture %in% colnames(phosphoProtResiduals))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))

H2phosphoRes = t(sapply(1:nrow(phosphoProtResiduals), function(i){
  cat(i,' ')
  y = phosphoProtResiduals[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2phosphoRes, file = "data/heritabilityphosphoRes.RData")

# test difference between real and permuted
H2phosphoResPermute = t(sapply(1:nrow(phosphoProtResiduals), function(i){
  cat(i,' ')
  y = phosphoProtResiduals[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = sample(strainNames))
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
H2phosphoResTest = wilcox.test(H2phosphoRes[,1], H2phosphoResPermute[,1], alternative = "two.sided")
pdfAndPng("graph/heritability/heritabilityHistogramsPhosphoRes", 8,8,expression({
  par(mfrow = c(2,1))
  hist(H2phosphoRes[,1], col = colors["phResQTL"], xlim = c(0,1),breaks = seq(0,1,0.04),
       xlab = "broad-sense heritability", ylab = "number of phospho-peptides", main = "")
  hist(H2phosphoResPermute[,1], col = colors["phResQTL"], xlim = c(0,1),breaks = seq(0,1,0.04),
       xlab = "broad-sense heritability of Permutations", ylab = "number of phospho-peptides", main = "")
}))
#####


# RNA prot Residuals #####
load("data/protRna.RData")
protRnaResiduals = protRnaResiduals[genesAvail,]
subMeta <- subset(meta,culture %in% colnames(protRnaResiduals))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))

H2protRna = t(sapply(1:nrow(protRnaResiduals), function(i){
  cat(i,' ')
  y = protRnaResiduals[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2protRna, file = "data/heritabilityProtRna.RData")

# test difference between real and permuted
H2protRnaPermute = t(sapply(1:nrow(protRnaResiduals), function(i){
  cat(i,' ')
  y = protRnaResiduals[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = sample(strainNames))
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
H2protRnaTest = wilcox.test(H2protRna[,1], H2protRnaPermute[,1], alternative = "two.sided")
pdfAndPng("graph/heritability/heritabilityHistogramsProtRna", 8,8,expression({
  par(mfrow = c(2,1))
  hist(H2protRna[,1], col = colors["phResQTL"], xlim = c(0,1),breaks = seq(0,1,0.04),
       xlab = "broad-sense heritability", ylab = "number of proteins", main = "")
  hist(H2protRnaPermute[,1], col = colors["phResQTL"], xlim = c(0,1),breaks = seq(0,1,0.04),
       xlab = "broad-sense heritability of Permutations", ylab = "number of proteins", main = "")
}))
#####


# plot heritabilies #####
load("data/heritabilityProtein.RData")
load("data/heritabilityExpression.RData")
load("data/heritabilityPhospholevel.RData")
load("data/heritabilityphosphoRes.RData")
load("data/heritabilityProtRna.RData")
# H2phospho, H2protein, H2expression, H2phosphoRes, H2protRna

#colors
library(RColorBrewer)
colors = brewer.pal(n = 12, name = "Paired")[c(2,3:6,8,10)]
names(colors) = c("eQTL",  "ptQTL", "pQTL", "phResQTL", "phQTL",  "RM",  "BY")

# histogram
cexlab = 1.6
cexaxis = 1.2
pdfAndPng("graph/heritability/heritabilityHistogramsWithResiduals", 8,12,Cairo = T,  expression({
  par(mfrow = c(5,1), mar=c(4,6,1,0)+.1)
  hist(H2expression[,1], col = colors["eQTL"], xlim = c(0,1),breaks = 20,
       cex.lab = cexlab, cex.axis = cexaxis, 
       xlab = "", ylab = "number of\ntranscripts", main = "")
  hist(H2protRna[,1], col = colors["ptQTL"], xlim = c(0,1), breaks = 20,
       cex.lab = cexlab, cex.axis = cexaxis,
       xlab = "", ylab = "number of\nRNA-prot residuals", main = "")
  hist(H2protein[,1], col = colors["pQTL"], xlim = c(0,1), breaks = 20,
       cex.lab = cexlab, cex.axis = cexaxis,
       xlab = "", ylab = "number of\nproteins", main = "")
  hist(H2phosphoRes[,1], col = colors["phResQTL"], xlim = c(0,1), breaks = 20,
       cex.lab = cexlab, cex.axis = cexaxis,
       xlab = "", ylab = "number of\nphospho-residuals", main = "")
  hist(H2phospho[,1], col = colors["phQTL"], xlim = c(0,1), breaks = 20,
       cex.lab = cexlab, cex.axis = cexaxis,
       xlab = "broad-sense heritability", ylab = "number of\nphospho-peptides", main = "")
}))

# violoin plot
library(ggplot2)
H2df = data.frame(rbind(H2expression, H2protRna, H2protein, H2phosphoRes, H2phospho),
                  source = c(rep("transcripts", nrow(H2expression)),
                             rep("protein-residuals", nrow(H2protRna)),
                             rep("proteins", nrow(H2protein)),
                             rep("phospho-residuals", nrow(H2phosphoRes)),
                             rep("phospho-peptides", nrow(H2phospho))))
H2df[,3] = factor(H2df[,3],levels(H2df[,3])[c(4,5,3,2,1)])
ggplot(H2df, aes(x = factor(source), y = H, fill=source)) + 
  scale_x_discrete(limits=as.character(unique(unique(H2df$source))))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  guides(fill=FALSE)+
  labs(y="broad-sense heritability", x = "")+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=20,face = "plain",  colour = "black", 
                                   angle=45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15, colour = "black",margin=margin(l=10)),
        axis.ticks =  element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(x = 4,units = "mm"),
        axis.line.y = element_line(size=0.5))+
  geom_violin(adjust = 1)+ # , draw_quantiles = c(0.25, 0.5, 0.75)
  stat_summary(fun.y = "median", geom = "point")+
  geom_boxplot(width = 0.2)+
  theme(panel.background = element_rect(fill = NA),
    panel.grid.major.x = element_blank())+
  #ylim(0,1)+
  scale_fill_manual(values=unname(colors[c(3,1,2,4,5)]))
ggsave(filename = "graph/heritability/heritabilityViolinplot_Box.pdf")
ggsave(filename = "graph/heritability/heritabilityViolinplot_Box.png")

png("graph/heritability/heritabilitySinaplot.png")
library(sinaplot)
dat = split(H2df$H, H2df$source)
dat = dat[c(2,3,1,4,5)]
sinaplot(dat, pch = 19, col =colors[1:5], las = 2 )
dev.off()

# layered density plot
pdfAndPng("graph/heritability/heritabilityDensityplot", width = 10, height = 10, expr = expression({
  plot(NA, xlim = c(0,1), ylim = c(0,3.5), xlab = "Broad-sense heritability", ylab = "Density", cex.lab = 1.3)
  kern = 1; linew = 4.5
  lines(density(H2protRna[,1], adjust = kern, cut = 0), col = colors["ptQTL"], lwd = linew)
  lines(density(H2protein[,1], adjust = kern, cut = 0), col = colors["pQTL"], lwd = linew)
  lines(density(H2phosphoRes[,1], adjust = kern, cut = 0), col = colors["phResQTL"], lwd = linew)
  lines(density(H2phospho[,1], adjust = kern, cut = 0), col = colors["phQTL"], lwd = linew)
  lines(density(H2expression[,1], adjust = kern, cut = 0), col = colors["eQTL"], lwd = linew)
  legend("topright", col = colors, lwd = linew,
         legend = c("RNA", "RNA-protein residuals", "protein", "protein-phospho residuals", "phospho"))
}))


# compare expression and protein
x = H2expression[,1]
SEx = H2expression[,2]/2
y = H2protein[,1]
SEy = H2protein[,2]/2
pdfAndPng("graph/heritability/heritabilityCorrProteinExpression", 8,6,expression({
  plot (x,y, xlim = c(0,1), ylim = c(0,1), pch = 19, cex = 0.5,
        xlab = "H2 expression", ylab = "H2 protein")
  arrows(x,y-SEy,x,y+SEy, code=3, length=0, angle = 90, col = "grey")
  arrows(x-SEx,y,x+SEx,y, code=3, length=0, angle = 90, col = "grey")
  points (x,y, xlim = c(0,1), ylim = c(0,1), pch = 19, cex = 0.5)
  abline(a = 0, b = 1)
}))

#####

# plot heritabilities against each other #####
# average heritabilities for phospho traits belonging to the same protein
load("data/heritabilityExpression.RData")
load("data/heritabilityProtRna.RData")
load("data/heritabilityProtein.RData")
load("data/heritabilityphosphoRes.RData")
load("data/heritabilityPhospholevel.RData")
# H2expression, H2protRna, H2protein, H2phosphoRes,  H2phospho 
peptides2Use = phospho2prot[phospho2prot[,2] %in% genesAvail,1]
protein2pep = split(phospho2prot[phospho2prot[,2] %in% genesAvail,1],
                    phospho2prot[phospho2prot[,2] %in% genesAvail,2])
protein2pep = protein2pep[genesAvail]
H2s = data.frame(gene = genesAvail,
                 transcript = H2expression[,1],
                 # transcriptSE = H2expression[,2],
                 pt = H2protRna[,1],
                 # ptSE = H2protRna[,2],
                 protein = H2protein[,1],
                 # proteinSE = H2protein[,2],
                 t(sapply(protein2pep, function(x){
                   sub = H2phosphoRes[peptides2Use %in% x,1]
                   return(c(phRes = mean(sub), phResSD = sd(sub)))
                 })),
                 t(sapply(protein2pep, function(x){
                   sub = H2phospho[peptides2Use %in% x,1]
                   return(c(phospho = mean(sub), phosphoSD = sd(sub)))
                 })))
png("graph/heritability/compareLayerHeritabilities.png", width=900, height=900, pointsize=20)
pairs(H2s[,c(2:5,7)], xlim=c(0, 1), ylim=c(0, 1),
      upper.panel=function(x,y, ...){ 
        points(x, y)
        abline(0,1)
      }, 
      lower.panel=function(x,y, ...){ 
        text(0.5,0.5,paste0("r = ",format(cor(x,y), digits = 2)), cex = 2)
      })
title(xlab = "broad-sense heritability", line = 3.5, cex.lab =1.2)
title(ylab = "broad-sense heritability", line = 2.8, cex.lab = 1.2)
dev.off()
# load("data/expressionLevel.RData")#use eBatchGeneLengthCorrected
# eBatchGeneLengthCorrected = eBatchGeneLengthCorrected[genesAvail,]
# load("data/protRna.RData")
# protRnaResiduals = protRnaResiduals[genesAvail,]
# load("data/proteinLevel.RData")# use proteinLevelBatchCorrected
# proteinLevelBatchCorrected = proteinLevelBatchCorrected[genesAvail,]
# load("data/phosphoNew/phosphoProt.RData") # phosphoProtResiduals
# 
# phosphoProtResiduals = phosphoProtResiduals[peptides2Use,]
# load("data/phosphoLevel.RData") # use phosphoLevelBatchCorrected
# peptides2Use = phospho2prot[phospho2prot[,2] %in% genesAvail,1]
# phosphoLevelBatchCorrected = phosphoLevelBatchCorrected[peptides2Use,]
#####


# compute heritabilities for all traits not only those in the overlap #####
#proteins
load("data/proteinLevel.RData")
#proteinLevelBatchCorrected = proteinLevelBatchCorrected[genesAvail,]
subMeta <- subset(meta,culture %in% colnames(proteinLevelBatchCorrected))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))

H2proteinFull = t(sapply(1:nrow(proteinLevelBatchCorrected), function(i){
  cat(i,' ')
  y = proteinLevelBatchCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
rownames(H2proteinFull) <- rownames(proteinLevelBatchCorrected)
save(H2proteinFull, file = "data/heritabilityProteinFull.RData")

#all transcripts
#eBatchGeneLengthCorrected = eBatchGeneLengthCorrected[genesAvail,]
load("data/expressionLevel.RData")
subMeta <- subset(meta,culture %in% colnames(eBatchGeneLengthCorrected))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))


H2expressionFull = t(sapply(1:nrow(eBatchGeneLengthCorrected), function(i){
  cat(i,' ')
  y = eBatchGeneLengthCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2expressionFull, file = "data/heritabilityExpressionFull.RData")

#all pt-traits
load("data/protRna.RData")
subMeta <- subset(meta,culture %in% colnames(protRnaResiduals))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))

H2ptFull = t(sapply(1:nrow(protRnaResiduals), function(i){
  cat(i,' ')
  y = protRnaResiduals[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2ptFull, file = "data/heritabilityPtFull.RData")

#all raw ph traits
load("data/phosphoLevel.RData")
subMeta <- subset(meta,culture %in% colnames(phosphoLevelBatchCorrected))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))
H2phosphoFull = t(sapply(1:nrow(phosphoLevelBatchCorrected), function(i){
  cat(i,' ')
  y = phosphoLevelBatchCorrected[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2phosphoFull, file = "data/heritabilityPhospholevelFull.RData")

#all ph res traits
load("data/phosphoNew/phosphoProt.RData")
subMeta <- subset(meta,culture %in% colnames(phosphoProtResiduals))
replicates = table(subMeta$strainX)
replicates = names(replicates[replicates > 2])
CultureByReplicates = sapply(replicates, function(i){ 
  subMeta[subMeta$strainX == i,"culture"] }, simplify = F)
strainNames = unlist(sapply(names(CultureByReplicates), function(x) 
  rep(x, length(CultureByReplicates[[x]]))))
H2phosphoResFull = t(sapply(1:nrow(phosphoProtResiduals), function(i){
  cat(i,' ')
  y = phosphoProtResiduals[i, unlist(CultureByReplicates)]
  y = scale(y)
  df = data.frame(y = y, strain = strainNames)
  results <- jackknife(1:32, theta, xdata=df)
  return(c(H = mean(results$jack.values), se = results$jack.se))
}))
save(H2phosphoResFull, file = "data/heritabilityphosphoResFull.RData")

###compute the proportion of traits with H2>=50%
load("data/heritabilityExpressionFull.RData")
load("data/heritabilityProteinFull.RData")
load("data/heritabilityPtFull.RData")
load("data/heritabilityPhospholevelFull.RData")
load("data/heritabilityphosphoResFull.RData")

sapply(list(H2expressionFull,H2ptFull,H2proteinFull,H2phosphoResFull,H2phosphoFull),FUN=function(mat){
  return(c(sum(mat[,1]>0.5),sum(mat[,1]>0.5)/nrow(mat)))
})
######