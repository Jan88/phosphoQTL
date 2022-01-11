library(xlsx)
library(RFQTL)
library(lme4)

# load data from Nogami et al. 2007 #####
# downloaded from the page of Gael Yvert http://www.ens-lyon.fr/LBMC/gisv/index.php/en/downloads 
# and the supplement here: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030031#
path = "data/Nogami2007/"
# phenotypes from parents are here
phenoBY = read.table(paste0(path, "RawMorphoData_PLGE2007/BY4716_morpho_9rep.dat"),
                     as.is = T, header = T, row.names=1)
phenoRM = read.table(paste0(path, "RawMorphoData_PLGE2007/YEF1946_morpho_9rep.dat"),
                     as.is = T, header = T, row.names=1)
# the offspring from the cross here:
bashcommand <- paste0("cat ", path, "RawMorphoData_PLGE2007/seg_morpho_3rep.dat", 
                      "| grep -v 'NA ' > ",
                      path, "RawMorphoData_PLGE2007/seg_morpho_3rep_cleaned.dat")
system(bashcommand)
phenoCross = read.table(paste0(path, "RawMorphoData_PLGE2007/seg_morpho_3rep_cleaned.dat"), 
                 as.is=T, header=T)
# does the order of the measured traits match between parent and offspring tables?
all(colnames(phenoCross[,-(1:2)]) == colnames(phenoBY)) & 
   all(colnames(phenoCross[,-(1:2)]) == colnames(phenoRM)) # yes.
# combine
pheno = cbind(
   sapply(unique(phenoCross$segname), function(seg){
      temp = phenoCross[phenoCross$segname == seg,-(1:2)]
      cor = cor(t(temp))
      if(any(cor[upper.tri(cor)] < 0.8))
         print(seg)
      colMeans(temp, na.rm = T)
   }), 
   "BY" = colMeans(phenoBY), 
   "RM" = colMeans(phenoRM))
# here is the information what the phenotypes stand for.
info = read.xlsx(paste0(path, "journal.pgen.0030031.st001_curated.XLS"),
                 sheetIndex=1, startRow=5, header=T,stringsAsFactors=F)
info$ID = gsub(x=info$ID, pattern="-",replacement=".")
#####


# subset to strains that are present both in our and their data #####
geno <- read.table("data/genotype_for_mapping.tsv",header=T)
markerPos = geno[,(1:3)]
genotype = geno[,-(1:3)]
ourStrains2theirStrains = sapply(colnames(genotype), function(x){
   if(x == "BY4716"){
      return("BY")
   } else if(x == "RM11.1a"){
      return("RM")
   } else{
      x = gsub(x=x, pattern=".", 
               replacement="_", fixed=T)
      x = gsub(x=x,pattern = "X", 
               replacement = "", fixed=T)
      return(x)
   }
})
theirStrains2ourStrains = setNames(names(ourStrains2theirStrains), ourStrains2theirStrains)
both = intersect(colnames(pheno), ourStrains2theirStrains)
pheno = pheno[,both]
save(pheno, info, file = paste0(path, "pheno.RData"))
# here is the mistake:

genotype = t(genotype[,theirStrains2ourStrains[both]])
save(genotype,markerPos, file = paste0(path, "Nogami_genotypes.RData"))
######


# mapping data preparation #####
sampleInfo = 1:nrow(genotype)
mappingData = preMap(genotype = genotype, phenotype=pheno, 
                     sampleInfo = sampleInfo, propVar=0.6)
save(mappingData, file = paste0(path, "Nogami_mappingData.RData"))
#####



# prepare data only for heritable traits ######
phenoReplicates = rbind(phenoCross[,-(1:2)], phenoBY, phenoRM)
phenoReplicatesStrains = c(phenoCross[,2],
                           rep("BY", nrow(phenoBY)), 
                           rep("RM", nrow(phenoRM)))
# the warning "boundary (singular) fit: see ?isSingular" can be ignored
H2 = sapply(1:ncol(phenoReplicates),function(i){
   df = data.frame(y = phenoReplicates[,i], strain = phenoReplicatesStrains)
   model = lmer(y~1+(1|strain), data=df)
   res = as.data.frame(VarCorr(model))
   if(!all(dim(res) == c(2,5)))
      return(NA)
   g = res$vcov[1]
   e = res$vcov[2]
   g/(g+e)
})
png("graph/Nogami2007/heritabilities.png")
hist(H2, las = 1)
dev.off()
heritableTraits = colnames(phenoReplicates)[H2 >= 0.5]
phenoHerit = cbind(
   sapply(unique(phenoCross$segname), function(seg){
      temp = phenoCross[phenoCross$segname == seg,heritableTraits]
      cor = cor(t(temp))
      if(any(cor[upper.tri(cor)] < 0.7))
         print(seg)
      colMeans(temp, na.rm = T)
   }), 
   "BY" = colMeans(phenoBY[,heritableTraits]), 
   "RM" = colMeans(phenoRM[,heritableTraits]))
phenoHerit = phenoHerit[,both]
save(phenoHerit, H2, heritableTraits, file = paste0(path, "phenoHerit.RData"))
sampleInfo = 1:nrow(genotype)
mappingDataHerit = preMap(genotype = genotype, phenotype=phenoHerit, 
                     sampleInfo = sampleInfo, propVar=0.6)
save(mappingDataHerit, file = paste0(path, "Nogami_mappingDataHerit.RData"))
#####

