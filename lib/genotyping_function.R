#' function to read VCF files
#' @author Mathieu Cleément-Ziza
#' 
#' @param genotype_vcf_file Path to the VCF file 
#' @param headerMax maximun number of line to read before the header
#' 
#' @return a dataFrame
#' 
read.vcf <- function(genotype_vcf_file, headerMax=200) {
  header <- grep("#CHROM",readLines(genotype_vcf_file,n=headerMax),value=T)
  if (length(header)==0) {
    stop("header not found")
  }
  header <- substring(header, 2) # remove the #
  vcf <- read.table(file=genotype_vcf_file,sep="\t",header=F,blank.lines.skip=TRUE,comment.char='#',stringsAsFactors=F)
  colnames(vcf) <- unlist(strsplit(header,"\t"))
  return(vcf)
}  



##############################################################################
#' @param chr chromosome index
#' @param pos genomic postion on the chromosome
#' @param chr.size vector with the size of each chromosome
#' @return return the absolute genomic postion (when putting Chr one after the other
##############################################################################
getGenomicPosition<-function(chr,pos,chromosome.size){
  if(is.na(chr)||is.na(pos)){return(NA)}
  if(is.character(chr)) {chr<-which(names(chromosome.size)==chr)}
  chromosome.size<-c(0,chromosome.size)
  sum(chromosome.size[1:chr],pos)
}



###########################################################
#' A function to plot genotype of segregants from a cross 
#' 
#' @param geno gentype data.frame. The vcf should be formated with the first 4 columns as  'chr','pos','ref','alternate'. The genotype of the strains should be in the remaining colums which name should be the strains. genotpye should be numbers or factors
#' @param chromosome names of the chromosome to consider. If NULL all the chr in chromosome length are plotted
#' @param chr_length numreric vector with the length of the chromosomes names of the vector, should be names of the chrs.
#' @param p1 parent1 genotype. default 0
#' @paran p2 parent1 genotype. default 1
#' @param title char. title of the plot. default, empty
#' @param col: colors of the genotype levels
#' @param colNA color for NA
#' @labels labels for the raw. if NULL. ames form the columns
plot_genotype <- function(geno,chromosome=NULL,chr_length,p1=0,p2=1,title='',col=c('#377EB8','#E41A1C'),colNA="black",labels=NULL,axis=TRUE,cex.y.axis=1,...){
  if(is.null(labels)){legend.strains <-colnames(geno)[5:ncol(geno)] } else {legend.strains<-labels}
  if(is.null(chromosome)){chromosome <- names(chr_length)}
  chr_length <- chr_length[chromosome]
  
  # coordinate of chr limimits for plotting 
  chr.lim <-sapply(1:length(chr_length),function(chr){
    c(getGenomicPosition(chr,1,chromosome.size=chr_length),getGenomicPosition(chr,chr_length[chr],chromosome.size=chr_length))
  })
  chr.lim <-as.data.frame( t(chr.lim))
  colnames(chr.lim)<-c('start','end')
  rownames(chr.lim)<-names(chr_length)
  
  geno[,1] <- as.character(geno[,1])
  geno <- geno[geno[,1]%in%names(chr_length),]
  
  g<-as.data.frame(geno[,5:ncol(geno)])
  
  plot.new()
  plot.window(xlim=c(-20,sum(chr_length)+20), 
              ylim=c(1-0.45,(ncol(g)+0.45)),
              xlab="chr coordinate",xaxs='i')
  
  
  if(!is.null(legend.strains)){
    axis(2,at=1:(ncol(geno)-4),labels=legend.strains,lwd=0,las=1,cex.axis=cex.y.axis, lines=-1)
  }
  if(axis){axis(1,at=seq(0,sum(chr_length),1000000),
       labels=paste(seq(0,sum(chr_length),1000000)/1000000,'Mb'),...)
  }
  x <- mapply(getGenomicPosition,
              MoreArgs=list(chromosome.size=chr_length),
              geno[,1], geno[,2])
  sapply(1:ncol(g),function(i){
    sel <-  is.na(g[,i])
    if(sum(sel,na.rm=T)>0){segments(x[sel],rep(i-0.45,sum(sel,na.rm=T)),x[sel],rep(i+0.45,sum(sel,na.rm=T)),col=colNA,xaxt='n',yaxt='n')}
    sel <-   g[,i]==p1
    if(sum(sel,na.rm=T)>0){segments(x[sel],rep(i-0.45,sum(sel,na.rm=T)),x[sel],rep(i+0.45,sum(sel,na.rm=T)),col=col[1],xaxt='n',yaxt='n')}
    sel <-   g[,i]==p2
    if(sum(sel,na.rm=T)>0){segments(x[sel],rep(i-0.45,sum(sel,na.rm=T)),x[sel],rep(i+0.45,sum(sel,na.rm=T)),col=col[2],xaxt='n',yaxt='n')}
  })
  segments(chr.lim[2:nrow(chr.lim),1],par("usr")[3],chr.lim[2:nrow(chr.lim),1],par("usr")[4],col='black',lwd=4)
}


export_vcf<-function(vcf,file){
  cat('##fileformat=VCFv4.0', 
      '##INFO=<ID=AB,Number=1,Type=Float,Description="Allele Balance for hets (ref/(ref+alt))">',
      '#############################################',
      '#####GENERTATED FOR INTERNAL USE ONLY #######',
      '#############################################',
      paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO',sep='\t'),
      sep='\n',
      file=file)
  out<-data.frame(vcf[,c(1,2)],'.',vcf[,c(3,4)],999,'PASS','AB=0.1')
  write.table(out,file=file,append=TRUE,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)   
}





find.suspect <- function(genotype, chromosome){
  id <- 1:length(genotype)
  out <-sapply(unique(chromosome),function(chr){ #for each chromsome
    id <- id[chromosome==chr]
    genotype <- genotype[chromosome==chr]
    
    #select for know genotype
    sel <- !is.na(genotype) 
    genotype <- genotype[sel] ; id<- id[sel] ;
    
    #sleect the suspicious
    if(length(genotype)<3) {return(NA)}
    before<-genotype[1:(length(genotype)-2)]
    actual <-genotype[2:(length(genotype)-1)]
    after <- genotype[3:length(genotype)]
    id <- id[2:(length(genotype)-1)]
    return(id[actual!=before & actual !=after])
  })
  out <- unlist(out)
  return(out[!is.na(out)])
}



fillTheGap <- function(geno, strain, max_distance_flanking_markers=10000){
  print(strain)
  chromosome <- unique(geno$CHR)
  out <-sapply(chromosome,function(chr){ #for each chromsome
    genotype <- geno[geno$CHR==chr,strain]
    postion<- geno[geno$CHR==chr,'POS']
    #select for know genotype
    sel <- which(!is.na(genotype))
    pos1<- sel[1:(length(sel)-1)]
    pos2<-sel[2:length(sel)]
    for(i in 1:length(pos1)){
      g1<-genotype[pos1[i]]
      g2<-genotype[pos2[i]]
      if(g1==g2 && postion[pos2[i]]-postion[pos1[i]]< max_distance_flanking_markers){
        genotype[pos1[i]:pos2[i]]<-g1
      }
    }
    if(-postion[1]+postion[pos1[1]] < max_distance_flanking_markers){
      genotype[1:pos1[1]]<-genotype[pos1[1]]
    }
    if(postion[length(genotype)]-postion[pos2[length(pos2)]] < max_distance_flanking_markers){
      genotype[pos2[length(pos2)]:length(genotype)]<-genotype[pos2[length(pos2)]]
    }
    #genotype[!genotype%in%c(genotype.p1,genotype.p2)]<- -1
    return(genotype)
  })
  return(unlist(out))
}



