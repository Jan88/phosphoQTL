   read.counts<-function(samples,gene_names,counts_path,verbose=FALSE, extension='.counts'){
   	a<-lapply(samples,function(s){
         if(verbose){cat(paste(counts_path,'/',s,extension,sep=''),'\n')}
         if(file.info(paste(counts_path,'/',s,extension ,sep=''))$size==0){out<-data.frame(NA)}else{
   		   out<-read.table(paste(counts_path,'/',s, extension, sep=''),stringsAsFactors=F,row.names=2)
         }
   		out<-out[gene_names,]
   		out[is.na(out)]<-0
   		return(out)
   	})
   	a<-do.call('cbind',a)
   	colnames(a)<-samples
   	rownames(a)<-gene_names
   	return(a)
   }

