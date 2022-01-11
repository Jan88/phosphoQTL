PlotHclustWithColClasses <- function(mat, colClasses, proportion=8, ...) {
        
        clust <- hclust(dist(t(mat)))
        dendo <- as.dendrogram(clust)
        h     <- max(clust$height)
        
        plot(dendo, axes=FALSE, leaflab="none", ylim=c(-h/proportion,h), ...)
        
        if (nrow(colClasses) == 1) {
                colClasses <- matrix(colClasses[, order.dendrogram(dendo)], nrow=1)
        } else {	
                colClasses <- colClasses[, order.dendrogram(dendo)]
        }
        
        for(i in 1:nrow(colClasses)) {
                rect(
                        xleft   = 1:ncol(colClasses)-0.5,
                        ybottom = rep(i*(-h/(proportion-1)/nrow(colClasses)), ncol(colClasses)),
                        xright  = 1:ncol(colClasses)+0.5,
                        ytop    = rep((i-1)*(-h/(proportion-1)/nrow(colClasses)), ncol(colClasses)),
                        col     = colClasses[i,], border="white", lwd=0.1
                )
        }
        
        axis(1, 
             at       = 1:ncol(mat), 
             labels   = colnames(mat)[order.dendrogram(dendo)], 
             cex.axis = 0.05 + 1/log10(ncol(mat)), 
             tick     = F, 
             lwd      = 0, 
             las      = 2, 
             line     = -1
        )

}
