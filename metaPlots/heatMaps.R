##
## Creates heatmaps on human dREG-HD sites.

require(bigWig)

dist <- 10000

## Load dREG-HD sites.
hd <- read.table("../dREG_HD/H-U_dREG_HD.bed"); hd <- hd[hd[,2]-dist > 0,]

writeHeatmap<- function(hMarkFile, name, hm_order= NULL, 
							cols= c("white","#00A63E"), dist= 10000, step=100,
							path="/local/storage/data/hg19/cd4/epiRoadmap_histone/") {
	## Load mark.
	hMark <- load.bigWig(paste(path, hMarkFile, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")

	## Get a matrix of counts.
	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(hd, dist, dist), step=step)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(hd), byrow=TRUE)+1)
	if(is.null(hm_order)) {
	  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
	}
	hmat <- hmat[hm_order,]

	## Average by rows of 10.
	

	## Write out a heatmap.
	library(pheatmap)
	bk <- seq(min(hmat), max(hmat), 0.01)
	hmcols <- colorRampPalette(cols)(length(bk)-1) # red

	png(paste(name,".png",sep=""), width=300, height = 800)     # width and height are in pixels

	pheatmap( hmat, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE )

	dev.off()
	
	return(hm_order)
}

ord <- writeHeatmap("H3K27ac.bw", "H3K27ac", cols=c("white","#00A63E"))
sup <- writeHeatmap("H3K4me1.bw", "H3K4me1", hm_order= ord, cols=c("white","#ff551c"))
sup <- writeHeatmap("H3K4me3.bw", "H3K4me3", hm_order= ord, cols=c("white","#fe0000"))

pth= "/local/storage/projects/NHP/AllData/All_Merge/"
sup <- writeHeatmap("H-U_plus.bw", "PROseq.plus", hm_order= ord, cols=c("white","#fe0000"), path=pth)
sup <- writeHeatmap("H-U_minus.bw", "PROseq.minus", hm_order= ord, cols=c("white","#0000fe"), path=pth)


