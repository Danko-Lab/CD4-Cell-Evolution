##
## Creates heatmaps on human dREG-HD sites.

require(bigWig)

dist <- 40000; step=50

## Load dREG-HD sites.
hd <- read.table("../dREG_HD/H-U_dREG_HD.bed"); hd <- hd[hd[,2]-dist > 0,]

## Load mark.
hMark <- load.bigWig("/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")

## Get a matrix of counts.
hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(hd, dist, dist), step=step)
hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(hd), byrow=TRUE)+1)
hmat <- hmat[order(rowSums(hmat), decreasing=TRUE),]

## 

## Write out a heatmap.
library(pheatmap)
bk <- seq(min(hmat), max(hmat), 0.01)
hmcols <- colorRampPalette(c("white","red"))(length(bk)-1)

png("H3K27ac.png", width=300, height = 800)     # width and height are in pixels

pheatmap( hmat, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE )

dev.off()

