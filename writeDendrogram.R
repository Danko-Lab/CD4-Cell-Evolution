## Get counts across species.
## 
require(cluster)

## Create a dendrogram.
## Using more efficient counting/ merging using bedops and bash scripts.						
source("readData.R")
#ca <- ca[grep("gc18", ca[,"annot_type"]),] ## Remove pause sites and TSS for this task...
ca <- ca[grep("protein_coding", ca[,"type"]),]
ca <- ca[(ca[,3]-ca[,2]) > 5000,]

yb.sig.pal <- function(n, scale=10) {
 ints<- c(0:(n-1))/(n-1)   ## Linear scale from 0:1 x N values.
 ints<- 1/(1+exp(scale*(0.5-ints)))## Transfer to sigmoidal scale.
 b<- min(ints)
 m<- 2*b/(n-1)
 ints<- ints+(m*c(0:(n-1)) -b)## Transform by linear function to fill colors out to maxes.
 
 ## Transfer to colorspace.
 # Yellow: 255, 255, 0
 # White:  255, 255, 255
 # Blue:   0, 0, 255
 YW <- ints[ints < 0.5] *2
 WB <- (ints[ints >= 0.5]-0.5) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
}

drawCor <- function(indx) {
	rpkm_df <- as.matrix(ca[,indx])/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
	for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

	cond <- Condition[indx]
	spec <- Species[indx]
	labs <- Labels[indx]

	cc <- cor(rpkm_df, method="spearman")
	#clu <- agnes(t(rpkm_df))

	pal2 <- c("#567C34", "#D14C29", "#567C34", "#69D270", "#7073C8", "#557571", "#CD9537", "#C0CB88")
	#pal2 <- c("#CE50CA", "#5D9D4C", "#D75631", "#5A8399", "#A18132", "#AD5687", "#7D71C7", "#AD4C4C")
	pal1 <- c("#B65BCB", "#C0513A", "#84CA54", "#92C2AF", "#4D4639", "#7B7EB5", "#BDA04D", "#B1517B")
    pal3 <- c("#E03CE9", "#17B92B", "#E6350D", "#6FD2F0", "#F9F77F", "#5B6C0C", "#68003D", "#310F08")
	
	## Print dendrogram and heatmap with latticeExtra.
	 library(latticeExtra)
	# hc1 <- agnes(1-cc, diss=TRUE, method="ward")
	# hc1 <- hclust(dist(t(rpkm_df), method = "canberra"))
	 hc1 <- hclust(dist(cc, method = "euclidean"),  method="centroid")
	 hc1 <- as.dendrogram(hc1)
	 ord.hc1 <- order.dendrogram(hc1)
	 hc2 <- reorder(hc1, cond[ord.hc1])
	 ord.hc2 <- order.dendrogram(hc2)
	 #region.colors <- trellis.par.get("superpose.polygon")$col

	 pl <- levelplot((cc)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="", #rev(cm.colors(100)),  # #c("white", "yellow", "blue") # c("#E9F231", "#B1EC2C", "#5DBDEF")
		 colorkey = list(space="left", labels=list(cex=1.5)), 
		 scales = list(x= list(rot=90, cex=1.5, labels=labs[ord.hc2]), y=list(draw=FALSE)), #scales = list(x = list(rot = 90)), 
		 legend = list(
			right = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "right", #lwd=2,
				 size = 7, size.add = 0.5, 
				 add = list(rect = list(col = "transparent", fill = pal3[c(7, 1, 8)][cond])),
				 type = "rectangle")), 
			top = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
				 size = 1, size.add = 0.5, 
				 add = list(rect = list(col = "transparent", fill = pal3[2:6][spec])),
				 type = "rectangle"))
				 ))
	 print(pl)
}

pdf("Agg.Hier.pdf")
 drawCor(indx.unt)
 drawCor(indx.good)
# drawCor(indx.all)
dev.off()



## Print correlation matrix w/ hexbin.
require(hexbin)
rpkm_df <- as.matrix(ca[,indx.all])/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

bin <- hexbin(x= log(rpkm_df[,2]+1,10), y=log(rpkm_df[,3]+1,10), 50)

indx <- rowSums(rpkm_df) > 0
 
png("scatterplots.png", width = 10000, height = 10000)
 hexplom(log(rpkm_df[indx,indx.unt-10]+min(rpkm_df[rpkm_df>0]),10), xbins=25, colramp = magent)
 #pairs(log(rpkm_df,10), pch=19)
dev.off()


