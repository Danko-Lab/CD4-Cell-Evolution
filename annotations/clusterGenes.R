##
## ClusterGenes.R -- Clusters genes.

load("fdr.RData")

## Get changed genes...
PVAL <- 0.01
change <- fdr_df$fdr_min < PVAL | fdr_df$fdr_min_pi < PVAL

## Clustering...
source("../lib/CCVgen.R")
CreateCoorelationMatrix(MATRIX=fc_t[change,], NAMES=rownames(fdr_df$name[change]), k=NULL, METH="ward", MODE="heatmap")


## Simplest clustering...
png("tmp.png")
 heatmap(fc_t[change,])
dev.off()

## Cluster changed genes...
hc <- hclust(dist(fc_t[change,], method = "euclidean"), method="single")
hc1 <- as.dendrogram(hc)
ord.hc1 <- order.dendrogram(hc1)


pl <- levelplot((fc_t[change,])[ord.hc1,], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="")
print(pl)

## Write plots...
 pl <- levelplot((cc)[ord.hc1, ord.hc1], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="", #rev(cm.colors(100)),  # #c("white", "yellow", "blue") # c("#E9F231", "#B1EC2C", "#5DBDEF")
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

