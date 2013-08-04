## Get counts across species.
## 

require(cluster)

## Create a dendrogram.
## Using more efficient counting/ merging using bedops and bash scripts.						
ca <- read.table("countall.tsv")
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mapSize", "mgi", "Jurkat ", "Human 1", "Human 3", "Chimp 3", "R. Macaque 1", "R. Macaque 2")
ca <- ca[!is.na(ca$V9),]
count_df <- as.matrix(ca[,c(10:15)])/(ca[,9]) ## Normalize counts ...

indx <- (ca$chromEnd-ca$chromStart)>50000

cor(count_df, method="spearman")
cor(count_df[indx,], method="spearman")
#clu <- agnes(t(count_df))

clu <- agnes((1-cor(count_df[indx,], method="spearman")), diss=TRUE, method="ward")

pdf("Agg.Hier.pdf")
  plot(clu)
dev.off()

