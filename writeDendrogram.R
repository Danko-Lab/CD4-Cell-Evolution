## Get countscross species.
## 
require(cluster)

## Create a dendrogram.
## Using more efficient counting/ merging using bedops and bash scripts.						
ca <- read.table("countall.tsv")
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", "Jurkat ", "Human 1", "Human 3", "Chimp 3", "R. Macaque 1", "R. Macaque 2")
ca <- ca[!is.na(ca[,10]),]
rpkm_df <- as.matrix(ca[,c(10:15)])/(ca[,9]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

cor(rpkm_df, method="spearman")
cor(rpkm_df, method="spearman")
#clu <- agnes(t(rpkm_df))

clu <- agnes((1-cor(rpkm_df, method="spearman")), diss=TRUE, method="ward")

pdf("Agg.Hier.pdf")
 plot(clu)
dev.off()

require(hexbin)
bin <- hexbin(x= log(rpkm_df[,2]+1,10), y=log(rpkm_df[,3]+1,10), 50)

indx <- rowSums(rpkm_df) > 0
 
png("scatterplots.png", width = 1000, height = 1000)
 hexplom(log(rpkm_df[indx,]+min(rpkm_df[rpkm_df>0]),10), xbins=25, colramp = magent)
 #pairs(log(rpkm_df,10), pch=19)
dev.off()

 across species.
## 
require(cluster)

## Create a dendrogram.
## Using more efficient counting/ merging using bedops and bash scripts.						
ca <- read.table("countall.tsv")
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", "Jurkat ", "Human 1", "Human 3", "Chimp 3", "R. Macaque 1", "R. Macaque 2")
ca <- ca[!is.na(ca[,10]),]
count_df <- as.matrix(ca[,c(10:15)])/(ca[,9]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(count_df)) count_df[,i] <- count_df[,i]/sum(count_df[,i])

indx <- (ca$mapSize)>50000

cor(count_df, method="spearman")
cor(count_df[indx,], method="spearman")
#clu <- agnes(t(count_df))

clu <- agnes((1-cor(count_df[indx,], method="spearman")), diss=TRUE, method="ward")

pdf("Agg.Hier.pdf")
 plot(clu)
dev.off()

require(hexbin)
bin <- hexbin(x= log(count_df[,2]+1,10), y=log(count_df[,3]+1,10), 50)

indx <- rowSums(count_df) > 0
 
png("scatterplots.png", width = 1000, height = 1000)
 hexplom(log(count_df[indx,]+min(count_df[count_df>0]),10), xbins=25, colramp = magent)
 #pairs(log(count_df,10), pch=19)
dev.off()

