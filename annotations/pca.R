##
## Apply the RUV methods of Gagnon-Bartsch & Speed 
##
## REF: Gagnon-Bartsch and Speed, Biostatistics (2012), 13, 3, pp. 539-552.

source("readData.R")

##
## Do PCA ...

# Get RPKM
rpkm_df <- as.matrix(ca[,indx.good[c(2:9,11:17)]]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- log(1000*(rpkm_df[,i]+1)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"]))
#for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"])

pca <- prcomp(rpkm_df[,1:8], center=FALSE, scale=FALSE) ## UNT
#pca <- prcomp(rpkm_df, scale=FALSE, center=FALSE) ## ALL

cols <- c(rep("red",3), rep("green",2), rep("blue", 3), rep("dark red", 3), rep("dark green", 2), rep("dark blue", 2), "black", "black")
pch <- c(rep(19,8), rep(6,7), 9, 24)

summary(pca) # Prints variance summary for all principal components.
#plot(pca$rotation[,1], pca$rotation[,2], pch=19, col=cols)
pairs(pca$rotation[,1:5], col=cols, pch=pch)

## 
## Correlate each PC w/ the expected effects of cell compositio on GE.


##
## Get a set of genes expected to be invariant ...
require(bigWig)
cd4 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd4/rnaseq/cd4.epigenome.rnaseq.bw")
cd14p1 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapPlusRep1.bigWig")
cd14p2 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapPlusRep2.bigWig")
cd14m1 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapMinusRep1.bigWig")
cd14m2 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapMinusRep2.bigWig")

## Count reads ... 
CD4c  <- bedQuery.bigWig(ca[,1:3], cd4, gapValue=0) #bed.region.bpQuery.bigWig(cd4, ca[,1:6])
CD14c <- bedQuery.bigWig(ca[,1:3], cd14p1, gapValue=0)+ #+bed.region.bpQuery.bigWig(cd14p2, ca[,1:6])+bed.region.bpQuery.bigWig(cd14m1, ca[,1:6])+bed.region.bpQuery.bigWig(cd14m2, ca[,1:6])
         bedQuery.bigWig(ca[,1:3], cd14p2, gapValue=0)+
         bedQuery.bigWig(ca[,1:3], cd14m1, gapValue=0)+
         bedQuery.bigWig(ca[,1:3], cd14m2, gapValue=0)

ctrl_m <- (log((CD4c+1)/sum(CD4c))-log((CD14c+1)/sum(CD14c))) < -5
ctrl_t <- (log((CD4c+1)/sum(CD4c))-log((CD14c+1)/sum(CD14c))) > 5

monocyte <- (log(1000* (CD14c+1)/sum(CD14c)* 1000/(ca[,3]-ca[,2]))) ## Poor approximation for RPKM.
cd4rnah  <- (log(1000* (CD4c+1)/sum(CD4c)* 1000/(ca[,3]-ca[,2])))

tst <- ctrl_m | ctrl_t

##
## Do correlation.
for(i in 1:5) {
 print(cor.test(pca$x[tst,i], (monocyte-cd4rnah)[tst], method="spearman"))
 boxplot(pca$x[ctrl_m,i], pca$x[ctrl_t,i])
}

## Actually include monocyte.
pca <- prcomp(cbind(rpkm_df, monocyte, cd4rnah), scale=FALSE, center=FALSE) ## ALL

