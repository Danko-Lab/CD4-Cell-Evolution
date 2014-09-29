##
## Compare to GE data taken from a number of distinct blood cells to ensure no contamination.

setwd("/usr/projects/GROseq/NHP/annotations")

source("readData.R")

require(bigWig)
cd4 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd4/rnaseq/cd4.epigenome.rnaseq.bw")
cd14p1 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapPlusRep1.bigWig")
cd14p2 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapPlusRep2.bigWig")
cd14m1 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapMinusRep1.bigWig")
cd14m2 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapMinusRep2.bigWig")
cd4M <- load.bigWig("/usr/data/GROseq.parser/hg19/cd4/rnaseq/cd4.Memory.epigenome.rnaseq.bw")
cd8 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd4/rnaseq/cd8.Naive.epigenome.rnaseq.bw")
pbmc <- load.bigWig("/usr/data/GROseq.parser/hg19/cd4/rnaseq/pbmc.epigenome.rnaseq.bw")


## Count reads ... 
CD4c  <- bedQuery.bigWig(ca[,1:3], cd4, gapValue=0) #bed.region.bpQuery.bigWig(cd4, ca[,1:6])
CD14c <- bedQuery.bigWig(ca[,1:3], cd14p1, gapValue=0)+ #+bed.region.bpQuery.bigWig(cd14p2, ca[,1:6])+bed.region.bpQuery.bigWig(cd14m1, ca[,1:6])+bed.region.bpQuery.bigWig(cd14m2, ca[,1:6])
         bedQuery.bigWig(ca[,1:3], cd14p2, gapValue=0)+
         bedQuery.bigWig(ca[,1:3], cd14m1, gapValue=0)+
         bedQuery.bigWig(ca[,1:3], cd14m2, gapValue=0)
CD4Mc <- bedQuery.bigWig(ca[,1:3], cd4M, gapValue=0)
CD8c  <- bedQuery.bigWig(ca[,1:3], cd8, gapValue=0)
PBMCc <- bedQuery.bigWig(ca[,1:3], pbmc, gapValue=0)

## Expression ratio chosen to be ~500 genes (ballpark).
ctrl_mon <- (log((CD4c+1)/sum(CD4c))-log((CD14c+1)/sum(CD14c))) < -5
ctrl_mem <- (log((CD4c+1)/sum(CD4c))-log((CD4Mc+1)/sum(CD4Mc))) < -4
ctrl_cd8 <- (log((CD4c+1)/sum(CD4c))-log((CD8c+1)/sum(CD8c))) < -3
ctrl_pbmc <- (log((CD4c+1)/sum(CD4c))-log((PBMCc+1)/sum(PBMCc))) < -5

## Correlate each species expression w/ diff cell types.
rpkm_df <- as.matrix(ca[,indx.good[c(2:9,11:17)]]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- log(1000*(rpkm_df[,i]+1)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"]))
#for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"])

Hu <- rowMeans(rpkm_df[,1:3])
Ch <- rowMeans(rpkm_df[,4:5])
Ma <- rowMeans(rpkm_df[,6:8])

## Boxplots ...
boxplot(Hu[ctrl_mon], Ch[ctrl_mon], Ma[ctrl_mon])
boxplot(Hu[ctrl_mem], Ch[ctrl_mem], Ma[ctrl_mem])
boxplot(Hu[ctrl_cd8], Ch[ctrl_cd8], Ma[ctrl_cd8])
boxplot(Hu[ctrl_pbmc], Ch[ctrl_pbmc], Ma[ctrl_pbmc])



 
