##
## Estimate cell population assuming mixed sample a lienar combination of RNA-seq data.
##

source("readData.R")
PC <- 0 ## Pseudocount

rpkm_df <- as.matrix(ca[,indx.good[c(2:9,11:17)]]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*(rpkm_df[,i]+PC)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"])

##
## Get reference files ...
source("../dataSanityChecks/compareBloodCells.R")

## Normalize
CD4c  <- 1000*(CD4c+PC)/sum(CD4c) *1000/(ca[,3]-ca[,2])
CD4Mc <- 1000*(CD4Mc+PC)/sum(CD4Mc) *1000/(ca[,3]-ca[,2])
CD14c <- 1000*(CD14c+PC)/sum(CD14c) *1000/(ca[,3]-ca[,2])
ct_df <- cbind(CD4c, CD4Mc, CD14c)

## Non-linear model.
df1 <- data.frame(H1= rpkm_df[,5], CD4= CD4c, CD4M= CD4Mc, CD14= CD14c)
#nlfit <- nls(H1~0+ I(A*CD4+C*CD14), data=df1, start=data.frame(A=0.25, C=0.25)) ## H1~I(A*CD4+B*CD4M+C*CD14)
#nlfit
glm(H1~0+CD4+CD4M+CD14, data=df1, family= quasipoisson)


#####################

## Does it work?
Hu <- rowMeans(rpkm_df[,1:3]); Hua <- rowMeans(adj_df[,1:3])
Ch <- rowMeans(rpkm_df[,4:5]); Cha <- rowMeans(adj_df[,4:5])
Ma <- rowMeans(rpkm_df[,6:8]); Maa <- rowMeans(adj_df[,6:8])

## Boxplots ...
library(vioplot)
vioplot(Hu[ctrl_mon], Ch[ctrl_mon], Ma[ctrl_mon], Hua[ctrl_mon], Cha[ctrl_mon], Maa[ctrl_mon])
vioplot(Hu[ctrl_mem], Ch[ctrl_mem], Ma[ctrl_mem], Hua[ctrl_mem], Cha[ctrl_mem], Maa[ctrl_mem])
vioplot(Hu[ctrl_cd8], Ch[ctrl_cd8], Ma[ctrl_cd8], Hua[ctrl_cd8], Cha[ctrl_cd8], Maa[ctrl_cd8])
vioplot(Hu[ctrl_pbmc], Ch[ctrl_pbmc], Ma[ctrl_pbmc], Hua[ctrl_pbmc], Cha[ctrl_pbmc], Maa[ctrl_pbmc])

