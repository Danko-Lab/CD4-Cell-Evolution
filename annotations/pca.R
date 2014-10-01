##
## Apply the RUV methods of Gagnon-Bartsch & Speed 
##
## REF: Gagnon-Bartsch and Speed, Biostatistics (2012), 13, 3, pp. 539-552.

source("readData.R")

##
## Do PCA ...
rpkm_df <- log(rpkm_df[,c(2:9,11:17)]+1e-7)

#pkm_df <- as.matrix(ca[,indx.good[c(2:9,11:17)]]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
#or(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- log(1000*(rpkm_df[,i]+0.25)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"]))


pca <- prcomp(rpkm_df[,1:8], center=FALSE, scale=FALSE) ## UNT
#pca <- prcomp(rpkm_df[,9:15], center=FALSE, scale=FALSE) ## PI
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
source("../dataSanityChecks/compareBloodCells.R")

#ctrl_m <- (log((CD4c+1)/sum(CD4c))-log((CD14c+1)/sum(CD14c))) < -5
#ctrl_t <- (log((CD4c+1)/sum(CD4c))-log((CD14c+1)/sum(CD14c))) > 5
#
#monocyte <- (log(1000* (CD14c+1)/sum(CD14c)* 1000/(ca[,3]-ca[,2]))) ## Poor approximation for RPKM.
#cd4rnah  <- (log(1000* (CD4c+1)/sum(CD4c)* 1000/(ca[,3]-ca[,2])))
#
#tst <- ctrl_m | ctrl_t

##
## Do correlation.
#for(i in 1:5) {
# print(cor.test(pca$x[tst,i], (monocyte-cd4rnah)[tst], method="spearman"))
# boxplot(pca$x[ctrl_m,i], pca$x[ctrl_t,i])
#}


##
## ADJUST using PCA.
i <- 3                     ## Factor out PCA3 (unt. only).  PCA3&&4 (all).
n <- NROW(pca$rotation)      ## Dimensionality
adj_df <- ((pca$x[,c(1:(i-1),(i+1):n)] %*% t(pca$rotation[,c(1:(i-1),(i+1):n)]))) 
#adj_df <- ((pca$x[,c(1:2,5:n)] %*% t(pca$rotation[,c(1:2,5:n)]))) ## 3 && 5 removed (both conditions...).


## Does it work?
Hu <- rowMeans(rpkm_df[,1:3]); Hua <- rowMeans(adj_df[,1:3])
Ch <- rowMeans(rpkm_df[,4:5]); Cha <- rowMeans(adj_df[,4:5])
Ma <- rowMeans(rpkm_df[,6:8]); Maa <- rowMeans(adj_df[,6:8])

## Boxplots ...
boxplot(Hu[ctrl_mon], Ch[ctrl_mon], Ma[ctrl_mon], Hua[ctrl_mon], Cha[ctrl_mon], Maa[ctrl_mon])
boxplot(Hu[ctrl_mem], Ch[ctrl_mem], Ma[ctrl_mem], Hua[ctrl_mem], Cha[ctrl_mem], Maa[ctrl_mem])
boxplot(Hu[ctrl_cd8], Ch[ctrl_cd8], Ma[ctrl_cd8], Hua[ctrl_cd8], Cha[ctrl_cd8], Maa[ctrl_cd8])
boxplot(Hu[ctrl_pbmc], Ch[ctrl_pbmc], Ma[ctrl_pbmc], Hua[ctrl_pbmc], Cha[ctrl_pbmc], Maa[ctrl_pbmc])


