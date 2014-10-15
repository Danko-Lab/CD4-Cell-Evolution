##
## Sanity checks for gene expression changes ... Look at the reported raw levels of gene expression.

source("readData.R")

rpkm_df <- as.matrix(ca[,indx.good[c(2:9,11:17)]]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- log(1000*(rpkm_df[,i]+0.25)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"]))

source("../lib/circplot.R")
snU <- c(rep("H",3), rep("C", 2), rep("M",3))

## Good ...
cd.circplot(rpkm_df[ca$name == "chr12_15771150_15943000", 1:8], snU)

## Not sure...
cd.circplot(rpkm_df[ca$name == "chrY_1523573_1606700", 1:8], snU)



