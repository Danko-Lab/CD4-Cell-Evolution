##
## Fold-changes from U to PI -- checking that all species are on the same scale.

#setwd("/usr/projects/GROseq/NHP/annotations")

load("../annotations/fdr.RData")
source("../lib/densScatterplot.R")

## Index extent of changes.
isresp <- fdr_df$U2PIFDR_H < 0.05 #& fdr_df$U2PIFDDR_C < 0.05 & fdr_df$U2PIFDR_M < 0.05 ## Clearly responding in all three.
summary(isresp)
cor(data.frame(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp]))

pdf("foldchange.scatterplot.pdf")

 densScatterplot(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp])
 densScatterplot(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_M[isresp])

 ## Plots ... this is the supplementary figure.
 pairs(data.frame(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp]))

dev.off()


