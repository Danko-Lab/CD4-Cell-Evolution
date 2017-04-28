##
## Fold-changes from U to PI -- checking that all species are on the same scale.
source("../lib/densScatterplot.R")
source("../lib/runLimmaQuantile.R")

#setwd("/local/storage/projects/NHP/annotations")
load("../annotations/fdr.RData")

## Index extent of changes.
isresp <- fdr_df$U2PIFDR_H < PVAL & fdr_df$U2PIFDR_C < PVAL & fdr_df$U2PIFDR_M < PVAL ## Clearly responding in any of the three.
summary(isresp)

## Take anything that's changed in human, not in rhesus or chimp.
MAXNEG <- 0.25
MAXFC  <- 0.75 #0.5
isrhespec <- fdr_df$U2PIFDR_M < PVAL & (fdr_df$U2PIFDR_C > MAXNEG & fdr_df$U2PIFDR_H > MAXNEG) & abs(fdr_df$U2PIFC_M) > 1 & abs(fdr_df$U2PIFC_C) < MAXFC & abs(fdr_df$U2PIFC_H) < MAXFC
isrheloss <- fdr_df$U2PIFDR_M > MAXNEG & (fdr_df$U2PIFDR_C < PVAL & fdr_df$U2PIFDR_H < PVAL) & ((fdr_df$U2PIFC_H>1 & fdr_df$U2PIFC_C>1)|(fdr_df$U2PIFC_H< -1 & fdr_df$U2PIFC_C< -1)) & abs(fdr_df$U2PIFC_M) < MAXFC
summary(isrhespec)
summary(isrheloss)

pdf("RMacaque.Differences_In_Induction.pdf")
 ## Sanity checks...
 hist((fdr_df$U2PIFC_M-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])))
 hist((fdr_df$U2PIFC_M-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")]))[isrhespec])

 ## Plots.
 plot(fdr_df$U2PIFC_M[isrhespec], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isrhespec], pch=19); abline(0,1)
 densScatterplot(fdr_df$U2PIFC_M[isrhespec], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isrhespec], xlab="R. Macaque", ylab="Mean Human-Chimp")
 abline(0,1)

 plot(fdr_df$U2PIFC_M[isresp] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isresp], pch=19, xlab="Mean Human-Chimp", ylab="R. Macaque")
 points(fdr_df$U2PIFC_M[isrhespec] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isrhespec], pch=19, col="red")
 points(fdr_df$U2PIFC_M[isrheloss] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isrheloss], pch=19, col="blue")
 abline(h=0); abline(v=0)

 ## Check for bias in direction.
 require(vioplot)
 vioplot(rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isresp & fdr_df$U2PIFC_M<0], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isresp & fdr_df$U2PIFC_M>0]); abline(h=0)
 vioplot(rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isrhespec & fdr_df$U2PIFC_M<0], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])[isrhespec & fdr_df$U2PIFC_M>0]); abline(h=0) ## There's some residual bias in direction.

 ## Add labels.
dev.off()

fdr_df$score <- fdr_df$U2PIFC_M - rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")])

## Write out REs for Zhong.
istss <- (fdr_df$annot_type=="dREG_ENH" | fdr_df$annot_type=="dREG_INGENE" | fdr_df$annot_type=="dREG_TSS")
#write.table(fdr_df[isrhespec & istss,], "RMacaque.IndChange.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

write.table(fdr_df[((isrhespec & fdr_df$U2PIFC_M < 0)) & istss,c(1:6)], "RMacaque.gainActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[((isrheloss & rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")]) < 0)) & istss,c(1:6)], "RMacaque.looseActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

write.table(fdr_df[((isrhespec & fdr_df$U2PIFC_M > 0)) & istss,c(1:6)], "RMacaque.gainSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[((isrheloss & rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_H")]) > 0)) & istss,c(1:6)], "RMacaque.looseSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#write.table(fdr_df[((isrhespec & fdr_df$U2PIFC_M < 0)) & istss,], "RMacaque.gainActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
#write.table(fdr_df[((isrhespec & fdr_df$U2PIFC_M > 0)) & istss,], "RMacaque.gainSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#################
