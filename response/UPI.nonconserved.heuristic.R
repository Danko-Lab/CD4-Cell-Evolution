##
## Fold-changes from U to PI -- checking that all species are on the same scale.
source("../lib/densScatterplot.R")
source("../lib/runLimmaQuantile.R")

#setwd("/local/storage/projects/NHP/annotations")
load("../annotations/fdr.RData")

## Index extent of changes.
isresp <- fdr_df$U2PIFDR_H < PVAL & fdr_df$U2PIFDR_C < PVAL & fdr_df$U2PIFDR_M < PVAL ## Clearly responding in any of the three.
summary(isresp)

pdf("foldchange.correlations.pdf")

 boxplot(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp])
 library(vioplot)
 vioplot(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp])

 cor(data.frame(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp]))

 ## Plots ... this is the supplementary figure.
 pairs(data.frame(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp]))

 cor.test(fdr_df$U2PIFC_H[isresp], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp])
 plot(fdr_df$U2PIFC_H[isresp], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp], pch=19); abline(0,1)
 densScatterplot(fdr_df$U2PIFC_H[isresp], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp], xlab="Human", ylab="Mean non-human primate"); abline(0,1)

dev.off()


## Take anything that's changed in human, not in rhesus or chimp.
MAXNEG <- 0.5
ishumspec <- fdr_df$U2PIFDR_H < PVAL & (fdr_df$U2PIFDR_C > MAXNEG & fdr_df$U2PIFDR_M > MAXNEG)
summary(ishumspec)

pdf("Human.Differences_In_Induction.pdf")
 ## Sanity checks...
 hist((fdr_df$U2PIFC_H-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])))
 hist((fdr_df$U2PIFC_H-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")]))[ishumspec])

 ## Plots.
 plot(fdr_df$U2PIFC_H[ishumspec], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumspec], pch=19); abline(0,1)
 densScatterplot(fdr_df$U2PIFC_H[ishumspec], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumspec], xlab="Human", ylab="Mean non-human primate")
 abline(0,1)

 plot(fdr_df$U2PIFC_H[isresp] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp], pch=19, xlab="Mean non-human primate", ylab="Human")
 points(fdr_df$U2PIFC_H[ishumspec] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumspec], pch=19, col="red")
 abline(h=0); abline(v=0)

 ## Add labels.
dev.off()

## Which genes are they?
fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 5 ,]
fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 5 & fdr_df$type == "protein_coding",]

fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 2 & fdr_df$type == "protein_coding" & abs(fdr_df$HumanFC_PI) > 1,]

## How frquently does the change in induction result in little-to-no differences in the PI condition?
pc<-NROW(fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 2 & fdr_df$type == "protein_coding",])
pc_dif_pi<-NROW(fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 2 & fdr_df$type == "protein_coding" & abs(fdr_df$HumanFC_PI) > 1,])

npc<-NROW(fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 2 & fdr_df$type != "protein_coding",])
npc_dif_pi<-NROW(fdr_df[ishumspec & abs(fdr_df$U2PIFC_H) > 2 & fdr_df$type != "protein_coding" & abs(fdr_df$HumanFC_PI) > 1,])

fisher.test(data.frame(c(pc, pc_dif_pi), c(npc, npc_dif_pi)))  ## No excess in protein-coding.


