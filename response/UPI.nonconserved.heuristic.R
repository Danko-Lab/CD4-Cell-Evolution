##
## Fold-changes from U to PI -- checking that all species are on the same scale.
source("../lib/densScatterplot.R")
source("../lib/runLimmaQuantile.R")

#setwd("/local/storage/projects/NHP/annotations")
load("../annotations/fdr.RData")

## Index extent of changes.
isresp <- fdr_df$U2PIFDR_H < PVAL & fdr_df$U2PIFDR_C < PVAL & fdr_df$U2PIFDR_M < PVAL ## Clearly responding in any of the three.
noresp <- fdr_df$U2PIFDR_H > 0.5 & fdr_df$U2PIFDR_C > 0.5 & fdr_df$U2PIFDR_M > 0.5 & abs(fdr_df$U2PIFC_H) < 0.25 & abs(fdr_df$U2PIFC_C) < 0.25 & abs(fdr_df$U2PIFC_M) < 0.25
summary(isresp)
summary(noresp)

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

istss <- fdr_df$annot_type=="dREG_ENH" | fdr_df$annot_type=="dREG_INGENE" | fdr_df$annot_type=="dREG_TSS"

## Take anything that's changed in human, not in rhesus or chimp.
MAXNEG <- 0.25#5
FCNEG  <- 0.5 #0.75
ishumspec <- fdr_df$U2PIFDR_H < PVAL & (fdr_df$U2PIFDR_C > MAXNEG & fdr_df$U2PIFDR_M > MAXNEG) & abs(fdr_df$U2PIFC_H) > 1 & abs(fdr_df$U2PIFC_C) < FCNEG & abs(fdr_df$U2PIFC_M) < FCNEG
ishumloss <- fdr_df$U2PIFDR_H > MAXNEG & (fdr_df$U2PIFDR_C < PVAL & fdr_df$U2PIFDR_M < PVAL) & ((fdr_df$U2PIFC_M>1 & fdr_df$U2PIFC_C>1)|(fdr_df$U2PIFC_M< -1 & fdr_df$U2PIFC_C< -1)) & abs(fdr_df$U2PIFC_H) < FCNEG
summary(ishumspec & istss)
summary(ishumloss & istss)

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
 points(fdr_df$U2PIFC_H[ishumloss] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumloss], pch=19, col="blue")
 abline(h=0); abline(v=0)

 ## Add labels.
dev.off()

fdr_df$score <- fdr_df$U2PIFC_H - rowMeans(fdr_df[,c("U2PIFC_M", "U2PIFC_C")])

## Write out REs for Zhong.
istss <- fdr_df$annot_type=="dREG_ENH" | fdr_df$annot_type=="dREG_INGENE" | fdr_df$annot_type=="dREG_TSS"
#write.table(fdr_df[ishumspec & istss,], "Human.IndChange.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[((ishumspec & fdr_df$U2PIFC_H < 0)) & istss,c(1:6)], "Human.gainActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[((ishumloss & rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")]) < 0)) & istss,c(1:6)], "Human.looseActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

write.table(fdr_df[((ishumspec & fdr_df$U2PIFC_H > 0)) & istss,c(1:6)], "Human.gainSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[((ishumloss & rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")]) > 0)) & istss,c(1:6)], "Human.looseSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Write out background sets.
fdr_df$score <- fdr_df$U2PIFC

write.table(fdr_df[isresp & fdr_df$U2PIFC_H < 0 & istss, c(1:6)], "all.conservedActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[isresp & fdr_df$U2PIFC_H > 0 & istss, c(1:6)], "all.conservedSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fdr_df[noresp & istss, c(1:6)], "all.noActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#write.table(fdr_df[(ishumspec & fdr_df$U2PIFC_H < 0) & istss,], "Human.gainActivation.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
#write.table(fdr_df[(ishumspec & fdr_df$U2PIFC_H > 0) & istss,], "Human.gainSuppression.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


## Write out all TSS that increase in human following PI...

write.table(fdr_df[(fdr_df$U2PIFDR_H < 0.05 & fdr_df$U2PIFC_H < 0 & (fdr_df$annot_type=="dREG_ENH" | fdr_df$annot_type=="dREG_INGENE" | fdr_df$annot_type=="dREG_TSS")),], "Human.Any.Change.tss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#################
## EXPERIMENTAL

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


