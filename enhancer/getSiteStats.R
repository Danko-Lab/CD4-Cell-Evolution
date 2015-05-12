## This analysis focuses on lineage-specific changes that are currently active in hg19.
load("../annotations/fdr.RData")

tss_aln <- fdr_df[grepl("dREG", ca$annot_type),]
tss <- read.table("counttss.tsv.tmp")
tss <- data.frame(tss[,1:20], tss_aln[match(tss$V4, tss_aln$name),31:42])

## Alignable fraction (V20) denotes a gap in either species.  Make sure gaps are in both.

## Classify as 'promoter'/ 'enhancer'
stab <- rowMax(tss[,17:18])
dist <- tss[,13]
class <- rep("tss", NROW(tss)) ## tss is then unclassified as a promoter or enhancer
class[stab < 0.1 & dist < 500]  <- "Prox_Stab" ## Clearly protein coding promoter
class[stab > 0.1  & dist > 10000] <- "Dist_UnSt" ## Clearly distal enhancer
class[stab < 0.001  & dist > 10000] <- "Dist_Stab" ## Clearly stable, but distal
summary(as.factor(class))
tss$V19 <- class

## Count types of elements ...
total <- summary(as.factor(tss$V19))
unmap <- summary(as.factor(tss$V19[tss$V20 == 0]))
chang <- summary(as.factor(tss$V19[tss$V20 > 0 & tss$fdr_min < 0.01]))

total
unmap
chang 

unmap/total
chang/total

## HOWEVER! Note that these numbers seem to include 1:many orthologies, not just gaps.  These really should be removed.
