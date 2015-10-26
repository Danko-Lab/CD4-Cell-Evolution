## This script identifies branch-specific changes in RE activities.
##
load("../annotations/fdr.RData")
source("../lib/normalizeSubsample.R")

require(boot)

tss_aln <- fdr_df[grepl("dREG", ca$annot_type),]
tss <- read.table("tss.tsv")
tss <- data.frame(tss, tss_aln[match(tss$V4, tss_aln$name),c(9,33:50)])

## Alignable fraction (V20) denotes a gap in either species.  Make sure gaps are in both.

## Classify as 'promoter'/ 'enhancer'
stab <- rowMax(tss[,17:18])
dist <- tss[,13]
class <- rep("tss", NROW(tss)) ## tss is then unclassified as a promoter or enhancer
class[stab < 0.1 & dist < 500]  <- "Prox_Stab" ## Clearly protein coding promoter
class[stab > 0.1  & dist > 10000] <- "Dist_UnSt" ## Clearly distal enhancer
class[stab < 0.1  & dist > 125000] <- "Dist_Stab" ## Clearly stable, but distal
summary(as.factor(class))
tss$V5 <- as.factor(class)

## Change unscored to 0
for(i in 7:12) { tss[is.na(tss[,i]),i] <- 0 }

## Change in basal T-cells.
              ##  1:1 ortholog,  mappable,             complete gain/ loss,                            gain/ loss in magnitude.
indx_hg19_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V7 > 0.7 & tss$V8 < 0.1 & tss$V9 < 0.1) | (tss$HumanFDR < 0.05 & tss$HumanFC > 0))
indx_hg19_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V7 < 0.1 & tss$V8 > 0.7 & tss$V9 > 0.7) | (tss$HumanFDR < 0.05 & tss$HumanFC < 0))

write.table(tss[indx_hg19_gain | indx_hg19_loss,], "hg19.gain.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

              ##     1:1 ortholog,  mappable,             complete gain/ loss,                            gain/ loss in magnitude.
indx_rheMac3_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V9 > 0.7 & tss$V8 < 0.1 & tss$V7 < 0.1) | (tss$MacaqueFDR < 0.05 & tss$MacaqueFC > 0))
indx_rheMac3_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V9 < 0.1 & tss$V8 > 0.7 & tss$V7 > 0.7) | (tss$MacaqueFDR < 0.05 & tss$MacaqueFC < 0))

write.table(tss[indx_rheMac3_gain | indx_rheMac3_loss,], "rheMac3.gain.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

              ##     1:1 ortholog,  mappable,             complete gain/ loss,                            gain/ loss in magnitude.
indx_panTro4_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V8 > 0.7 & tss$V9 < 0.1 & tss$V7 < 0.1) | (tss$ChimpFDR < 0.05 & tss$ChimpFC > 0))
indx_panTro4_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V8 < 0.1 & tss$V9 > 0.7 & tss$V7 > 0.7) | (tss$ChimpFDR < 0.05 & tss$ChimpFC < 0))

write.table(tss[indx_panTro4_gain | indx_panTro4_loss,], "panTro4.gain.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



## Data playtime!

## 
indx_hg19_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V7 > 0.7 & tss$V8 < 0.1 & tss$V9 < 0.1))
indx_hg19_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V7 < 0.1 & tss$V8 > 0.7 & tss$V9 > 0.5))

sum(indx_hg19_gain) ## Differences in gain/ loss rates.  I'd guess due to incomplete power in rm3 dREG sites?!
sum(indx_hg19_loss)


indx_rm3_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$MacaqueFDR < 0.05 & tss$MacaqueFC > 0))
indx_rm3_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$MacaqueFDR < 0.05 & tss$MacaqueFC < 0))


## Change in inducibility 
              ##  1:1 ortholog,  mappable,             gain/ loss in inducability.
#indx_hg19_gaini<- tss$V20 == 0 & !is.na(tss$mapSize) & 
