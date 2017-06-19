## This script identifies branch-specific changes in RE activities.
##
load("../annotations/fdr.RData")
source("../lib/normalizeSubsample.R")

highth <- 0.30
lowth  <- 0.05

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
indx_hg19_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V7 > highth & tss$V8 < lowth & tss$V9 < lowth) | (tss$HumanFDR < PVAL & tss$HumanFC > 0))
indx_hg19_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V7 < lowth & tss$V8 > highth & tss$V9 > highth) | (tss$HumanFDR < PVAL & tss$HumanFC < 0))

sum(indx_hg19_gain)
sum(indx_hg19_loss)

write.table(tss[indx_hg19_gain | indx_hg19_loss,], "hg19.gain.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[indx_hg19_gain,], "hg19.gain.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[indx_hg19_loss,], "hg19.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

write.table(tss[indx_hg19_gain | indx_hg19_loss,1:3], "hg19.gain.loss.insight.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[(indx_hg19_gain | indx_hg19_loss) & abs(tss$HumanFC) > 5^(1/2) & tss$HumanFDR < 0.01,1:3], "hg19.gl.fold-GT5.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[(indx_hg19_gain | indx_hg19_loss) & tss$HumanFDR < 0.01 & tss$HumanFDR_PI < 0.01,1:3], "hg19.gl-UPI-HC.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


              ##     1:1 ortholog,  mappable,             complete gain/ loss,                            gain/ loss in magnitude.
indx_rheMac3_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V9 > highth & tss$V8 < lowth & tss$V7 < lowth) | (tss$MacaqueFDR < PVAL & tss$MacaqueFC > 0))
indx_rheMac3_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V9 < lowth & tss$V8 > highth & tss$V7 > highth) | (tss$MacaqueFDR < PVAL & tss$MacaqueFC < 0))

sum(indx_rheMac3_gain)
sum(indx_rheMac3_loss)

write.table(tss[indx_rheMac3_gain | indx_rheMac3_loss,], "rheMac3.gain.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[indx_rheMac3_gain,], "rheMac3.gain.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[indx_rheMac3_loss,], "rheMac3.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

              ##     1:1 ortholog,  mappable,             complete gain/ loss,                            gain/ loss in magnitude.
indx_panTro4_gain <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V8 > highth & tss$V9 < lowth & tss$V7 < lowth) | (tss$ChimpFDR < PVAL & tss$ChimpFC > 0))
indx_panTro4_loss <- tss$V20 == 0 & !is.na(tss$mapSize) & ((tss$V8 < lowth & tss$V9 > highth & tss$V7 > highth) | (tss$ChimpFDR < PVAL & tss$ChimpFC < 0))

sum(indx_panTro4_gain)
sum(indx_panTro4_loss)

write.table(tss[indx_panTro4_gain | indx_panTro4_loss,], "panTro4.gain.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[indx_panTro4_gain,], "panTro4.gain.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(tss[indx_panTro4_loss,], "panTro4.loss.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Conserved in all species.
indx <- tss$V20 == 0 & !is.na(tss$mapSize) & (tss$V7 > highth & tss$V8 > highth & tss$V9 > highth) & (tss$HumanFDR > 0.5 & tss$ChimpFDR > 0.5 & tss$MacaqueFDR > 0.5) & (abs(tss$HumanFC) < 0.25 & abs(tss$ChimpFC) < 0.25 & abs(tss$MacaqueFC) < 0.25)
sum(indx) 

write.table(tss[indx,], "all.conserved.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


## Validation in humans.
require(bigWig)
source("../lib/avg.metaprofile.R")
random_sites <- read.table("random-sites.bed.gz")
makePlot <- function(bed, mark, bwpath= "/local/storage/data/hg19/cd4/epiRoadmap_histone/", halfWindow= 5000, step= 25, ...) {
  bw <- load.bigWig(paste(bwpath, mark, ".bw", sep=""))
  mp <- avg.metaprofile.bigWig(center.bed(bed[,1:3], halfWindow, halfWindow), bw, step=step, ...)
  plot(mp)

  bed.region.bpQuery.bigWig(bw, bed[,1:3]) * 1000/(bed[,3]-bed[,2])
}

pdf("dREG-Changes.pdf")

a <- makePlot(tss[indx,], "H3K27ac", name="H3K27ac")
b <- makePlot(tss[indx_hg19_gain & (tss$V8 > lowth | tss$V9 > lowth),], "H3K27ac", name="H3K27ac gain")
c <- makePlot(tss[indx_hg19_loss & tss$V7 > lowth,], "H3K27ac", name="H3K27ac loss") ## These include sites that are decreases.
cc<- makePlot(tss[indx_hg19_gain & tss$V8 < lowth & tss$V9 < lowth,], "H3K27ac", name="H3K27ac complete gain") ## tss$HumanFDR > PVAL
d <- makePlot(tss[indx_hg19_loss & tss$V7 < lowth,], "H3K27ac", name="H3K27ac complete loss") ## tss$HumanFDR > PVAL
e <- makePlot(random_sites, "H3K27ac", name="H3K27ac random")
boxplot(list(conserved= a, gain= b, loss= c, complete.gain= cc, complete.loss= d, random= e), ylab="Reads per kilobase", main="H3K27ac", outline=FALSE)

a <- makePlot(tss[indx,], "H3K27me3", name="H3K27me3")
b <- makePlot(tss[indx_hg19_gain & (tss$V8 > lowth | tss$V9 > lowth),], "H3K27me3", name="H3K27me3 gain")
c <- makePlot(tss[indx_hg19_loss & tss$V7 > lowth,], "H3K27me3", name="H3K27me3 loss")
cc<- makePlot(tss[indx_hg19_gain & tss$V8 < lowth & tss$V9 < lowth,], "H3K27me3", name="H3K27me3 complete gain")
d <- makePlot(tss[indx_hg19_loss & tss$V7 < lowth,], "H3K27me3", name="H3K27me3 complete loss")
e <- makePlot(random_sites, "H3K27me3", name="H3K27me3 random")
boxplot(list(conserved= a, gain= b, loss= c, complete.gain= cc, complete.loss= d, random= e), ylab="Reads per kilobase", main="H3K27me3", outline=FALSE)

a <- makePlot(tss[indx,], "H3K4me3", name="H3K4me3")
b <- makePlot(tss[indx_hg19_gain & (tss$V8 > lowth | tss$V9 > lowth),], "H3K4me3", name="H3K4me3 gain")
c <- makePlot(tss[indx_hg19_loss & tss$V7 > lowth,], "H3K4me3", name="H3K4me3 loss")
cc<- makePlot(tss[indx_hg19_gain & tss$V8 < lowth & tss$V9 < lowth,], "H3K4me3", name="H3K4me3 complete gain")
d <- makePlot(tss[indx_hg19_loss & tss$V7 < lowth,], "H3K4me3", name="H3K4me3 complete loss")
e <- makePlot(random_sites, "H3K4me3", name="H3K4me3 random")
boxplot(list(conserved= a, gain= b, loss= c, complete.gain= cc, complete.loss= d, random= e), ylab="Reads per kilobase", main="H3K4me3", outline=FALSE)

a <- makePlot(tss[indx,], "H3K4me1", name="H3K4me1")
b <- makePlot(tss[indx_hg19_gain & (tss$V8 > lowth | tss$V9 > lowth),], "H3K4me1", name="H3K4me1 gain")
c <- makePlot(tss[indx_hg19_loss & tss$V7 > lowth,], "H3K4me1", name="H3K4me1 loss")
cc<- makePlot(tss[indx_hg19_gain & tss$V8 < lowth & tss$V9 < lowth,], "H3K4me1", name="H3K4me1 complete gain")
d <- makePlot(tss[indx_hg19_loss & tss$V7 < lowth,], "H3K4me1", name="H3K4me1 complete loss")
e <- makePlot(random_sites, "H3K4me1", name="H3K4me1 random")
boxplot(list(conserved= a, gain= b, loss= c, complete.gain= cc, complete.loss= d, random= e), ylab="Reads per kilobase", main="H3K4me1", outline=FALSE)

a <- makePlot(tss[indx,], "MeDIP-Seq", name="MeDIP-Seq")
b <- makePlot(tss[indx_hg19_gain & (tss$V8 > lowth | tss$V9 > lowth),], "MeDIP-Seq", name="MeDIP-Seq gain")
c <- makePlot(tss[indx_hg19_loss & tss$V7 > lowth,], "MeDIP-Seq", name="MeDIP-Seq loss")
cc<- makePlot(tss[indx_hg19_gain & tss$V8 < lowth & tss$V9 < lowth,], "MeDIP-Seq", name="MeDIP-Seq complete gain")
d <- makePlot(tss[indx_hg19_loss & tss$V7 < lowth,], "MeDIP-Seq", name="MeDIP-Seq complete loss")
e <- makePlot(random_sites, "MeDIP-Seq", name="MeDIP-Seq random")
boxplot(list(conserved= a, gain= b, loss= c, complete.gain= cc, complete.loss= d, random= e), ylab="Reads per kilobase", main="MeDIP-seq", outline=FALSE)

dev.off()

makeHeatmap <- function(bed, path, halfWindow=5000, step=25) {
bw <- load.bigWig(paste("/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw", sep="")) 
hm <- bed.step.bpQuery.bigWig(bw, center.bed(tss[indx,1:3], 5000, 5000), step=25)
hm_mat <- t(matrix(unlist(hm), NROW(hm[[1]])))
}

#2# Data playtime!

