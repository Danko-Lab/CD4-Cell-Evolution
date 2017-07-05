## This analysis focuses on lineage-specific changes that are currently active in hg19.
load("../annotations/fdr.RData")
source("../lib/normalizeSubsample.R")

stp <- read.table("../tss_caller/proximity.strand.tsv")

require(boot)
require(pheatmap)

lnTH= 0.05
hcTH= 0.30

tss_aln <- fdr_df[grepl("dREG", ca$annot_type),]
tss <- read.table("tss.tsv")
tss <- data.frame(tss, tss_aln[match(tss$V4, tss_aln$name),c(9,33:50)])

summary(tss$V13 - abs(stp$V10)) ## Sanity check ... should all be 0.  Yep.
tss$V13 <- stp$V10

## Alignable fraction (V20) denotes a gap in either species.  Make sure gaps are in both.lignable fraction (V20) denotes a gap in either species.  Make sure gaps are in both.

## Classify as 'promoter'/ 'enhancer'
stab <- rowMax(tss[,17:18])
dist <- tss[,13]
class <- rep("tss", NROW(tss)) ## tss is then unclassified as a promoter or enhancer
class[dist < 100]  <- "Prox_Stab" ## stab < 0.1
class[dist > 5000] <- "Dist_UnSt" ## stab > 0.1
#class[stab < 0.1 & dist > 1e6] <- "Dist_Stab" ## stab < 0.1
summary(as.factor(class))
tss$V5 <- factor(class, levels= c("Dist_Stab", "Dist_UnSt", "Prox_Stab", "tss"))

## Scatterplots
## Change unscored to 0.  Untested (reflects possible unmappable regions) to 1.
for(i in 7:12) { tss[is.na(tss[,i]),i] <- 0 }

## Create a barplot.
require(ggplot2)
library(reshape2)

## Density scatterplots at promoters and enhancers
source("../lib/densScatterplot.R")



## Classify as 'promoter'/ 'enhancer'
stab <- rowMax(tss[,17:18])
dist <- tss[,13]
class <- rep("tss", NROW(tss)) ## tss is then unclassified as a promoter or enhancer
class[dist < 100]  <- "Prox_Stab" ## stab < 0.1
class[dist > 5000] <- "Dist_UnSt" ## stab > 0.1
#class[stab < 0.1 & dist > 1e6] <- "Dist_Stab" ## stab < 0.1
summary(as.factor(class))
tss$V5 <- factor(class, levels= c("Dist_Stab", "Dist_UnSt", "Prox_Stab", "tss"))

## Scatterplots
## Change unscored to 0.  Untested (reflects possible unmappable regions) to 1.
for(i in 7:12) { tss[is.na(tss[,i]),i] <- 0 }

## Create a barplot.
require(ggplot2)
library(reshape2)

## Density scatterplots at promoters and enhancers
source("../lib/densScatterplot.R")

##########################################################
## Questions ... :

##########################################################
## Do SE change less frequently?
cmpFracConserved <- function(tss1, tss2, i=2, i2=i) { ## defaults to enhancer (i=2) ... promoters (i=3)
 total1 <- summary(as.factor(tss1$V5))
 one2m1 <- summary(as.factor(tss1$V5[tss1$V20 > 0])) ## Possible 1:many orthology
 unmap1 <- summary(as.factor(tss1$V5[tss1$V20 == 0 & is.na(tss1$mapSize)])) ## INDEL
 lccng1 <- summary(as.factor(tss1$V5[(tss1$V7 < lnTH | tss1$V8 < lnTH | tss1$V9 < lnTH) & (tss1$V7 > hcTH | tss1$V8 > hcTH | tss1$V9 > hcTH) & tss1$V20 == 0 & !is.na(tss1$mapSize)])) # Complete change.
 chang1 <- summary(as.factor(tss1$V5[tss1$V20 == 0 & !is.na(tss1$mapSize) & tss1$fdr_min < PVAL & (tss1$V7 > hcTH | tss1$V8 > hcTH | tss1$V9 > hcTH) & (tss1$V7 > lnTH & tss1$V8 > lnTH & tss1$V9 > lnTH)])) # Partial change.
 allcng1<- summary(as.factor(tss1$V5[(tss1$V20 == 0 & !is.na(tss1$mapSize) & tss1$fdr_min < PVAL)]))

 tot1 <- total1[i]-one2m1[i]
 con1 <- total1[i]-one2m1[i]-unmap1[i]-chang1[i]-lccng1[i]
 print(paste("1: ", con1/tot1))

 total2 <- summary(as.factor(tss2$V5))
 one2m2 <- summary(as.factor(tss2$V5[tss2$V20 > 0])) ## Possible 1:many orthology
 unmap2 <- summary(as.factor(tss2$V5[tss2$V20 == 0 & is.na(tss2$mapSize)])) ## INDEL
 lccng2 <- summary(as.factor(tss2$V5[(tss2$V7 < lnTH | tss2$V8 < lnTH | tss2$V9 < lnTH) & (tss2$V7 > hcTH | tss2$V8 > hcTH | tss2$V9 > hcTH) & tss2$V20 == 0 & !is.na(tss2$mapSize)])) # Complete change.
 chang2 <- summary(as.factor(tss2$V5[tss2$V20 == 0 & !is.na(tss2$mapSize) & tss2$fdr_min < PVAL & (tss2$V7 > hcTH | tss2$V8 > hcTH | tss2$V9 > hcTH) & (tss2$V7 > lnTH & tss2$V8 > lnTH & tss2$V9 > lnTH)])) # Partial change.
 allcng2<- summary(as.factor(tss2$V5[(tss2$V20 == 0 & !is.na(tss2$mapSize) & tss2$fdr_min < PVAL)]))

 tot2 <- total2[i2]-one2m2[i2]
 con2 <- total2[i2]-one2m2[i2]-unmap2[i2]-chang2[i2]-lccng2[i2]
 print(paste("2: ", con2/tot2))

 ft <- fisher.test(data.frame(c(tot1, con1), c(tot2, con2)))
 print(ft)
 print(paste("Numerical p-value ==>", ft$p.value))
}

fracConserved <- function(tss, ret="C") { ## defaults to ... 
 total <- NROW(tss)
 one2m <- sum(tss$V20 > 0) ## Possible 1:many orthology
 unmap <- sum(tss$V20 == 0 & is.na(tss$mapSize)) ## INDEL
 lccng <- sum((tss$V7 < lnTH | tss$V8 < lnTH | tss$V9 < lnTH) & (tss$V7 > hcTH | tss$V8 > hcTH | tss$V9 > hcTH) & tss$V20 == 0 & !is.na(tss$mapSize)) # Complete change.
 chang <- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < PVAL & (tss$V7 > hcTH | tss$V8 > hcTH | tss$V9 > hcTH) & (tss$V7 > lnTH & tss$V8 > lnTH & tss$V9 > lnTH)) # 'High-confidence' partial change.
 allcng<- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < PVAL)

 if(ret=="C") {
  tot <- total-one2m
  con <- total-one2m-unmap-chang-lccng

  return(con/tot)
 }
 if(ret=="GL") {
  tot <- total-one2m
  lcc <- lccng

  return(lcc/tot)
 }
 if(ret=="CNG") {
  tot <- total-one2m
  cng <- allcng ## A santiy check to make sure this shows the same trend as lccng.

  return(cng/tot)
 }
 if(ret=="UNM") {
  tot <- total-one2m
  unm <- unmap #allcng ## A santiy check to make sure this shows the same trend as lccng.

  return(unm/tot)
 }

}

#vect <- as.numeric(cut2(log(tss$V13), g=100)); summary(as.factor(vect))
require(Hmisc)

xaxis <- c(seq(-6, 6, 0.1))#, 0, seq(0.5, 6, 0.1))
cvals <- rep(0, NROW(tss)); cvals[tss$V13 > 0] <-  log(tss$V13[tss$V13 > 0], 10); cvals[tss$V13 < 0] <-  -1* log(-1*tss$V13[tss$V13 < 0], 10)
vect <- as.numeric(cut2(cvals, cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
dist <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,])})
#glc  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "GL")})
#cng  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "CNG")})
#unm  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "UNM")})

getCex <- function(n) { y=0.00538888*n+0.1; y[y>3] <- 3; y[y<0.01] <- 0.01; y*1.5 }
n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})

size_key <- seq(0, 200,by= 20)

## Last point is [>=xaxis]. Each other point is [<point].
idx <- 1:NROW(xaxis) ## As below, index 3 encodes [1000,10^3.02)

pdf("ChangeOverDistance-sanity-check-direction.pdf")
 plot(xaxis[idx], dist[idx], type="p", xlab="Log-10 distance from TSS [bp]", ylab="Fraction conserved", xlim=c(-6, 6), ylim=c(0.3, 0.8), pch=19, cex=getCex(n[idx]), cex.axis= 1.5, cex.lab=1.5)
dev.off()

