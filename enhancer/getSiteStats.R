## This analysis focuses on lineage-specific changes that are currently active in hg19.
load("../annotations/fdr.RData")
source("../lib/normalizeSubsample.R")

tss_aln <- fdr_df[grepl("dREG", ca$annot_type),]
tss <- read.table("tss.tsv")
tss <- data.frame(tss, tss_aln[match(tss$V4, tss_aln$name),c(9,31:42)])

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

## Scatterplots
pdf("dist.stability.summaries.pdf")

source("../lib/densScatterplot.R")
#densScatterplot(dist, log(stab,10), xlim=c(0,50000))

plot(dist, stab, log="y", xlab="Distance to nearest gene [bp]", ylab="Predicted instability score")
abline(v=500, col="red", lty="dotted");abline(v=10000, col="red", lty="dotted");
abline(h= 0.001, col="red", lty="dotted");abline(h= 0.1, col="red", lty="dotted")
abline(v=250000, col="red", lty="dotted")

plot(dist, stab, log="y", xlim=c(0,30000), xlab="Distance to nearest gene [bp]", ylab="Predicted instability score")

plot(dist, stab, log="y", xlim=c(0,30000), xlab="Distance to nearest gene [bp]", ylab="Predicted instability score")
abline(v=500, col="red", lty="dotted");abline(v=10000, col="red", lty="dotted");
abline(h= 0.001, col="red", lty="dotted");abline(h= 0.1, col="red", lty="dotted")

hist(log(dist+1, 10))
hist(log(stab, 10))

dev.off()

## Change unscored to 0
for(i in 7:12) { tss[is.na(tss[,i]),i] <- 0 }

##################################################
## Look at species...
## Count types of elements ... and create discriptive plots.
total <- summary(as.factor(tss$V5))
one2m <- summary(as.factor(tss$V5[tss$V20 > 0])) ## Possible 1:many orthology
unmap <- summary(as.factor(tss$V5[tss$V20 == 0 & is.na(tss$mapSize)])) ## INDEL
lccng <- summary(as.factor(tss$V5[(tss$V7 < 0.1 | tss$V8 < 0.1 | tss$V9 < 0.1) & tss$V20 == 0 & !is.na(tss$mapSize)])) # 'Low-confidence'
chang <- summary(as.factor(tss$V5[tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05 & (tss$V7 > 0.7 & tss$V8 > 0.7 & tss$V9 > 0.7)])) # 'High-confidence'

chang_U2PI <- summary(as.factor(tss$V5[tss$V20 > 0 & !is.na(tss$mapSize) & tss$U2PIFDR_H < 0.05]))
lccng_U2PI <- summary(as.factor(tss$V5[tss$V20 > 0 & !is.na(tss$mapSize) & ((tss$V7 < 0.1 & tss$V10 > 0.7) | (tss$V7 > 0.7 & tss$V10 < 0.1))]))

total
one2m
unmap
chang
lccng 
chang_U2PI
lccng_U2PI

one2m/total ## Seems strange that there's such an enrichment in genes here ...
unmap/(total - one2m)
chang/(total - one2m - unmap) ## ADDED unmap to these sites b/c need to factor out for paper.
lccng/(total - one2m - unmap)

chang_U2PI/(total - one2m - unmap)
lccng_U2PI/(total - one2m - unmap)

## Compare changes at superenhancers.
SEtotal <- summary(as.factor(tss$V5[tss$V19 == 1]))
SEone2m <- summary(as.factor(tss$V5[tss$V20 > 0 & tss$V19 == 1])) ## Possible 1:many orthology
SEunmap <- summary(as.factor(tss$V5[tss$V20 == 0 & is.na(tss$mapSize) & tss$V19 == 1])) ## INDEL
SElccng <- summary(as.factor(tss$V5[(tss$V7 < 0.1 | tss$V8 < 0.1 | tss$V9 < 0.1) & tss$V20 == 0 & !is.na(tss$mapSize) & tss$V19 == 1])) # 'Low-confidence'
SEchang <- summary(as.factor(tss$V5[tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05 & (tss$V7 > 0.7 & tss$V8 > 0.7 & tss$V9 > 0.7) & tss$V19 == 1])) # 'High-confidence'

SEone2m/SEtotal ## Seems strange that there's such an enrichment in genes here ...
SEunmap/(SEtotal - SEone2m)
SEchang/(SEtotal - SEone2m - SEunmap) ## ADDED unmap to these sites b/c need to factor out for paper.
SElccng/(SEtotal - SEone2m - SEunmap)

## Create a barplot.
require(ggplot2)
library(reshape2)

dat <- data.frame(type= names(unmap), gap= unmap/(total - one2m), change= chang/(total - one2m), lc_change= lccng/(total - one2m), same=(total - one2m - unmap - chang - lccng)/(total - one2m))
dat <- melt(dat)

SEdat <- data.frame(type= names(SEunmap), gap= SEunmap/(SEtotal - SEone2m), change= SEchang/(SEtotal - SEone2m), lc_change= SElccng/(SEtotal - SEone2m), same=(SEtotal - SEone2m - SEunmap - SEchang - SElccng)/(SEtotal - SEone2m))
SEdat <- melt(SEdat)

a <-  ggplot(dat, aes(x=variable, y=value)) + xlab("") + ylab("Fraction of Transcripts Changed") +
        geom_bar(aes(fill=type), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")

b <-  ggplot(dat, aes(x=type, y=value)) + xlab("") + ylab("Fraction of Transcripts Changed") +
        geom_bar(aes(fill=variable), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")

c <-  qplot(x= factor(type), y= value, stat="identity", data= dat, geom="bar", fill=factor(variable))
SEc<-  qplot(x= factor(type), y= value, stat="identity", data= SEdat, geom="bar", fill=factor(variable))


pdf("enh.type.change.barplot.pdf")
        a
        b
        c
dev.off()

## Density scatterplots at promoters and enhancers
source("../lib/densScatterplot.R")

#pdf("dreg.scatterplots.pdf")
# max_hc <- rowMax(tss[,c(7,8)])
# densScatterplot(tss$V7, tss$V8, xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
# densScatterplot(tss$V7[tss$V19 == "Prox_Stab" & max_hc>0.7], tss$V8[tss$V19 == "Prox_Stab" & max_hc>0.7], xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
# densScatterplot(tss$V7[tss$V19 == "Dist_UnSt" & max_hc>0.7], tss$V8[tss$V19 == "Dist_UnSt" & max_hc>0.7], xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
#
# max_hm <- rowMax(tss[,c(7,9)])
# densScatterplot(tss$V7[max_hm>0.7], tss$V9[max_hm>0.7], xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")
# densScatterplot(tss$V7[tss$V19 == "Prox_Stab" & max_hm>0.7], tss$V9[tss$V19 == "Prox_Stab" & max_hm>0.7], xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")
# densScatterplot(tss$V7[tss$V19 == "Dist_UnSt" & max_hm>0.7], tss$V9[tss$V19 == "Dist_UnSt" & max_hm>0.7], xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")
#
# max_up <- rowMax(tss[,c(7,10)])
# densScatterplot(tss$V7, tss$V10, xlab="Human Unt.", ylab="Human PI", main="dREG Scores")
# densScatterplot(tss$V7[tss$V19 == "Prox_Stab" & max_up>0.7], tss$V10[tss$V19 == "Prox_Stab" & max_up>0.7], xlab="Human Unt.", ylab="Human PI", main="dREG Scores")
# densScatterplot(tss$V7[tss$V19 == "Dist_UnSt" & max_up>0.7], tss$V10[tss$V19 == "Dist_UnSt" & max_up>0.7], xlab="Human Unt.", ylab="Human PI", main="dREG Scores")
#dev.off()


##########################################################
## Questions ... :

##########################################################
## Do SE change less frequently?
cmpFracConserved <- function(tss1, tss2, i=2) { ## defaults to enhancer (i=2) ... promoters (i=3)
 total1 <- summary(as.factor(tss1$V5))
 one2m1 <- summary(as.factor(tss1$V5[tss1$V20 > 0])) ## Possible 1:many orthology
 unmap1 <- summary(as.factor(tss1$V5[tss1$V20 == 0 & is.na(tss1$mapSize)])) ## INDEL
 lccng1 <- summary(as.factor(tss1$V5[(tss1$V7 < 0.1 | tss1$V8 < 0.1 | tss1$V9 < 0.1) & tss1$V20 == 0 & !is.na(tss1$mapSize)])) # 'Low-confidence'
 chang1 <- summary(as.factor(tss1$V5[tss1$V20 == 0 & !is.na(tss1$mapSize) & tss1$fdr_min < 0.05 & (tss1$V7 > 0.7 & tss1$V8 > 0.7 & tss1$V9 > 0.7)])) # 'High-confidence'
 allcng1<- summary(as.factor(tss1$V5[(tss1$V20 == 0 & !is.na(tss1$mapSize) & tss1$fdr_min < 0.05)]))

 tot1 <- total1[i]-one2m1[i]
 con1 <- total1[i]-one2m1[i]-unmap1[i]-chang1[i]-lccng1[i]
 print(paste("1: ", con1/tot1))

 total2 <- summary(as.factor(tss2$V5))
 one2m2 <- summary(as.factor(tss2$V5[tss2$V20 > 0])) ## Possible 1:many orthology
 unmap2 <- summary(as.factor(tss2$V5[tss2$V20 == 0 & is.na(tss2$mapSize)])) ## INDEL
 lccng2 <- summary(as.factor(tss2$V5[(tss2$V7 < 0.1 | tss2$V8 < 0.1 | tss2$V9 < 0.1) & tss2$V20 == 0 & !is.na(tss2$mapSize)])) # 'Low-confidence'
 chang2 <- summary(as.factor(tss2$V5[tss2$V20 == 0 & !is.na(tss2$mapSize) & tss2$fdr_min < 0.05 & (tss2$V7 > 0.7 & tss2$V8 > 0.7 & tss2$V9 > 0.7)])) # 'High-confidence'
 allcng2<- summary(as.factor(tss2$V5[(tss2$V20 == 0 & !is.na(tss2$mapSize) & tss2$fdr_min < 0.05)]))

 tot2 <- total2[i]-one2m2[i]
 con2 <- total2[i]-one2m2[i]-unmap2[i]-chang2[i]-lccng2[i]
 print(paste("2: ", con2/tot2))

 ft <- fisher.test(data.frame(c(tot1, con1), c(tot2, con2)))
 print(ft)
 print(paste("Numerical p-value ==>", ft$p.value))
}

fracConserved <- function(tss, ret="C") { ## defaults to ... 
 total <- NROW(tss)
 one2m <- sum(tss$V20 > 0) ## Possible 1:many orthology
 unmap <- sum(tss$V20 == 0 & is.na(tss$mapSize)) ## INDEL
 lccng <- sum((tss$V7 < 0.1 | tss$V8 < 0.1 | tss$V9 < 0.1) & tss$V20 == 0 & !is.na(tss$mapSize)) # 'Low-confidence'
 chang <- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05 & (tss$V7 > 0.7 & tss$V8 > 0.7 & tss$V9 > 0.7)) # 'High-confidence'
 allcng<- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05)

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
  cng <- chang #allcng ## A santiy check to make sure this shows the same trend as lccng.

  return(cng/tot)
 }
}

cmpFracConserved(tss, tss[tss$V19==1,], i=2) ## Enhancers
cmpFracConserved(tss, tss[tss$V19==1,], i=3) ## Promoters

## Is this more so than expected based on distance bias alone?
sum(tss$V19==1 & tss$V5=="Dist_UnSt") ## Get the number of elements to use for selecting a reasonable number for re-sampling.
sum(tss$V19==0 & tss$V5=="Dist_UnSt")

ns <- norm.subsample(log(tss[tss$V19==1 & tss$V5=="Dist_UnSt","V13"]), log(tss[tss$V19==0 & tss$V5=="Dist_UnSt","V13"]), nBins=25, nsamp=1000, plot.cdf=TRUE)
cmpFracConserved(tss[tss$V19==1 & tss$V5=="Dist_UnSt",][ns$s1,], tss[tss$V19==0 & tss$V5=="Dist_UnSt",][ns$s2,], i=2) ## Enhancers

## Bootstrap the subsampling to create error bars.
cons_se  <- double()
cons_nose<- double()
for(i in 1:1000) {
 ns <- norm.subsample(log(tss[tss$V19==1 & tss$V5=="Dist_UnSt","V13"]), log(tss[tss$V19==0 & tss$V5=="Dist_UnSt","V13"]), nsamp=1000)
 cons_se  <- c(cons_se,  fracConserved(tss[tss$V19==1 & tss$V5=="Dist_UnSt",][ns$s1,]))
 cons_nose<- c(cons_nose,fracConserved(tss[tss$V19==0 & tss$V5=="Dist_UnSt",][ns$s2,]))
}

## Plot the data as a barplot.
pdf("SEBarplot.pdf")
 source("../lib/barplot.R")
 bars <- c(mean(cons_se), mean(cons_nose))
 errs <- c(sqrt(var(cons_se)), sqrt(var(cons_nose)))
 names<- c("SE", "No SE")
 cd.barplot(bars, errs, names, fill=TRUE)
dev.off()

########################################################
## Evolutionary conservation->distance from genes.
fracConserved(tss[tss$V5 == "Dist_UnSt",])
fracConserved(tss[tss$V5 == "Prox_Stab",])

#vect <- as.numeric(cut2(log(tss$V13), g=100)); summary(as.factor(vect))
require(Hmisc)

xaxis <- c(0, seq(3, 6, 0.02))
vect <- as.numeric(cut2(log(tss$V13, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
dist <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,])})
glc  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "GL")})
cng  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "CNG")})

getCex <- function(n) { y=0.0138888*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})

## Last point is [>=xaxis]. Each other point is [<point].
idx <- 3:NROW(xaxis) ## As below, index 3 encodes [1000,10^3.02)

pdf("ChangeOverDistance.pdf")
 plot(10^xaxis[idx], dist[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction conserved", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
 plot(10^xaxis[idx], glc[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction gain/ loss", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
 plot(10^xaxis[idx], cng[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction change", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
dev.off()

##################################################################
## Evolutionary conservation compared to max{expression across species}.
## No correlation.  But to really test intensity, have to remove those sites that are intragenic.
## Still no correlation...
indx_match <- match(tss$V4, tss_aln$name)

exp_max <- rowSums(rpkm_df[,2:9])[indx_match] # rowMax(rpkm_df[,2:9])[indx_match]
type_max <- ca$annot_type[indx_match]
exp_vect <- as.numeric(cut2(exp_max, g=100))
exp_vect[is.na(exp_vect)] <- 0
conserved <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x & type_max == "dREG_ENH",])})
gainloss <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x & type_max == "dREG_ENH",], "GL")})
change   <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x & type_max == "dREG_ENH",], "CNG")})

cor.test(conserved, c(1:NROW(conserved))) ## p= (for manuscript)

pdf("ChangeOverActivity.pdf")
 plot(conserved, type="b")
 plot(gainloss,  type="b")
 plot(change,    type="b")
dev.off()

###############################################################
## Do looped enhancers change more slowly.  
## Yep.
loop <- read.table("tss.tsv.loop")
fracConserved(tss[rowSums(loop[,5:6]) >  0 & tss$V5 == "Dist_UnSt",])  ## Conservation of looped REs.
fracConserved(tss[rowSums(loop[,5:6]) == 0 & tss$V5 == "Dist_UnSt",])  ## Conservation of unlooped. 

cmpFracConserved(tss[rowSums(loop[,5:6]) >  0,], tss[rowSums(loop[,5:6]) == 0,], i=2) ## Enhancers
cmpFracConserved(tss[rowSums(loop[,5:6]) >  0,], tss[rowSums(loop[,5:6]) == 0,], i=3) ## Promoters

## Is this more so than expected based on distance bias alone?
sum(rowSums(loop[,5:6]) >  0 & tss$V5=="Dist_UnSt") ## Get the number of elements to use for selecting a reasonable number for re-sampling.
sum(rowSums(loop[,5:6]) == 0 & tss$V5=="Dist_UnSt")

ns <- norm.subsample(log(tss[rowSums(loop[,5:6]) >  0 & tss$V5=="Dist_UnSt","V13"]), log(tss[rowSums(loop[,5:6]) == 0 & tss$V5=="Dist_UnSt","V13"]), nBins=25, nsamp=1000, plot.cdf=TRUE)
cmpFracConserved(tss[rowSums(loop[,5:6]) >  0,][ns$s1,], tss[rowSums(loop[,5:6]) == 0,][ns$s2,], i=2)

## Bootstrap the subsampling to create error bars.
cons_looped <- double()
cons_nonloop<- double()
for(i in 1:1000) {
 ns <- norm.subsample(log(tss[rowSums(loop[,5:6]) >  0 & tss$V5=="Dist_UnSt","V13"]), log(tss[rowSums(loop[,5:6]) == 0 & tss$V5=="Dist_UnSt","V13"]), nsamp=1000)
 cons_looped <- c(cons_looped, fracConserved(tss[rowSums(loop[,5:6]) >  0 & tss$V5=="Dist_UnSt",][ns$s1,]))
 cons_nonloop<- c(cons_nonloop,fracConserved(tss[rowSums(loop[,5:6]) == 0 & tss$V5=="Dist_UnSt",][ns$s2,]))
}

## Plot the data as a barplot.
pdf("LoopBarplot.pdf")
 source("../lib/barplot.R")
 bars <- c(mean(cons_looped), mean(cons_nonloop))
 errs <- c(sqrt(var(cons_looped)), sqrt(var(cons_nonloop)))
 names<- c("Looped", "No loop")
 cd.barplot(bars, errs, names, fill=TRUE)
# drawBars(bars, errs, names) ## From the eRNA regression code in the dREG paper.  Does not work as well out of the box.
dev.off()

####################################################
## How much are SE over-represented in the loop set?
sum(tss$V19==1 & rowSums(loop[,5:6])>0)/ sum(tss$V19==1) ## 49% of SE loop.
sum(tss$V19==0 & rowSums(loop[,5:6])>0)/ sum(tss$V19==0) ## 16% of all enhancers.

########################################################################
## Compare looped and SEs in scatterplot. 

xaxis <- c(0, seq(3, 6, 0.2))

vect <- as.numeric(cut2(log(tss$V13, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
dist <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,])})
dist_se <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x & tss$V19 == 1,])})
dist_loop <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x & rowSums(loop[,5:6]) >  0,])})

## Size points based on 'n'
getCex <- function(n) { y=0.0138888*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})
n_se <- sapply(1:max(vect), function(x) {NROW(tss[vect == x & tss$V19 == 1,])})
n_loop <- sapply(1:max(vect), function(x) {NROW(tss[vect == x & rowSums(loop[,5:6]) >  0,])})


## Create barplot
xaxis <- c(0,seq(3, 6, 1))

vect <- as.numeric(cut2(log(tss$V13, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
b_dist <- boot(data= tss, R=1000, statistic= function(a, i) {sapply(1:max(vect), function(x) {fracConserved(a[i,][vect == x,])})})
b_dist_se <- boot(data= tss[tss$V19 == 1,], R=1000, statistic= function(a, i) {sapply(1:max(vect), function(x) {fracConserved(a[i,][vect[tss$V19 == 1] == x,])})})
b_dist_loop <- boot(data= tss[rowSums(loop[,5:6]) >  0,], R=1000, statistic= function(a, i) {sapply(1:max(vect), function(x) {fracConserved(a[i,][vect[rowSums(loop[,5:6]) > 0] == x,])})})

idx <- 3:5 #summary(cut2(log(tss$V13, 10), cuts=xaxis)) ## We want idx: 3 ([3.00,4.00)) - 5 ([5.00,6.00))
names<- paste(rep(c("[1-10)", "[10-100)", "[100-1000)"),3), c(rep("all", 3), rep("se", 3), rep("loop", 3)))
bars <- c(b_dist$t0[idx], b_dist_se$t0[idx], b_dist_loop$t0[idx])
serr <- c(sapply(1:NROW(b_dist$t0) , function(x) {sd(b_dist$t[,x], na.rm=TRUE)})[idx], 
	sapply(1:NROW(b_dist_se$t0) , function(x) {sd(b_dist_se$t[,x], na.rm=TRUE)})[idx],
	sapply(1:NROW(b_dist_loop$t0) , function(x) {sd(b_dist_loop$t[,x], na.rm=TRUE)})[idx])

ord <- c(seq(1,9,3), seq(1,9,3)+1, seq(1,9,3)+2)

pdf("Distance_SE_loop.pdf")
 ## Scatterplot
 plot(10^xaxis[1:NROW(dist)], dist, type="p", xlab="Distance from TSS [bp]", ylab="Fraction conserved", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n), ylim=c(0,0.8))
 points(10^xaxis[1:NROW(dist)], dist_loop, col="dark red", type="p", pch=19, cex=getCex(n_loop))
 points(10^xaxis[1:NROW(dist)], dist_se, col="dark green", type="p", pch=19, cex=getCex(n_se))

 ## Barplot
 cd.barplot(bars[ord], serr[ord], names[ord], fill=TRUE, order=FALSE)
dev.off()




