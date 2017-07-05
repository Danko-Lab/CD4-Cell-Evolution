## This analysis focuses on lineage-specific changes that are currently active in hg19.
load("../annotations/fdr.RData")
source("../lib/normalizeSubsample.R")

require(boot)
require(pheatmap)

lnTH= 0.05
hcTH= 0.30

tss_aln <- fdr_df[grepl("dREG", ca$annot_type),]
tss <- read.table("tss.tsv")
tss <- data.frame(tss, tss_aln[match(tss$V4, tss_aln$name),c(9,33:50)])

## Alignable fraction (V20) denotes a gap in either species.  Make sure gaps are in both.

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

## Change unscored to 0.  Untested (reflects possible unmappable regions) to 1.
for(i in 7:12) { tss[is.na(tss[,i]),i] <- 0 }

## H-C
lccng <- summary(as.factor(tss$V5[(tss$V7 < lnTH | tss$V8 < lnTH) & tss$V20 == 0 & !is.na(tss$mapSize)])) # 'Low-confidence'
chang <- summary(as.factor(tss$V5[tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < PVAL & (tss$V7 > hcTH & tss$V8 > hcTH & tss$V9 > hcTH)])) # 'High-confidence'



##################################################
## Look at species...
## Count types of elements ... and create discriptive plots.
total <- summary(as.factor(tss$V5))
one2m <- summary(as.factor(tss$V5[tss$V20 > 0])) ## Possible 1:many orthology
unmap <- summary(as.factor(tss$V5[tss$V20 == 0 & is.na(tss$mapSize)])) ## INDEL
lccng <- summary(as.factor(tss$V5[(tss$V7 < lnTH | tss$V8 < lnTH | tss$V9 < lnTH) & (tss$V7 > hcTH | tss$V8 > hcTH | tss$V9 > hcTH) & tss$V20 == 0 & !is.na(tss$mapSize)])) # 'Low-confidence'
chang <- summary(as.factor(tss$V5[tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < PVAL & (tss$V7 > hcTH | tss$V8 > hcTH | tss$V9 > hcTH) & (tss$V7 > lnTH & tss$V8 > lnTH & tss$V9 > lnTH)])) # 'High-confidence'
hctot <- summary(as.factor(tss$V5[tss$V20 == 0 & tss$fdr_min < PVAL & !is.na(tss$mapSize)])) # 'HC all'

lccng_U2PI <- summary(as.factor(tss$V5[((tss$V7 < lnTH & tss$V10 > hcTH) | (tss$V7 > hcTH & tss$V10 < lnTH)) & !is.na(tss$mapSize)]))
hccng_U2PI <- summary(as.factor(tss$V5[(tss$U2PIFDR_H < PVAL) & !((tss$V7 < lnTH & tss$V10 > hcTH) | (tss$V7 > hcTH & tss$V10 < lnTH)) & !is.na(tss$mapSize)]))
ALLchang_U2PI <- summary(as.factor(tss$V5[tss$U2PIFDR_H < PVAL & !is.na(tss$mapSize)]))

## Create a barplot.
require(ggplot2)
library(reshape2)

dat <- data.frame(type= names(unmap), gap= unmap/(total - one2m), change= chang/(total - one2m), lc_change= lccng/(total - one2m), same=(total - one2m - unmap - chang - lccng)/(total - one2m), hctot=hctot/(total - one2m) )
dat <- melt(dat)
dat

U2PIdat <- data.frame( type= names(unmap), gap= rep(0,4), change= hccng_U2PI/(total - one2m), lc_change= lccng_U2PI/(total - one2m), same=(total - one2m - unmap - hccng_U2PI - lccng_U2PI)/(total - one2m), hctot=ALLchang_U2PI/(total - one2m) )
U2PIdat <- melt(U2PIdat)
U2PIdat

a <-  ggplot(dat, aes(x=variable, y=value)) + xlab("") + ylab("Fraction of Transcripts Changed") +
        geom_bar(aes(fill=type), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")

b <-  ggplot(dat, aes(x=type, y=value)) + xlab("") + ylab("Fraction of Transcripts Changed") +
        geom_bar(aes(fill=variable), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")

aU2PI <-  ggplot(U2PIdat, aes(x=variable, y=value)) + xlab("") + ylab("Fraction of Transcripts Changed") +
        geom_bar(aes(fill=type), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")

bU2PI <-  ggplot(U2PIdat, aes(x=type, y=value)) + xlab("") + ylab("Fraction of Transcripts Changed") +
        geom_bar(aes(fill=variable), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")

pdf("enh.type.change.barplot.pdf")
        a
        b
	aU2PI
	bU2PI
dev.off()

## Create a color scale for complete changes.
pval_du <- tss$fdr_min[( (tss$V7 < lnTH | tss$V8 < lnTH | tss$V9 < lnTH) & (tss$V7 > hcTH | tss$V8 > hcTH | tss$V9 > hcTH) & tss$V20 == 0 & !is.na(tss$mapSize) & tss$V5 == "Dist_UnSt" )]
pval_ps <- tss$fdr_min[( (tss$V7 < lnTH | tss$V8 < lnTH | tss$V9 < lnTH) & (tss$V7 > hcTH | tss$V8 > hcTH | tss$V9 > hcTH) & tss$V20 == 0 & !is.na(tss$mapSize) & tss$V5 == "Prox_Stab" )]

pval_du_upi <- tss$fdr_min[((tss$V7 < lnTH & tss$V10 > hcTH) | (tss$V7 > hcTH & tss$V10 < lnTH)) & !is.na(tss$mapSize) & tss$V5 == "Dist_UnSt" ]
pval_ps_upi <- tss$fdr_min[((tss$V7 < lnTH & tss$V10 > hcTH) | (tss$V7 > hcTH & tss$V10 < lnTH)) & !is.na(tss$mapSize) & tss$V5 == "Prox_Stab" ]

breaks     <- rev(c(seq(0.01, 0.10, 0.005), seq(0.10, 0.20, 0.05)))
colorscale <- colorRampPalette(c("#EEFFF4","#00BFC4"))(length(breaks)-1)

pdf("pv.colorscale.lcchanges.pdf")

pheatmap(rev(sort(pval_du)), cluster_rows = FALSE, cluster_cols = FALSE, col= colorscale, breaks = breaks, legend=TRUE, legend_breaks= c(0.01, 0.05, 0.1, 0.2), legend_labels= c(0.01, 0.05, 0.1, 0.2), show_rownames=FALSE, show_colnames=FALSE)
pheatmap(rev(sort(pval_ps)), cluster_rows = FALSE, cluster_cols = FALSE, col= colorscale, breaks = breaks, legend=TRUE, legend_breaks= c(0.01, 0.05, 0.1, 0.2), legend_labels= c(0.01, 0.05, 0.1, 0.2), show_rownames=FALSE, show_colnames=FALSE)

pheatmap(rev(sort(pval_du_upi)), cluster_rows = FALSE, cluster_cols = FALSE, col= colorscale, breaks = breaks, legend=TRUE, legend_breaks= quantile(breaks), legend_labels= signif(quantile(breaks)), show_rownames=FALSE, show_colnames=FALSE)
pheatmap(rev(sort(pval_ps_upi)), cluster_rows = FALSE, cluster_cols = FALSE, col= colorscale, breaks = breaks, legend=TRUE, legend_breaks= quantile(breaks), legend_labels= signif(quantile(breaks)), show_rownames=FALSE, show_colnames=FALSE)

dev.off()

## Summary stats.
sum(pval_du < 0.2)/ NROW(pval_du)  ## 64%
sum(pval_du < 0.1)/ NROW(pval_du)  ## 54%
sum(pval_du < 0.05)/ NROW(pval_du) ## 46%
sum(pval_du < 0.01)/ NROW(pval_du) ## 32%

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
unm  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "UNM")})

getCex <- function(n) { y=0.0138888*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})

size_key <- seq(0, 200,by= 20)

## Last point is [>=xaxis]. Each other point is [<point].
idx <- 3:NROW(xaxis) ## As below, index 3 encodes [1000,10^3.02)

## Compute for Supplementary Table 2
cor.test(dist[idx], idx)
corr(cbind(dist[idx], idx), w= n[idx]/sum(n[idx]))

pdf("ChangeOverDistance.pdf")
 plot(10^xaxis[idx], dist[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction conserved", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
 plot(10^xaxis[idx], glc[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction gain/ loss", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
 plot(10^xaxis[idx], cng[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction change", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
 plot(10^xaxis[idx], unm[idx], type="p", xlab="Distance from TSS [bp]", ylab="Fraction unmappable", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))

 plot(size_key, rep(1, NROW(size_key)), cex=getCex(size_key))
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
gaps     <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x & type_max == "dREG_ENH",], "UNM")})


cor.test(conserved, c(1:NROW(conserved))) ## p= (for manuscript)

pdf("ChangeOverActivity.pdf")
 plot(conserved, type="b")
 plot(gainloss,  type="b")
 plot(change,    type="b")
 plot(gaps,      type="b")
dev.off()

###############################################################
## Do looped enhancers change more slowly.  
## Yep.
loop <- read.table("tss.tsv.loop")

fracConserved(tss[rowSums(loop[,5:6]) >  0 & tss$V5 == "Dist_UnSt",])  ## Conservation of looped REs.
fracConserved(tss[rowSums(loop[,5:6]) == 0 & tss$V5 == "Dist_UnSt",])  ## Conservation of unlooped. 

fracConserved(tss[rowSums(loop[,5:6]) >  0 & tss$V5 == "Prox_Stab",])  ## Conservation of looped REs.
fracConserved(tss[rowSums(loop[,5:6]) == 0 & tss$V5 == "Prox_Stab",])  ## Conservation of unlooped. 

cmpFracConserved(tss[rowSums(loop[,5:6]) >  0,], tss[rowSums(loop[,5:6]) == 0,], i=2) ## Enhancers
cmpFracConserved(tss[rowSums(loop[,5:6]) >  0,], tss[rowSums(loop[,5:6]) == 0,], i=3) ## Promoters

## Reproducible with promoter capture Hi-C data?
cmpFracConserved(tss[loop[,8] >  0,], tss[loop[,8] == 0,], i=2) ## Enhancers
cmpFracConserved(tss[rowSums(loop[,7:8]) >  0,], tss[rowSums(loop[,7:8]) == 0,], i=2) ## Enhancers
cmpFracConserved(tss[loop[,8] >  0,], tss[loop[,8] == 0,], i=3) ## Promoters
cmpFracConserved(tss[rowSums(loop[,7:8]) >  0,], tss[rowSums(loop[,7:8]) == 0,], i=3) ## Promoters

## Compare looped enhancers to promoters.
cmpFracConserved(tss[rowSums(loop[,5:6]) >  0,], tss, i=2, i2=3)
b_du <- boot(data= tss[tss$V5=="Dist_UnSt",], R=1000, statistic= function(a, i) {fracConserved(a[i,])})
b_du_loop <- boot(data= tss[tss$V5=="Dist_UnSt" & rowSums(loop[,5:6]) >  0,], R=1000, statistic= function(a, i) {fracConserved(a[i,])})
b_ps <- boot(data= tss[tss$V5=="Prox_Stab",], R=1000, statistic= function(a, i) {fracConserved(a[i,])})
b_ps_loop <- boot(data= tss[tss$V5=="Prox_Stab" & rowSums(loop[,5:6]) >  0,], R=1000, statistic= function(a, i) {fracConserved(a[i,])})

pdf("LoopEnhPromoter.pdf")
 source("../lib/barplot.R")
 bars <- c(b_du$t0, b_du_loop$t0, b_ps$t0, b_ps_loop$t0)
 errs <- c(sqrt(var(b_du$t)), sqrt(var(b_du_loop$t)), sqrt(var(b_ps$t)), sqrt(var(b_ps_loop$t)))
 names<- c("All Enhancer", "Looped Enhancer", "All Promoter", "Looped Promoter")
 cd.barplot(bars, errs, names, fill=TRUE)
# drawBars(bars, errs, names) ## From the eRNA regression code in the dREG paper.  Does not work as well out of the box.
dev.off()

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

##########################################################
## Compare to strength of loop interaction.
## 
require(Hmisc)

usethese <- loop[,7]==0 & tss$V5 == "Dist_UnSt" ## Loop[,7] == 0 --> Excludes any TREs that overlap Hi-C capture anchors.

xaxis <- c(0, seq(5, max(loop[,9]), 2))
vect <- as.numeric(cut2(loop$V9, cuts=xaxis)); vect[is.na(vect)] <- 0; summary(as.factor(vect))
strn <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x & usethese,])})
dist <- sapply(1:max(vect), function(x) {median(tss[vect == x & loop[,7]==0 & usethese,13])})

getCex <- function(n) { y=0.005*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})

size_key <- seq(0, 1000,by= 50)

## Last point is [>=xaxis]. Each other point is [<point].
idx <- !is.nan(strn)

cor.test(xaxis[idx], strn[idx], method="spearman")
corr(cbind(xaxis[idx], strn[idx]), w = n[idx]/sum(n[idx]))

## ADD A LOOP SWAP HERE..
## Evaluate significance of this weighted correlation using a label swap.
loop_swap <- boot(data= vect, R=1000, statistic= function(a,i) {
        strn_b <- sapply(1:max(vect), function(x) {fracConserved(tss[a[i] == x & usethese,])})  ## Conservation of unlooped. 
	n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})
        indx <- !is.na(strn_b)
        corr(cbind(xaxis, strn_b)[indx,], w = n[indx]/sum(n[indx]))
}) ## Here, using boot to swap labels.
sum(loop_swap$t >= loop_swap$t0)/NROW(loop_swap$t)

pdf("dHi-C.conservation.pdf")
 plot(xaxis[idx], strn[idx], type="p", xlab="Loop Strength [dHi-C]", ylab="Fraction conserved", pch=19, cex=getCex(n[idx]), xlim=c(0,50), ylim=c(0.3, 0.8))
 plot(loop_swap)
 plot(size_key, rep(1, NROW(size_key)), cex=getCex(size_key))
dev.off()

## Use a weighted regression to factor out effect of distance.
df <- data.frame(loop= xaxis[idx], cons= strn[idx], dist= dist[idx])

mod <- glm(cons~loop+dist, data=df, weights=n[idx])
mod

glm(cons~loop, data=df, weights=n[idx])$aic
glm(cons~dist, data=df, weights=n[idx])$aic
glm(cons~loop+dist+loop:dist, data=df, weights=n[idx])$aic

plot(xaxis[idx], dist[idx])


#################################################################
## Do more loops make a promtoer more likely to be conserved?

#getCex <- function(n) { y=0.075*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
getCex <- function(n) { y=0.01*n+0.3; y[y>3] <- 3; y[y<0.3] <- 0.3; y } 
n <- summary(as.factor(rowSums(loop[,5:6])))
nloops <- 0:10

## Correlation between promoter conservation and the number of chromatin loops to a promoter ...
loop_cons <- sapply(nloops, function(i) {fracConserved(tss[rowSums(loop[,5:6]) == i & tss$V5=="Prox_Stab",])})  ## Conservation of unlooped. 
n_s <- sapply(nloops, function(i) { sum(rowSums(loop[,5:6]) == i & tss$V5=="Prox_Stab") })
loop_cons_df <- data.frame(nloops= nloops, conservation= loop_cons, n= n_s, t0=rep(NA), sd=rep(NA))
cor.test(nloops, loop_cons)
indx <- !is.na(loop_cons)
corr(loop_cons_df[indx,1:2], w = n[indx]/sum(n[indx]))

## Weighted regression
fit_line <- lm(conservation~nloops, data= loop_cons_df[indx,1:2], w = n[indx]/sum(n[indx]))
fit_line

## Evaluate significance of this weighted correlation using a label swap.
## Note that this is not quire right, because boot is going to sample with replacement.  
## Real label-swap samples without replacemenet.
loop_swap <- boot(data= rowSums(loop[,5:6]), R=1000, statistic= function(a,i) {
	loop_cons <- sapply(nloops, function(x) {fracConserved(tss[a[i] == x & tss$V5=="Prox_Stab",])})  ## Conservation of unlooped. 
	loop_cons_df <- data.frame(nloops= nloops, conservation= loop_cons, n= summary(as.factor(rowSums(loop[,5:6]))), t0=rep(NA), sd=rep(NA))
	indx <- !is.na(loop_cons)
	corr(loop_cons_df[indx,1:2], w = n[indx]/sum(n[indx]))
}) ## Here, using boot to swap labels.
sum(loop_swap$t >= loop_swap$t0)/NROW(loop_swap$t)

## Use bootstrap to get a confidence for each conservation fraction.
for(x in nloops) {
  data <- tss[tss$V5=="Prox_Stab" & rowSums(loop[,5:6]) == x,]
  if(NROW(data) > 0) {
   bt <- boot(data= tss[tss$V5=="Prox_Stab" & rowSums(loop[,5:6]) == x,], R=1000, statistic= function(a, i) {fracConserved(a[i,])})
   loop_cons_df$t0[x+1] <- bt$t0
   loop_cons_df$sd[x+1] <- sd(bt$t, na.rm=TRUE)
  }
}

use <- !is.na(loop_cons_df$t0)
pdf("NumberOfLoops.Promoter.pdf")
  plot(nloops, loop_cons, pch=19, cex=3*getCex(n_s), xlab= "Number of loops to promoter", ylab= "Fraction conserved", ylim=c(0.5,1.0), xlim=c(-1,12))
  abline(fit_line)
  cd.barplot(loop_cons_df$t0[use], loop_cons_df$sd[use], as.character(nloops)[use], fill=TRUE, order=FALSE)
  plot(loop_swap)
dev.off()

cmpFracConserved(tss[rowSums(loop[,5:6]) >= 0 & rowSums(loop[,5:6]) <= 0,], tss[rowSums(loop[,5:6]) >= 3,], i=3) ## Promoters
cmpFracConserved(tss[rowSums(loop[,5:6]) >= 0 & rowSums(loop[,5:6]) <= 1,], tss[rowSums(loop[,5:6]) >= 3,], i=3) ## Promoters
cmpFracConserved(tss[rowSums(loop[,5:6]) >= 0 & rowSums(loop[,5:6]) <= 2,], tss[rowSums(loop[,5:6]) >= 3,], i=3) ## Promoters

## Do differences in expression explain loop correlation?
expr <- rowSums(rpkm_df[,2:9])
loop_expr <- sapply(nloops, function(i) {median(expr[rowSums(loop[,5:6]) == i & tss$V5=="Prox_Stab"])})  ## Conservation of unlooped. 
plot(loop_expr)

###################################################
## What about enhancers?

loop_cons_enh <- sapply(nloops, function(i) {fracConserved(tss[rowSums(loop[,5:6]) == i & tss$V5=="Dist_UnSt",])})  ## Conservation of unlooped. 
n_s <- sapply(nloops, function(i) { sum(rowSums(loop[,5:6]) == i & tss$V5=="Dist_UnSt") })
loop_cons_enh_df <- data.frame(nloops= nloops, conservation= loop_cons_enh, n= n_s, t0=rep(NA), sd=rep(NA))
cor.test(nloops, loop_cons_enh)
indx <- !is.na(loop_cons_enh)
corr(loop_cons_enh_df[indx,1:2], w = n[indx]/sum(n[indx]))

pdf("NumberOfLoops.Enhancer.pdf")
  plot(nloops, loop_cons_enh, pch=19, cex=3*getCex(n_s), xlab= "Number of loops to enhancer", ylab= "Fraction conserved")
dev.off()

cmpFracConserved(tss[(rowSums(loop[,5:6]) ==2 | rowSums(loop[,5:6]) ==3) & tss$V5=="Dist_UnSt",], tss[rowSums(loop[,5:6]) >= 4 & tss$V5=="Dist_UnSt",], i=2) ## Enhancers
cmpFracConserved(tss[(rowSums(loop[,5:6]) ==2 | rowSums(loop[,5:6]) ==3),], tss[rowSums(loop[,5:6]) >= 4,], i=2) ## Enhancers
cmpFracConserved(tss[(rowSums(loop[,5:6]) ==2 | rowSums(loop[,5:6]) ==3),], tss[rowSums(loop[,5:6]) >= 4,], i=3) ## Promoters

cmpFracConserved(tss[(rowSums(loop[,5:6]) ==1 | rowSums(loop[,5:6]) ==2 | rowSums(loop[,5:6]) ==3),], tss[rowSums(loop[,5:6]) >= 4,], i=2) ## Enhancers


####################################################
## How much are SE over-represented in the loop set?
sum(tss$V19==1 & rowSums(loop[,5:6])>0)/ sum(tss$V19==1) ## 49% of SE loop.
sum(tss$V19==0 & rowSums(loop[,5:6])>0)/ sum(tss$V19==0) ## 16% of all enhancers.

sum(tss$V19==1 & rowSums(loop[,5:6])>0)/ sum(rowSums(loop[,5:6])>0) ## 17% of loops part of a SE.
sum(tss$V19==1 & rowSums(loop[,5:6])==0)/ sum(rowSums(loop[,5:6])==0) ## 5% of non-looped.



########################################################################
## Compare looped and SEs in scatterplot. 

xaxis <- c(0, seq(3, 6, 0.2))

vect <- as.numeric(cut2(log(tss$V13, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
dist <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,])})
dist_se <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x & tss$V19 == 1,])})
dist_loop <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x & rowSums(loop[,5:6]) >  0,])})

## Size points based on 'n'
#getCex <- function(n) { y=0.0138888*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
n <- sapply(1:max(vect), function(x) {NROW(tss[vect == x,])})
n_se <- sapply(1:max(vect), function(x) {NROW(tss[vect == x & tss$V19 == 1,])})
n_loop <- sapply(1:max(vect), function(x) {NROW(tss[vect == x & rowSums(loop[,5:6]) >  0,])})


## Create barplot
xaxis <- c(0,seq(3, 6, 1))

vect <- as.numeric(cut2(log(tss$V13, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
for(i in 1:6) print(summary(tss$V13[vect == i]))

n_bootstrap_iter=10
b_dist <- boot(data= tss, R=n_bootstrap_iter, statistic= function(a, i) {sapply(1:max(vect), function(x) {fracConserved(a[i,][vect[i] == x,])})}) ## NOTE: Right order here?
b_dist_se <- boot(data= tss[tss$V19 == 1,], R=n_bootstrap_iter, statistic= function(a, i) {sapply(1:max(vect), function(x) {fracConserved(a[i,][vect[tss$V19 == 1][i] == x,])})})
b_dist_loop <- boot(data= tss[rowSums(loop[,5:6]) >  0,], R=n_bootstrap_iter, statistic= function(a, i) {sapply(1:max(vect), function(x) {fracConserved(a[i,][vect[rowSums(loop[,5:6]) > 0][i] == x,])})})

idx <- 1:5#3:5 #summary(cut2(log(tss$V13, 10), cuts=xaxis)) ## We want idx: 3 ([3.00,4.00)) - 5 ([5.00,6.00))
names<- paste(rep(c("[0]", "(0,1)", "[1-10)", "[10-100)", "[100-1000)"),3), c(rep("all", NROW(idx)), rep("se", NROW(idx)), rep("loop", NROW(idx))))
bars <- c(b_dist$t0[idx], b_dist_se$t0[idx], b_dist_loop$t0[idx])
serr <- c(sapply(1:NROW(b_dist$t0) , function(x) {sd(b_dist$t[,x], na.rm=TRUE)})[idx], 
	sapply(1:NROW(b_dist_se$t0) , function(x) {sd(b_dist_se$t[,x], na.rm=TRUE)})[idx],
	sapply(1:NROW(b_dist_loop$t0) , function(x) {sd(b_dist_loop$t[,x], na.rm=TRUE)})[idx])

ord <- c(seq(1,NROW(bars),5), seq(1,NROW(bars),5)+1, seq(1,NROW(bars),5)+2, seq(1,NROW(bars),5)+3, seq(1,NROW(bars),5)+4)

pdf("Distance_SE_loop.pdf")
 source("../lib/barplot.R")

 ## Barplot
 cd.barplot(bars[ord], serr[ord], names[ord], fill=TRUE, order=FALSE)
 cd.barplot(bars, serr, names, fill=TRUE, order=FALSE)
dev.off()




