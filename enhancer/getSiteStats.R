## This analysis focuses on lineage-specific changes that are currently active in hg19.
load("../annotations/fdr.RData")

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

## Look at species...
## Count types of elements ...
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

## Do SE change less frequently?
testSE <- function(i) {
 SEtot <- SEtotal[i]-SEone2m[i]
 SEcon <- SEtotal[i]-SEone2m[i]-SEunmap[i]-SEchang[i]-SElccng[i]
 print(paste("SE:", SEcon/SEtot))

 tot <- total[i]-one2m[i]
 con <- total[i]-one2m[i]-unmap[i]-chang[i]-lccng[i]
 print(paste("ALL:", con/tot))

 print(fisher.test(data.frame(c(SEtot, SEcon), c(tot, con))))
}

testSE(2) # enhancers
testSE(3) # promoters

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

for(i in 7:12) tss[is.na(tss[,i]),i] <- 0 ## Assign 'NA' scores to 0.
source("../lib/densScatterplot.R")

pdf("dreg.scatterplots.pdf")
 max_hc <- rowMax(tss[,c(7,8)])
 densScatterplot(tss$V7, tss$V8, xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
 densScatterplot(tss$V7[tss$V19 == "Prox_Stab" & max_hc>0.7], tss$V8[tss$V19 == "Prox_Stab" & max_hc>0.7], xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
 densScatterplot(tss$V7[tss$V19 == "Dist_UnSt" & max_hc>0.7], tss$V8[tss$V19 == "Dist_UnSt" & max_hc>0.7], xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")

 max_hm <- rowMax(tss[,c(7,9)])
 densScatterplot(tss$V7[max_hm>0.7], tss$V9[max_hm>0.7], xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")
 densScatterplot(tss$V7[tss$V19 == "Prox_Stab" & max_hm>0.7], tss$V9[tss$V19 == "Prox_Stab" & max_hm>0.7], xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")
 densScatterplot(tss$V7[tss$V19 == "Dist_UnSt" & max_hm>0.7], tss$V9[tss$V19 == "Dist_UnSt" & max_hm>0.7], xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")

 max_up <- rowMax(tss[,c(7,10)])
 densScatterplot(tss$V7, tss$V10, xlab="Human Unt.", ylab="Human PI", main="dREG Scores")
 densScatterplot(tss$V7[tss$V19 == "Prox_Stab" & max_up>0.7], tss$V10[tss$V19 == "Prox_Stab" & max_up>0.7], xlab="Human Unt.", ylab="Human PI", main="dREG Scores")
 densScatterplot(tss$V7[tss$V19 == "Dist_UnSt" & max_up>0.7], tss$V10[tss$V19 == "Dist_UnSt" & max_up>0.7], xlab="Human Unt.", ylab="Human PI", main="dREG Scores")
dev.off()

## Compare distance to conservation...
fact <- cut(tss$V13,c(-Inf,1,5,10,20,50,100,Inf)*1000)
sapply(fact, function(i) { sum(tss$HumanFDR[fact == i] < 0.01)/sum(fact == i) })
boxplot()


########################################################
## Evolutionary conservation->distance from genes.

fracConserved <- function(tss, ret="C") { ## defaults to ... 
 total <- NROW(tss)
 one2m <- sum(tss$V20 > 0) ## Possible 1:many orthology
 unmap <- sum(tss$V20 == 0 & is.na(tss$mapSize)) ## INDEL
 lccng <- sum((tss$V7 < 0.1 | tss$V8 < 0.1 | tss$V9 < 0.1) & tss$V20 == 0 & !is.na(tss$mapSize)) # 'Low-confidence'
 chang <- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05 & (tss$V7 > 0.7 & tss$V8 > 0.7 & tss$V9 > 0.7)) # 'High-confidence'

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
  cng <- chang

  return(cng/tot)
 }
}
fracConserved(tss[tss$V5 == "Dist_UnSt",])
fracConserved(tss[tss$V5 == "Prox_Stab",])

#vect <- as.numeric(cut2(log(tss$V13), g=100)); summary(as.factor(vect))
xaxis <- c(0, seq(3, 6, 0.02))
vect <- as.numeric(cut2(log(tss$V13, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
dist <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,])})
glc  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "GL")})
cng  <- sapply(1:max(vect), function(x) {fracConserved(tss[vect == x,], "CNG")})

pdf("ChangeOverDistance.pdf")
 plot(10^xaxis[1:NROW(dist)], dist, type="p", xlab="Distance from TSS [bp]", ylab="Fraction conserved", xlim=10^c(3, 6), log="x", pch=19, cex=1.5)
 plot(10^xaxis[1:NROW(dist)], glc, type="p", xlab="Distance from TSS [bp]", ylab="Fraction gain/ loss", xlim=10^c(3, 6), log="x", pch=19, cex=1.5)
 plot(10^xaxis[1:NROW(dist)], cng, type="p", xlab="Distance from TSS [bp]", ylab="Fraction change", xlim=10^c(3, 6), log="x", pch=19, cex=1.5)
dev.off()

## Now based on expression level (max across species).
exp_max <- rowMax(rpkm_df[,2:9])[match(tss$V4, tss_aln$name)]
exp_vect <- as.numeric(cut2(exp_max, g=100))
exp_vect[is.na(exp_vect)] <- 0
conserved <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x,])})
gainloss <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x,], "GL")})
change   <- sapply(1:max(exp_vect), function(x) {fracConserved(tss[exp_vect == x,], "CNG")})

plot(conserved, type="b")
plot(gainloss,  type="b")
plot(change,    type="b")



