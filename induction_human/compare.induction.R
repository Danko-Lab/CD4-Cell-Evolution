## Compares induction in the simple species-specific hg19 analysis


## Read files.
hs <- read.table("results/human-changed.genes.tsv"); hs$V4 <- paste(hs$V1,hs$V2,hs$V3,hs$V4,sep="_")
pt <- read.table("results/chimp-changed.genes.tsv"); pt$V4 <- paste(pt$V1,pt$V2,pt$V3,pt$V4,sep="_")
rm <- read.table("results/rhesus-changed.genes.tsv");rm$V4 <- paste(rm$V1,rm$V2,rm$V3,rm$V4,sep="_")

## Get same order.
pt <- pt[match(hs$V4, pt$V4),]
rm <- rm[match(hs$V4, rm$V4),]

## Remove NAs.
skip <-is.na(hs$V9) | is.na(pt$V9) | is.na(rm$V9)
sum(skip); sum(!skip)

hs<- hs[!skip,]; pt<- pt[!skip,]; rm<- rm[!skip,]

## Get correlations.
cc<- cor(cbind(hs$V9, pt$V9, rm$V9))
cc
pairs(cbind(hs$V9, pt$V9, rm$V9))

plot(c(0, 6, 25+25-6), cc[,1])

## Get number of genes.
PVAL <- 0.01
cngOne <- hs$V13 < PVAL | pt$V13 < PVAL | rm$V13 < PVAL; cngOne <- cngOne & !is.na(cngOne)
cngAll <- hs$V13 < PVAL & pt$V13 < PVAL & rm$V13 < PVAL; cngAll <- cngAll & !is.na(cngAll)

sum(cngOne); sum(cngAll)

NROW(unique(hs$V7[cngAll])) ## Yep. 2,953.  A few that change directions.  Not many.

cor(cbind(hs$V9[cngAll], pt$V9[cngAll], rm$V9[cngAll]))
pairs(cbind(hs$V9[cngAll], pt$V9[cngAll], rm$V9[cngAll]))

pdf("CompareFoldChanges.pdf")
 source("../lib/densScatterplot.R")
 densScatterplot(hs$V9[cngAll], pt$V9[cngAll], xlab="Log2 fold-change Human", ylab="Log2 fold-change Chimp")
 densScatterplot(hs$V9[cngAll], rm$V9[cngAll], xlab="Log2 fold-change Human", ylab="Log2 fold-change Rhesus")
 densScatterplot(pt$V9[cngAll], rm$V9[cngAll], xlab="Log2 fold-change Chimp", ylab="Log2 fold-change Rhesus")
dev.off() 

cor(cbind(hs$V9[cngOne], pt$V9[cngOne], rm$V9[cngOne]))
pairs(cbind(hs$V9[cngOne], pt$V9[cngOne], rm$V9[cngOne]))


