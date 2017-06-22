## Compares induction in the simple species-specific hg19 analysis


## Read files.
hs <- read.table("results/human-changed.genes.tsv"); hs$V4 <- paste(hs$V1,hs$V2,hs$V3,hs$V4,sep="_")
pt <- read.table("results/chimp-changed.genes.tsv"); pt$V4 <- paste(pt$V1,pt$V2,pt$V3,pt$V4,sep="_")
rm <- read.table("results/rhesus-changed.genes.tsv");rm$V4 <- paste(rm$V1,rm$V2,rm$V3,rm$V4,sep="_")

## Get same order.
pt <- pt[match(hs$V4, pt$V4),]
rm <- rm[match(hs$V4, rm$V4),]

## Remove NAs.
skip <-is.na(hs$V15) | is.na(pt$V15) | is.na(rm$V15)
sum(skip); sum(!skip)

hs<- hs[!skip,]; pt<- pt[!skip,]; rm<- rm[!skip,]

## Get correlations.
cc<- cor(cbind(hs$V15, pt$V15, rm$V15))
cc
pairs(cbind(hs$V15, pt$V15, rm$V15))

plot(c(0, 12, 25+25-12), cc[,1])

## Get number of genes.
PVAL <- 0.01
cngOne <- hs$V19 < PVAL | pt$V19 < PVAL | rm$V19 < PVAL; cngOne <- cngOne & !is.na(cngOne)
cngAll <- hs$V19 < PVAL & pt$V19 < PVAL & rm$V19 < PVAL; cngAll <- cngAll & !is.na(cngAll)

sum(cngOne); sum(cngAll)

NROW(unique(hs$V7[cngAll])) ## Yep. 3,157.  A few that change directions.  Not many.

cor(cbind(hs$V15[cngAll], pt$V15[cngAll], rm$V15[cngAll]))
pairs(cbind(hs$V15[cngAll], pt$V15[cngAll], rm$V15[cngAll]))

pdf("CompareFoldChanges.pdf")
 source("../lib/densScatterplot.R")
 densScatterplot(hs$V15[cngAll], pt$V15[cngAll], xlab="Log2 fold-change Human", ylab="Log2 fold-change Chimp")
 densScatterplot(hs$V15[cngAll], rm$V15[cngAll], xlab="Log2 fold-change Human", ylab="Log2 fold-change Rhesus")
 densScatterplot(pt$V15[cngAll], rm$V15[cngAll], xlab="Log2 fold-change Chimp", ylab="Log2 fold-change Rhesus")
dev.off() 

cor(cbind(hs$V15[cngOne], pt$V15[cngOne], rm$V15[cngOne]))
pairs(cbind(hs$V15[cngOne], pt$V15[cngOne], rm$V15[cngOne]))


