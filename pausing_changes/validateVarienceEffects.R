#
# Designed to validate the causes of varience
#
# Specifically: 
#  * Does read density affect varience?
#  * Does absolute number of reads affect varience?
#
source("../lib/densScatterplot.R")

indxH <- c(11:13) #21:23 #12:13 -> Remove H1 (outlier for pause)
#indxC <- c(15:16) #25:26
#indxM <- c(17:19) #28:29

#ps <- read.table("../annotations/countpause.tsv") ## Read pause sites.
gb <- read.table("../annotations/countall.tsv") ## Read gene bodies.
gb <- gb[!is.na(gb[,10]),]

rowVar <- function(mat) {
 sapply(1:NROW(mat), function(x) {var(as.double(mat[x,]))})
}

## Choose a range of lengths, and plot varience between reps as a function of read depth.
genelength <- gb[,3]-gb[,2]
hist(genelength)
useL <- genelength > 55000
exprL <- rowSums(gb[useL,indxH])
varsL <- rowVar(gb[useL,indxH])
fcvL  <- rowVar(gb[useL,indxH]/exprL/NROW(indxH)) ## Varience in fold-change from the mean.
fcL   <- gb[useL,indxH[1]]/gb[useL,indxH[2]]

## Choose a range of expression, and plot varience between reps as a function of gene length.
expr <- log(rowSums(gb[,indxH])+1)
hist(expr)
useE <- expr > 7.5 & expr < 8
lensE <- genelength[useE]
varsE <- rowVar(gb[useE,indxH])
fcvE  <- rowVar(gb[useE,indxH]/expr[useE]/NROW(indxH)) ## Varience in fold-change from the mean.
fcE   <- gb[useE,indxH[1]]/gb[useE,indxH[2]]

pdf("VariencePropertiesWithReadCountsAndLength.pdf")
par(mfrow=c(1,2))

## Plot raw Fold-Change.
cor.test(exprL, log(fcL), method="spearman")
densScatterplot(log(exprL+1), log(fcL), main="constant length, varying expression, fold change", xlab="expression", ylab="fold-change")

cor.test(lensE, fcE, method="spearman")
densScatterplot(lensE, log(fcE+1), main="constant expr, varying length, fold change", xlab="length", ylab="fold-change")


## Varience in fold-change.
cor.test(exprL, fcvL, method="spearman")
densScatterplot(log(exprL+1), log(fcvL+1), main="constant length, varying expression, fold change varience", xlab="expression", ylab="var in fold-change")

cor.test(lensE, fcvE, method="spearman")
densScatterplot(lensE, log(fcvE+1), main="constant expr, varying length, fold change varience", xlab="length", ylab="var in fold-change")


## Plot raw varience.
cor.test(exprL, varsL, method="spearman")
densScatterplot(log(exprL+1), log(varsL+1), main="constant length, varying expression, absolute varience", xlab="expression", ylab="varience")

cor.test(lensE, varsE, method="spearman")
densScatterplot(lensE, log(varsE+1), main="constant expr, varying length", xlab="length", ylab="varience")

dev.off()


