##
## Fold-changes from U to PI -- checking that all species are on the same scale.
source("../lib/densScatterplot.R")
source("../lib/runLimmaQuantile.R")

#setwd("/local/storage/projects/NHP/annotations")
load("../annotations/fdr.RData")

## Index extent of changes.
isresp <- fdr_df$U2PIFDR_H < PVAL & fdr_df$U2PIFDR_C < PVAL & fdr_df$U2PIFDR_M < PVAL ## Clearly responding in any of the three.
summary(isresp)

pdf("foldchange.correlations.pdf")

 boxplot(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp])
 library(vioplot)
 vioplot(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp])

 cor(data.frame(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp]))

 ## Plots ... this is the supplementary figure.
 pairs(data.frame(fdr_df$U2PIFC_H[isresp], fdr_df$U2PIFC_C[isresp], fdr_df$U2PIFC_M[isresp]))

 cor.test(fdr_df$U2PIFC_H[isresp], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp])
 plot(fdr_df$U2PIFC_H[isresp], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp], pch=19); abline(0,1)
 densScatterplot(fdr_df$U2PIFC_H[isresp], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp], xlab="Human", ylab="Mean non-human primate"); abline(0,1)

dev.off()

## Now fit a new model using LIMMA that includes the interaction terms.
condition <- species
condition[species == "H"] <- "Human"
condition[species != "H"] <- "NHP"

count.dat <- counts

condA <- "conditionNHP.treatmentU"
condB <- "treatmentU"

## use LIMMA
require(limma)
design <- model.matrix(~0+condition+treatment+condition*treatment)
colnames(design)[4] <- "conditionNHP.treatmentU"

counts.log.dat=as.matrix(log2(count.dat+1)) ## CGD: Cast appears to be required in Limma 3.18.13
counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')
dat=counts.log.norm.dat
fit=lmFit(dat,design)

#################################################################
## Not sure of the best way to write contrast matrix ############
contrast.matrix <- makeContrasts(contrasts= paste(condA,"-",condB), levels=design) ## Contrast btwn. interaction term and treatment term.
#contrast.matrix <- makeContrasts(contrasts= paste(condA), levels=design) ## This might be more what we want!? Sig. interaction term (either direction).

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res=decideTests(fit2,p.value=PVAL)
tab<-topTable(fit2, adjust = "BH", number=nrow(fit2), sort.by='none')

## Select genes that are differentially regulated in human.
ishumspec <- ((tab$adj.P.Val < PVAL))# & 
#	(abs(fdr_df$U2PIFC_H-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])) > 1) & 
#	(abs(rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])) < 1))
summary(ishumspec)

pdf("Human.Differences_In_Induction.pdf")
 ## Sanity checks...
 hist((fdr_df$U2PIFC_H-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])))
 hist((fdr_df$U2PIFC_H-rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")]))[ishumspec])

 ## Plots.
 plot(fdr_df$U2PIFC_H[ishumspec], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumspec], pch=19); abline(0,1)
 densScatterplot(fdr_df$U2PIFC_H[ishumspec], rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumspec], xlab="Human", ylab="Mean non-human primate")
 abline(0,1)

 plot(fdr_df$U2PIFC_H[isresp] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[isresp], pch=19, xlab="Mean non-human primate", ylab="Human")
 points(fdr_df$U2PIFC_H[ishumspec] ~ rowMeans(fdr_df[,c("U2PIFC_C", "U2PIFC_M")])[ishumspec], pch=19, col="red")
 abline(h=0); abline(v=0)

 ## Add labels.
dev.off()

## Which genes are they?


