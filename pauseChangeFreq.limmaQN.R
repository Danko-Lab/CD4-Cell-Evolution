##
## pauseChanges.limmaQN.R -- 
## GOAL: Compare official changes in pause sites to changes in protein-coding genes.
##
##

## Read pause sites.
ps <- read.table("../annotations/countpause.tsv")
ps <- ps[!is.na(ps[,10]) & rowSums(ps[,c(12:14,16:20)-1])>0,]

## Read gene bodies.
gb <- read.table("../annotations/countall.tsv")
gb <- gb[gb$V7 == "protein_coding" & !is.na(gb[,10]) & rowSums(gb[,c(12:14,16:20)-1])>0,]

###
## Fancy-shmancy selection of gb to be similar in coutns to ps.
source("../lib/normalizeSubsample.R")
gb_Sums <- rowSums(gb[,c(12:14,16:20)-1])
ps_Sums <- rowSums(ps[,c(12:14,16:20)-1])

idx <- norm.subsample(log(gb_Sums+1), log(ps_Sums+1), nsamp=20000, plot.cdf=TRUE)

## Quantile normalize each independently.
source("../lib/runLimmaQuantile.R")
conditions <- c("Human", "Human", "Human", rep("NHP",5))#"Chimp", "Chimp", "RM", "RM", "RM")

## Use Limma to ID changes.
hg <- runLimmaQuantile(gb[idx$s1,c(12:14,16:20)-1], conditions, gb[idx$s1,1:8], condA="Human", condB="NHP", q.cut=0.05)
hp <- runLimmaQuantile(ps[idx$s2,c(12:14,16:20)-1], conditions, ps[idx$s2,1:8], condA="Human", condB="NHP", q.cut=0.05)

## Use ...
pv <- 0.01
sum(hg$tab$P.Value < pv) ## Do NOT use adjusted p-values here.
sum(hp$tab$P.Value < pv)

boxplot(hg$tab$logFC, hp$tab$logFC, ylim=c(-3,3), names=c("body", "pause"))

## Plot out actual numbers...
source("../lib/densScatterplot.R")
cor.test(log(hg$tab$mean_Human), log(hg$tab$mean_NHP), method="spearman")
cor.test(log(hp$tab$mean_Human), log(hp$tab$mean_NHP), method="spearman")

par(mfrow=c(1,2))
densScatterplot(log(hg$tab$mean_Human), log(hg$tab$mean_NHP), main="body", xlab="Human", ylab="Mean NHP")
densScatterplot(log(hp$tab$mean_Human), log(hp$tab$mean_NHP), main="pause", xlab="Human", ylab="Mean NHP")


