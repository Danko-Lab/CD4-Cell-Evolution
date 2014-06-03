##
## pauseChanges.limmaQN.R -- 
## GOAL: Compare official changes in pause sites to changes in protein-coding genes.
##
##
indxH <- c(11:13) #21:23 #12:13 -> Remove H1 (outlier for pause)
indxC <- c(15:16) #25:26
indxM <- c(17:19) #28:29
indxHCM <- c(indxH, indxC, indxM)

ps <- read.table("../annotations/countpause.tsv") ## Read pause sites.
gb <- read.table("../annotations/countall.tsv") ## Read gene bodies.

makeNormPlots <- function(indx1, indx2, name1, name2, pv=0.01, nsamp=50000) { 
 indx <- c(indx1, indx2)
 ps <- ps[!is.na(ps[,10]),]# & rowSums(ps[,indx])>0,]
 gb <- gb[gb$V7 == "protein_coding" & !is.na(gb[,10]),]# & rowSums(gb[,indx])>0,]

 ###
 ## Fancy-shmancy selection of gb to be similar in coutns to ps.
 source("../lib/normalizeSubsample.R")
 gb_Sums <- rowSums(gb[,indx1])
 ps_Sums <- rowSums(ps[,indx1])

 idx <- norm.subsample(log(gb_Sums+1), log(ps_Sums+1), nsamp=nsamp, plot.cdf=TRUE)
 
 ## Quantile normalize each independently.
 source("../lib/runLimmaQuantile.R")
 conditions <- c(rep(name1, NROW(indx1)), rep(name2, NROW(indx2)))

 ## Use Limma to ID changes.
 hg <- runLimmaQuantile(gb[idx$s1,indx], conditions, gb[idx$s1,1:8], condA=name1, condB=name2, q.cut=0.05)
 hp <- runLimmaQuantile(ps[idx$s2,indx], conditions, ps[idx$s2,1:8], condA=name1, condB=name2, q.cut=0.05)
 
 ## Use ...
 print("Number of Changes:")
 print(paste("body:", sum(hg$tab$P.Value < pv), "pause:", sum(hp$tab$P.Value < pv))) ## Must not use adjusted p-value
 print(fisher.test(data.frame(c(sum(hg$tab$P.Value < pv), nsamp), c(sum(hp$tab$P.Value < pv), nsamp))))

 boxplot(hg$tab$logFC, hp$tab$logFC, ylim=c(-3,3), names=c("body", "pause"))
 
 ## Plot out actual numbers...
 source("../lib/densScatterplot.R")
 print("Correlations: (Spearman)")
 print(paste("Body:",cor(log(hg$tab[,12]), log(hg$tab[,13]), method="spearman")))
 print(paste("Pause:",cor(log(hp$tab[,12]), log(hp$tab[,13]), method="spearman")))

 par(mfrow=c(1,2))
 densScatterplot(log(hg$tab[,12]), log(hg$tab[,13]), main="body", xlab=name1, ylab=name2)
 densScatterplot(log(hp$tab[,12]), log(hp$tab[,13]), main="pause", xlab=name1, ylab=name2)
}

makeNormPlots(indxH, indxC, "Human", "Chimp")
makeNormPlots(indxC, indxH, "Chimp", "Human")


makeNormPlots(indxH, indxM, "Human", "Macaque")
makeNormPlots(indxM, indxH, "Macaque", "Human")

makeNormPlots(indxC, indxM, "Chimp", "Macaque")
makeNormPlots(indxM, indxC, "Macaque", "Chimp")



makeNormPlots(indxH, c(indxC, indxM), "Human", "NHP")

