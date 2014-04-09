##
## pauseChanges.limmaQN.R -- 
## GOAL: Compare official changes in pause sites to changes in protein-coding genes.
##
##
indxH <- c(11:13) #21:23 #12:13 -> Remove H1 (outlier for pause)
indxC <- c(15:16) #25:26
indxM <- c(17:19) #28:29
indxHCM <- c(indxH, indxC, indxM)

## Read pause sites.
ps <- read.table("../annotations/countpause.tsv")
ps <- ps[!is.na(ps[,10]),]

## Read gene bodies.
gb <- read.table("../annotations/countall.tsv")
gb <- gb[gb$V7 == "protein_coding" & !is.na(gb[,10]),]

makePBPlot <- function(indx1, indx2, name1, name2, hv=3, lv=hv) {
 indx <- c(indx1, indx2)

 ## Quantile normalize each independently.
 source("../lib/runLimmaQuantile.R")
 conditions <- c(rep(name1, NROW(indx1)), rep(name2, NROW(indx2)))
 
 ## Use Limma to ID changes.
 hg <- runLimmaQuantile(gb[,indx], conditions, gb[,1:8], condA=name1, condB=name2, q.cut=0.05)
 hp <- runLimmaQuantile(ps[,indx], conditions, ps[,1:8], condA=name1, condB=name2, q.cut=0.05)
 
 ## Get ps and gb into the same indices, for the same genes ...
 body <-  hg$tab[hg$tab[,4] %in% hp$tab[,4],] 
 pause <- hp$tab[match(as.character(body[,4]), as.character(hp$tab[,4])),] 
 stopifnot(sum(body[,4] == as.character(pause[,4])) == NROW(pause)) ## SANTIY CHECK.
 
 source("../lib/densScatterplot.R")
 #boxplot(body$logFC, pause$logFC, ylim=c(-3,3))
 cc <- cor(body$logFC, pause$logFC, method="pearson")
 print(paste("Correlation: ", cc, "Fraction of Varience: ", cc*cc))
 #plot(body$logFC, pause$logFC)
 densScatterplot(body$logFC, pause$logFC, xlab="Body", ylab="Pause", main=paste(name1,"vs",name2))
 abline(v=c(hv,-hv), h=c(hv,-hv))
 
 ## Human changes in pausing >3
 ## # Pause / # Body
 pc <- sum(abs(pause$logFC)>hv & abs(body$logFC)<lv)
 bc <- sum(abs(body$logFC)>hv & abs(pause$logFC)<lv)
 pv <- fisher.test(data.frame(c(pc, NROW(pause)), c(bc, NROW(body))))$p.value
 print(paste("Fraction Pause_Q(2,8)/ Body_Q(4,6):", pc,"/",bc,"=", pc/ bc, "P=", pv))
}

makePBPlot(indxH, indxC, "Human", "Chimp", hv=5, lv=3)
makePBPlot(indxH, indxM, "Human", "Macaque", hv=5, lv=3)
makePBPlot(indxC, indxM, "Chimp", "Macaque", hv=5, lv=3)



