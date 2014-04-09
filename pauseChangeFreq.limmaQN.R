##
## pauseChanges.limmaQN.R -- 
## GOAL: Compare official changes in pause sites to changes in protein-coding genes.
##
##

## Read pause sites.
ps <- read.table("../annotations/countpause.tsv")
ps <- ps[!is.na(ps[,10]),]

## Read gene bodies.
gb <- read.table("../annotations/countall.tsv")
gb <- gb[gb$V7 == "protein_coding" & !is.na(gb[,10]),]

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
 

q("no")

## Get ps and gb into the same indices, for the same genes ...
body <-  hg$tab[hg$tab[,4] %in% hp$tab[,4],] 
pause <- hp$tab[match(as.character(body[,4]), as.character(hp$tab[,4])),] 
stopifnot(sum(body[,4] == as.character(pause[,4])) == NROW(pause)) ## SANTIY CHECK.

source("../lib/densScatterplot.R")
#boxplot(body$logFC, pause$logFC, ylim=c(-3,3))
cor.test(body$logFC, pause$logFC, method="spearman")
#plot(body$logFC, pause$logFC)
densScatterplot(body$logFC, pause$logFC, xlab="Body", ylab="Pause")
abline(v=c(3,-3), h=c(3,-3))

## Human changes in pausing >3
hv <- 3; lv <- 0.5
## # Pause / # Body
sum(abs(pause$logFC)>hv & abs(body$logFC)<lv)/ sum(abs(body$logFC)>hv & abs(pause$logFC)<lv)




