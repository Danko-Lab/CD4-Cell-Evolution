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


## Quantile normalize each independently.
source("../lib/runLimmaQuantile.R")
conditions <- c("Human", "Human", "Human", rep("NHP",5))#"Chimp", "Chimp", "RM", "RM", "RM")

## Use Limma to ID changes.
hg <- runLimmaQuantile(gb[,c(12:14,16:20)-1], conditions, gb[,1:8], condA="Human", condB="NHP")
hp <- runLimmaQuantile(ps[,c(12:14,16:20)-1], conditions, ps[,1:8], condA="Human", condB="NHP")



