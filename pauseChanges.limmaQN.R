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
hg <- runLimmaQuantile(gb[,c(12:14,16:20)-1], conditions, gb[,1:8], condA="Human", condB="NHP", q.cut=0.05)
hp <- runLimmaQuantile(ps[,c(12:14,16:20)-1], conditions, ps[,1:8], condA="Human", condB="NHP", q.cut=0.05)

#### NEW
## Get ps and gb into the same ...

body <-  hg$tab[hg$tab[,4] %in% hp$tab[,4],] # gb[gb[,4] %in% ps[,4],]
pause <- hp$tab[match(as.character(body[,4]), as.character(hp$tab[,4])),] # ps[match(as.character(body[,4]), as.character(ps[,4])),]
stopifnot(sum(body[,4] == as.character(pause[,4])) == NROW(pause)) ## SANTIY CHECK.



indx <- pause[,"Human 1"]>0 & pause[,"Chimp 4"]>0 & pause[,"R. Macaque 3"]>0 ## Unnecesarily conservative.
PI <- (((pause[,c(11:NCOL(pause))]+1)/pause[,"mapSize"]) / ((body[,c(11:NCOL(body))]+1)/body[,"mapSize"]))[indx,]


