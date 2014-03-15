##
## Compare dREG scores across species.
##
##

source("../lib/densScatterplot.R")

a <- read.table("HCM-U-PI.dREG-tss-clusters.tsv")
for(i in 7:12) a[is.na(a[,i]),i] <- 0 ## Assign 'NA' scores to 0.


densScatterplot(a$V7, a$V8, xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
densScatterplot(a$V7, a$V9, xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores") 

densScatterplot(a$V7, a$V10, xlab="Human Unt.", ylab="Human P/I", main="dREG Scores") 

