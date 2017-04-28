##
## Compare dREG scores across species.
##
##

source("../lib/densScatterplot.R")

a <- read.table("HCM-U-PI.dREG-tss-clusters.tsv")
for(i in 7:12) a[is.na(a[,i]),i] <- 0 ## Assign 'NA' scores to 0.

siteth<- 0.25
upper <- 0.4
lower <- 0.05

printSummaryData <- function(v1, v2) {
 b1   <-sum(v1 < lower & v2 > upper)
 b2   <-sum(v2 < lower & v1 > upper)
 inall<-sum(v1 > upper & v2 > upper)
 sum(v1 < lower & v2 < lower)
 all<- sum(v1 > siteth & v2 > siteth) #sum(b1, b2, inall)

 print(paste(b1, b2, inall, b1/(all+b1+b2), b2/(all+b1+b2), (b1+b2)/(all+b1+b2)))
}

pdf("dREG.DensityScatterplot.pdf")
densScatterplot(a$V7, a$V8, xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
abline(h=upper); abline(h=lower); abline(v=upper); abline(v=lower)
printSummaryData(a$V7, a$V8)

densScatterplot(a$V7, a$V9, xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores") 
abline(h=upper); abline(h=lower); abline(v=upper); abline(v=lower)
printSummaryData(a$V7, a$V9)

densScatterplot(a$V8, a$V9, xlab="Chimpanzee Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores")
abline(h=upper); abline(h=lower); abline(v=upper); abline(v=lower)
printSummaryData(a$V8, a$V9)

densScatterplot(a$V7, a$V10, xlab="Human Unt.", ylab="Human P/I", main="dREG Scores") 
abline(h=upper); abline(h=lower); abline(v=upper); abline(v=lower)
printSummaryData(a$V7, a$V10)
printSummaryData(a$V8, a$V11)
printSummaryData(a$V9, a$V12)

dev.off()
