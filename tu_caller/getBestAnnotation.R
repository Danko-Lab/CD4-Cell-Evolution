#
# Gets the longest isoform of an annotation, such that
# the 5' end has a dREG site, and the 3' end has substantial
#  expression.

require(bigWig)

a <- read.table("annot/gc.tss.tsv.gz")

## get 3' ends.
startWindow <- 1000 #bp (How close does the TSS have to be?!)
endWindow <- 1000 #bp (What's the window size to measure expression in?!)
bed <- a[,1:6]
bed[bed[,6] == "+",2] <- bed[bed[,6] == "+",3] - endWindow
bed[bed[,6] == "-",3] <- bed[bed[,6] == "-",2] + endWindow

## Get data from a bigWig.
## Returns vector of counts ...
getCounts <- function(bwPlus, bwMinus, bed6=bed) {
  bwp <- load.bigWig(bwPlus)
  bwm <- load.bigWig(bwMinus)

 bed6.region.bpQuery.bigWig(bw.plus= bwp, bw.minus= bwm, bed6= bed6)
}

Hc <- getCounts("/local/storage/projects/NHP/AllData/All_Merge/H-U_plus.bw", "/local/storage/projects/NHP/AllData/All_Merge/H-U_minus.bw") ## Human
Cc <- getCounts("/local/storage/projects/NHP/AllData/All_Merge/C-U_plus.hg19.bw", "/local/storage/projects/NHP/AllData/All_Merge/C-U_minus.hg19.bw") ## Chimp
Mc <- getCounts("/local/storage/projects/NHP/AllData/All_Merge/M-U_plus.hg19.bw", "/local/storage/projects/NHP/AllData/All_Merge/M-U_minus.hg19.bw") ## Rhesus

## Compute the threshold number of reads in this window.
e_rpkb <- 0.04*1e8/10751533*1000/endWindow ## 0.04 taken from Core Waterfall Lis, 10751533 was their library size.  1e8 is (roughly) the NHP library sizes.
bgrc <- ppois(1:100, (e_rpkb), lower.tail=FALSE)*NROW(bed)
nrth <- min(which(bgrc < 0.01)) ## Minimum number of reads with meets the Bonferroni 0.01 bg corrected significance level (n=6).

## Get those with significant expression in the 3' end.
expr3p <- abs(Hc)>nrth | abs(Cc)>nrth | abs(Mc)>nrth

## Get those with a TSS at 5' end.
tss5p <- a$V14<startWindow

## TU that can be added back!
sum(tss5p & expr3p) # 837 (4-17-2014)

## Now get the longest incarnation of a TU that can be added back.
addBack <- a[tss5p & expr3p, ]
addBack$length <- addBack[,3]-addBack[,2]
addBack <- addBack[order(addBack$length, decreasing=TRUE),]

source("/local/storage/projects/NHP/lib/CountUnique.R")
addBack <- addBack[is.firstUnique(addBack),]
NROW(addBack)
write.table(addBack, "annot/addBack.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)



