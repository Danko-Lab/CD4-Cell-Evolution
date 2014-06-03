#
# Collect systematica differences between pause sites at 
# orthologous promoters.  Distinguish using PCA.

dist <- 50  ## Differences within an individual pause site.
dist <- 200 ## Changes in which pause wite is used.

## Get the 3' patterns at TSS.  Fiter for complete mappability.
ps <- read.table("../annotations/pausesites.tsv")
cen <- ps[,3]-ps[,2]
rps <- ps
rps[,2] <- cen-dist
rps[,3] <- cen+dist
#ps.n <- unique(ps[,4])
#sapply

require(bigWig)

bwp_f <- dir(path="../AllData/", pattern="-U.bed.gz_plus.hg19.bw", full.names=TRUE)
bwm_f <- dir(path="../AllData/", pattern="-U.bed.gz_minus.hg19.bw", full.names=TRUE)

bwp <- lapply(bwp_f, load.bigWig)
bwm <- lapply(bwm_f, load.bigWig)

getPsPc <- function(x) {
 ns <- NROW(bwp)

 ## For each gene, construct a matrix of ... reads x sample.
 ## https://github.com/andrelmartins/bigWig/blob/master/bigWig/R/basepair.query.R
 mat <- lapply(1:ns, function(i) {
  bed6.step.bpQuery.bigWig(bwp[[i]], bwm[[i]], rps, 1)
 }) ## Change this to one matrix per condition.

 ## PCA. Filter for pause sites that don't change in quantity.
 

 ## Get pause sites whose PCA is significantly different across species.  
 ##  Quantify withing-species divergence; -> between species divergence.

}


## Quantify frequency of shifts.  Look at types of shifts.  

## Correlate shifts with specific patterns of TXN.
