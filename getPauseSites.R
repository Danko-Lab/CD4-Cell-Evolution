#
# getPauseSites.R -- Returns position(argmax(counts in 100bp in the specified interval)).
#
# Usage: R --no-save --args bed_file out_file big_wig_plus big_wig_minus < getPauseSites.R
#
# TEST: cat tsssites.tsv | grep "TSS" > tmp.tmp
#       R --no-save --args tmp.tmp tmp.out.tmp ../AllData/H1-U_plus.bw ../AllData/H1-U_minus.bw < getPauseSites.R

require(bigWig)

size=100
overlap=5
nwindows <- size/overlap

## Process command arguments.
args <- commandArgs(trailingOnly=TRUE)
bed_file <- args[1]
out_file <- args[2]
bw_plus_file <- args[3]
bw_minus_file <- args[4]

## Load the bed and bigWigs. 
bed <- read.table(bed_file)
bw_plus <- load.bigWig(bw_plus_file)
bw_minus<- load.bigWig(bw_minus_file)

## Get the start position of 100bp window with the max counts inside each dREG site.
startPlus <- integer()
startMinus <- integer()

for(i in 1:NROW(bed)) {
  pad <- overlap- ((bed[i,3]-bed[i,2]) %% overlap) ## Padding, to make the window a multiple of the step size.
  bases  <- seq(bed[i,2], bed[i,3]+pad, overlap)

  countsp <- abs(step.bpQuery.bigWig(bw_plus, chrom=bed[i,1], start= bed[i,2], end= bed[i,3]+pad, step=overlap))
  countsp_SUM <- sapply(1:(length(countsp)-nwindows), function(x) {sum(countsp[x:(x+nwindows)])})
  stopifnot(NROW(bases)-1-nwindows == NROW(countsp_SUM)) ## SANTIY CHECK

  countsm <- abs(step.bpQuery.bigWig(bw_minus, chrom=bed[i,1], start= bed[i,2], end= bed[i,3]+pad, step=overlap))
  countsm_SUM <- sapply(1:(length(countsm)-nwindows), function(x) {sum(countsm[x:(x+nwindows)])})
  stopifnot(NROW(bases)-1-nwindows == NROW(countsm_SUM)) ## SANTIY CHECK

  startPlus  <- c(startPlus, bases[1]+(which.max(countsp_SUM)-1)*overlap)
  startMinus <- c(startMinus, bases[1]+(which.max(countsm_SUM)-1)*overlap)
}

data <- rbind(data.frame(chrom= bed[,1], start= startPlus, end= startPlus+size, name= bed[,4], score= bed[,5], strand= "+"),
		data.frame(chrom= bed[,1], start= startMinus, end= startMinus+size, name= bed[,4], score= bed[,5], strand= "-"))

options("scipen"=100, "digits"=4)
write.table(data, pipe(paste("sort-bed - >", out_file)), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
