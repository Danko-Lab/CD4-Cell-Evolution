#
# getCounts.R -- Uses bigWig package to get counts for all genes.
#
# Usage: R --no-save --args bed_file out_file big_wig_plus big_wig_minus < getCounts.R
#

require(bigWig)

## Process command arguments.
args <- commandArgs(trailingOnly=TRUE)
bed_file <- args[1]
out_file <- args[2]
bw_plus_file <- args[3]
bw_minus_file <- args[4]

## Get counts using bigWig
bed <- read.table(bed_file)
bw_plus <- load.bigWig(bw_plus_file)
bw_minus<- load.bigWig(bw_minus_file)

counts <- bed6.region.bpQuery.bigWig(bw_plus, bw_minus, bed, abs.value=TRUE) #bedQuery.bigWig(bed, bw_plus, bw_minus, gapValue=0)
#counts <- abs(counts)  ## Reverse strand for minus.

bed[,5] <- counts
write.table(bed, out_file, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
