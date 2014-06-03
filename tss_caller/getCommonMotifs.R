#
# getCommonMotifs.R -- Writes a file with commonly occuring motifs. 
#
# Usage: R --no-save --args in_bed_file min_occurences out_motif_name_file < getCommonMotifs.R
#


## Process command arguments.
args <- commandArgs(trailingOnly=TRUE)
inf <- args[1]
thr <- as.integer(args[2])
outf<- args[3]

ff <- read.table(inf)
commonMotifs <- names(summary(ff$V4))[grep("(Other)", names(summary(ff$V4))[summary(ff$V4) > thr], invert=TRUE)]
write.table(commonMotifs, outf, row.names=FALSE, col.names=FALSE, quote=FALSE)

