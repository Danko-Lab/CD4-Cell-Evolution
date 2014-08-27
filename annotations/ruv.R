##
## Apply the RUV methods of Gagnon-Bartsch & Speed 
##
## REF: Gagnon-Bartsch and Speed, Biostatistics (2012), 13, 3, pp. 539-552.

source("readData.R")


##
## Get a set of genes expected to be invariant ...
require(bigWig)
cd4 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd4/rnaseq/cd4.epigenome.rnaseq.bw")
cd14p1 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapPlusRep1.bigWig")
cd14p2 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapPlusRep2.bigWig")
cd14m1 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapMinusRep1.bigWig")
cd14m2 <- load.bigWig("/usr/data/GROseq.parser/hg19/cd14/rnaseq/wgEncodeCshlLongRnaSeqMonocd14CellPapMinusRep2.bigWig")

## Count reads ... 
CD4c  <- bed.region.bpQuery.bigWig(cd4, ca[,1:6])
CD14c <- bed.region.bpQuery.bigWig(cd14p1, ca[,1:6])+bed.region.bpQuery.bigWig(cd14p2, ca[,1:6])+bed.region.bpQuery.bigWig(cd14m1, ca[,1:6])+bed.region.bpQuery.bigWig(cd14m2, ca[,1:6])

##
## Now apply the RUV method...
require(ruv.DE)

Y <- as.matrix(ca[,indx.good[c(2:9,11:17)]])
#X <- as.integer(as.factor(rbind(c("H","H","H","C","C","M","M","M",    "H","H","H","C","C","M","M"), c(rep("U",8), rep("PI",7)))))

X <- rbind(c(1,1,1,2,2,3,3,3,    1,1,1,2,2,3,3), c(rep(4,8), rep(5,7)))
ctrl <- ca$type == "protein_coding" ## Assume no changes in protein-coding...

norm_vals <- RUVinv(Y, X, ctrl)




