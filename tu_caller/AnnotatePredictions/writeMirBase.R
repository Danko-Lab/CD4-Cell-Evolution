##
## Converts the mirBase.gff file into a bed format!
##

mb <- read.table("hsa.gff3", comment.char="#")
mb14bed <- data.frame(chrom=paste("chr",mb$V1, sep=""), chromStart=mb$V4, chromEnd=mb$V5,
                name=paste(mb$V9,mb$V10,sep=""), score=rep(1,NROW(mb)), strand=mb$V7)
write.table(mb14bed, "hsa.mirbase20.hg19.bed", row.names=F, col.names=F, quote=F, sep="\t")

