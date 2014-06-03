Do <- read.table("GENCODE/divergent.Overlap.bed")
Dn <- read.table("GENCODE/divergent.NoOverlap.bed")
Dn <- Dn[!(Dn$V4 %in% Do$V4),]
write.table(Dn, "GENCODE/divergent.NoOverlap.bed", row.names=F, quote=F, col.names=F, sep="\t")
