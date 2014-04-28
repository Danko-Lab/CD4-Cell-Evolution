Do <- read.table("TMP/divergent.Overlap.bed")
Dn <- read.table("TMP/divergent.NoOverlap.bed")
Dn <- Dn[!(Dn$V4 %in% Do$V4),]
write.table(Dn, "TMP/divergent.NoOverlap.bed", row.names=F, quote=F, col.names=F, sep="\t")
