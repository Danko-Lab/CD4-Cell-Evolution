#
# gene_enhancer_correlation.R -- correlates gene expression with enhancer initiation levels.
#

as <- rbind(read.table("tmp/H.change-U.genepc.as"), read.table("tmp/C.change-U.genepc.as"), read.table("tmp/M.change-U.genepc.as"))

as <- as[!is.na(as$V7),]

print(cor.test(as$V5, as$V7))
print(cor.test(as$V5, as$V7, method="spearman"))

plot(as$V5, as$V7, xlab= "Gene Expression", ylab="UAS Expression", pch=19)

pdf("gene_uas_correlation.pdf")
 plot(as$V5, as$V7, xlab= "Gene Expression", ylab="UAS Expression", pch=19)

 source("../lib/densScatterplot.R")
 densScatterplot(as$V5, as$V7, xlab= "Gene Expression", ylab="UAS Expression")
 abline(a=0, b=1)
dev.off()
