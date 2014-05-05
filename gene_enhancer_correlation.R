#
# gene_enhancer_correlation.R -- correlates gene expression with enhancer initiation levels.
#

as <- rbind(read.table("tmp/H.change-U.tsspc.enh"), read.table("tmp/C.change-U.tsspc.enh"), read.table("tmp/M.change-U.tsspc.enh"))

print(cor.test(as$V5, as$V7))
print(cor.test(as$V5, as$V7, method="spearman"))

plot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)

pdf("gene_enhancer_correlation.pdf")
 plot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)

 source("../lib/densScatterplot.R")
 densScatterplot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers")
dev.off()
