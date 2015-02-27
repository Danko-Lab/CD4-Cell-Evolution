##
##
## compareH2_H2.2_PI.R -- Compares two independent draws of human 2.  
##                        Differences: 
##                           (1) Completely separate days for the prep and PRO-seq.
##                           (2) H2 processed immediately, H2.2 processed after 24h at 4C.
##                           (3) ...

epsilon <- 1e-10

gc18 <- read.table("../annotations/gencode.v18.transcript.tsv")
gc18 <- gc18[(gc18$V3-gc18$V2) > 5000,] ## Remove short transcripts.

require(bigWig)

h2_p <- load.bigWig("../Alignments_2ndPrep/H2-PI.bed.gz_plus.bw")
h2_m <- load.bigWig("../Alignments_2ndPrep/H2-PI.bed.gz_minus.bw")

h2.2_p <- load.bigWig("../Alignments_3rdPrep/H2.2-PI.bed.gz_plus.bw")
h2.2_m <- load.bigWig("../Alignments_3rdPrep/H2.2-PI.bed.gz_minus.bw")

h1u_p <- load.bigWig("../../CD4/Alignments/H1-U_plus.bw")
h1u_m <- load.bigWig("../../CD4/Alignments/H1-U_minus.bw")
h1pi_p <- load.bigWig("../../CD4/Alignments/H1-PI_plus.bw")
h1pi_m <- load.bigWig("../../CD4/Alignments/H1-PI_minus.bw")

h2 <- bed6.region.bpQuery.bigWig(h2_p, h2_m, gc18)/22464429*1000/(gc18$V3-gc18$V2) ## RPKM
h2.2 <- bed6.region.bpQuery.bigWig(h2.2_p, h2.2_m, gc18)/29758309*1000/(gc18$V3-gc18$V2) ## RPKM
h1u <- bed6.region.bpQuery.bigWig(h1u_p, h1u_m, gc18)/38036167*1000/(gc18$V3-gc18$V2) ## RPKM
h1pi<- bed6.region.bpQuery.bigWig(h1pi_p, h1pi_m, gc18)/39753494*1000/(gc18$V3-gc18$V2) ## RPKM

## Correlations.
cor.test(abs(h2), abs(h2.2), method="spearman") ## 0.9456079
cor.test(log(abs(h2)+epsilon), log(abs(h2.2)+epsilon), method="pearson") ## 0.9166077
cor.test(log(abs(h2[h2>0 & h2.2>0])), log(abs(h2.2[h2>0 & h2.2>0])), method="pearson") ## 0.9290175

## All correlations, including with H1.
cor(cbind(log(abs(h2)+epsilon), log(abs(h2.2)+epsilon), log(abs(h1u)+epsilon), log(abs(h1pi)+epsilon)))
cor(cbind(log(abs(h2)+epsilon), log(abs(h2.2)+epsilon), log(abs(h1u)+epsilon), log(abs(h1pi)+epsilon)), method="spearman")

## Write out some plots.
plot(log(abs(h2)+epsilon), log(abs(h2.2)+epsilon), pch=19)
plot(log(abs(h2[h2>0 & h2.2>0])), log(abs(h2.2[h2>0 & h2.2>0])), pch=19)
abline(0,1, col="blue")

densScatterplot(log(abs(h2)+epsilon), log(abs(h2.2)+epsilon))

## Sanity check some longer genes...
df <- cbind(gc18[,c(1:3.6,8)], lh2=log(abs(h2)+epsilon), lh2.2=log(abs(h2.2)+epsilon), lh1u=log(abs(h1u)+epsilon), lh1pi=log(abs(h1pi)+epsilon))
df[grep("NFKB1", gc18$V8),]

## Write out PDF ...
pdf("Processing.Before.After.pdf")
 densScatterplot(log(abs(h2)[h2>0 & h2.2>0]), log(abs(h2.2)[h2>0 & h2.2>0]), xlab="Overnight @4C", ylab="Process immediately")
dev.off()

