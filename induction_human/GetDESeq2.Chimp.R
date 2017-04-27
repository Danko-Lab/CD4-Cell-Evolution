## MA Plot
require(bigWig)

refGene <- read.table("tuSelecter/final_tus.txt", header=TRUE)  #refGene.bed.gz")

refGene <- refGene[grep("random|Un|hap", refGene$TXCHROM, invert=TRUE),]
refGene <- refGene[(refGene$TXEND-refGene$TXSTART)>1000,]

bodies <- refGene
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+500
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-500

bodies$TXEND[(bodies$TXEND - bodies$TXSTART) > 60000 & bodies$TXSTRAND == "+"] <- bodies$TXSTART[(bodies$TXEND - bodies$TXSTART) > 60000 & bodies$TXSTRAND == "+"]+60000
bodies$TXSTART[(bodies$TXEND - bodies$TXSTART) > 60000 & bodies$TXSTRAND == "-"] <- bodies$TXEND[(bodies$TXEND - bodies$TXSTART) > 60000 & bodies$TXSTRAND == "-"]-60000

getCounts <- function(plus, minus, path, intervals= bodies) {
  pl <- load.bigWig(paste(path, plus, sep=""))
  mn <- load.bigWig(paste(path, minus, sep=""))
  counts <- bed6.region.bpQuery.bigWig(pl, mn, intervals, abs.value = TRUE) ## Get counts
  counts #* (1000/(bodies$V3-bodies$V2)) * (1e6/ (pl$mean*pl$basesCovered+mn$mean*mn$basesCovered)) ## Normalize to RPKM
}

raw_counts <- cbind(
chimp_1_U= getCounts("C3-U.bed.gz_plus.hg19.bw", "C3-U.bed.gz_minus.hg19.bw", "../AllData/"),
chimp_2_U= getCounts("C4-U.bed.gz_plus.hg19.bw", "C4-U.bed.gz_minus.hg19.bw", "../AllData/"),
chimp_3_U= getCounts("C5-U_plus.hg19.bw", "C5-U_minus.hg19.bw", "../AllData/"),
chimp_1_PI= getCounts("C3-PI.bed.gz_plus.hg19.bw", "C3-PI.bed.gz_minus.hg19.bw", "../AllData/"),
chimp_2_PI= getCounts("C4-PI.bed.gz_plus.hg19.bw", "C4-PI.bed.gz_minus.hg19.bw", "../AllData/"),
chimp_3_PI= getCounts("C5-PI_plus.hg19.bw", "C5-PI_minus.hg19.bw", "../AllData/")
)
print(cor(raw_counts, method="spearman"))

library("DESeq2")
colData <- data.frame(Condition= c(rep("U",3), rep("PI",3)), row.names=colnames(raw_counts))

## Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData= raw_counts, colData= colData, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="U") ## Set the reference condition as the primary tumor.

dds <- DESeq(dds)
res <- results(dds)

print(paste("Number of changes: ", sum(res$padj < 0.01, na.rm=TRUE))) ## Number of transcripts.
print(paste("Number of unique genes: ", NROW(unique(refGene$GENENAME[res$padj < 0.01])))) ## Number of genes.

## Add specific genes to the plot.
addlab <- function(gene_ID, deRes, genes, ...) {
 idx<-sapply(gene_ID, function(gene_ID) {
  ig <- which(genes[,7] == gene_ID)
  io <- ig[which.min(deRes$padj[ig])]
  if(NROW(io)>0) {
        text(deRes$baseMean[io], deRes$log2FoldChange[io], labels= gene_ID, cex= 1, pos= 3, ...)
        points(deRes$baseMean[io], deRes$log2FoldChange[io], col="blue", cex=1.5)
        io
  }
 })
 print(idx)
 idx <- unlist(idx); idx <- idx[!is.null(idx)]
 return(data.frame(Gene= genes[idx,7], AveExpr= deRes$baseMean[idx], logFC= deRes$log2FoldChange[idx], adj.P.Val= deRes$padj[idx]))
}

## MA plot.
pdf("results/chimp.MAplot.pdf")

plotMA(res, ylim=c(-11,11), xlim=c(1e-1, 2e5), alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, U/ PI")

plotMA(res, ylim=c(-11,11), xlim=c(1e-1, 2e5), alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, U/ PI")
addlab(c("CDK5", "GIMAP4", "VMAC", "ZKSCAN4", "MRS2", "ETS1"), res, bodies) ## Downreg.
addlab(c("IL2", "IL21", "IL2RA", "TNF", "IFNG", "NFKB1", "EGR3", "NR4A1", "C5orf52"), res, bodies) ## Upreg.

plot(0, 0, ylim=c(-11,11), xlim=c(1e-1, 2e5), log="x")
addlab(c("CDK5", "GIMAP4", "VMAC", "ZKSCAN4", "MRS2", "ETS1"), res, bodies) ## Downreg.
addlab(c("IL2", "IL21", "IL2RA", "TNF", "IFNG", "NFKB1", "EGR3", "NR4A1", "C5orf52"), res, bodies) ## Upreg.

dev.off()

## Write out genes.
PX <- cbind(bodies, res)
PX <- PX[order(PX$padj),]
write.table(PX, "results/chimp-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

q("no")

######################
## Now get REs.
hdU <- read.table("../dREG_HD/C-U_dREG_HD.bed")
hdPI<- read.table("../dREG_HD/C-PI_dREG_HD.bed")

hd <- read.table("../dREG_HD/dREG_HD.merge.HCM.UPI.chimp.bed") #rbind(hdU, hdPI)
hd$V2 <- hd$V2-250; hd$V2[hd$V2 < 0] = 0
hd$V3 <- hd$V3+250

hd <- hd[grep("random|Un", hd$V1, invert=TRUE),]

getCountsE <- function(plus, minus, path, intervals= hd) {
  pl <- load.bigWig(paste(path, plus, sep=""))
  mn <- load.bigWig(paste(path, minus, sep=""))
  counts <- bed.region.bpQuery.bigWig(pl, intervals, abs.value = TRUE)+bed.region.bpQuery.bigWig(mn, intervals, abs.value = TRUE) ## Get counts
  counts #* (1000/(bodies$V3-bodies$V2)) * (1e6/ (pl$mean*pl$basesCovered+mn$mean*mn$basesCovered)) ## Normalize to RPKM
}

###################################################################################
### NOW FOCUS ON dREG-HD
raw_counts <- cbind(
chimp_1_U= getCountsE("C3-U_plus.bw", "C3-U_minus.bw", "../Alignments_2ndPrep/"),
chimp_2_U= getCountsE("C4-U_plus.bw", "C4-U_minus.bw", "../Alignments_3rdPrep/"),
chimp_3_U= getCountsE("C5-U_plus.bw", "C5-U_minus.bw", "../Alignments_4thPrep/"),
chimp_1_PI= getCountsE("C3-PI_plus.bw", "C3-PI_minus.bw", "../Alignments_2ndPrep/"),
chimp_2_PI= getCountsE("C4-PI_plus.bw", "C4-PI_minus.bw", "../Alignments_3rdPrep/"),
chimp_3_PI= getCountsE("C5-PI_plus.bw", "C5-PI_minus.bw", "../Alignments_4thPrep/")
)

print(cor(raw_counts, method="spearman"))

colData <- data.frame(Condition= c(rep("U",3), rep("PI",3)), row.names=colnames(raw_counts))

## Create DESeq2 object.
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData= raw_counts, colData= colData, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="U") ## Set the reference condition as the primary tumor.

dds <- DESeq(dds)
res <- results(dds)

print(paste("Number of changes: ", sum(res$padj < 0.01, na.rm=TRUE))) ## Number of transcripts.
#print(paste("Number of unique genes: ", NROW(unique(refGene$V7[res$padj < 0.01])))) ## Number of genes.

## Write out dREG-HD TREs for motif analysis.
hd <- read.table("../dREG_HD/dREG_HD.merge.HCM.UPI.chimp.bed") #rbind(hdU, hdPI); 
hd <- hd[grep("random|Un", hd$V1, invert=TRUE),]

PE <- cbind(hd, res)
PE <- PE[order(PE$padj),]
write.table(PE, "results/chimp-changed.TREs.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


