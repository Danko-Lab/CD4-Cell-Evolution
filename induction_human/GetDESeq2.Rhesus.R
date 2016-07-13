## MA Plot
require(bigWig)

refGene <- read.table("refGene.bed.gz")
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene <- refGene[(refGene$V3-refGene$V2)>1000,]

bodies <- refGene
bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+500
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-500

bodies$V3[(bodies$V3 - bodies$V2) > 60000 & bodies$V6 == "+"] <- bodies$V2[(bodies$V3 - bodies$V2) > 60000 & bodies$V6 == "+"]+60000
bodies$V2[(bodies$V3 - bodies$V2) > 60000 & bodies$V6 == "-"] <- bodies$V3[(bodies$V3 - bodies$V2) > 60000 & bodies$V6 == "-"]-60000

getCounts <- function(plus, minus, path, intervals= bodies) {
  pl <- load.bigWig(paste(path, plus, sep=""))
  mn <- load.bigWig(paste(path, minus, sep=""))
  counts <- bed6.region.bpQuery.bigWig(pl, mn, intervals, abs.value = TRUE) ## Get counts
  counts #* (1000/(bodies$V3-bodies$V2)) * (1e6/ (pl$mean*pl$basesCovered+mn$mean*mn$basesCovered)) ## Normalize to RPKM
}

raw_counts <- cbind(
macaque_1_U= getCounts("M1-U.bed.gz_plus.hg19.bw", "M1-U.bed.gz_minus.hg19.bw", "../AllData/"),
macaque_2_U= getCounts("M2-U_plus.hg19.bw", "M2-U_minus.hg19.bw", "../AllData/"),
macaque_3_U= getCounts("M3-U_plus.hg19.bw", "M3-U_minus.hg19.bw", "../AllData/"),
#macaque_1_PI= getCounts("C3-PI.bed.gz_plus.hg19.bw", "C3-PI.bed.gz_minus.hg19.bw", "../AllData/"),
macaque_2_PI= getCounts("M2-PI_plus.hg19.bw", "M2-PI_minus.hg19.bw", "../AllData/"),
macaque_3_PI= getCounts("M3-PI_plus.hg19.bw", "M3-PI_minus.hg19.bw", "../AllData/")
)
print(cor(raw_counts, method="spearman"))

library("DESeq2")
colData <- data.frame(Condition= c(rep("U",3), rep("PI",2)), row.names=colnames(raw_counts))

## Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData= raw_counts, colData= colData, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="U") ## Set the reference condition as the primary tumor.

dds <- DESeq(dds)
res <- results(dds)

print(paste("Number of changes: ", sum(res$padj < 0.01, na.rm=TRUE))) ## Number of transcripts.
print(paste("Number of unique genes: ", NROW(unique(refGene$V7[res$padj < 0.01])))) ## Number of genes.

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
pdf("results/rhesus_macaque.MAplot.pdf")

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
write.table(PX, "results/rhesus-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

q("no")

hdU <- read.table("../dREG_HD/M-U_dREG_HD.bed")
hdPI<- read.table("../dREG_HD/M-PI_dREG_HD.bed")

hd <- read.table("../dREG_HD/dREG_HD.merge.HCM.UPI.rhesus.bed") #rbind(hdU, hdPI)
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
rhesus_1_U= getCountsE("M1-U.bed.gz_plus.bw", "M1-U.bed.gz_minus.bw", "../Alignments_1stPrep/"),
rhesus_2_U= getCountsE("M2-U_plus.bw", "M2-U_minus.bw", "../Alignments_2ndPrep/"),
rhesus_3_U= getCountsE("M3-U_plus.bw", "M3-U_minus.bw", "../Alignments_3rdPrep/"),
rhesus_2_PI= getCountsE("M2-PI_plus.bw", "M2-PI_minus.bw", "../Alignments_2ndPrep/"),
rhesus_3_PI= getCountsE("M3-PI_plus.bw", "M3-PI_minus.bw", "../Alignments_3rdPrep/")
)

print(cor(raw_counts, method="spearman"))

colData <- data.frame(Condition= c(rep("U",3), rep("PI",2)), row.names=colnames(raw_counts))

## Create DESeq2 object.
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData= raw_counts, colData= colData, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="U") ## Set the reference condition as the primary tumor.

dds <- DESeq(dds)
res <- results(dds)

print(paste("Number of changes: ", sum(res$padj < 0.01, na.rm=TRUE))) ## Number of transcripts.
#print(paste("Number of unique genes: ", NROW(unique(refGene$V7[res$padj < 0.01])))) ## Number of genes.

## Write out dREG-HD TREs for motif analysis.
hd <- read.table("../dREG_HD/dREG_HD.merge.HCM.UPI.rhesus.bed") #rbind(hdU, hdPI); 
hd <- hd[grep("random|Un", hd$V1, invert=TRUE),]

PE <- cbind(rbind(hd), res)
PE <- PE[order(PE$padj),]
write.table(PE, "results/rhesus-changed.TREs.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

