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
human_1_U= getCounts("H1-U_plus.bw", "H1-U_minus.bw", "../AllData/"),
human_2_U= getCounts("H2-U.bed.gz_plus.bw", "H2-U.bed.gz_minus.bw", "../AllData/"),
human_3_U= getCounts("H3-U.bed.gz_plus.bw", "H3-U.bed.gz_minus.bw", "../AllData/"),
human_1_PI= getCounts("H1-PI_plus.bw", "H1-PI_minus.bw", "../AllData/"),
human_2_PI= getCounts("H2-PI_plus.bw", "H2-PI_minus.bw", "../AllData/"),
human_3_PI= getCounts("H3-PI.bed.gz_plus.bw", "H3-PI.bed.gz_minus.bw", "../AllData/")
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

## Exploratory...
cbind(bodies, res)[res$log2FoldChange > 9 & !is.na(res$log2FoldChange),]
cbind(bodies, res)[res$log2FoldChange < -5 & !is.na(res$log2FoldChange) & res$baseMean > 1e2,]

## MA plot.
pdf("results/human.MAplot.pdf")

plotMA(res, ylim=c(-11,11), xlim=c(1e-1, 2e5), alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, U/ PI")

plotMA(res, ylim=c(-11,11), xlim=c(1e-1, 2e5), alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, U/ PI")
addlab(c("CDK5", "GIMAP4", "VMAC", "ZKSCAN4", "MRS2"), res, bodies) ## Downreg.
addlab(c("IL2", "IL21", "IL2RA", "TNF", "IFNG", "NFKB1", "EGR3", "NR4A1", "C5orf52"), res, bodies) ## Upreg.

plot(0, 0, ylim=c(-11,11), xlim=c(1e-1, 2e5), log="x")
addlab(c("CDK5", "GIMAP4", "VMAC", "ZKSCAN4", "MRS2"), res, bodies) ## Downreg.
addlab(c("IL2", "IL21", "IL2RA", "TNF", "IFNG", "NFKB1", "EGR3", "NR4A1", "C5orf52"), res, bodies) ## Upreg.

dev.off()

## Write out genes.
PX <- cbind(bodies, res)
PX <- PX[order(PX$padj),]
write.table(PX, "results/human-changed.genes.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


######################
## Now get REs.
dREG<- read.table("../tss_caller/Human-U.dREG.bed.gz")

hdU <- read.table("../dREG_HD/H-U_dREG_HD.bed")
hdPI<- read.table("../dREG_HD/H-PI_dREG_HD.bed")

hd <- rbind(hdU, hdPI)
hd$V2 <- hd$V2-250; hd$V2[hd$V2 < 0] = 0
hd$V3 <- hd$V3+250

getCountsE <- function(plus, minus, path, intervals= hd) {
  pl <- load.bigWig(paste(path, plus, sep=""))
  mn <- load.bigWig(paste(path, minus, sep=""))
  counts <- bed.region.bpQuery.bigWig(pl, intervals, abs.value = TRUE)+bed.region.bpQuery.bigWig(mn, intervals, abs.value = TRUE) ## Get counts
  counts #* (1000/(bodies$V3-bodies$V2)) * (1e6/ (pl$mean*pl$basesCovered+mn$mean*mn$basesCovered)) ## Normalize to RPKM
}

################################################################################
## First focus on dREG ... just to get number of changed sites for the paper.
counts_dreg <- cbind(
human_1_U= getCountsE("H1-U_plus.bw", "H1-U_minus.bw", "../AllData/", intervals= dREG),
human_2_U= getCountsE("H2-U.bed.gz_plus.bw", "H2-U.bed.gz_minus.bw", "../AllData/", intervals= dREG),
human_3_U= getCountsE("H3-U.bed.gz_plus.bw", "H3-U.bed.gz_minus.bw", "../AllData/", intervals= dREG),
human_1_PI= getCountsE("H1-PI_plus.bw", "H1-PI_minus.bw", "../AllData/", intervals= dREG),
human_2_PI= getCountsE("H2-PI_plus.bw", "H2-PI_minus.bw", "../AllData/", intervals= dREG),
human_3_PI= getCountsE("H3-PI.bed.gz_plus.bw", "H3-PI.bed.gz_minus.bw", "../AllData/", intervals= dREG)
)


print(cor(counts_dreg, method="spearman"))

library("DESeq2")
colData <- data.frame(Condition= c(rep("U",3), rep("PI",3)), row.names=colnames(counts_dreg))

## Number of dREG sites that change.
dds <- DESeqDataSetFromMatrix(countData= counts_dreg, colData= colData, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="U") ## Set the reference condition as the primary tumor.

dds <- DESeq(dds)
res <- results(dds)

print(paste("Number of changes: ", sum(res$padj < 0.01, na.rm=TRUE))) ## Number of transcripts.

###################################################################################
### NOW FOCUS ON dREG-HD
raw_counts <- cbind(
human_1_U= getCountsE("H1-U_plus.bw", "H1-U_minus.bw", "../AllData/"),
human_2_U= getCountsE("H2-U.bed.gz_plus.bw", "H2-U.bed.gz_minus.bw", "../AllData/"),
human_3_U= getCountsE("H3-U.bed.gz_plus.bw", "H3-U.bed.gz_minus.bw", "../AllData/"),
human_1_PI= getCountsE("H1-PI_plus.bw", "H1-PI_minus.bw", "../AllData/"),
human_2_PI= getCountsE("H2-PI_plus.bw", "H2-PI_minus.bw", "../AllData/"),
human_3_PI= getCountsE("H3-PI.bed.gz_plus.bw", "H3-PI.bed.gz_minus.bw", "../AllData/")
)

print(cor(raw_counts, method="spearman"))

## Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData= raw_counts, colData= colData, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="U") ## Set the reference condition as the primary tumor.

dds <- DESeq(dds)
res <- results(dds)

print(paste("Number of changes: ", sum(res$padj < 0.01, na.rm=TRUE))) ## Number of transcripts.
#print(paste("Number of unique genes: ", NROW(unique(refGene$V7[res$padj < 0.01])))) ## Number of genes.

## Write out dREG-HD TREs for motif analysis.
PE <- cbind(rbind(hdU, hdPI), res)
PE <- PE[order(PE$padj),]
write.table(PE, "results/human-changed.TREs.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


