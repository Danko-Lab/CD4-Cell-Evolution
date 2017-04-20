################################################################################
## ID changes in GE with DESeq2 (v 1.2.8)
require(DESeq2)

PVAL <- 0.01
FOLD <- 1
#FOLD <- log(10,2)

## Read count data.
source("readData.R")

## Use only untreated, and get switch to RPKM
#ca <- ca[,c(1:10,indx.unt)]

counts <- ca[,c(12:20)]
countsPI <- ca[,c(22:30)]
genes  <- ca[,c(1:10)]

## Build experimental design matrix
sampleID <- c("Human 1", "Human 2", "Human 4", "Chimp 3", "Chimp 4", "Chimp 5", "R. Macaque 2", "R. Macaque 3", "R. Macaque 4")
prepDay  <- c(1, 3, 5, 2, 3, 4, 2, 3, 5)
subject  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
species  <- c("H", "H", "H", "C", "C", "C", "M", "M", "M")
Design <- data.frame(sampleID, prepDay, subject, species)

## Fit the model spearately for each species comparison, treating 'outgroups' as 'other'.
fitModel <- function(species, cnts) {
  Design <- data.frame(sampleID, prepDay, subject, species)

  ## Estimate dispersions.
  dds <- DESeqDataSetFromMatrix(countData= cnts, colData= Design, design= ~ species)
  dds <- DESeq(dds, betaPrior=FALSE)
  dds <- replaceOutliersWithTrimmedMean(dds) ## NAs otherwise.
  resultsNames(dds)
  res <- results(dds, name="species_SS_vs_SO")
  ss <- data.frame(genes, res)

  return(ss)
}

## Treating alternatie primates as a single 'group'.
species <- c("SS", "SS", "SS", "SO", "SO", "SO", "SO", "SO", "SO")
hs <- fitModel(species, counts, Design)
hs_pi <- fitModel(species, countsPI, Design)
head(hs[order(hs$padj), ])

species <- c("SO", "SO", "SO", "SS", "SS", "SS", "SO", "SO", "SO")
cs <- fitModel(species, counts, Design)
cs_pi <- fitModel(species, countsPI, Design)
head(cs[order(cs$padj), ])

species <- c("SO", "SO", "SO", "SO", "SO", "SO", "SS", "SS", "SS")
ms <- fitModel(species, counts)
ms_pi <- fitModel(species, countsPI)
head(ms[order(ms$padj), ])

## Fit model U->PI
species  <- c("H", "H", "H", "C", "C", "C", "M", "M", "M")
fitModel_U2PI <- function(condition, cnts, indx) {
  Design <- data.frame(sampleID= rep(sampleID, 2), prepDay= rep(prepDay,2), subject= rep(subject,2), species= rep(species,2))[indx,]; Design <- cbind(Design, condition)
  cnts <- cnts[,indx]

  print(Design)

  ## Estimate dispersions.
  dds <- DESeqDataSetFromMatrix(countData= cnts, colData= Design, design= ~ condition)
  dds <- DESeq(dds, betaPrior=FALSE)
  dds <- replaceOutliersWithTrimmedMean(dds) ## NAs otherwise.
  resultsNames(dds)
  res <- results(dds, name="condition_U_vs_PI")
  ss <- data.frame(genes, res)

  return(ss)
}

## Now look at U2PI conditions.
conditions <- c("U", "U", "U", "PI", "PI", "PI")
hs_tfc <- fitModel_U2PI(conditions, cbind(counts, countsPI), c(1:3,10:12))
cs_tfc <- fitModel_U2PI(conditions, cbind(counts, countsPI), c(4:6,13:15))
ms_tfc <- fitModel_U2PI(conditions, cbind(counts, countsPI), c(7:9,16:18))
tfc    <- fitModel_U2PI(c(rep("U", 9), rep("PI", 9)), cbind(counts, countsPI), c(1:18))

## Append tables.  
fdr_t <- data.frame(HumanFDR= hs$padj, ChimpFDR= cs$padj, MacaqueFDR= ms$padj, HumanFDR_PI= hs_pi$padj, ChimpFDR_PI= cs_pi$padj, MacaqueFDR_PI= ms_pi$padj, 
			U2PIFDR_H= hs_tfc$padj, U2PIFDR_C= cs_tfc$padj, U2PIFDR_M= ms_tfc$padj, U2PIFDR= tfc$padj)
fc_t  <- data.frame(HumanFC= hs$log2FoldChange, ChimpFC= cs$log2FoldChange, MacaqueFC= ms$log2FoldChange, HumanFC_PI= hs_pi$log2FoldChange, ChimpFC_PI= cs_pi$log2FoldChange, MacaqueFC= ms_pi$log2FoldChange,
			U2PIFC_H= hs_tfc$log2FoldChange, U2PIFC_C= cs_tfc$log2FoldChange, U2PIFC_M= ms_tfc$log2FoldChange, U2PIFC= tfc$log2FoldChange)

fdr_t[is.na(fdr_t)] <- 1
fc_t[is.na(fc_t)] <- 0

fdr_min<- rowMin(fdr_t[,1:3]) # sapply(c(1:NROW(fdr_t)), function(x) {min(fdr_t[x,])} )
fdr_min_pi <- rowMin(fdr_t[,4:6])

fc_min <- sapply(c(1:NROW(fdr_t)), function(x) {
        mm<-which.min(fdr_t[x,]);
        if(NROW(mm) == 0) {return(NA)};
        return(fc_t[x,which.min(fdr_t[x,1:3])])})
fc_min_pi <- sapply(c(1:NROW(fdr_t)), function(x) {
        mm<-which.min(fdr_t[x,]);
        if(NROW(mm) == 0) {return(NA)};
        return(fc_t[x,which.min(fdr_t[x,4:6])])})


fdr_df <- (data.frame(ca, fdr_t, fdr_min, fdr_min_pi, fc_t, fc_min, fc_min_pi))

## Count absolute numbers of changes.
ntxn_units <- NROW(unique(ca[,"mgi"])) # HUMAN

NROW(unique(ca[fdr_t[,1] < PVAL & abs(hs$log2FoldChange) > FOLD,"mgi"]))/ ntxn_units # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & abs(cs$log2FoldChange) > FOLD,"mgi"]))/ ntxn_units # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & abs(ms$log2FoldChange) > FOLD,"mgi"]))/ ntxn_units # RHESUS

## Up/ down-regulation.
NROW(unique(ca[fdr_t[,1] < PVAL & hs$log2FoldChange > FOLD,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,1] < PVAL & hs$log2FoldChange < -FOLD,"mgi"])) # HUMAN

NROW(unique(ca[fdr_t[,2] < PVAL & cs$log2FoldChange > FOLD,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,2] < PVAL & cs$log2FoldChange < -FOLD,"mgi"])) # CHIMP

NROW(unique(ca[fdr_t[,3] < PVAL & ms$log2FoldChange > FOLD,"mgi"])) # RHESUS
NROW(unique(ca[fdr_t[,3] < PVAL & ms$log2FoldChange < -FOLD,"mgi"])) # RHESUS


NROW(unique(ca[fdr_t[,1] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # RHESUS

## Write a MA-plot for GC18.
png("img/hs.MA.png")
 plot(hs$baseMean, hs$log2FoldChange, pch=19, log="x")
 points(hs$baseMean[hs$padj<PVAL], hs$log2FoldChange[hs$padj<PVAL], col="red", pch=19)
# points(hs$baseMean[hs$padj<PVAL & abs(hs$log2FoldChange) > FOLD], hs$log2FoldChange[hs$padj<PVAL & abs(hs$log2FoldChange) > FOLD], col="red", pch=19)
dev.off()

## Write table of genes to use in GO.
writeExprChange <- function(ss, name, indx, indx_pi) {
 ## For GO.
 write.table(unique(ca[fdr_t[,indx] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".U.mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(unique(ca[fdr_t[,indx_pi] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".PI.mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(unique(ca[fdr_t[,indx] < PVAL & fdr_t[,indx_pi] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".U+PI.mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

 ## Tables ...
 write.table(cbind(genes[,1:4], fc_t[,indx],genes[,6:10], fdr_t, fc_t)[fdr_t[,indx]<PVAL, ], paste("chage_expr/",name,".change-U.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes[,1:4], fc_t[,indx], genes[,6:10], fdr_t, fc_t)[fdr_t[,indx_pi]<PVAL, ], paste("chage_expr/",name,".change-PI.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes, fdr_t, fc_t)[fdr_t[,indx] < PVAL & fdr_t[,indx_pi]<PVAL, ], paste("chage_expr/",name,".change-U+PI.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

 ## Tables ... all genes (no p-value).
 write.table(cbind(genes[,1:4], fc_t[,indx],genes[,6:10], fdr_t, fc_t), paste("chage_expr/",name,".change-U.all.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes[,1:4], fc_t[,indx_pi],genes[,6:10], fdr_t, fc_t), paste("chage_expr/",name,".change-PI.all.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

writeExprChange(ss, "H", 1,4)
writeExprChange(ss, "C", 2,5)
writeExprChange(ss, "M", 3,6)

isExpr <- rowMax(rpkm_df[,2:9]) > 1e-4
write.table(unique(ca[isExpr & ca[,"annot_type"] == "gc18","mgi"]), "chage_expr/exprbg.mgi.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

## Compute frequency of changes for 'expressed' genes in each class.
save.image("fdr.RData")

