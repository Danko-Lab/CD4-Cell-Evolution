################################################################################
## ID changes in GE with DESeq2 (v 1.2.8)
require(DESeq2)

PVAL <- 0.01
FOLD <- 3

## Read count data.
source("readData.R")

## Use only untreated, and get switch to RPKM
ca <- ca[,c(1:10,indx.unt)]

rpkm_df <- as.matrix(ca[,c(11:NCOL(ca))])/(ca[,"mapSize"]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

counts <- ca[,c(11:19)]
genes  <- ca[,c(1:10)]

## Build experimental design matrix
sampleID <- c("Jurkat", "Human 1", "Human 2", "Human 3", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3")
prepDay  <- c(1, 1, 4, 2, 3, 4, 2, 3, 4)
subject  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
species  <- c("zJ", "H", "H", "H", "C", "C", "M", "M", "M")
Design <- data.frame(sampleID, prepDay, subject, species)

#awd <- DESeqDataSetFromMatrix(countData= counts, colData= Design, design= ~ species)
#awd <- DESeq(awd, betaPrior=FALSE)
##awd <- DESeq(awd, betaPrior=FALSE, test="LRT", reduced = ~ 1)
#resultsNames(awd)

## Fit the model spearately for each species comparison, treating 'outgroups' as 'other'.
fitModel <- function(species) {
  Design <- data.frame(sampleID, prepDay, subject, species)

  ## Estimate dispersions.
  dds <- DESeqDataSetFromMatrix(countData= counts, colData= Design, design= ~ species)
  dds <- DESeq(dds, betaPrior=FALSE)
  dds <- replaceOutliersWithTrimmedMean(dds) ## NAs otherwise.
  resultsNames(dds)
  res <- results(dds, name="species_SS_vs_SO")
  ss <- data.frame(genes, res)

  return(ss)
}

## Treating alternatie primates as a single 'group'.
species <- c("zJ", "SS", "SS", "SS", "SO", "SO", "SO", "SO", "SO")
hs <- fitModel(species)
head(hs[order(hs$padj), ])

species <- c("zJ", "SO", "SO", "SO", "SS", "SS", "SO", "SO", "SO")
cs <- fitModel(species)
head(cs[order(cs$padj), ])

species <- c("zJ", "SO", "SO", "SO", "SO", "SO", "SS", "SS", "SS")
ms <- fitModel(species)
head(ms[order(ms$padj), ])

## Append tables.  
fdr_t <- data.frame(HumanFDR= hs$padj, ChimpFDR= cs$padj, MacaqueFDR= ms$padj)
fc_t  <- data.frame(HumanFC= hs$log2FoldChange, ChimpFC= cs$log2FoldChange, MacaqueFC= ms$log2FoldChange)
fdr_min<- sapply(c(1:NROW(fdr_t)), function(x) {min(fdr_t[x,])} )
fc_min <- sapply(c(1:NROW(fdr_t)), function(x) {fc_t[x,which.min(fdr_t[x,])]})
fdr_df <- (data.frame(ca, fdr_t, fdr_min, fc_min))

## Count absolute numbers of changes.
NROW(unique(ca[fdr_t[,1] < PVAL & abs(hs$log2FoldChange) > FOLD,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & abs(cs$log2FoldChange) > FOLD,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & abs(ms$log2FoldChange) > FOLD,"mgi"])) # RHESUS

NROW(unique(ca[fdr_t[,1] < PVAL & ca[,"annot_type"] == "gc18" & abs(hs$log2FoldChange) > FOLD,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & ca[,"annot_type"] == "gc18" & abs(cs$log2FoldChange) > FOLD,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & ca[,"annot_type"] == "gc18" & abs(ms$log2FoldChange) > FOLD,"mgi"])) # RHESUS

## Compute frequency of changes for 'expressed' genes in each class.
save.image("fdr.RData")

