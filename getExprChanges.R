################################################################################
## ID changes in GE with edgeR (v.3.2.4).
require(edgeR)

PVAL <- 0.05

## Read count data.
ca <- read.table("countall.tsv")
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", "Jurkat", "Human 1", "Human 3", "Chimp 3", "R. Macaque 1", "R. Macaque 2")
ca <- ca[!is.na(ca[,10]),]
rpkm_df <- as.matrix(ca[,c(10:15)])/(ca[,9]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

dge <- DGEList(counts=ca[,c(10:15)], genes=ca[,1:8])

## Build experimental design matrix
sampleID <- c("Jurkat", "Human 1", "Human 3", "Chimp 3", "R. Macaque 1", "R. Macaque 2")
prepDay  <- c(1, 1, 2, 3, 2, 3)
species  <- c("J", "H", "H", "C", "M", "M")
data.frame(sampleID, prepDay, species)

design <- model.matrix(~prepDay+species)
rownames(design) <- colnames(dge)

## Estimate dispersions.
dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE)
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge,design)

## Fit neg. binom. GLM.  Compute p-values using LRT.
fit <- glmFit(dge,design)
#lrt <- glmLRT(fit,coef=2)

HumSpec <- c(0,0,1,-0.5,-0.5)
hs <- glmLRT(fit,contrast=HumSpec)
topTags(hs)

ChpSpec <- c(0,0,-0.334,-0.334,-0.334)
cs <- glmLRT(fit,contrast=ChpSpec)
topTags(cs)

MacSpec <- c(0,0,-0.5,-0.5,1)
ms <- glmLRT(fit,contrast=MacSpec)
topTags(ms)

## Append tables.  
if(sum(ca$name == lrt$genes$name)/NROW(ca) !=1) print("ERROR!  Sanity check failed")
fdr_t <- data.frame(HumanFDR= p.adjust(hs$table$PValue), 
			ChimpFDR= p.adjust(cs$table$PValue),
			MacaqueFDR= p.adjust(ms$table$PValue))
fdr_min<- sapply(c(1:NROW(fdr_t)), function(x) {min(fdr_t[x,])} )
fdr_df <- (data.frame(ca, fdr_t, fdr_min))

## Compute frequency of changes for 'expressed' genes in each class.
save.image("fdr.RData")

