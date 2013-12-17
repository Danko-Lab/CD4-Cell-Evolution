################################################################################
## ID changes in GE with edgeR (v.3.2.4).
require(edgeR)

PVAL <- 0.05

## Read count data.
ca <- read.table("countall.tsv")
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", 
						"Jurkat", "Human 1", "Human 2", "Human 3", "Chimp 2", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3",
						"PI Jurkat ", "PI Human 1", "PI Human 2", "PI Human 3", "PI Chimp 2", "PI Chimp 3", "PI Chimp 4", "PI R. Macaque 1", "PI R. Macaque 2", "PI R. Macaque 3", 
						"K562", "GM12878", "IMR90")
indx.unt <- c(1:9,10:13,15:19,30:32)## ONLY UNTREATED, GOOD Remove C2-U
ca <- ca[!is.na(ca[,10]),indx.unt]

rpkm_df <- as.matrix(ca[,c(10:NCOL(ca))])/(ca[,9]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

dge <- DGEList(counts=ca[,c(10:18)], genes=ca[,1:9])

## Build experimental design matrix
sampleID <- c("Jurkat", "Human 1", "Human 2", "Human 3", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3")
prepDay  <- c(1, 1, 4, 2, 3, 4, 2, 3, 4)
species  <- c("J", "H", "H", "H", "C", "C", "M", "M", "M")
data.frame(sampleID, prepDay, species)

design <- model.matrix(~species)
rownames(design) <- colnames(dge)

## Estimate dispersions.
dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE)
dge <- estimateGLMTrendedDisp(dge,design)

## Fit neg. binom. GLM.  Compute p-values using LRT.
fit <- glmFit(dge,design)
#lrt <- glmLRT(fit,coef=2)

HumSpec <- c(0,1,0,-0.5)
hs <- glmLRT(fit,contrast=HumSpec)
topTags(hs)

ChpSpec <- c(0,-0.5,0,-0.5)
cs <- glmLRT(fit,contrast=ChpSpec)
topTags(cs)

MacSpec <- c(0,-0.5,0,1)
ms <- glmLRT(fit,contrast=MacSpec)
topTags(ms)

## Append tables.  
if(sum(ca$name == fit$genes$name)/NROW(ca) !=1) print("ERROR!  Sanity check failed")
fdr_t <- data.frame(HumanFDR= p.adjust(hs$table$PValue), 
			ChimpFDR= p.adjust(cs$table$PValue),
			MacaqueFDR= p.adjust(ms$table$PValue))
fdr_min<- sapply(c(1:NROW(fdr_t)), function(x) {min(fdr_t[x,])} )
fdr_df <- (data.frame(ca, fdr_t, fdr_min))

## Count absolute numbers of changes.
NROW(unique(ca[fdr_t[,1] < PVAL,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL,"mgi"])) # RHESUS

## Compute frequency of changes for 'expressed' genes in each class.
save.image("fdr.RData")

