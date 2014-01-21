################################################################################
## ID changes in GE with edgeR (v.3.2.4).
require(edgeR)

PVAL <- 0.05

## Read count data.
source("readData.R")

## Use only untreated, and get switch to RPKM
ca <- ca[!is.na(ca[,11]),c(1:10,indx.unt)]

rpkm_df <- as.matrix(ca[,c(11:NCOL(ca))])/(ca[,"mapSize"]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

dge <- DGEList(counts=ca[,c(11:19)], genes=ca[,1:10])

## Build experimental design matrix
sampleID <- c("Jurkat", "Human 1", "Human 2", "Human 3", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3")
prepDay  <- c(1, 1, 4, 2, 3, 4, 2, 3, 4)
subject  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
species  <- c("J", "H", "H", "H", "C", "C", "M", "M", "M")
data.frame(sampleID, prepDay, subject, species)

fitModel <- function(species) {
  design <- model.matrix(~species) #+prepDay)
  rownames(design) <- colnames(dge)

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE)
  dge <- estimateGLMTrendedDisp(dge,design)

  ## Fit neg. binom. GLM.  Compute p-values using LRT.
  fit <- glmFit(dge,design)

  SpeSpec <- makeContrasts(speciesSS-speciesSO, levels=design)
  ss <- glmLRT(fit,contrast=SpeSpec)
  return(ss)
}

## Treating alternatie primates as a single 'group'.
species <- c("J", "SS", "SS", "SS", "SO", "SO", "SO", "SO", "SO")
hs <- fitModel(species)
topTags(hs)

species <- c("J", "SO", "SO", "SO", "SS", "SS", "SO", "SO", "SO")
cs <- fitModel(species)
topTags(cs)

species <- c("J", "SO", "SO", "SO", "SO", "SO", "SS", "SS", "SS")
ms <- fitModel(species)
topTags(ms)

## Append tables.  
#if(sum(ca$name == fit$genes$name)/NROW(ca) !=1) print("ERROR!  Sanity check failed")
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

