################################################################################
## ID changes in GE with Limma Voom
require(limma)
require(edgeR)

PVAL <- 0.01
FOLD <- 0 #log(10,2)

## Read count data.
source("readData.R")

## Use only untreated, and get switch to RPKM
ca <- ca[,c(1:10,indx.unt)]

rpkm_df <- as.matrix(ca[,c(11:NCOL(ca))])/(ca[,"mapSize"]) #/ (colSums(ca[,c(10:15)])) ## Normalize counts ... RPKM
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*rpkm_df[,i]/sum(rpkm_df[,i])

counts <- ca[,c(12:19)]
genes  <- ca[,c(1:10)]

## Build experimental design matrix
sampleID <- c("Human 1", "Human 2", "Human 3", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3")
prepDay  <- c(1, 4, 2, 3, 4, 2, 3, 4)
subject  <- c(2, 3, 4, 5, 6, 7, 8, 9)
species  <- c("H", "H", "H", "C", "C", "M", "M", "M")
Design <- data.frame(sampleID, prepDay, subject, species)

colnames(counts) <- sampleID

dge <- DGEList(counts=counts, genes=genes)
dge <- calcNormFactors(dge)

## Fit the model spearately for each species comparison, treating 'outgroups' as 'other'.
fitModel <- function(species) {
  design <- model.matrix(~species, Design)
  rownames(design) <- colnames(dge)

  voom_obj <- voom(dge, design, plot=FALSE)

#  plotMDS(dge,top=500,labels=species,gene.selection="common")

  ## Estimate dispersions.
  fit <- lmFit(voom_obj, design)
  fit <- eBayes(fit)

#  plotMA(fit, array=2, status=as.factor(decideTests(fit)[,2]), col=c("black", "red", "blue"))
#  plotMA(fit, array=2, status=(p.adjust(fit$p.value[,2], method="fdr")<PVAL), col=c("black", "red"))

  return(fit)
}

## Treating alternatie primates as a single 'group'.
species <- c("SS", "SS", "SS", "SO", "SO", "SO", "SO", "SO")
Design <- data.frame(sampleID, prepDay, subject, species)
hs <- fitModel(species)
topTable(hs, coef=2, n=16, sort="p")

species <- c("SO", "SO", "SO", "SS", "SS", "SO", "SO", "SO")
Design <- data.frame(sampleID, prepDay, subject, species)
cs <- fitModel(species)
topTable(cs, coef=2, n=16, sort="p")

species <- c("SO", "SO", "SO", "SO", "SO", "SS", "SS", "SS")
Design <- data.frame(sampleID, prepDay, subject, species)
ms <- fitModel(species)
topTable(ms, coef=2, n=16, sort="p")

## Append tables.  
fdr_t <- matrix(p.adjust(c(hs$p.value[,2], cs$p.value[,2], ms$p.value[,2]), method="fdr"), ncol=3)
colnames(fdr_t) <- c("HumanFDR", "ChimpFDR", "MacaqueFDR")

fc_t  <- data.frame(HumanFC= hs$coefficients[,2], ChimpFC= cs$coefficients[,2], MacaqueFC= ms$coefficients[,2])

fdr_min<- sapply(c(1:NROW(fdr_t)), function(x) {min(fdr_t[x,])} )
fc_min <- sapply(c(1:NROW(fdr_t)), function(x) {
	mm<-which.min(fdr_t[x,]); 
	if(NROW(mm) == 0) {return(NA)}; 
	return(fc_t[x,which.min(fdr_t[x,])])})

fdr_df <- (data.frame(ca, fdr_t, fdr_min, fc_t, fc_min))

## Count absolute numbers of changes.
NROW(unique(ca[fdr_t[,1] < PVAL,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL,"mgi"])) # RHESUS

NROW(unique(ca[fdr_t[,1] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # RHESUS

## Write table of genes to use in GO.
writeExprChange <- function(ss, name, indx) {
 write.table(unique(ca[fdr_t[,indx] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes, fdr_t, fc_t)[fdr_t[,indx]<PVAL, ], paste("chage_expr/",name,".change.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

writeExprChange(ss, "H", 1)
writeExprChange(ss, "C", 2)
writeExprChange(ss, "M", 3)

isExpr <- rowMeans(rpkm_df[,2:9]) > 1e-4
write.table(unique(ca[isExpr & ca[,"annot_type"] == "gc18","mgi"]), "chage_expr/exprbg.mgi.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

## Compute frequency of changes for 'expressed' genes in each class.
save.image("fdr.RData")

