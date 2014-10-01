################################################################################
## ID changes in GE with Limma Voom
PVAL <- 0.05
FOLD <- 0 #log(10,2)

## Read count data.
source("readData.R")

## Use only untreated, and get switch to RPKM
ca <- ca[,c(1:10,indx.good)]

counts <- ca[,c(12:19,21:27)]
genes  <- ca[,c(1:10)]

## Build experimental design matrix
sampleID <- c("Human 1", "Human 2", "Human 3", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3", 
		"PI Human 1", "PI Human 2", "PI Human 3", "PI Chimp 3", "PI Chimp 4", "PI R. Macaque 2", "PI R. Macaque 3")
prepDay  <- c(1,4,2,3,4,2,3,4,    1,4,2,3,4,3,4)
subject  <- c(2,3,4,5,6,7,8,9,    2,3,4,5,6,8,9)
species  <- c("H","H","H","C","C","M","M","M",    "H","H","H","C","C","M","M")
treatment<- c(rep("U",8), rep("PI",7))
Design   <- data.frame(sampleID, prepDay, subject, species, treatment)
colnames(counts) <- sampleID

lib.size <- c(H1U= 37609124, H2U= 38413904, H3U= 20356701, C2U= 11613548, C3U= 22784541, C4U= 43539929, M1U= 11395858, M2U= 21926503, M3U= 35784519,
                H1PI= 39231203, H2PI= 51228809, H3PI= 4774117, C2PI= 3040431, C3PI= 23702213, C4PI= 35190892, M1PI= 127687, M2PI= 15569212, M3PI= 32421627)[c(1:3,5:12,14:15,17:18)]

## Indices for different combinations of data.
dindx_u  <- c(1:8)
dindx_pi <- c(9:15)
dindx_h  <- c(1:3,9:11)
dindx_c  <- c(4:5,12:13)
dindx_m  <- c(6:8,14:15)

## Fit the model spearately for each species comparison, treating 'outgroups' as 'other'.
fitModel <- function(Design, counts, lib.size) {
  dge <- DGEList(counts=counts, genes=genes)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~condition, Design)
  rownames(design) <- colnames(dge)

  voom_obj <- voom(dge, design, normalize.method = "quantile", lib.size= lib.size, plot=FALSE)

#  plotMDS(dge,top=500,labels=species,gene.selection="common")

  ## Estimate dispersions.
  fit <- lmFit(voom_obj, design)
  fit <- eBayes(fit)

#  plotMA(fit, array=2, status=as.factor(decideTests(fit)[,2]), col=c("black", "red", "blue"))
#  plotMA(fit, array=2, status=(p.adjust(fit$p.value[,2], method="fdr")<PVAL), col=c("black", "red"))
#  abline(h=0, col="blue")

  print(topTable(fit, coef=2, n=16, sort="p"))
  return(fit)
}

## Treating alternatie primates as a single 'group' in either the U or PI condition separately.
condition <- species
condition[species == "H"] <- "SS"
condition[species != "H"] <- "SO"
Design <- data.frame(sampleID, prepDay, subject, condition)
hs    <- fitModel(Design[dindx_u,], counts[,dindx_u], lib.size[dindx_u])
hs_pi <- fitModel(Design[dindx_pi,], counts[,dindx_pi], lib.size[dindx_pi])

condition <- species
condition[species == "C"] <- "SS"
condition[species != "C"] <- "SO"
Design <- data.frame(sampleID, prepDay, subject, condition)
cs <- fitModel(Design[dindx_u,], counts[,dindx_u], lib.size[dindx_u])
cs_pi <- fitModel(Design[dindx_pi,], counts[,dindx_pi], lib.size[dindx_pi])

condition <- species
condition[species == "M"] <- "SS"
condition[species != "M"] <- "SO"
Design <- data.frame(sampleID, prepDay, subject, condition)
ms <- fitModel(Design[dindx_u,], counts[,dindx_u], lib.size[dindx_u])
ms_pi <- fitModel(Design[dindx_pi,], counts[,dindx_pi], lib.size[dindx_pi])

modelCondition <- function() {
  Design <- data.frame(sampleID, prepDay, subject, species, treatment)

  dge <- DGEList(counts=counts, genes=genes)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~species+treatment, Design)
  rownames(design) <- colnames(dge)

  voom_obj <- voom(dge, design, plot=FALSE)
  fit <- lmFit(voom_obj, design)
  fit <- eBayes(fit)

  plotMA(fit, array=2, status=(p.adjust(fit$p.value[,2], method="fdr")<PVAL), col=c("black", "red"))
  abline(h=0, col="blue")

  cor.test(fit$coefficients[,4], tfc$coefficients[,2]) ## Computed fold-changes ... almost identical
}

condition <- treatment
Design <- data.frame(sampleID, prepDay, subject, condition)
hs_tfc <- fitModel(Design[dindx_h,], counts[,dindx_h], lib.size[dindx_h])
cs_tfc <- fitModel(Design[dindx_c,], counts[,dindx_c], lib.size[dindx_c])
ms_tfc <- fitModel(Design[dindx_m,], counts[,dindx_m], lib.size[dindx_m])
tfc    <- fitModel(Design, counts, lib.size)

## Append tables.  
#fdr_t <- matrix(p.adjust(c(hs$p.value[,2], cs$p.value[,2], ms$p.value[,2], hs_pi$p.value[,2], cs_pi$p.value[,2], ms_pi$p.value[,2], 
#				hs_tfc$p.value[,2], cs_tfc$p.value[,2], ms_tfc$p.value[,2], tfc$p.value[,2]), method="fdr"), ncol=10)

fdr_t <- cbind(p.adjust(hs$p.value[,2], method="fdr"), p.adjust(cs$p.value[,2], method="fdr"), p.adjust(ms$p.value[,2], method="fdr"), 
		p.adjust(hs_pi$p.value[,2], method="fdr"), p.adjust(cs_pi$p.value[,2], method="fdr"), p.adjust(ms_pi$p.value[,2], method="fdr"), 
                p.adjust(hs_tfc$p.value[,2], method="fdr"), p.adjust(cs_tfc$p.value[,2], method="fdr"), p.adjust(ms_tfc$p.value[,2], method="fdr"), 
		p.adjust(tfc$p.value[,2], method="fdr")) ## Each dataset should have a 1% FDR (within the dataset).

colnames(fdr_t) <- c("HumanFDR", "ChimpFDR", "MacaqueFDR", "HumanFDR_PI", "ChimpFDR_PI", "MacaqueFDR_PI", "U2PI_H", "U2PI_C", "U2PI_M", "U2PI")

fc_t  <- data.frame(HumanFC= hs$coefficients[,2], ChimpFC= cs$coefficients[,2], MacaqueFC= ms$coefficients[,2], 
			HumanFC_PI= hs_pi$coefficients[,2], ChimpFC_PI= cs_pi$coefficients[,2], MacaqueFC_PI= ms_pi$coefficients[,2], 
			U2PI_H= hs_tfc$coefficients[,2], U2PI_C= cs_tfc$coefficients[,2], U2PI_M= ms_tfc$coefficients[,2], U2PI= tfc$coefficients[,2])

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
NROW(unique(ca[fdr_t[,1] < PVAL,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL,"mgi"])) # RHESUS
NROW(unique(ca[fdr_t[,4] < PVAL,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,5] < PVAL,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,6] < PVAL,"mgi"])) # RHESUS

NROW(unique(ca[fdr_t[,1] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # RHESUS
NROW(unique(ca[fdr_t[,4] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,5] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,6] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # RHESUS

NROW(unique(ca[fdr_t[,1] < PVAL & fdr_t[,4] < PVAL,"mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & fdr_t[,5] < PVAL,"mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & fdr_t[,6] < PVAL,"mgi"])) # MACAQUE

NROW(unique(ca[fdr_t[,1] < PVAL & fdr_t[,4] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # HUMAN
NROW(unique(ca[fdr_t[,2] < PVAL & fdr_t[,5] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # CHIMP
NROW(unique(ca[fdr_t[,3] < PVAL & fdr_t[,6] < PVAL & ca[,"annot_type"] == "gc18","mgi"])) # MACAQUE

## Write table of genes to use in GO.
writeExprChange <- function(ss, name, indx, indx_pi) {
 write.table(unique(ca[fdr_t[,indx] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".U.mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes, fdr_t, fc_t)[fdr_t[,indx]<PVAL, ], paste("chage_expr/",name,".change-U.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

 write.table(unique(ca[fdr_t[,indx_pi] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".PI.mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes, fdr_t, fc_t)[fdr_t[,indx_pi]<PVAL, ], paste("chage_expr/",name,".change-PI.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

 write.table(unique(ca[fdr_t[,indx] < PVAL & fdr_t[,indx_pi] < PVAL & ca[,"annot_type"] == "gc18","mgi"]), paste("chage_expr/",name,".U+PI.mgi.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
 write.table(cbind(genes, fdr_t, fc_t)[fdr_t[,indx] < PVAL & fdr_t[,indx_pi]<PVAL, ], paste("chage_expr/",name,".change-U+PI.tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

writeExprChange(ss, "H", 1,4)
writeExprChange(ss, "C", 2,5)
writeExprChange(ss, "M", 3,6)

isExpr <- rowMax(rpkm_df[,2:9]) > 1e-4
write.table(unique(ca[isExpr & ca[,"annot_type"] == "gc18","mgi"]), "chage_expr/exprbg.mgi.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

## Compute frequency of changes for 'expressed' genes in each class.
save.image("fdr.RData")

