################################################################################
## ID changes in GE with Limma Voom
PVAL <- 0.05
FOLD <- 0 #log(10,2)

## Read count data.
source("readData.R")

## Use only untreated, and get switch to RPKM
ca <- ca[,c(1:10,indx.good)]

counts <- ca[,c(12:19,21:27)] #rpkm_df[,c(2:9,11:17)] #ca[,c(12:19,21:27)]
genes  <- ca[,c(1:10)]
species  <- c("H","H","H","C","C","M","M","M",    "H","H","H","C","C","M","M")
treatment<- c(rep("U",8), rep("PI",7))

## Indices for different combinations of data.
dindx_u  <- c(1:8)
dindx_pi <- c(9:15)
dindx_h  <- c(1:3,9:11)
dindx_c  <- c(4:5,12:13)
dindx_m  <- c(6:8,14:15)

## Import function from lib.
source("../lib/runLimmaQuantile.R")

## Treating alternatie primates as a single 'group' in either the U or PI condition separately.
condition <- species
condition[species == "H"] <- "Human"
condition[species != "H"] <- "NHP"
pdf("MAPlotHuman.pdf")
hs    <- runLimmaQuantile(counts[,dindx_u], condition[dindx_u], genes, condA="Human", condB="NHP", q.cut=PVAL, lfc=FOLD, plotMA=TRUE, remove.pc=3)
dev.off()
hs_pi <- runLimmaQuantile(counts[,dindx_pi], condition[dindx_pi], genes, condA="Human", condB="NHP", q.cut=PVAL, lfc=FOLD, remove.pc=3)

condition <- species
condition[species == "C"] <- "Chimp"
condition[species != "C"] <- "PRIMATES"
cs <- runLimmaQuantile(counts[,dindx_u], condition[dindx_u], genes, condA="Chimp", condB="PRIMATES", q.cut=PVAL, lfc=FOLD, remove.pc=3)
cs_pi <- runLimmaQuantile(counts[,dindx_pi], condition[dindx_pi], genes, condA="Chimp", condB="PRIMATES", q.cut=PVAL, lfc=FOLD, remove.pc=3)

condition <- species
condition[species == "M"] <- "RM"
condition[species != "M"] <- "PRIMATES"
ms <- runLimmaQuantile(counts[,dindx_u], condition[dindx_u], genes, condA="RM", condB="PRIMATES", q.cut=PVAL, lfc=FOLD, remove.pc=3) 
ms_pi <- runLimmaQuantile(counts[,dindx_pi], condition[dindx_pi], genes, condA="RM", condB="PRIMATES", q.cut=PVAL, lfc=FOLD, remove.pc=3)

condition <- treatment
hs_tfc <- runLimmaQuantile(counts[,dindx_h], condition[dindx_h], genes, condA="U", condB="PI", q.cut=PVAL, lfc=FOLD) 
cs_tfc <- runLimmaQuantile(counts[,dindx_c], condition[dindx_c], genes, condA="U", condB="PI", q.cut=PVAL, lfc=FOLD) 
ms_tfc <- runLimmaQuantile(counts[,dindx_m], condition[dindx_m], genes, condA="U", condB="PI", q.cut=PVAL, lfc=FOLD)
tfc    <- runLimmaQuantile(counts, condition, genes, condA="U", condB="PI", q.cut=PVAL, lfc=FOLD)

## Append tables.  
#fdr_t <- matrix(p.adjust(c(hs$p.value[,2], cs$p.value[,2], ms$p.value[,2], hs_pi$p.value[,2], cs_pi$p.value[,2], ms_pi$p.value[,2], 
#				hs_tfc$p.value[,2], cs_tfc$p.value[,2], ms_tfc$p.value[,2], tfc$p.value[,2]), method="fdr"), ncol=10)

fdr_t <- cbind(hs$tab$adj.P.Val, cs$tab$adj.P.Val, ms$tab$adj.P.Val, 
		hs_pi$tab$adj.P.Val, cs_pi$tab$adj.P.Val, ms_pi$tab$adj.P.Val, 
                hs_tfc$tab$adj.P.Val, cs_tfc$tab$adj.P.Val, ms_tfc$tab$adj.P.Val, 
		tfc$tab$adj.P.Val) ## Each dataset should have a 1% FDR (within the dataset).

colnames(fdr_t) <- c("HumanFDR", "ChimpFDR", "MacaqueFDR", "HumanFDR_PI", "ChimpFDR_PI", "MacaqueFDR_PI", "U2PIFDR_H", "U2PIFDDR_C", "U2PIFDR_M", "U2PIFDR")

fc_t  <- cbind(hs$tab$logFC, cs$tab$logFC, ms$tab$logFC,  
                hs_pi$tab$logFC, cs_pi$tab$logFC, ms_pi$tab$logFC,  
                hs_tfc$tab$logFC, cs_tfc$tab$logFC, ms_tfc$tab$logFC,
                tfc$tab$logFC) 

colnames(fc_t) <- c("HumanFC", "ChimpFC", "MacaqueFC", "HumanFC_PI", "ChimpFC_PI", "MacaqueFC_PI", "U2PIFC_H", "U2PIFC_C", "U2PIFC_M", "U2PIFC")

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

