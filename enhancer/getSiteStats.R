## This analysis focuses on lineage-specific changes that are currently active in hg19.
load("../annotations/fdr.RData")

## Currently, using dREG_TSS == promoters (proximal); dREG_ENH == enhancers (distal).
#cat ../annotations/tsssites.tsv | grep "dREG_TSS" -c
#cat ../annotations/tsssites.tsv | grep "dREG_ENH" -c



## Fraction of changes that are ... 
getFracChanges <- function(prefix) {
	dat_indx <- which(ca$type == "dREG_TSS")
	dat <- ca[dat_indx,]
	dat_all <- NROW(dat)

	# unmappable
        write.table(dat[,c(1:3)], file= pipe(" sort-bed - > tmp.bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	dat_gap <- as.double(system("bedmap --count --fraction-ref 0.5 tmp.bed ../annotations/gapMerge | grep -c -v '^0$'", intern=TRUE)) 
	print(dat_gap)

	# complete gain/loss
	as.double(system("cat ../annotations/tsssites.tsv | grep 'dREG_TSS' | sed s/NA/0/g | awk '($7>0.7 && $8<0.3 && $9<0.3)' | grep '' -c", intern=TRUE))	

	# changed in activity
	dat_chg <- sum(fdr_df$HumanFDR[dat_indx] < PVAL)
	print(dat_chg)

	# unchanged
}


