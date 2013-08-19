load("fdr.RData")

summary(fdr_df$fdr_min < 0.05) ## ~12k transcripts that change expression.
changeExpr <- fdr_df$fdr_min < 0.05

## Comput RPKM
isExpr <- rowMeans(rpkm_df) > 1e-4

sum(changeExpr & isExpr)/ sum(isExpr)
summary(fdr_df$type[changeExpr & isExpr])/summary(fdr_df$type[isExpr])

fractionChanged <- list(
	lincRNA= sum(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)/sum((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr),
	protein_coding=  sum(changeExpr & fdr_df$type=="protein_coding" & isExpr)/sum(fdr_df$type=="protein_coding" & isExpr),
	antisense=  sum(changeExpr & fdr_df$type=="antisense" & isExpr)/sum(fdr_df$type=="antisense" & isExpr),
	pseudogene=  sum(changeExpr & fdr_df$type=="pseudogene" & isExpr)/sum(fdr_df$type=="pseudogene" & isExpr),
	all_expressed= sum(changeExpr & isExpr)/sum(isExpr)
)

pval <- list(
	lincRNA= fisher.test(data.frame(c(sum(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr), sum((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,
	protein_coding=   fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="protein_coding" & isExpr), sum(fdr_df$type=="protein_coding" & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,
	antisense=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="antisense" & isExpr), sum(fdr_df$type=="antisense" & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,
	pseudogene=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="pseudogene" & isExpr), sum(fdr_df$type=="pseudogene" & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value
)

print(paste("non-coding RNAs:",fractionChanged$lincRNA, "p = ", pval$lincRNA))
print(paste("Protein coding:", fractionChanged$protein_coding, "p = ", pval$protein_coding))
print(paste("Antisense:", fractionChanged$antisense, "p = ", pval$antisense))
print(paste("Pseudogene:",fractionChanged$pseudogene, "p = ", pval$pseudogene))
print(paste("All expressed:", sum(changeExpr & isExpr)/sum(isExpr) ))

require(ggplot2)
library(reshape2)
m_f <- melt(fractionChanged)
names(m_f) <- c("Changed", "Type") 
# http://learnr.wordpress.com/2009/03/17/ggplot2-barplots/
# http://stackoverflow.com/questions/10352894/barplot-using-ggplot2
ggplot(m_f, aes(Type, Changed), panel.background = element_blank()) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Set1")

a <- ggplot(melt(fractionChanged))

pdf("FractionChagned.pdf")
 barplot(unlist(fractionChanged))
dev.off()


##############################################################################
## Compare frequency of changes in expression for transcripts inside re-arrangements.

ca <- read.table("countall.tsv")
gap <- read.table("genes.inGap")[!is.na(ca[,10]),]
ca <- ca[!is.na(ca[,10]),]

summary(gap$V2[fdr_df$ChimpFDR< 0.05])/ summary(gap$V2)
summary(gap$V4[fdr_df$MacaqueFDR< 0.05])/ summary(gap$V4)

summary(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr])/ summary(gap$V2[isExpr])
summary(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr])/ summary(gap$V4[isExpr])

fisher.test(data.frame(c(sum(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr] == "nonSyn"), 
						sum(gap$V2[isExpr] == "nonSyn")), 
						c(sum(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr] == "NONE"), 
						sum(gap$V2[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr] == "nonSyn"), 
						sum(gap$V4[isExpr] == "nonSyn")), 
						c(sum(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr] == "NONE"), 
						sum(gap$V4[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr] == "inv"), 
						sum(gap$V2[isExpr] == "inv")), 
						c(sum(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr] == "NONE"), 
						sum(gap$V2[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr] == "inv"), 
						sum(gap$V4[isExpr] == "inv")), 
						c(sum(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr] == "NONE"), 
						sum(gap$V4[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr] == "syn"), 
						sum(gap$V2[isExpr] == "syn")), 
						c(sum(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr] == "NONE"), 
						sum(gap$V2[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr] == "syn"), 
						sum(gap$V4[isExpr] == "syn")), 
						c(sum(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr] == "NONE"), 
						sum(gap$V4[isExpr] == "NONE"))))$p.value

