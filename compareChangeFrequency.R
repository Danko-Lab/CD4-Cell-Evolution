
##############################################################################
## Compare frequency of changes in expression for genes, and non-coding RNAs.
load("fdr.RData")

summary(fdr_df$fdr_min < 0.05) ## ~12k transcripts that change expression.
changeExpr <- fdr_df$fdr_min < 0.05

## Comput RPKM
isExpr <- rowMeans(rpkm_df) > 1e-4

sum(changeExpr & isExpr)/ sum(isExpr)
summary(fdr_df$type[changeExpr & isExpr])/summary(fdr_df$type[isExpr])

fractionChanged <- list(
	"LincRNA"= sum(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)/sum((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr),
	"Protein Coding"=  sum(changeExpr & fdr_df$type=="protein_coding" & isExpr)/sum(fdr_df$type=="protein_coding" & isExpr),
	"Antisense"=  sum(changeExpr & fdr_df$type=="antisense" & isExpr)/sum(fdr_df$type=="antisense" & isExpr),
	"Pseudogene"=  sum(changeExpr & fdr_df$type=="pseudogene" & isExpr)/sum(fdr_df$type=="pseudogene" & isExpr),
	"All Expressed"= sum(changeExpr & isExpr)/sum(isExpr)
)

pval <- list(
	"LincRNA"= fisher.test(data.frame(c(sum(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr), sum((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,
	"Protein Coding"=   fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="protein_coding" & isExpr), sum(fdr_df$type=="protein_coding" & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,
	"Antisense"=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="antisense" & isExpr), sum(fdr_df$type=="antisense" & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,
	"Pseudogene"=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="pseudogene" & isExpr), sum(fdr_df$type=="pseudogene" & isExpr)), c(sum(changeExpr & isExpr), sum(isExpr))))$p.value
)

#print(paste("Non-coding RNAs:",fractionChanged$"LincRNA", "p = ", pval$"LincRNA"))
#print(paste("Protein coding:", fractionChanged$"Protein Coding", "p = ", pval$"Protein Coding"))
#print(paste("Antisense:", fractionChanged$"Antisense", "p = ", pval$"Antisense"))
#print(paste("Pseudogene:",fractionChanged$"Pseudogene", "p = ", pval$"Pseudogene"))
#print(paste("All expressed:", sum(changeExpr & isExpr)/sum(isExpr) ))

require(ggplot2)
library(reshape2)
data_df <- data.frame(Type= names(fractionChanged), Changed= as.double(unlist(fractionChanged)))
data_df$Type <- factor(data_df$Type, levels= data_df$Type[c(1,3,4,2,5)])
data_df

# http://learnr.wordpress.com/2009/03/17/ggplot2-barplots/
# http://stackoverflow.com/questions/10352894/barplot-using-ggplot2
a <- ggplot(data_df, aes(Type, Changed)) + xlab("") + ylab("Fraction of Transcripts Changed") +
	geom_bar(aes(fill=Type), stat="identity") + 
	scale_fill_brewer(palette = "Set1")

#a	 


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

						

data_df <- data.frame(Species= factor(c(rep("Chimpanzee", 4), rep("Rhesus Macaque", 4))), 
				Type= factor(rep(c("Inversion", "None", "Different Chrom", "Same Chrom"), 2), levels=c("Different Chrom", "Same Chrom", "Inversion", "None")), 
				Value= as.double(c(c(summary(gap$V2[fdr_df$ChimpFDR< 0.05 & isExpr])/ summary(gap$V2[isExpr])), 
					c(summary(gap$V4[fdr_df$MacaqueFDR< 0.05 & isExpr])/ summary(gap$V4[isExpr])))))

b <-  ggplot(data_df, aes(x=Species, y=Value)) + xlab("") + ylab("Fraction of Transcripts Changed") + 
	geom_bar(aes(fill=Type), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")


pdf("FractionChagned.pdf")
font_size <- 16

theme_update(axis.text.x = element_text(angle = 90, hjust = 1, size=font_size), 
     axis.text.y = element_text(size=font_size), 
	 axis.title = element_text(size=font_size),
	 panel.border = element_rect(colour = "black", fill= "transparent", size=2),
	 panel.grid.major = element_line(colour = "grey90"),
     panel.grid.minor = element_blank(), 
	 panel.background = element_blank(),
     axis.ticks = element_line(size=1), 
	 legend.position = "none")

a

theme_update(axis.text.x = element_text(angle = 0, hjust = 0.5, size=font_size), 
     axis.text.y = element_text(size=font_size), 
	 axis.title = element_text(size=font_size),
	 panel.border = element_rect(colour = "black", fill= "transparent", size=2),
	 panel.grid.major = element_line(colour = "grey90"),
     panel.grid.minor = element_blank(), 
	 panel.background = element_blank(),
     axis.ticks = element_line(size=1), 
	 legend.position = "right", legend.text = element_text(size=font_size), legend.title = element_text(size=font_size))

b
	 
dev.off()

