##
##

load("fdr.RData")

##############################################################################
## Compare frequency of changes in expression for genes, and non-coding RNAs.

summary(fdr_df$fdr_min < PVAL) ## ~12k transcripts that change expression.
changeExpr <- fdr_df$fdr_min < PVAL & !is.na(fdr_df$fdr_min) & abs(fdr_df$fc_min) > FOLD #1

## Comput RPKM
isExpr <- rowMax(rpkm_df[,2:9]) > 5e-4
sum(isExpr)/NROW(isExpr)

sum(changeExpr & isExpr)/ sum(isExpr)
summary(fdr_df$type[changeExpr & isExpr])/summary(fdr_df$type[isExpr])

fractionChanged <- list(
	"eRNA"= NROW(unique(fdr_df$name[changeExpr & (fdr_df$type=="PromEnh") & isExpr]))/NROW(unique(fdr_df$name[(fdr_df$type=="PromEnh") & isExpr])),
	"LincRNA"= NROW(unique(fdr_df$name[(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)]))/ NROW(unique(fdr_df$name[((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)])),
	"Protein Coding"=  NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)]))/NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])),
	"Antisense"=  NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="antisense" & isExpr)]))/NROW(unique(fdr_df$name[(fdr_df$type=="antisense" & isExpr)])),
	"Pseudogene"=  NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)]))/NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)])),
	"All Expressed"= NROW(unique(fdr_df$name[changeExpr & isExpr]))/NROW(unique(fdr_df$name[isExpr]))
)

pval <- list(
        "eRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & (fdr_df$type=="PromEnh") & isExpr])), NROW(unique(fdr_df$name[(fdr_df$type=="PromEnh") & isExpr]))),  
                                c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
#                                c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

	"LincRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)])), NROW(unique(fdr_df$name[((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)]))), 
                                c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
#				c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

	"Protein Coding"=   fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)]))), 
                                c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
#				c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

	"Antisense"=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="antisense" & isExpr), sum(fdr_df$type=="antisense" & isExpr)), 
                                c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
#				c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,

	"Pseudogene"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)]))), 
                                c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
#                                c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

	"All Expressed"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr]))),
                                c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value
#                                c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value
)


pval.pg <- list(
        "eRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & (fdr_df$type=="PromEnh") & isExpr])), NROW(unique(fdr_df$name[(fdr_df$type=="PromEnh") & isExpr]))),
                               c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)]))) ))$p.value,

        "LincRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)])), NROW(unique(fdr_df$name[((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)]))),
                               c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)])))  ))$p.value,

        "Protein Coding"=   fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)]))),
                               c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)])))  ))$p.value,

        "Antisense"=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="antisense" & isExpr), sum(fdr_df$type=="antisense" & isExpr)),
                               c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)])))   ))$p.value,

        "Pseudogene"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)]))),
                               c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)])))  ))$p.value,

        "All Expressed"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr]))),
                               c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)]))) ))$p.value
)



require(ggplot2)
library(reshape2)
data_df <- data.frame(Type= names(fractionChanged), Changed= as.double(unlist(fractionChanged)), Pval= as.double(unlist(pval)), Pval.pg= as.double(unlist(pval.pg)))
data_df$Type <- factor(data_df$Type, levels= data_df$Type[c(1,2,4,5,3,6)])
data_df

# http://learnr.wordpress.com/2009/03/17/ggplot2-barplots/
# http://stackoverflow.com/questions/10352894/barplot-using-ggplot2
a <- ggplot(data_df, aes(Type, Changed)) + xlab("") + ylab("Fraction of Transcripts Changed") +
	geom_bar(aes(fill=Type), stat="identity") + 
	scale_fill_brewer(palette = "Set1")

#a	 


##############################################################################
## Compare frequency of changes in expression for transcripts inside re-arrangements.

## Re-read data, and restrict search to known genes.
ca <- read.table("countall.tsv") ## HRM.  Re-reading here is ugly.  Correct, however.  No use in changing it.
gap <- read.table("genes.inGap")[!is.na(ca[,10]),]
ca <- ca[!is.na(ca[,10]),]
fdr_df <- fdr_df[1:NROW(ca),]
rpkm_df <- rpkm_df[1:NROW(ca),]
isExpr <- isExpr[1:NROW(ca)]
changeExpr <- changeExpr[1:NROW(ca)]

## Shouldn't have to do this?!
fdr_df[is.na(fdr_df)] <- 1

chimpchange <- fdr_df$ChimpFDR < PVAL & abs(fdr_df$ChimpFC) > FOLD
macaquechange <- fdr_df$MacaqueFDR < PVAL & abs(fdr_df$MacaqueFC) > FOLD

summary(gap$V2[chimpchange])/ summary(gap$V2)
summary(gap$V4[macaquechange])/ summary(gap$V4)

summary(gap$V2[chimpchange & isExpr])/ summary(gap$V2[isExpr])
summary(gap$V4[macaquechange & isExpr])/ summary(gap$V4[isExpr])

fisher.test(data.frame(c(sum(gap$V2[chimpchange & isExpr] == "nonSyn"), 
						sum(gap$V2[isExpr] == "nonSyn")), 
						c(sum(gap$V2[chimpchange & isExpr] == "NONE"), 
						sum(gap$V2[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V4[macaquechange & isExpr] == "nonSyn"), 
						sum(gap$V4[isExpr] == "nonSyn")), 
						c(sum(gap$V4[macaquechange & isExpr] == "NONE"), 
						sum(gap$V4[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V2[chimpchange & isExpr] == "inv"), 
						sum(gap$V2[isExpr] == "inv")), 
						c(sum(gap$V2[chimpchange & isExpr] == "NONE"), 
						sum(gap$V2[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V4[macaquechange & isExpr] == "inv"), 
						sum(gap$V4[isExpr] == "inv")), 
						c(sum(gap$V4[macaquechange & isExpr] == "NONE"), 
						sum(gap$V4[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V2[chimpchange & isExpr] == "syn"), 
						sum(gap$V2[isExpr] == "syn")), 
						c(sum(gap$V2[chimpchange & isExpr] == "NONE"), 
						sum(gap$V2[isExpr] == "NONE"))))$p.value

fisher.test(data.frame(c(sum(gap$V4[macaquechange & isExpr] == "syn"), 
						sum(gap$V4[isExpr] == "syn")), 
						c(sum(gap$V4[macaquechange & isExpr] == "NONE"), 
						sum(gap$V4[isExpr] == "NONE"))))$p.value

						

data_df <- data.frame(Species= factor(c(rep("Chimpanzee", 4), rep("Rhesus Macaque", 4))), 
				Type= factor(rep(c("Inversion", "None", "Different Chrom", "Same Chrom"), 2), levels=c("Different Chrom", "Same Chrom", "Inversion", "None")), 
				Value= as.double(c(c(summary(gap$V2[chimpchange & isExpr])/ summary(gap$V2[isExpr])), 
					c(summary(gap$V4[macaquechange & isExpr])/ summary(gap$V4[isExpr])))))

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



