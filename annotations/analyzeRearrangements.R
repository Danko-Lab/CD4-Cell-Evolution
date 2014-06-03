##
##
require(ggplot2)
library(reshape2)

load("fdr.RData")
PVAL <- 0.01

summary(fdr_df$fdr_min < PVAL) ## ~12k transcripts that change expression.
changeExpr <- fdr_df$fdr_min < PVAL & !is.na(fdr_df$fdr_min) & abs(fdr_df$fc_min) > FOLD #1

## Comput RPKM
isExpr <- rowMax(rpkm_df[,2:9]) > 5e-5
sum(isExpr)/NROW(isExpr)

sum(changeExpr & isExpr)/ sum(isExpr)


##############################################################################
## Compare frequency of changes in expression for transcripts inside re-arrangements.

## Re-read data, and restrict search to known genes.
use_indx <- which(ca$type == "protein_coding") # which(ca$annot_type == "gc18") 

gap <- gap[use_indx,]
fdr_df <- fdr_df[use_indx,]
rpkm_df <- rpkm_df[use_indx,]
isExpr <- isExpr[use_indx] #rep(TRUE, NROW(use_indx)) #isExpr[use_indx]
changeExpr <- changeExpr[use_indx]

ca <- ca[use_indx,]

## Shouldn't have to do this?!
#fdr_df[is.na(fdr_df)] <- 1

chimpchange <- fdr_df$ChimpFDR < PVAL & abs(fdr_df$ChimpFC) > FOLD
macaquechange <- fdr_df$MacaqueFDR < PVAL & abs(fdr_df$MacaqueFC) > FOLD

chimpData <- list( nonSyn=c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "nonSyn"])),   
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "nonSyn"]))),

		   syn=c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "syn"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "syn"]))), 

		    inv=c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "inv"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "inv"]))), 

		 none=c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "NONE"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "NONE"]))) 
		)
chimpData


macaqueData <- list( nonSyn=c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "nonSyn"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "nonSyn"]))),

                   syn=c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "syn"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "syn"]))),

                    inv=c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "inv"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "inv"]))),

                 none=c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "NONE"])),
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "NONE"])))   
                )
macaqueData


fisher.test(data.frame(c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "nonSyn"])), 
						NROW(unique(ca$mgi[isExpr & gap$V2 == "nonSyn"]))), 
						c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "NONE"])), 
						NROW(unique(ca$mgi[isExpr & gap$V2 == "NONE"]))) ))$p.value

fisher.test(data.frame(c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "nonSyn"])),   
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "nonSyn"]))),    
                                                c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "NONE"])),     
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "NONE"]))) ))$p.value

fisher.test(data.frame(c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "inv"])),   
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "inv"]))),    
                                                c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "NONE"])),     
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "NONE"]))) ))$p.value

fisher.test(data.frame(c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "inv"])),      
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "inv"]))),
                                                c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "NONE"])),     
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "NONE"]))) ))$p.value

fisher.test(data.frame(c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "syn"])),   
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "syn"]))),    
                                                c(NROW(unique(ca$mgi[chimpchange & isExpr & gap$V2 == "NONE"])),     
                                                NROW(unique(ca$mgi[isExpr & gap$V2 == "NONE"]))) ))$p.value

fisher.test(data.frame(c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "syn"])),      
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "syn"]))),
                                                c(NROW(unique(ca$mgi[macaquechange & isExpr & gap$V4 == "NONE"])),     
                                                NROW(unique(ca$mgi[isExpr & gap$V4 == "NONE"]))) ))$p.value



data_df <- data.frame(Species= factor(c(rep("Chimpanzee", 4), rep("Rhesus Macaque", 4))), 
				Type= factor(rep(c("Inversion", "None", "Different Chrom", "Same Chrom"), 2), levels=c("Different Chrom", "Same Chrom", "Inversion", "None")), 
				Value= as.double(c(c(chimpData$inv[1]/chimpData$inv[2], chimpData$none[1]/chimpData$none[2], chimpData$nonSyn[1]/chimpData$nonSyn[2], chimpData$syn[1]/chimpData$syn[2]), 
					c(macaqueData$inv[1]/macaqueData$inv[2], macaqueData$none[1]/macaqueData$none[2], macaqueData$nonSyn[1]/macaqueData$nonSyn[2], macaqueData$syn[1]/macaqueData$syn[2]) )))
data_df


b <-  ggplot(data_df, aes(x=Species, y=Value)) + xlab("") + ylab("Fraction of Transcripts Changed") + 
	geom_bar(aes(fill=Type), stat="identity",position=position_dodge(), colour="black") + scale_fill_brewer(palette = "Set3")


pdf("rearrangements.pdf")
font_size <- 16

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



