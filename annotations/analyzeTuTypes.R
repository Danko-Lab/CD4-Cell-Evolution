##
##
load("fdr.RData")
PVAL <- 0.01
EXPR <- 5e-5

toLogical <- function(indx, expLength) {
  var <- rep(FALSE, expLength)
  var[indx] <- TRUE
  return(var)
}

## Define TU types.
  indx.eRNA <- grepl("dREG_ENH", fdr_df$type)
  indx.lincRNA <- grepl("lincRNA|processed_transcript|sense_intronic|sense_overlapping", fdr_df$type)
  indx.unannot <- grepl("INTERGENIC|GENE_BadMatch|AS_BadMatch", fdr_df$type)
  indx.pseudogene.rep <- grepl("pseudogene|GERST_PG|PSEUDOGENE+REP", fdr_df$type)
  indx.protein_coding <- grepl("protein_coding", fdr_df$type)
  indx.antisense <- grepl("antisense", fdr_df$type)
  indx.uas <- grepl("ups_antisense", fdr_df$type)
  indx.srna<- grepl("sRNA", fdr_df$type)

  ## Adjust data to normalize for expression levels ...
  source("../lib/normalizeSubsample.R")
  indx.good.df <- c(12:19,21:27)
  #head(fdr_df[,indx.good.df])
  data_sums <- log(rowSums(fdr_df[,indx.good.df])+1)
  lnorm <- list(indx.eRNA= data_sums[indx.eRNA], 
		indx.lincRNA= data_sums[indx.lincRNA], 
		indx.unannot= data_sums[indx.unannot], 
		indx.pseudogene.rep= data_sums[indx.pseudogene.rep], 
		indx.protein_coding= data_sums[indx.protein_coding],
		indx.antisense= data_sums[indx.antisense],
		indx.uas= data_sums[indx.uas],
		indx.srna= data_sums[indx.srna])
  ns <- norm.subsample.n(lnorm, plot.cdf=TRUE, dist=data_sums[indx.eRNA])
  indx.eRNA <- toLogical(which(indx.eRNA)[ns[[1]]])
  indx.lincRNA <- toLogical(which(indx.lincRNA)[ns[[2]]])
  indx.unannot <- toLogical(which(indx.unannot)[ns[[3]]])
  indx.pseudogene.rep <- toLogical(which(indx.pseudogene.rep)[ns[[4]]])
  indx.protein_coding <- toLogical(which(indx.protein_coding)[ns[[5]]])
  indx.antisense <- toLogical(which(indx.antisense)[ns[[6]]])
  indx.uas <- toLogical(which(indx.uas)[ns[[7]]])
  indx.srna <- toLogical(which(indx.srna)[ns[[8]]])

## UNDO THIS ONE!
  indx.pseudogene.rep <- grepl("pseudogene|GERST_PG|PSEUDOGENE+REP", fdr_df$type)


getDifferences <- function(changeExpr, isExpr) {
  sum(changeExpr & isExpr)/ sum(isExpr)
  print(summary(fdr_df$type[changeExpr & isExpr])/summary(fdr_df$type[isExpr]))

  fractionChanged <- list(
    "eRNA"= NROW(unique(fdr_df$name[changeExpr & indx.eRNA & isExpr]))/NROW(unique(fdr_df$name[indx.eRNA & isExpr])),
    "LincRNA"= NROW(unique(fdr_df$name[changeExpr & indx.lincRNA & isExpr]))/ NROW(unique(fdr_df$name[indx.lincRNA & isExpr])),
    "Unannot"= NROW(unique(fdr_df$name[changeExpr & indx.unannot & isExpr]))/ NROW(unique(fdr_df$name[indx.unannot & isExpr])),
    "sRNA"= NROW(unique(fdr_df$name[changeExpr & indx.srna & isExpr]))/ NROW(unique(fdr_df$name[indx.srna & isExpr])),
    "Protein Coding"=  NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr]))/NROW(unique(fdr_df$name[indx.protein_coding & isExpr])),
    "Antisense"=  NROW(unique(fdr_df$name[changeExpr & indx.antisense & isExpr]))/NROW(unique(fdr_df$name[indx.antisense & isExpr])),
    "ups_Antis"= NROW(unique(fdr_df$name[changeExpr & indx.uas & isExpr]))/NROW(unique(fdr_df$name[indx.uas & isExpr])),
    "Pseudogene"=  NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr]))/NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr])),
    "All Expressed"= NROW(unique(fdr_df$name[changeExpr & isExpr]))/NROW(unique(fdr_df$name[isExpr]))
  )

  pval <- list(
    "eRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.eRNA & isExpr])), NROW(unique(fdr_df$name[indx.eRNA & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "LincRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.lincRNA & isExpr])), NROW(unique(fdr_df$name[indx.lincRNA & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "Unannot"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.unannot & isExpr])), NROW(unique(fdr_df$name[indx.unannot & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "sRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.srna & isExpr])), NROW(unique(fdr_df$name[indx.srna & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "Protein Coding"=   fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "Antisense"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.antisense & isExpr])), NROW(unique(fdr_df$name[indx.antisense & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "ups_Antis"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.uas & isExpr])), NROW(unique(fdr_df$name[indx.uas & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "Pseudogene"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))), 
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value,

    "All Expressed"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr])))))$p.value
  )

  pval.pg <- list(
    "eRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.eRNA & isExpr])), NROW(unique(fdr_df$name[indx.eRNA & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,

    "LincRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.lincRNA & isExpr])), NROW(unique(fdr_df$name[indx.lincRNA & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,

    "Unannot"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.unannot & isExpr])), NROW(unique(fdr_df$name[indx.unannot & isExpr]))),       
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,  

    "sRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.srna & isExpr])), NROW(unique(fdr_df$name[indx.srna & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,  

    "Protein Coding"=   fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.protein_coding & isExpr])), NROW(unique(fdr_df$name[indx.protein_coding & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,

    "Antisense"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.antisense & isExpr])), NROW(unique(fdr_df$name[indx.antisense & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,

    "ups_Antis"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.uas & isExpr])), NROW(unique(fdr_df$name[indx.uas & isExpr]))),                
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,   

    "Pseudogene"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value,

    "All Expressed"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr]))),
      c(NROW(unique(fdr_df$name[changeExpr & indx.pseudogene.rep & isExpr])), NROW(unique(fdr_df$name[indx.pseudogene.rep & isExpr]))) ))$p.value
  )

  require(ggplot2)
  library(reshape2)
  data_df <- data.frame(Type= names(fractionChanged), Changed= as.double(unlist(fractionChanged)), Pval= as.double(unlist(pval)), Pval.pg= as.double(unlist(pval.pg)))
#  data_df$Type <- factor(data_df$Type, levels= data_df$Type[c(1,2,4,5,3,6)])

  # http://learnr.wordpress.com/2009/03/17/ggplot2-barplots/
  # http://stackoverflow.com/questions/10352894/barplot-using-ggplot2
  #a <- ggplot(data_df, aes(Type, Changed)) + xlab("") + ylab("Fraction of Transcripts Changed") +
  #  geom_bar(aes(fill=Type), stat="identity") + 
  #  scale_fill_brewer(palette = "Set1")

  return(data_df)
}

##############################################################################
## Compare frequency of changes in expression for genes, and non-coding RNAs.

summary(fdr_df$fdr_min < PVAL) ## ~12k transcripts that change expression.
changeExpr <- fdr_df$fdr_min < PVAL & !is.na(fdr_df$fdr_min) & abs(fdr_df$fc_min) > FOLD #1

## Comput RPKM
isExpr <- rowMax(rpkm_df[,2:9]) > EXPR
sum(isExpr)/NROW(isExpr)

isExprPI<-rowMax(rpkm_df[,11:17]) > EXPR
sum(isExprPI)/NROW(isExprPI)

getDifferences(changeExpr, isExpr)

## Untreated
hu <- getDifferences(fdr_t[,1]<PVAL, isExpr); hu # HUMAN
cu <- getDifferences(fdr_t[,2]<PVAL, isExpr); cu # CHIMP
mu <- getDifferences(fdr_t[,3]<PVAL, isExpr); mu # MACAQUE

## P/I
hp <- getDifferences(fdr_t[,4]<PVAL, isExprPI); hp # HUMAN
cp <- getDifferences(fdr_t[,5]<PVAL, isExprPI); cp # CHIMP
mp <- getDifferences(fdr_t[,6]<PVAL, isExprPI); mp # MACAQUE


## Now actually writes out barplots.
require(ggplot2)
library(reshape2)

pdf("tuTypeChanges.pdf")
font_size <- 16

  a <- ggplot(hu, aes(Type, Changed)) + xlab("") + ylab("Fraction of Transcripts Changed") +
    geom_bar(aes(fill=Type), stat="identity") +
    scale_fill_brewer(palette = "Set1")

  b <- ggplot(hu, aes(Type, Changed)) + xlab("") + ylab("Fraction of Transcripts Changed") +
    geom_bar(aes(fill=Type), stat="identity") +
    scale_fill_brewer(palette = "Set1")


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

b
	 
dev.off()



