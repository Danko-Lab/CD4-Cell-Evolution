##
##
load("fdr.RData")
PVAL <- 0.01

getDifferences <- function(changeExpr, isExpr) {
  sum(changeExpr & isExpr)/ sum(isExpr)
  summary(fdr_df$type[changeExpr & isExpr])/summary(fdr_df$type[isExpr])

  fractionChanged <- list(
    "eRNA"= NROW(unique(fdr_df$name[changeExpr & (fdr_df$type=="dREG_ENH") & isExpr]))/NROW(unique(fdr_df$name[(fdr_df$type=="dREG_ENH") & isExpr])),
    "LincRNA"= NROW(unique(fdr_df$name[(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)]))/ NROW(unique(fdr_df$name[((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)])),
    "Protein Coding"=  NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)]))/NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])),
    "Antisense"=  NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="antisense" & isExpr)]))/NROW(unique(fdr_df$name[(fdr_df$type=="antisense" & isExpr)])),
    "Pseudogene"=  NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)]))/NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)])),
    "All Expressed"= NROW(unique(fdr_df$name[changeExpr & isExpr]))/NROW(unique(fdr_df$name[isExpr]))
  )

  pval <- list(
    "eRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & (fdr_df$type=="dREG_ENH") & isExpr])), NROW(unique(fdr_df$name[(fdr_df$type=="dREG_ENH") & isExpr]))),
      c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
      # c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

    "LincRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & (fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)])), NROW(unique(fdr_df$name[((fdr_df$type=="processed_transcript" | fdr_df$type=="lincRNA") & isExpr)]))), 
      c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
      #	c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

    "Protein Coding"=   fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)]))), 
      c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
      #	c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

    "Antisense"=  fisher.test(data.frame(c(sum(changeExpr & fdr_df$type=="antisense" & isExpr), sum(fdr_df$type=="antisense" & isExpr)), 
      c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
      # c(sum(changeExpr & isExpr), sum(isExpr))))$p.value,

    "Pseudogene"=  fisher.test(data.frame(c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="GERST_PG" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="GERST_PG" & isExpr)]))), 
      c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value,
      # c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value,

    "All Expressed"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr]))),
      c(NROW(unique(fdr_df$name[(changeExpr & fdr_df$type=="protein_coding" & isExpr)])), NROW(unique(fdr_df$name[(fdr_df$type=="protein_coding" & isExpr)])))))$p.value
      # c(NROW(unique(fdr_df$name[changeExpr & isExpr])), NROW(unique(fdr_df$name[isExpr])))))$p.value
  )

  pval.pg <- list(
    "eRNA"= fisher.test(data.frame(c(NROW(unique(fdr_df$name[changeExpr & (fdr_df$type=="dREG_ENH") & isExpr])), NROW(unique(fdr_df$name[(fdr_df$type=="dREG_ENH") & isExpr]))),
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
isExpr <- rowMax(rpkm_df[,2:9]) > 5e-4
sum(isExpr)/NROW(isExpr)

isExprPI<-rowMax(rpkm_df[,11:17]) > 5e-4
sum(isExprPI)/NROW(isExprPI)

getDifferences(changeExpr, isExpr)

## Untreated
hu <- getDifferences(fdr_t[,1]<PVAL, isExpr) # HUMAN
cu <- getDifferences(fdr_t[,2]<PVAL, isExpr) # CHIMP
mu <- getDifferences(fdr_t[,3]<PVAL, isExpr) # MACAQUE

## P/I
hp <- getDifferences(fdr_t[,4]<PVAL, isExprPI) # HUMAN
cp <- getDifferences(fdr_t[,5]<PVAL, isExprPI) # CHIMP
mp <- getDifferences(fdr_t[,6]<PVAL, isExprPI) # MACAQUE

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



