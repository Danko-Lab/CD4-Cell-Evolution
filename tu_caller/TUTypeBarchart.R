# Barchart of biotypes.

tu <- read.table("annot/TU.final.bed.gz")

tu$V7[grepl("lincRNA|processed_transcript|sense_overlapping|sense_intronic", as.character(tu$V7))] <- "lincRNA"
tu$V7[grepl("pseudogene|PSEUDOGENE|GERST_PG", as.character(tu$V7))] <- "PSEUDOGENE+REP"
tu$V7[grepl("protein_coding|gene", as.character(tu$V7))] <- "protein_coding"
tu$V7[grepl("INTERGENIC|BadMatch", as.character(tu$V7))] <- "INTERGENIC"

tu$V7 <- as.factor(as.character(tu$V7))

summary(tu$V7)

  require(ggplot2)
  library(reshape2)

a <- ggplot(tu$V7) + ylab("Fraction of Transcripts Changed") +
    geom_bar(aes(fill=Type), stat="identity") +
    scale_fill_brewer(palette = "Set1")

