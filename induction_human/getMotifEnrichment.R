## R script to get motif enrichment using rtfbsdb
enh <- read.table("results/human-changed.TREs.tsv")

PVAL <- 0.01
FOLD <- 1

enh.up <- enh[enh$V9 < PVAL & enh$V5 > FOLD,]
enh.dn <- enh[enh$V9 < PVAL & enh$V5 < (-1*FOLD),]
enh.unc<- enh[enh$V9 > 0.3 & abs(enh$V5) < 0.5, ]

library(rtfbsdb)

file.twoBit_path       <- "/local/storage/data/hg19/hg19.2bit";

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db);

tf_up <- tfbs.enrichmentTest( tfs,
        file.twoBit_path,
        enh.up,
        enh.unc,
        gc.correction=TRUE,
        ncores = 21);

tfbs.reportEnrichment(tfs, tf_up, file.pdf="tfbs.comp.full.pdf", sig.only=F, report.title="TEST FULL");

save.image(file="tfbs.comp.full.rdata");

