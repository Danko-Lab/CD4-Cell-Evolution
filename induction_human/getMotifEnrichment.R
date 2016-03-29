## R script to get motif enrichment using rtfbsdb

## CREATE MOTIF tfs object.

#db <- CisBP.extdata("Homo_sapiens");
#tfs <- tfbs.createFromCisBP(db);
#
## Clustering... and expression
#tfs <- tfbs.clusterMotifs(tfs, method="apcluster", ncores=30)
#tfs <- tfbs.getExpression(tfs, file.twoBit_path, gencode.gtf, h.PROseq.plus, h.PROseq.minus, ncore=10);
#tfs <- tfbs.selectByGeneExp(tfs)
#
#save.image(file="APCluster.rdata")

library(rtfbsdb)
load("APCluster.rdata")

## Now get enriched motifs...
PVAL <- 0.01
FOLD <- 3
mTH  <- 7.5

## Do human
file.twoBit_path        <- "/local/storage/data/hg19/hg19.2bit";
gencode.gtf             <- "/local/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz"
h.PROseq.plus           <- "/local/storage/projects/NHP/AllData/All_Merge/H-U_plus.bw"
h.PROseq.minus          <- "/local/storage/projects/NHP/AllData/All_Merge/H-U_minus.bw"


enh <- read.table("results/human-changed.TREs.tsv"); enh$V4 <- 1:NROW(enh$V4)
enh <- enh[(enh$V3 - enh$V2)>0,] ## CURRENTLY BUGGED.
enh.up <- enh[enh$V9 < PVAL & enh$V5 > FOLD & !is.na(enh$V9),]
enh.dn <- enh[enh$V9 < PVAL & enh$V5 < (-1*FOLD) & !is.na(enh$V9),]
enh.unc<- enh[enh$V9 > 0.1 & abs(enh$V5) < 1 & !is.na(enh$V9), ]


tf_up <- tfbs.enrichmentTest(
        tfbs= tfs,
        file.twoBit= file.twoBit_path,
        positive.bed= enh.up,
        negative.bed= enh.unc,
        gc.correction=TRUE,
        use.cluster=TRUE,
        threshold = mTH,
        ncores = 21);
tfbs.reportEnrichment(tfs, tf_up, file.pdf="Human.TF-up.full.pdf", sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);

tf_dn <- tfbs.enrichmentTest(
        tfbs= tfs,
        file.twoBit= file.twoBit_path,
        positive.bed= enh.dn,
        negative.bed= enh.unc,
        gc.correction=TRUE,
        use.cluster=TRUE,
        threshold = mTH,
        ncores = 21);
tfbs.reportEnrichment(tfs, tf_dn, file.pdf="Human.TF-dn.full.pdf", sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);

##
## Now do the non-human primates.

## Do Chimp in Chimp genome coordinates.
chimp.twoBit_path        <- "/local/storage/data/2bit/panTro4.2bit";
enh <- read.table("results/chimp-changed.TREs.tsv"); enh$V4 <- 1:NROW(enh$V4)
enh <- enh[(enh$V3 - enh$V2)>0,] ## CURRENTLY BUGGED.
enh.up <- enh[enh$V9 < PVAL & enh$V5 > FOLD & !is.na(enh$V9),]
enh.dn <- enh[enh$V9 < PVAL & enh$V5 < (-1*FOLD) & !is.na(enh$V9),]
enh.unc<- enh[enh$V9 > 0.1 & abs(enh$V5) < 1 & !is.na(enh$V9), ]

tf_up <- tfbs.enrichmentTest(
        tfbs= tfs,
        file.twoBit= chimp.twoBit_path,
        positive.bed= enh.up,
        negative.bed= enh.unc,
        gc.correction=TRUE,
        use.cluster=TRUE,
        threshold = mTH,
        ncores = 21);
tfbs.reportEnrichment(tfs, tf_up, file.pdf="Chimp.TF-up.full.pdf", sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);

tf_dn <- tfbs.enrichmentTest(
        tfbs= tfs,
        file.twoBit= chimp.twoBit_path,
        positive.bed= enh.dn,
        negative.bed= enh.unc,
        gc.correction=TRUE,
        use.cluster=TRUE,
        threshold = mTH,
        ncores = 21);
tfbs.reportEnrichment(tfs, tf_dn, file.pdf="Chimp.TF-dn.full.pdf", sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);

## Do Rhesus in Rhesus genome coordinates.
rhesus.twoBit_path        <- "/local/storage/data/2bit/rheMac3.2bit";
enh <- read.table("results/rhesus-changed.TREs.tsv"); enh$V4 <- 1:NROW(enh$V4)
enh <- enh[(enh$V3 - enh$V2)>0,] ## CURRENTLY BUGGED.
enh.up <- enh[enh$V9 < PVAL & enh$V5 > FOLD & !is.na(enh$V9),]
enh.dn <- enh[enh$V9 < PVAL & enh$V5 < (-1*FOLD) & !is.na(enh$V9),]
enh.unc<- enh[enh$V9 > 0.1 & abs(enh$V5) < 1 & !is.na(enh$V9), ]

tf_up <- tfbs.enrichmentTest(
        tfbs= tfs,
        file.twoBit= rhesus.twoBit_path,
        positive.bed= enh.up,
        negative.bed= enh.unc,
        gc.correction=TRUE,
        use.cluster=TRUE,
        threshold = mTH,
        ncores = 21);
tfbs.reportEnrichment(tfs, tf_up, file.pdf="Rhesus.TF-up.full.pdf", sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);

tf_dn <- tfbs.enrichmentTest(
        tfbs= tfs,
        file.twoBit= rhesus.twoBit_path,
        positive.bed= enh.dn,
        negative.bed= enh.unc,
        gc.correction=TRUE,
        use.cluster=TRUE,
        threshold = mTH,
        ncores = 21);
tfbs.reportEnrichment(tfs, tf_dn, file.pdf="Rhesus.TF-dn.full.pdf", sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);


