## Find the location of TFBS in humans.

## CREATE MOTIF tfs object.
file.twoBit_path        <- "/local/storage/data/hg19/hg19.2bit";
gencode.gtf             <- "/local/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz"
h.PROseq.plus           <- "/local/storage/projects/NHP/AllData/All_Merge/H-U_plus.bw"
h.PROseq.minus          <- "/local/storage/projects/NHP/AllData/All_Merge/H-U_minus.bw"

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

## Now scan the human genome for all TFBS.
tfbs.scanTFsite(tfs, file.twoBit_path, return.type="writedb", file.prefix= "H.motifs.th8", ncores=20, threshold= 8, threshold.type="score")


