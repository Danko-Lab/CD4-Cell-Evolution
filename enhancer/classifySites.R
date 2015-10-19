## Classifies sites as either (1) Promoter == proximal, stable.
##            :OR:            (2) Enhancer == distal (>10kb), unstable.

## Get the location of sites.
tss <- read.table("../tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")

## Reposition to -110bp from the highest site.
## TODO: READ COMBINED DATA.  

## Classify sites as stable or unstable ... based on hg19.
require(stabilityHMM)
require(twoBit)
require(rqhmm)

hg19 <- "~/storage/data/hg19/hg19.2bit"
seqLen= 1000

seqs <- collectSequences(hg19, tss, seq.length = seqLen)
mData <- prepareData(seqs)
pl_Scores <- unstableScore(mData)
#modelPathSummary(mData)


## Then classify sites as proximal/ distal based on RefGene annotations (don't want enhancer TU in the annotations).


## Sanity check.
refGene <- read.table("~/tmp/refGene.hg19.chr21.bed.gz")
hg19 <- "~/storage/data/hg19/hg19.2bit"
seqLen= 1000

seqs <- collectSequences(hg19, refGene, seq.length = seqLen)
mData <- prepareData(seqs)
gene_Scores <- unstableScore(mData)

hist(gene_Scores)

## Get outliers.
cbind(refGene, gene_Scores)[gene_Scores > 0.8,]

## Checked by hand, they look like the real deal.
