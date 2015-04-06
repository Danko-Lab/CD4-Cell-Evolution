## Classifies sites as either (1) Promoter == proximal, stable.
##            :OR:            (2) Enhancer == distal (>10kb), unstable.

## Get the location of sites.
intersection_call <- paste("zcat refGene.hg19.bed.gz | awk 'BEGIN{OFS=","\t","} {print $1,$6==","+","?$2:$3-1,$6==","+","?$2+1:$3,$4,$5,$6}' | sort-bed - | bedtools closest -d -b stdin -a ", sep='"')
tss <- read.table(pipe(paste(intersection_call, "../tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")))
tss <- tss[grep("chrUn|random", tss$V1, invert=TRUE),]

## Reposition to -110bp from the highest site.
## TODO: READ COMBINED DATA. 
source("../lib/center_dREG_site.R") 
tss_centers <- center_dREG_site(bed=tss, bw_plus_file="../AllData/All_Merge/H-U_plus.bw", bw_minus_file="../AllData/All_Merge/H-U_minus.bw")
tss_centers <- center.bed(tss_centers, upstreamWindow=0, downstreamWindow=0)
tss_centers_plus <- tss_centers
tss_centers_plus[,6] <- rep("+")
tss_centers_minus <- tss_centers
tss_centers_minus[,6] <- rep("-")

## Classify sites as stable or unstable ... based on hg19.
require(stabilityHMM)
require(twoBit)
require(rqhmm)

hg19 <- "/local/storage/data/hg19/hg19.2bit"
seqLen= 400
dist <- 1000
Umin <- 0.001
Smax <- 0.1

seqs <- collectSequences(hg19, tss_centers_plus, seq.length = seqLen)
mData <- prepareData(seqs)
pl_Scores <- unstableScore(mData)

seqs <- collectSequences(hg19, tss_centers_minus, seq.length = seqLen)
mData <- prepareData(seqs)
mn_Scores <- unstableScore(mData)
#modelPathSummary(mData)

## Then classify sites as proximal/ distal based on RefGene annotations (don't want enhancer TU in the annotations).
stab_tss <- cbind(tss, tss_centers[,1:3], pl_Scores, mn_Scores)


## summary stats. 
mins <- sapply(1:NROW(pl_Scores), function(x) {min(pl_Scores[x], mn_Scores[x])})
cor.test(mins, stab_tss[,7], method="spearman")
