## Classifies sites as either (1) Promoter == proximal, stable.
##            :OR:            (2) Enhancer == distal (>10kb), unstable.

## Get the location of sites.
tss <- read.table("HCM-U-PI.dREG-tss-clusters.dist.tsv")
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
<<<<<<< HEAD

hg19 <- "~/storage/data/hg19/hg19.2bit"
seqLen= 1000

seqs <- collectSequences(hg19, tss, seq.length = seqLen)
mData <- prepareData(seqs)
pl_Scores <- unstableScore(mData)
#modelPathSummary(mData)

=======

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
>>>>>>> d5c0702fe42f1d98ce3524d8d0a2a20b91cb71a9

## Then classify sites as proximal/ distal based on RefGene annotations (don't want enhancer TU in the annotations).
stab_tss <- cbind(tss, tss_centers[,1:3], pl_Scores, mn_Scores)


## summary stats. 
mins <- sapply(1:NROW(pl_Scores), function(x) {min(pl_Scores[x], mn_Scores[x])})
cor.test(mins, stab_tss[,7], method="spearman")

## classes ...
stable_prox <- (stab_tss$V13 < dist & mins < Smax) & !is.na(mins)
unstable_dist <- stab_tss$V13 > dist & mins > Umin & !is.na(mins)

sum(stable_prox, na.rm=TRUE)
sum(unstable_dist, na.rm=TRUE)

require(vioplot)
vioplot(stab_tss[stable_prox & !is.na(stab_tss[,7]),7], stab_tss[unstable_dist & !is.na(stab_tss[,7]),7], names=c("Stable-Proximal", "Unstable-Distal"))

## Remove NAs.
for(x in 7:13){ stab_tss[is.na(stab_tss[,x]),x] <- 0 }

## Of sites in human, complete gains.
sum(stab_tss[stable_prox,7] > high & stab_tss[stable_prox, 8] < low & stab_tss[stable_prox, 9] < low)/ sum(stab_tss[stable_prox, 7] > high)
sum(stab_tss[unstable_dist,7] > high & stab_tss[unstable_dist, 8] < low & stab_tss[unstable_dist, 9] < low)/ sum(stab_tss[unstable_dist, 7] > high)

## Changes.
sum(stab_tss[stable_prox,7] > 0.7, na.rm=TRUE)/ sum(!is.na(stab_tss[stable_prox,7]))
sum(stab_tss[unstable_dist,7] > 0.7, na.rm=TRUE)/ sum(!is.na(stab_tss[unstable_dist,7]))


## Scatterplots of Hu-Mu, Hu-Hpi.
plot(stab_tss[stable_prox,7]-stab_tss[stable_prox,9], stab_tss[stable_prox,7]-stab_tss[stable_prox,10], xlab="species")
plot(stab_tss[unstable_dist,7]-stab_tss[unstable_dist,9], stab_tss[unstable_dist,7]-stab_tss[unstable_dist,10], xlab="species")
source("../lib/densScatterplot.R")
densScatterplot(stab_tss[stable_prox,7]-stab_tss[stable_prox,9], stab_tss[stable_prox,7]-stab_tss[stable_prox,10], xlab="species")
densScatterplot(stab_tss[unstable_dist,7]-stab_tss[unstable_dist,9], stab_tss[unstable_dist,7]-stab_tss[unstable_dist,10], xlab="species")



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
