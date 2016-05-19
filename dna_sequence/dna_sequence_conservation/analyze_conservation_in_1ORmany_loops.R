## Compare DNA sequence conservation in different classes of RE.
require(bigWig)
require(vioplot)
require(Hmisc)
require(stats)

source("../../lib/getOverlap.R")
source("../../lib/densScatterplot.R")

## Read in dREG-HD peaks.
peaksize <- 5 # 100 ## control peak size.

bed <- read.table("dREG-HD.2compare.bed"); #bed <- bed[bed[,1]=="chr1",]
bed_data <- bed
bed_data <- center.bed(bed, peaksize, peaksize)

bed1 <- read.table("oneTRE.dREG-HD.2compare.bed");
bed_data1 <- bed1
bed_data1 <- center.bed(bed1, peaksize, peaksize)

bed2 <- read.table("moreTREs.dREG-HD.2compare.bed");
bed_data2 <- bed2
bed_data2 <- center.bed(bed2, peaksize, peaksize)

## Read in MAFs
#con <- load.bigWig("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw")
con <- load.bigWig("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw")
word <- "100way"#"Primate"

## Get mean conservation for each h_bed.
mean1_con <- bed.region.bpQuery.bigWig(con, bed_data1[,1:3], op="avg")# op="max")
mean2_con <- bed.region.bpQuery.bigWig(con, bed_data2[,1:3], op="avg")# op="max")
mean_con <- bed.region.bpQuery.bigWig(con, bed_data[,1:3], op="avg")# op="max")


## General.
boxplot(mean1_con, mean2_con)

pdf(paste("OneVMany.DNAsequence.",word,".phyloP.Conservation.pdf", sep=""))
 ld<- 1; xlim_s=c(-0.75,3) #xlim_s=c(-0.75, 0.5) #c(-0.4, 0.5) #PRIMATE

 plot(cd.cdf(mean1_con), col="dark red", lwd=ld, xlim=xlim_s, type="l")
 plot(cd.cdf(mean2_con), col="dark green", add=TRUE, lwd=ld, type="l")
 plot(cd.cdf(mean_con),  col="dark gray", add=TRUE, lwd=ld, type="l")
 
dev.off()


wilcox.test(mean1_con, mean2_con)
ks.test(mean1_con, mean2_con)

########################
## Add loop information.

## DONE!
