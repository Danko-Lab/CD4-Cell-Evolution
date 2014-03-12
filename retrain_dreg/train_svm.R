## train_svm.R -- Trains an SVM to detect regulatory regions based on DNAse-1 data.
##

require(featureDetector)

## Read PRO-seq data.
ps_plus_path  <- "H-U_plus.bw"
ps_minus_path <- "H-U_minus.bw"

## Get positive regions.
dnase <- read.table("dnase1.peaks_peaks.narrowPeak")
extra_enrich_bed <- read.table("GencodeMerge.IntersectOpStrand.bed")
allow_bed <- read.table("CD4.chromHMM.Ernst2010.hg19.Prom.Enh.bed")

## Train the SVM.
inf_positions <- lapply(c(1:NROW(ps_plus_path)), function(x) {get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE)}) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))
asvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, dnase, n_train=75000, n_eval=0, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

remove(inf_positions, ps_plus_path, ps_minus_path, extra_enrich_bed, allow_bed)
save.image("cd4.dnase1.asvm.RData")

