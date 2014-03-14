## train_svm.R -- Trains an DBN to detect regulatory regions based on DNAse-1 data.
##

require(featureDetector)

## Read PRO-seq data.
ps_plus_path  <- list("H-U_plus.bw", "H-PI_plus.bw", "C-U_plus.bw", "C-PI_plus.bw", "M-U_plus.bw", "M-PI_plus.bw")
ps_minus_path <- list("H-U_minus.bw", "H-PI_minus.bw", "C-U_minus.bw", "C-PI_minus.bw", "M-U_minus.bw", "M-PI_minus.bw")

## Get positive regions.
dnase <- lapply(rep("dnase1.peaks_peaks.narrowPeak",6), function(x) { read.table(x)})
extra_enrich_bed <- lapply(rep("GencodeMerge.IntersectOpStrand.bed",6), function(x) {read.table(x)})
allow_bed <- lapply(rep("CD4.chromHMM.Ernst2010.hg19.Prom.Enh.bed",6), function(x) {read.table(x)})

## Train the SVM.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE)
print(paste("Number of inf. positions: ", NROW(inf_positions)))

## Create the genomc data model and the dbn.
gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))
adbn <- dbn(layer_sizes= c(360,300,300,500), batch_size=100, cd_n=1, momentum_decay= 0.9, weight_cost= 1e-5, learning_rate=0.1)

## Learn features basd on all species.
adbn <- regulatory_dbn(gdm, adbn, ps_plus_path, ps_minus_path, inf_positions, dnase, n_train=30000, n_eval=0, extra_enrich_bed= extra_enrich_bed, allow= allow_bed, training_mode="pretrain")

## Refine based on H-U through backprop.
adbn <- regulatory_dbn(gdm, adbn, ps_plus_path[[1]], ps_minus_path[[1]], inf_positions[[1]], dnase[[1]], n_train=150000, n_eval=0, extra_enrich_bed= extra_enrich_bed[[1]], allow= allow_bed[[1]], training_mode="refine")

remove(inf_positions, ps_plus_path, ps_minus_path, extra_enrich_bed, allow_bed)
save.image("cd4.dnase1.adbn.RData")

