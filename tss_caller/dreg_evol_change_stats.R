##
## Basic evolutinoary change stats to add to the text... 

hcm <- read.table("HCM-U-PI.dREG-tss-clusters.tsv")

## Overall conservation stats.
th <- 0.8; th_l <- 0.4

sum(hcm[,7]>th & hcm[,8]>th & hcm[,9]>th)/ sum(hcm[,7]>th | hcm[,8]>th | hcm[,9]>th)
sum((hcm[,7]>th_l & hcm[,8]>th_l & hcm[,9]>th_l) & (hcm[,7]>th | hcm[,8]>th | hcm[,9]>th))/ sum(hcm[,7]>th | hcm[,8]>th | hcm[,9]>th)

## Compare conservation between species.

