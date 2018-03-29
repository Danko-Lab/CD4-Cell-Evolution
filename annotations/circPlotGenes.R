##
## Sanity checks for gene expression changes ... Look at the reported raw levels of gene expression.

load("fdr.RData")

source("readData.R")

rpkm_df <- as.matrix(ca[,c(2:6,7:10,12:17,19:20,24:25)+10]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- log(1000*(rpkm_df[,i]+0.25)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"]), 2)

source("../lib/circplot.R")
snU <- c(rep("H-U",3), rep("C-U", 3), rep("M-U",3), rep("H-PI", 3), rep("C-PI", 3), rep("M-PI", 2), rep("Rodent-U",2))

## Write out all human
#hc <- fdr_df[fdr_df$HumanFDR < 0.05,][order(fdr_df$HumanFDR[fdr_df$HumanFDR < 0.05]),]
#hcn<- as.character$name(fdr_df[fdr_df$HumanFDR < 0.05][order(fdr_df$HumanFDR[fdr_df$HumanFDR < 0.05])])
#pdf("circPlot.hc.pdf")
#for(n in hcn) {
#  cd.circplot(rpkm_df[ca$name == n, ], snU)
#}
#dev.off()

## Specific examples...

pdf("CCR7.circplot.pdf")

cd.circplot(rpkm_df[ca$name == "chr17_38679750_38722150", ], snU) ## CCR7
cd.circplot(rpkm_df[ca$name == "chr17_38679750_38722150", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -7)) ## CCR7

dev.off()

pdf("SGPL1.circplot.pdf")

cd.circplot(rpkm_df[ca$name == "chr10_72575600_72671450", ], snU) ## CCR7

dev.off()

pdf("MasterRegulators.circplot.pdf")

cd.circplot(rpkm_df[ca$name == "chr17_45810350_45839600", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-6, -14)) ## TBX21
cd.circplot(rpkm_df[ca$name == "chr10_8091100_8177100", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-6, -12)) ## GATA3
cd.circplot(rpkm_df[ca$name == "chrX_49093550_49121250", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-4, -12)) ## FOXP3
cd.circplot(rpkm_df[ca$name == "chr1_198588600_198749550", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10))  # PTPRC
cd.circplot(rpkm_df[ca$name == "chr3_187431800_187468650", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) # BCL6
cd.circplot(rpkm_df[ca$name == "chr2_191766950_192016300", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-6, -12)) # STAT4
cd.circplot(rpkm_df[ca$name == "chr12_57472950_57506350", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) # STAT6
cd.circplot(rpkm_df[ca$name == "chr17_40453500_40541350", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) # STAT3
cd.circplot(rpkm_df[ca$name == "chr1_151770200_151812800", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-8, -18)) # RORC
cd.circplot(rpkm_df[ca$name == "chr11_47344050_47416950", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -14)) #  PU1/ SPI1, Wrong TSS

dev.off()

pdf("ETS1_ELF1_NRF1.circplot.pdf")
cd.circplot(rpkm_df[ca$name == "chr11_128261550_128395950",], snU) # ETS1
#cd.circplot(rpkm_df[ca$name == "chr13_41366600_41593550",], snU) # ELF1
cd.circplot(rpkm_df[ca$name == "chr7_129251550_129420450",], snU) # NRF1
dev.off()

#cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
#cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
#cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
#cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
#cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 


## Neat.
cd.circplot(rpkm_df[ca$name == "chr4_123528550_123542300", ], snU) ## IL21


## Note chimp is lower in FOXP3 (Treg) and higher in SPI1 (TH9) than other species.
## Possibly swap of similar cells? http://www.nature.com/icb/journal/v88/n6/fig_tab/icb201073f1.html
##
## Note, however, that there's virtually no IL9 production.

## Good ...
cd.circplot(rpkm_df[ca$name == "chr12_15771150_15943000", ], snU)
cd.circplot(rpkm_df[ca$name == "chr13_77989450_77990050", ], snU)
cd.circplot(rpkm_df[ca$name == "chr12_6308800_6352600", ], snU) ## CD9.  Quite clear.
cd.circplot(rpkm_df[ca$name == "chr12_2173350_2176350", ], snU)
cd.circplot(rpkm_df[ca$name == "chr5_140005300_140013300", ], snU) ## CD14.  Might be biased, but looks right on all accounts.

pdf("ANTXR2.pdf")
cd.circplot(rpkm_df[ca$mgi == "ANTXR2", ], snU) ## ANTXR2.  Quite clear.
dev.off()

cd.circplot(rpkm_df[ca$mgi == "TBX21", ], snU)
## Untested

## CD4
cd.circplot(rpkm_df[ca$mgi == "CD4",], snU)

## IL2RA
cd.circplot(rpkm_df[ca$mgi == "IL2RA", ][1,], snU)

## MALT1
cd.circplot(rpkm_df[ca$mgi == "MALT1", ], snU) 

## Candidates for Post-Transcriptional Regulation
#pdf("/home/cgd24/work/evol_rnaseq/PrimaryTXN.pdf")
#cd.circplot(rpkm_df[ca$name == "chr7_148844000_148892549",1:8], snU[1:8])
#cd.circplot(rpkm_df[ca$name == "chr1_47799400_47884950",1:8], snU[1:8])
#cd.circplot(rpkm_df[ca$name == "chr2_232570750_232585500",1:8], snU[1:8])
#dev.off()


