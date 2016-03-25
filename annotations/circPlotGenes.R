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
q("no")

pdf("CCR7.circplot.pdf")

cd.circplot(rpkm_df[ca$name == "chr17_38672950_38722300", ], snU) ## CCR7
cd.circplot(rpkm_df[ca$name == "chr17_38672950_38722300", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -7)) ## CCR7

dev.off()

pdf("MasterRegulators.circplot.pdf")

cd.circplot(rpkm_df[ca$name == "chr17_45810350_45842800", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-6, -14)) ## TBX21
cd.circplot(rpkm_df[ca$name == "chr10_8095400_8179250", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-6, -12)) ## GATA3
cd.circplot(rpkm_df[ca$name == "chrX_49093550_49130150", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) ## FOXP3, note problem with TSS
cd.circplot(rpkm_df[ca$name == "chr1_198588550_198748950", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10))  # PTPRC
cd.circplot(rpkm_df[ca$name == "chr3_187431250_187468300", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) # BCL6
cd.circplot(rpkm_df[ca$name == "chr2_191783100_192016300", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-6, -12)) # STAT4
cd.circplot(rpkm_df[ca$name == "chr12_57413700_57506300", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) # STAT6, problems with TTS.
cd.circplot(rpkm_df[ca$name == "chr17_40453350_40541300", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -10)) # STAT3
cd.circplot(rpkm_df[ca$name == "chr1_151769750_151788550", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-8, -16)) # RORC, very poor match to the annotation.
cd.circplot(rpkm_df[ca$name == "chr11_47363250_47416450", c(1:9,18:19)], snU[c(1:9,18:19)], lims= c(-2, -12)) #  PU1/ SPI1, Wrong TSS

dev.off()

pdf("ETS1_ELF1_NRF1.circplot.pdf")
cd.circplot(rpkm_df[ca$name == "chr11_128221800_128395950",], snU) # ETS1
#cd.circplot(rpkm_df[ca$name == "chr13_41366700_41593550",], snU) # ELF1
cd.circplot(rpkm_df[ca$name == "chr7_129251550_129420450",], snU) # NRF1
dev.off()

cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 
cd.circplot(rpkm_df[ca$name == "", c(1:9,18:19)], snU[c(1:9,18:19)]) # 


## Neat.
cd.circplot(rpkm_df[ca$name == "chr4_123527900_123542600", ], snU) ## IL21
cd.circplot(rpkm_df[ca$name == "chr4_123499950_123506250", ], snU) ## IL21 enhancer +.
cd.circplot(rpkm_df[ca$name == "chr4_123498550_123499600", ], snU) ## IL21 enhancer -.


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
cd.circplot(rpkm_df[ca$name == "ANTXR2", ], snU) ## ANTXR2.  Quite clear.
dev.off()



## Not a clean case on the browser, but looks good from availiable data.
cd.circplot(rpkm_df[ca$name == "chr14_104327650_104338350", ], snU)
cd.circplot(rpkm_df[ca$name == "chr12_8230900_8235450", ], snU)
cd.circplot(rpkm_df[ca$name == "chr6_34393250_34394100", ], snU)

## Not sure...
cd.circplot(rpkm_df[ca$name == "chrY_1523573_1606700", ], snU) ## H and M clearly both outlier species.  This one's the setup.
cd.circplot(rpkm_df[ca$name == "chr6_41301800_41302500", ], snU) ## H and C are pretty similar... H higher, but not by much.
cd.circplot(rpkm_df[ca$name == "chr6_36879950_36880650", ], snU) ## Again, H and C pretty similar ... Again H clearly higher.

## Bad...
cd.circplot(rpkm_df[ca$name == "chr17_38776650_38777400", ], snU) ## Extra chimp might help?!


cd.circplot(rpkm_df[ca$mgi == "TBX21", ], snU)
## Untested

## MYC
cd.circplot(rpkm_df[ca$name == "chr8_128746500_128768850",], snU)
cd.circplot(rpkm_df[ca$name == "chr8_128806600_129211750",], snU)

## CD4
cd.circplot(rpkm_df[ca$name == "chr12_6895750_6959750",], snU)

## IL2RA
cd.circplot(rpkm_df[ca$name == "chr10_6041300_6104350", ][1,], snU)

## IL2RA-enhancer
cd.circplot(rpkm_df[ca$name == "chr10_6079150_6108800", ], snU) 

## MALT1
cd.circplot(rpkm_df[ca$name == "chr18_56338250_56435300", ], snU) 


cd.circplot(rpkm_df[ca$name == "chr15_67356550_67513750",],snU)

## CXCR4
cd.circplot(rpkm_df[ca$name == "chr2_136862150_136897250"], snU)

## Candidates for Post-Transcriptional Regulation
pdf("/home/cgd24/work/evol_rnaseq/PrimaryTXN.pdf")
cd.circplot(rpkm_df[ca$name == "chr7_148844000_148892549",1:8], snU[1:8])
cd.circplot(rpkm_df[ca$name == "chr1_47799400_47884950",1:8], snU[1:8])
cd.circplot(rpkm_df[ca$name == "chr2_232570750_232585500",1:8], snU[1:8])
dev.off()


