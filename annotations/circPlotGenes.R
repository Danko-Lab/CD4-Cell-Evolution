##
## Sanity checks for gene expression changes ... Look at the reported raw levels of gene expression.

source("readData.R")

rpkm_df <- as.matrix(ca[,indx.good[c(2:9,11:17)]]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- log(1000*(rpkm_df[,i]+0.25)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"]))

source("../lib/circplot.R")
snU <- c(rep("H",3), rep("C", 2), rep("M",3))

## Good ...
cd.circplot(rpkm_df[ca$name == "chr12_15771150_15943000", 1:8], snU)
cd.circplot(rpkm_df[ca$name == "chr13_77989450_77990050", 1:8], snU)
cd.circplot(rpkm_df[ca$name == "chr12_6308800_6352600", 1:8], snU) ## CD9.  Quite clear.
cd.circplot(rpkm_df[ca$name == "chr12_2173350_2176350", 1:8], snU)
cd.circplot(rpkm_df[ca$name == "chr5_140005300_140013300", 1:8], snU) ## CD14.  Might be biased, but looks right on all accounts.


## Not a clean case on the browser, but looks good from availiable data.
cd.circplot(rpkm_df[ca$name == "chr14_104327650_104338350", 1:8], snU)
cd.circplot(rpkm_df[ca$name == "chr12_8230900_8235450", 1:8], snU)
cd.circplot(rpkm_df[ca$name == "chr6_34393250_34394100", 1:8], snU)

## Not sure...
cd.circplot(rpkm_df[ca$name == "chrY_1523573_1606700", 1:8], snU) ## H and M clearly both outlier species.  This one's the setup.
cd.circplot(rpkm_df[ca$name == "chr6_41301800_41302500", 1:8], snU) ## H and C are pretty similar... H higher, but not by much.
cd.circplot(rpkm_df[ca$name == "chr6_36879950_36880650", 1:8], snU) ## Again, H and C pretty similar ... Again H clearly higher.

## Bad...
cd.circplot(rpkm_df[ca$name == "chr17_38776650_38777400", 1:8], snU) ## Extra chimp might help?!

## Untested

