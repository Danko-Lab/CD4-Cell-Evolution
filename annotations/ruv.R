##
## Apply the RUV methods of Gagnon-Bartsch & Speed 
##
## REF: Gagnon-Bartsch and Speed, Biostatistics (2012), 13, 3, pp. 539-552.

source("readData.R")
require(ruv.DE)

Y <- as.matrix(ca[,indx.good[c(2:9,11:17)]])
#X <- as.integer(as.factor(rbind(c("H","H","H","C","C","M","M","M",    "H","H","H","C","C","M","M"), c(rep("U",8), rep("PI",7)))))

X <- rbind(c(1,1,1,2,2,3,3,3,    1,1,1,2,2,3,3), c(rep(4,8), rep(5,7)))
ctrl <- ca$type == "protein_coding" ## Assume no changes in protein-coding...

norm_vals <- RUVinv(Y, X, ctrl)




