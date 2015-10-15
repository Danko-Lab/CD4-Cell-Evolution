#
# gene_enhancer_correlation.R -- correlates gene expression with enhancer initiation levels.
#

as <- rbind(read.table("tmp/H.change-U.tsspc.enh"), read.table("tmp/C.change-U.tsspc.enh"), read.table("tmp/M.change-U.tsspc.enh"))

#as <- as[abs(as$V7)>1,]

print(cor.test(as$V5, as$V7))
print(cor.test(as$V5, as$V7, method="spearman"))

plot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)

pdf("gene_enhancer_correlation.pdf")
 plot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)

 source("../lib/densScatterplot.R")
 densScatterplot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers")
dev.off()

####################################################################################################
## Not sure why promoters are not being picked up in some cases using bedmap/ bedops?!  Repeat in R.
source("../lib/getOverlap.R")
source("../lib/densScatterplot.R")

getEnhNear <- function(prefix="H", column=21, dist=50000, post_pro= ".change-U.tsv", post_enh= ".change-U.all.tsv") {
 ## Get TSS of changed, annotated protein coding genes.
 HC <- read.table(paste("../annotations/chage_expr/",prefix,post_pro, sep=""))
 genes <- HC[HC$V7 == "protein_coding",]
 tss   <- genes
 tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]-250; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
 tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]+251; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1

 ## dREG sites.
 tres <- read.table(paste("../annotations/chage_expr/",prefix, post_enh, sep=""))
 tres <- tres[grep("dREG", tres$V10),]

 ## nearby ... 
 nearby <- tss
 nearby[,2] <- tss[,2]-dist
 nearby[,3] <- tss[,3]+dist

 ## Find out which loops intersect.
 enh_pro_change <- NULL
 for(i in c(1:NROW(tss))) {

  ## Get all nearby REs, excluding overlap with the TSS.
  indxtss <- getOverlap(tss[i,], tres)
  indx <- getOverlap(nearby[i,], tres)
  truth_ <- rep(FALSE, NROW(tres)); truth_[indx] <- TRUE; truth_[indxtss] <- FALSE
  indx <- which(truth_)

  ## If REs are left over, include them.
  if(NROW(indx)>0) {
        enh_pro_change <- rbind(enh_pro_change, data.frame(enh=mean(tres[indx, column]), pro=tss[i,column]))
  }
 }
 return(enh_pro_change)
}

enh_pro_change <- rbind(getEnhNear("H", 21),
                        getEnhNear("C", 22),
                        getEnhNear("M", 23))
cor.test(enh_pro_change$pro, enh_pro_change$enh)
plot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)
densScatterplot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Sum Enhancers")



