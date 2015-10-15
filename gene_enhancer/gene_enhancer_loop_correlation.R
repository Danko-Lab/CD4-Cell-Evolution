#
# gene_enhancer_loop_correlation.R -- correlates gene expression with enhancer initiation levels.
#


## Returns indices in BED2 that intersect BED1.
source("../lib/getOverlap.R")

##Prepare loops.
dist<- 1000 # 15000
loops <- read.table("/local/storage/data/hg19/cd4/chiapet_h3k4me2/H3K4me2_interact_hg19.bed.gz")
loopdist <- function(i) { ## Get the actual distance in the detected loop interaction.
 loop1 <- sapply(strsplit(as.character(loops$V4), split=";"), function(x) {x[[i]]})
 sapply(strsplit(loop1, split="-"), function(x) {as.double(x[[2]])-as.double(x[[1]])})
}
loops1<- loops; loops1[,2] <- loops[,2]+as.integer(loopdist(1)/2)-dist; loops1[,3] <- loops[,2]+as.integer(loopdist(1)/2)+dist
loops2<- loops; loops2[,2] <- loops[,3]+as.integer(loopdist(2)/2)-dist; loops2[,3] <- loops[,3]+as.integer(loopdist(2)/2)+dist
#loops1<- loops; loops1[,3] <- loops[,2]+loopdist(1) #+dist
#loops2<- loops; loops2[,2] <- loops[,3]-loopdist(2) #-dist
hist(log(loops[,3]-loops[,2], 10), 20) ## Summary stats on loop distances.

############
## Get enhancers interacting with changed promoters...
getEnhInt <- function(prefix="H", column=21, post_pro= ".change-U.tsv", post_enh= ".change-U.all.tsv") {
 ## Get TSS of changed, annotated protein coding genes.
 HC <- read.table(paste("../annotations/chage_expr/",prefix,post_pro, sep=""))
 genes <- HC[HC$V7 == "protein_coding",]
 tss   <- genes
 tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]-250; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
 tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]+251; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1
 
 ## dREG sites.
 tres <- read.table(paste("../annotations/chage_expr/",prefix, post_enh, sep=""))
 tres <- tres[grep("dREG", tres$V10),]
 
 ## Find out which loops intersect.
 enh_pro_change <- NULL
 for(i in c(1:NROW(tss))) {
  #print(i)
  indx1 <- getOverlap(tss[i,], loops1)
  indx2 <- getOverlap(tss[i,], loops2)
  indxTRE <- integer(0)
  if(NROW(indx1)>0) {
 	indxTRE <- c(indxTRE, getOverlap(loops2[indx1,], tres))
  }
  if(NROW(indx2)>0) {
 	indxTRE <- c(indxTRE, getOverlap(loops1[indx2,], tres))
  }
  if(NROW(indxTRE)>0) {
  	enh_pro_change <- rbind(enh_pro_change, data.frame(enh=sum(tres[indxTRE, column]), pro=tss[i,column]))
  }
 }
 return(enh_pro_change)
}

################################
## Actually do correlation.

enh_pro_change <- rbind(getEnhInt("H", 21, post_enh= ".change-U.tsv"), 
			getEnhInt("C", 22, post_enh= ".change-U.tsv"), 
			getEnhInt("M", 23, post_enh= ".change-U.tsv")) #Conditional on both promoter and enhancer change.
cor.test(enh_pro_change$pro, enh_pro_change$enh)
plot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19); abline(0,1)

enh_pro_change <- rbind(getEnhInt("H", 21), 
                        getEnhInt("C", 22), 
                        getEnhInt("M", 23))
cor.test(enh_pro_change$pro, enh_pro_change$enh)
plot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)
## Add violin plot.
require(vioplot)
vioplot(enh_pro_change$enh[enh_pro_change$pro < 0], enh_pro_change$enh[enh_pro_change$pro > 0]); abline(h=0)

pdf("gene_enhancer_correlation.pdf")
 plot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)

 source("../lib/densScatterplot.R")
 densScatterplot(as$V5, as$V7, xlab= "Gene Expression", ylab="Sum Enhancers")
dev.off()

####################################################################
## Find all genes that make multiple loops.
 ## dREG sites.
 tres <- read.table(paste("../annotations/chage_expr/H.change-U.all.tsv", sep=""))
 tres <- tres[grep("dREG", tres$V10),]

 ## Get TSS of changed, annotated protein coding genes.
 HC <- read.table(paste("../annotations/chage_expr/",prefix,post_pro, sep=""))
 genes <- HC[HC$V7 == "protein_coding",]
 tss   <- genes
 tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]-250; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
 tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]+251; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1

 

####################################################################
## Look at sig. enhancer changes.  What happens to target genes?

getGeneChanges <- function(prefix="H", column=21, post_pro= ".change-U.all.tsv", post_enh= ".change-U.tsv") {
 ## dREG sites.
 tres <- read.table(paste("../annotations/chage_expr/",prefix, post_enh, sep=""))
 tres <- tres[grep("dREG", tres$V10),]

 ## Get TSS of changed, annotated protein coding genes.
 HC <- read.table(paste("../annotations/chage_expr/",prefix,post_pro, sep=""))
 genes <- HC[HC$V7 == "protein_coding",]
 tss   <- genes
 tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]-250; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
 tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]+251; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1

 ## Find out which loops intersect.
 enh_pro_change <- NULL
 for(i in c(1:NROW(tres))) {
  #print(i)
  indx1 <- getOverlap(tres[i,], loops1)
  indx2 <- getOverlap(tres[i,], loops2)
  indxGENE <- integer(0)
  if(NROW(indx1)>0) {
        indxGENE <- c(indxGENE, getOverlap(loops2[indx1,], tss))
  }
  if(NROW(indx2)>0) {
        indxGENE <- c(indxGENE, getOverlap(loops1[indx2,], tss))
  }
  if(NROW(indxGENE)>0) {
        enh_pro_change <- rbind(enh_pro_change, data.frame(pro=sum(tss[indxGENE, column]), enh=tres[i,column]))
  }
 }
 return(enh_pro_change)
}

enh_pro_change <- rbind(getGeneChanges("H", 21),
                        getGeneChanges("C", 22),
                        getGeneChanges("M", 23))
cor.test(enh_pro_change$pro, enh_pro_change$enh)
plot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Sum Enhancers", pch=19)

## What fraction of the time are target promoters unchanged.
sum(abs(enh_pro_change$pro) < 1)/NROW(enh_pro_change)

##
enh_pro_change_gene <- rbind(getGeneChanges("H", 21, post_pro= ".change-U.tsv"),
                        getGeneChanges("C", 22, post_pro= ".change-U.tsv"),
                        getGeneChanges("M", 23, post_pro= ".change-U.tsv"))
cor.test(enh_pro_change_gene$pro, enh_pro_change_gene$enh)
plot(enh_pro_change_gene$pro, enh_pro_change_gene$enh, xlab= "Sum Gene Expression", ylab="Enhancers", pch=19)


