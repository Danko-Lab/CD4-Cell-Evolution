#
# enhancer_conservation_at_conserved_genes.R -- 
# If gene expression is conserved, are enhancers that make one loop interaction more conserved than those which make several?
#

## Returns indices in BED2 that intersect BED1.
source("../lib/getOverlap.R")

##Prepare loops.
dist<- 2500 # 1000
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

## Load fdr.RData
load("../annotations/fdr.RData")

############
## Get enhancers interacting with conserved promoters...
indxTREs <- list(); for(i in 1:50) {indxTREs[[i]] <- integer()}

 ## Get TSS of conserved, annotated protein coding genes.
 indx <- fdr_df$fdr_min > 0.1 & ca$type == "protein_coding" & abs(fdr_df$fc_min) < 1
 tss   <- ca[indx,1:8]
 tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]-250; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
 tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]+251; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1
 
 ## dREG sites.
 tre_aln <- fdr_df[grepl("dREG", ca$annot_type),]
 tre <- read.table("tss.tsv")
 tres <- data.frame(tre, tre_aln[match(tre$V4, tre_aln$name),c(9,33:50,7:8,10)])
 tres$V7[is.na(tres$V7)] = 0; tres$V8[is.na(tres$V8)] = 0; tres$V9[is.na(tres$V9)] = 0

 ## Find out which loops intersect.
 for(i in c(1:NROW(tss))) {
  if((i %% 100) == 0) print(i)

  ## If it finds a TSS ...
  indx1 <- getOverlap(tss[i,], loops1)
  indx2 <- getOverlap(tss[i,], loops2)
  nloops <- NROW(indx1) + NROW(indx2)

  ## ... then compare with loops on the other end.
  ## Compare loops ...
  if(NROW(indx1)>0) {
    indxTREs[[nloops]] <- c(indxTREs[[nloops]], getOverlap(loops2[indx1,], tres))
  }
  if(NROW(indx2)>0) {
    indxTREs[[nloops]] <- c(indxTREs[[nloops]], getOverlap(loops1[indx2,], tres))
  }

 }

save.image("enhancer_conservation_at_conserved_genes.RData")

#########
## Do enhancers that loop to exactly 1 TRE change less frequently?

doesChange <- function(re) {

 total <- NROW(re) 
 one2m <- sum(re$V20 > 0 & !is.na(re$mapSize)) ## Possible 1:many orthology
 unmap <- sum(re$V20 == 0 & is.na(re$mapSize)) ## INDEL
 lccng <- sum((re$V7 < 0.1 | re$V8 < 0.1 | re$V9 < 0.1) & (re$V7 > 0.7 | re$V8 > 0.7 | re$V9 > 0.7) & re$V20 == 0 & !is.na(re$mapSize)) # 'Low-confidence'
 chang <- sum(re$V20 == 0 & !is.na(re$mapSize) & re$fdr_min < 0.05 & (re$V7 > 0.7 | re$V8 > 0.7 | re$V9 > 0.7) & (re$V7 > 0.1 & re$V8 > 0.1 & re$V9 > 0.1)) # 'High-confidence'
 allcng<- sum(re$V20 == 0 & !is.na(re$mapSize) & re$fdr_min < 0.05)

 tot <- total-one2m
 con <- total-one2m-unmap-chang-lccng

 print(paste(tot, con))

 return(con/tot)
}

dc <- sapply(1:10, function(i) {doesChange(tres[unique(sort(indxTREs[[i]])),])})

data.frame(nloops= 1:10, conservation=dc)
plot(1:10, dc, xlab="Number of loops", ylab="Conservation")





