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
indx1TRE <- integer(0)
indxGT1TRE<-integer(0)
loopsToEnh<-data.frame(nloops=integer(), nenh=integer())

 ## Get TSS of conserved, annotated protein coding genes.
 indx <- fdr_df$fdr_min > 0.15 & ca$type == "protein_coding" & abs(fdr_df$fc_min) < 0.75
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
  ## Compare 1 loop ...
  enhLoops <- integer(0)
  if(nloops == 1) {
    if(NROW(indx1)>0) {
      indx1TRE <- c(indx1TRE, getOverlap(loops2[indx1,], tres))
      enhLoops <- c(enhLoops, getOverlap(loops2[indx1,], tres))
    }
    if(NROW(indx2)>0) {
      indx1TRE <- c(indx1TRE, getOverlap(loops1[indx2,], tres))
      enhLoops <- c(enhLoops, getOverlap(loops1[indx2,], tres))
    }
  } 

  ## ... to more than 2 loops.
  if(nloops > 2) {
    if(NROW(indx1)>0) {
      indxGT1TRE <- c(indxGT1TRE, getOverlap(loops2[indx1,], tres))
      enhLoops <- c(enhLoops, getOverlap(loops2[indx1,], tres))
    }
    if(NROW(indx2)>0) {
      indxGT1TRE <- c(indxGT1TRE, getOverlap(loops1[indx2,], tres))
      enhLoops <- c(enhLoops, getOverlap(loops1[indx2,], tres))
    }
  }

 loopsToEnh <- rbind(loopsToEnh, data.frame(nloops, NROW(enhLoops)))
 # if(nloops > 0) print(paste(nloops, NROW(enhLoops)))

 }

indx1TRE <- unique(sort(indx1TRE)); NROW(indx1TRE)
indxGT1TRE <- unique(sort(indxGT1TRE)); NROW(indxGT1TRE)

plot(loopsToEnh)

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

doesChange(tres[indx1TRE,])
doesChange(tres[indxGT1TRE,])

## Write out location files for DNA sequence analysis.
write.table(tres[indx1TRE,1:4],   "loopsToOneTRE.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(tres[indxGT1TRE,1:4], "loopsToMoreTRE.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


