#
# enhancer_conservation_at_conserved_genes.R -- 
# If gene expression is conserved, are enhancers that make one loop interaction more conserved than those which make several?
#

require(Hmisc)
require(boot)

lnTH= 0.05
hcTH= 0.30

## Returns indices in BED2 that intersect BED1.
source("../lib/getOverlap.R")

##Prepare loops.
dist<- 0 # 1000
loops <- read.table("/local/storage/data/hg19/cd4/chiapet_h3k4me2/H3K4me2_interact_hg19.bed.gz")
loopdist <- function(i) { ## Get the actual distance in the detected loop interaction.
 loop1 <- sapply(strsplit(as.character(loops$V4), split=";"), function(x) {x[[i]]})
 sapply(strsplit(loop1, split="-"), function(x) {as.double(x[[2]])-as.double(x[[1]])})
}
#loops1<- loops; loops1[,2] <- loops[,2]+as.integer(loopdist(1)/2)-dist; loops1[,3] <- loops[,2]+as.integer(loopdist(1)/2)+dist
#loops2<- loops; loops2[,2] <- loops[,3]+as.integer(loopdist(2)/2)-dist; loops2[,3] <- loops[,3]+as.integer(loopdist(2)/2)+dist
loops1<- loops; loops1[,3] <- loops[,2]+loopdist(1) +dist
loops2<- loops; loops2[,2] <- loops[,3]-loopdist(2) -dist
hist(log(loops[,3]-loops[,2], 10), 20) ## Summary stats on loop distances.

## Load fdr.RData
load("../annotations/fdr.RData")

############
## Get enhancers interacting with conserved promoters...
indxTREs <- list(); for(i in 1:50) {indxTREs[[i]] <- integer()}

 ## Get TSS of protein coding genes.
 indx <- ca$type == "protein_coding" #fdr_df$fdr_min > 0.1 & ca$type == "protein_coding" & abs(fdr_df$fc_min) < 1 ## To focus on conserved ... or not?
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
 one2m <- sum(re$V20 > 0) ## Possible 1:many orthology
 unmap <- sum(re$V20 == 0 & is.na(re$mapSize)) ## INDEL
 lccng <- sum((re$V7 < lnTH | re$V8 < lnTH | re$V9 < lnTH) & (re$V7 > hcTH | re$V8 > hcTH | re$V9 > hcTH) & re$V20 == 0 & !is.na(re$mapSize))  # 'Low-confidence'
 chang <- sum(re$V20 == 0 & !is.na(re$mapSize) & re$fdr_min < PVAL & (re$V7 > hcTH | re$V8 > hcTH | re$V9 > hcTH) & (re$V7 > lnTH & re$V8 > lnTH & re$V9 > lnTH)) # 'High-confidence'
 allcng<- sum(re$V20 == 0 & !is.na(re$mapSize) & re$fdr_min < PVAL)

 tot <- total-one2m
 con <- total-one2m-unmap-chang-lccng

# print(paste(tot, con))

 return(con/tot)
}

nloops <- 1:8
dc <- sapply(nloops, function(i) {doesChange(tres[unique(sort(indxTREs[[i]])),])})

data.frame(nloops= nloops, conservation=dc)
plot(nloops, dc, xlab="Number of loops", ylab="Conservation")

#########
## Scatterplot, with points sized by number of points.

getCex <- function(n) { y=0.035*n+0.05; y[y>3] <- 3; y[y<0.05] <- 0.05; y }
n <- sapply(nloops, function(x) {NROW(unique(sort(indxTREs[[x]])))})

size_key <- seq(0, 100,by= 10)

########
## Get the best fit line.

fit_line <- lm(dc~nloops, w = n/sum(n))
fit_line

pdf("DistalEnhancer.ChangeOverDistance.pdf")
 plot(nloops, dc, type="p", xlab="Number of Loops (proximal end)", ylab="% Conserved", pch=19, cex=3*getCex(n), xlim=c(-1,8), ylim=c(0.4, 0.8))
 abline(fit_line)
 plot(size_key, rep(1, NROW(size_key)), cex=3*getCex(size_key))
dev.off()

## Correlations...
require(boot)

cor.test(nloops, dc)
corr(as.matrix(cbind(nloops, dc)), w = n/sum(n))

## Bootstrap to test significance of negative correlation.
# 1. combine all data points.
all_data_points <- unlist(indxTREs)
n_non_unique    <- sapply(nloops, function(i) {NROW(unique(indxTREs[[i]]))})

boot_data <- sapply(1:1000, function(qwi) {
 # 2. split back up into hte same sizes at random.
 indxTREs_boot <- list()
 for(i in nloops) {
        indxTREs_boot[[i]] <- sample(all_data_points, n_non_unique[i], replace=FALSE)
 }

 # 3. Get the corr.
 dc_boot <- sapply(nloops, function(i) {doesChange(tres[unique(indxTREs_boot[[i]]),])})
 nelem_boot <- sapply(nloops, function(i) {NROW(unique(sort(indxTREs_boot[[i]])))})
 corr(data.frame(nloops, dc_boot), w=nelem_boot/sum(nelem_boot))
})
hist(boot_data, 100)
real_cor <- corr(as.matrix(cbind(nloops, dc)), w = n/sum(n))
abline(v=real_cor, lty="dotted", lwd=2, col="red")
sum(boot_data <= -1*abs(real_cor) | boot_data >= abs(real_cor))/ NROW(boot_data) ## Actual statistical test.


