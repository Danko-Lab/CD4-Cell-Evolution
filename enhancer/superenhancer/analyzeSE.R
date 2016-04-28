## Analyze enhancers/ promtoers based on cluster density.

load("../../annotations/fdr.RData")

require(boot)

clusterdensity <- read.table("dREG.cluster.tsv.gz")

tss_aln <- fdr_df[grepl("dREG", ca$annot_type),]
tss <- read.table("tss.nNearby.tsv")
tss <- data.frame(tss, tss_aln[match(tss$V4, tss_aln$name),c(9,33:50)])

## Classify as 'promoter'/ 'enhancer'
stab <- rowMax(tss[,17:18])
dist <- tss[,13]
class <- rep("tss", NROW(tss)) ## tss is then unclassified as a promoter or enhancer
class[stab < 0.1 & dist < 500]  <- "Prox_Stab" ## Clearly protein coding promoter; 50
class[stab > 0.1  & dist > 10000] <- "Dist_UnSt" ## Clearly distal enhancer; 1000
class[stab < 0.1  & dist > 125000] <- "Dist_Stab" ## Clearly stable, but distal
summary(as.factor(class))
tss$V5 <- as.factor(class)

for(i in 7:12) { tss[is.na(tss[,i]),i] <- 0 }


## Compare conservation of looped/ distal/ etc. enhancers
fracConserved <- function(tss, ret="C") { ## defaults to ... 
 total <- NROW(tss)
 one2m <- sum(tss$V20 > 0) ## Possible 1:many orthology
 unmap <- sum(tss$V20 == 0 & is.na(tss$mapSize)) ## INDEL
 lccng <- sum((tss$V7 < 0.1 | tss$V8 < 0.1 | tss$V9 < 0.1) & (tss$V7 > 0.7 | tss$V8 > 0.7 | tss$V9 > 0.7) & tss$V20 == 0 & !is.na(tss$mapSize)) # 'Low-confidence'
 chang <- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05 & (tss$V7 > 0.7 | tss$V8 > 0.7 | tss$V9 > 0.7) & (tss$V7 > 0.1 & tss$V8 > 0.1 & tss$V9 > 0.1)) # 'High-confidence'
 allcng<- sum(tss$V20 == 0 & !is.na(tss$mapSize) & tss$fdr_min < 0.05)

 if(ret=="C") {
  tot <- total-one2m
  con <- total-one2m-unmap-chang-lccng

  return(con/tot)
 }
 if(ret=="GL") {
  tot <- total-one2m
  lcc <- lccng

  return(lcc/tot)
 }
 if(ret=="CNG") {
  tot <- total-one2m
  cng <- chang #allcng ## A santiy check to make sure this shows the same trend as lccng.

  return(cng/tot)
 }
 if(ret=="UNM") {
  tot <- total-one2m
  unm <- unmap #allcng ## A santiy check to make sure this shows the same trend as lccng.

  return(unm/tot)
 }

}

ndreg= 1:20
prom_cluster_cons <- sapply(ndreg, function(i) {fracConserved(tss[tss$V21 == i & tss$V5=="Prox_Stab",])})  ## Conservation of unlooped. 
prom_n_s <- sapply(ndreg, function(i) { sum(tss$V21 == i & tss$V5=="Prox_Stab") })

enh_cluster_cons <- sapply(ndreg, function(i) {fracConserved(tss[tss$V21 == i & tss$V5=="Dist_UnSt",])})  ## Conservation of unlooped.
enh_n_s <- sapply(ndreg, function(i) { sum(tss$V21 == i & tss$V5=="Dist_UnSt") })

data.frame(prom_cluster_cons, prom_n_s, enh_cluster_cons, enh_n_s)

plot(prom_cluster_cons)
plot(enh_cluster_cons)



