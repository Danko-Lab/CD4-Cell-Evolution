## Compare DNA sequence conservation in different classes of RE.
require(bigWig)
require(vioplot)
require(Hmisc)
require(stats)

source("../../lib/getOverlap.R")
source("../../lib/densScatterplot.R")

## Read in dREG-HD peaks.
peaksize <- 5 # 100 ## control peak size.

bed <- read.table("dREG-HD.2compare.bed"); #bed <- bed[bed[,1]=="chr1",]
bed_data <- bed
bed_data <- center.bed(bed, peaksize, peaksize)

bed_null1 <-  read.table("dREG-HD.NULL-MODEL.bed"); #bed_null1 <- bed_null1[bed_null1[,1] == "chr1",]
bed_null1 <- bed_null1[grep("Un|random|hap", bed_null1$V1, invert=TRUE),]; 
bed_null1 <- center.bed(bed_null1, peaksize, peaksize)

bed_null2 <- read.table("dREG-HD.NULL-MODEL.bed.nogap"); #bed_null2 <- bed_null2[bed_null2[,1] == "chr1",]
bed_null2 <- bed_null2[grep("Un|random|hap", bed_null2$V1, invert=TRUE),]
bed_null2 <- center.bed(bed_null2, peaksize, peaksize)

## Read in MAFs
#con <- load.bigWig("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw")
con <- load.bigWig("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw")
word <- "100way"#"Primate"

## Get mean conservation for each h_bed.
mean_con <- bed.region.bpQuery.bigWig(con, bed_data[,1:3], op="avg")# op="max")
mean_con_null1 <- bed.region.bpQuery.bigWig(con, bed_null1[,1:3], op="avg")# op="max")
mean_con_null2 <- bed.region.bpQuery.bigWig(con, bed_null2[,1:3], op="avg")# op="max")

################################################
## Break it down based on distance, etc.
cor.test(mean_con, bed$V4)
idx <- rowSums(bed_data[,c(5:6)])>0; cor.test(mean_con[idx], bed$V4[idx])

## LOESS smoothing.
x_val <- seq(-1500000, 1500000, 1000)
spn <- 0.1

ls.data <-data.frame(dist= bed_data$V4, con= mean_con)
ls.fit <- loess(con~dist, ls.data, span=spn)
y_pred<- predict(ls.fit, newdata= x_val)

ls.loop.fit <- loess(con~dist, ls.data[idx,], span=spn)
y_loop<- predict(ls.loop.fit, newdata= x_val)

null1.data <-data.frame(dist= bed_null1$V4, null1= mean_con_null1)
null1.fit <- loess(null1~dist, null1.data, span=spn)
y_null1<-predict(null1.fit, newdata=x_val)

null2.data <-data.frame(dist= bed_null2$V4, null2= mean_con_null2)
null2.fit <- loess(null2~dist, null2.data, span=spn)
y_null2<-predict(null2.fit, newdata=x_val)

## Plot NULLs
pdf(paste(word,".PhyloPByDist.pdf", sep=""))
ylims=c(-3,8) #c(-2, 0.75)

plot(bed_data$V4, mean_con, xlim=c(-500000,500000), ylim=ylims, main="Mean Primate PhyloP per TFBS", xlab="Distance from nearest TSS", ylab="Phylo P")
points(x_val, y_pred, type="l", col="dark red", lwd=3)
abline(h=0, lty="dotted", col="gray")

plot(bed_data$V4[idx], mean_con[idx], xlim=c(-500000,500000), ylim=ylims, main="Mean Primate PhyloP per LOOPED TFBS", xlab="Distance from nearest TSS", ylab="Phylo P")
points(x_val, y_pred, type="l", col="dark gray", lwd=3)
points(x_val, y_loop, type="l", col="dark red", lwd=3)
abline(h=0, lty="dotted", col="gray")

plot(bed_null2$V4, mean_con_null2, xlim=c(-500000,500000), ylim=ylims, main="Mean Primate PhyloP, NULL2", xlab="Distance from nearest TSS", ylab="Phylo P")
points(x_val, y_null2, type="l", col="dark red", lwd=3)
abline(h=0, lty="dotted", col="gray")

plot(bed_null1$V4, mean_con_null1, xlim=c(-500000,500000), ylim=ylims, main="Mean Primate PhyloP, NULL1", xlab="Distance from nearest TSS", ylab="Phylo P")
points(x_val, y_null1, type="l", col="dark red", lwd=3)
abline(h=0, lty="dotted", col="gray")

dev.off()

par(mfrow=c(1,3))
plot(bed_data$V4, mean_con, xlim=c(0,1500000), ylim=c(-2,8))
plot(bed_null1$V4, mean_con_null1, xlim=c(0,1500000), ylim=c(-2,8))
plot(bed_null2$V4, mean_con_null2, xlim=c(0,1500000), ylim=c(-2,8))

densScatterplot(bed_data$V4, mean_con, xlim=c(0,1500000), ylim=c(-2,8))
densScatterplot(bed_null1$V4, mean_con_null1, xlim=c(0,1500000), ylim=c(-2,8))
densScatterplot(bed_null2$V4, mean_con_null2, xlim=c(0,1500000), ylim=c(-2,8))

par(mfrow=c(1,1))


plot(bed$V4+1, mean_con, log="x")
points((bed$V4+1)[idx], mean_con[idx], col="red")

cor.test(log(bed$V4+1), mean_con)
cor.test(log(bed$V4+1)[idx], mean_con[idx])

## Look at signal V. length.
cor.test(bed[,3]-bed[,2], mean_con)
plot(bed[,3]-bed[,2], mean_con)

median(mean_con[bed$V4 < 10])
median(mean_con[bed$V4 > 10 & bed$V4 < 1000])
median(mean_con[bed$V4 > 1000 & bed$V4 < 10000])
median(mean_con[bed$V4 > 10000 & bed$V4 < 100000])
median(mean_con[bed$V4 > 100000])

vioplot(mean_con[bed$V4 < 10], mean_con[bed$V4 > 10 & bed$V4 < 10000], mean_con[bed$V4 > 10000 & bed$V4 < 100000], mean_con[bed$V4 > 100000])

## Scatterplot, like for conservation of function.
xaxis <- c(0, seq(3, 6, 0.02))
vect <- as.numeric(cut2(log(bed$V4, 10), cuts=xaxis)); vect[is.na(vect)] <- 1; summary(as.factor(vect))
dist <- sapply(1:max(vect), function(x) {median(mean_con[vect == x])})

getCex <- function(n) { y=0.0138888*n+0.1; y[y>3] <- 3; y[y<0.1] <- 0.1; y }
n <- sapply(1:max(vect), function(x) {NROW(bed[vect == x,])})

## Last point is [>=xaxis]. Each other point is [<point].
idx <- 3:NROW(xaxis) ## As below, index 3 encodes [1000,10^3.02)

pdf("ChangeOverDistance.pdf")
 plot(10^xaxis[idx], dist[idx], type="p", xlab="Distance from TSS [bp]", ylab="PhyloP", xlim=10^c(3, 6), log="x", pch=19, cex=getCex(n[idx]))
dev.off()

########################
## Superenhancer.
vioplot(mean_con[bed[,7]>0], mean_con[bed$V7==0])

summary(mean_con[bed[,7]>0])
summary(mean_con[bed$V7==0])

pdf(paste("DNASequence.",word,".phyloP.Conservation.pdf", sep=""))
 ld<- 5; xlim_s=c(-0.75,3) #xlim_s=c(-0.75, 0.5) #c(-0.4, 0.5) #PRIMATE

 plot(ecdf(mean_con_null1), col="gray", lwd=ld, xlim=xlim_s)
 plot(ecdf(mean_con_null2), col="dark gray", add=TRUE, lwd=ld)
 plot(ecdf(mean_con[bed_data[,7]>0]), col="#00A63E", add=TRUE, lwd=ld)
 plot(ecdf(mean_con[bed_data[,7]==0 & rowSums(bed_data[,c(5:6)])==0]), col="black", add=TRUE, lwd=ld)
 plot(ecdf(mean_con[rowSums(bed_data[,c(5:6)])>0]), col="#b70000", add=TRUE, lwd=ld)

## Plot smaller version for figs.
# npts <- 5000
# plot(ecdf(mean_con[sample(which(bed_data[,7]>0), npts)]), col="#00A63E", xlim=xlim_s, lwd=ld)
# plot(ecdf(mean_con[sample(which(bed_data[,7]==0 & rowSums(bed_data[,c(5:6)])==0), npts)]), col="black", add=TRUE, lwd=ld)
# plot(ecdf(mean_con[sample(which(rowSums(bed_data[,c(5:6)])>0), npts)]), col="#b70000", add=TRUE, lwd=ld)
# plot(ecdf(mean_con_null1[sample(1:NROW(mean_con_null1), npts)]), col="gray", add=TRUE, lwd=ld)
# plot(ecdf(mean_con_null2[sample(1:NROW(mean_con_null2), npts)]), col="dark gray", add=TRUE, lwd=ld)


## Cut at 0.75
 th <- 0.75
 abline(v=th)

 sum(mean_con[bed_data[,7]==0 & rowSums(bed_data[,c(5:6)])==0] > th)/NROW(mean_con[bed_data[,7]==0 & rowSums(bed_data[,c(5:6)])==0])
 sum(mean_con[bed_data[,7]>0] > th)/NROW(mean_con[bed_data[,7]>0])
 sum(mean_con[rowSums(bed_data[,c(5:6)])>0] > th)/NROW(mean_con[rowSums(bed_data[,c(5:6)])>0])
 sum(mean_con_null1 > th)/ NROW(mean_con_null1)
 sum(mean_con_null2 > th)/ NROW(mean_con_null2)

 ## Then add barplots...
 require(boot)
 b_dreg <- boot(data= mean_con[bed_data[,7]==0 & rowSums(bed_data[,c(5:6)])==0], R=1000, statistic= function(a, i) {sum(a[i] > th)/NROW(i) })
 b_se   <- boot(data= mean_con[bed_data[,7]>0], R=1000, statistic= function(a, i) {sum(a[i] > th)/NROW(i) })
 b_loop <- boot(data= mean_con[rowSums(bed_data[,c(5:6)])>0], R=1000, statistic= function(a, i) {sum(a[i] > th)/NROW(i) })
 b_null1   <- boot(data= mean_con_null1, R=1000, statistic= function(a, i) {sum(a[i] > th)/NROW(i) })
 b_null2   <- boot(data= mean_con_null2, R=1000, statistic= function(a, i) {sum(a[i] > th)/NROW(i) })

 source("../../lib/barplot.R")
 bars <- c(b_dreg$t0, b_se$t0, b_loop$t0, b_null1$t0, b_null2$t0)
 errs <- c(sqrt(var(b_dreg$t)), sqrt(var(b_se$t)), sqrt(var(b_loop$t)), sqrt(var(b_null1$t)), sqrt(var(b_null2$t)))
 names<- c("Enhancer, noloop/ SE", "Super Enhancer", "Looped Enhancer", "NULL1", "NULL2")
 cd.barplot(bars, errs, names, fill=TRUE)

## Now seperate by distance. 

 ## 1-10k
 indx <- abs(bed$V4) < 10000
 plot(ecdf(mean_con[bed[,7]>0 & indx]), col="#00A63E", xlim=xlim_s, lwd=ld)
 plot(ecdf(mean_con[bed[,7]==0 & rowSums(bed[,c(5:6)])==0 & indx]), col="black", add=TRUE, lwd=ld)
 plot(ecdf(mean_con[rowSums(bed[,c(5:6)])>0 & indx]), col="#b70000", add=TRUE, lwd=ld)
 plot(ecdf(mean_con_null1[abs(bed_null1$V4) < 10000]), col="gray", add=TRUE, lwd=ld)
 plot(ecdf(mean_con_null2[abs(bed_null2$V4) < 10000]), col="dark gray", add=TRUE, lwd=ld)


 ## 10-100k
 indx <- abs(bed$V4) < 100000 & abs(bed$V4) > 10000
 plot(ecdf(mean_con[bed[,7]>0 & indx]), col="#00A63E", xlim=xlim_s, lwd=ld)
 plot(ecdf(mean_con[bed[,7]==0 & rowSums(bed[,c(5:6)])==0 & indx]), col="black", add=TRUE, lwd=ld)
 plot(ecdf(mean_con[rowSums(bed[,c(5:6)])>0 & indx]), col="#b70000", add=TRUE, lwd=ld)
 plot(ecdf(mean_con_null1[abs(bed_null1$V4) > 10000 & abs(bed_null1$V4) < 100000]), col="gray", add=TRUE, lwd=ld)
 plot(ecdf(mean_con_null2[abs(bed_null2$V4) > 10000 & abs(bed_null2$V4) < 100000]), col="dark gray", add=TRUE, lwd=ld)

 ## 100-1,000k
 indx <- abs(bed$V4) < 1000000 & abs(bed$V4) > 100000
 plot(ecdf(mean_con[bed[,7]==0 & rowSums(bed[,c(5:6)])==0  & indx]), col="black", xlim=xlim_s, lwd=ld)
 plot(ecdf(mean_con[bed[,7]>0 & indx]), col="#00A63E", add=TRUE, lwd=ld)
 plot(ecdf(mean_con[rowSums(bed[,c(5:6)])>0 & indx]), col="#b70000", add=TRUE, lwd=ld)
 plot(ecdf(mean_con_null1[abs(bed_null1$V4) > 100000 & abs(bed_null1$V4) < 1000000]), col="gray", add=TRUE, lwd=ld)
 plot(ecdf(mean_con_null2[abs(bed_null2$V4) > 100000 & abs(bed_null2$V4) < 1000000]), col="dark gray", add=TRUE, lwd=ld)

 ## Difference based on distance...
 ld<- 3; xlim_s=c(-1, 2) #c(-0.4, 0.5) #PRIMATE
 indx <- bed$V4 < 10000
 plot(ecdf(mean_con[bed[,7]==0 & rowSums(bed[,c(5:6)])==0 & indx]), col="black", xlim=xlim_s, lwd=ld)
 indx <- bed$V4 < 100000 & bed$V4 > 10000
 plot(ecdf(mean_con[bed[,7]==0 & rowSums(bed[,c(5:6)])==0 & indx]), col="dark gray", add=TRUE, lwd=ld)
 indx <- bed$V4 < 1000000 & bed$V4 > 100000
 plot(ecdf(mean_con[bed[,7]==0 & rowSums(bed[,c(5:6)])==0  & indx]), col="gray", lwd=ld, add=TRUE)

dev.off()

########################
## Add loop information.
vioplot(mean_con[rowSums(bed[,c(5:6)])>0], mean_con[rowSums(bed[,c(5:6)])==0])

summary(mean_con[rowSums(bed[,c(5:6)])>0])
summary(mean_con[rowSums(bed[,c(5:6)])==0])

plot(ecdf(mean_con[rowSums(bed[,c(5:6)])>0]), col="red", xlim=xlim_s)
plot(ecdf(mean_con[rowSums(bed[,c(5:6)])==0]), col="black", add=TRUE)

wt<- wilcox.test(mean_con[rowSums(bed[,c(5:6)])>0], mean_con[rowSums(bed[,c(5:6)])==0])
wt; wt$p.value

ks<- ks.test(mean_con[rowSums(bed[,c(5:6)])>0], mean_con[rowSums(bed[,c(5:6)])==0])
ks; ks$p.value

wt<- wilcox.test(mean_con[bed[,7]>0], mean_con[bed[,7]==0])
wt; wt$p.value

ks<- ks.test(mean_con[bed[,7]>0], mean_con[bed[,7]==0])
ks; ks$p.value



##########################################################################################################

######################################################
## Specifically get loops with gene at the other end.
loops <- read.table("/local/storage/data/hg19/cd4/chiapet_h3k4me2/H3K4me2_interact_hg19.bed.gz")
loopdist <- function(i) { ## Get the actual distance in the detected loop interaction.
 loop1 <- sapply(strsplit(as.character(loops$V4), split=";"), function(x) {x[[i]]})
 sapply(strsplit(loop1, split="-"), function(x) {as.double(x[[2]])-as.double(x[[1]])})
}
loops1<- loops; loops1[,3] <- loops[,2]+loopdist(1) #+dist
loops2<- loops; loops2[,2] <- loops[,3]-loopdist(2) #-dist

tss   <- read.table("../../tss_caller/refGene.hg19.bed.gz")
tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1

## Loop through the loops.  Each time there's a match at one or more TFBS, increment the numbers.
nloops <- rep(0, NROW(bed))
ngenesl <- rep(0, NROW(bed))
ngenes  <- rep(0, NROW(bed))

for(i in c(1:NROW(loops))) {
 ## Check for the tfbs | genes on either end of each loop.
 indx1 <- getOverlap(loops1[i,], bed)
 indx2 <- getOverlap(loops2[i,], bed)
 tss1 <- getOverlap(loops1[i,], tss)
 tss2 <- getOverlap(loops2[i,], tss)

 ## Record the number of loops for each tfbs.
 nloops[indx1] <- nloops[indx1]+1
 nloops[indx2] <- nloops[indx2]+1

 ## Record the number of genes.
 if(NROW(tss1)>0) {
  ngenesl[indx2] <- ngenesl[indx2] + 1
  ngenes[indx1] <- ngenes[indx1] + 1
 }
 if(NROW(tss2)>0) {
  ngenesl[indx1] <- ngenesl[indx1] + 1
  ngenes[indx2] <- ngenes[indx2] + 1
 }

}

## Plot results.
plot(ecdf(mean_con[nloops>0]), col="red", xlim=c(-1,6)) ## Should match those above.
plot(ecdf(mean_con[nloops==0]), col="black", add=TRUE)

## Not any convincing additional conservation by looping to a gene.
## No loop/ Loop but no gene/ Loop and gene. 
plot(ecdf(mean_con[nloops>0 & rowSums(data.frame(ngenes, ngenesl))==0]), col="red", xlim=c(-1,6)) ## Should match those above.
plot(ecdf(mean_con[nloops==0]), col="black", add=TRUE)
plot(ecdf(mean_con[rowSums(data.frame(ngenes, ngenesl))>0]), col="dark blue", add=TRUE) 

## No loop/ Loop but no gene/ Loop and gene AND no gene on current side.
plot(ecdf(mean_con[nloops>0 & ngenesl==0 & ngenes==0]), col="red", xlim=c(-1,6)) ## Should match those above.
plot(ecdf(mean_con[nloops==0 & ngenes==0]), col="black", add=TRUE)
plot(ecdf(mean_con[ngenesl>0 & ngenes==0]), col="dark blue", add=TRUE)            


## DONE!


