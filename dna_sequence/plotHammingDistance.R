require(bigWig)
source("../lib/avg.metaprofile.R")

center.bed.one.off <- function(bed, halfWindow) {

  ## No strand, easy passthrough.
  if(dim(bed)[2] < 6) return(center.bed(bed, halfWindow, halfWindow))

  ## Get plus.
  plc <- center.bed(bed[bed[,6]=="+",], halfWindow, halfWindow)

  ## Correct one-off & get minus.
#  bed[,2]<-bed[,2]+1; bed[,3]<-bed[,3]+1
  plm <- center.bed(bed[bed[,6]=="-",], halfWindow, halfWindow)

  return(rbind(plc, plm))
#  return(plc)
}

exact.bed <- function(bed, halfWindow) {
  bed <- bed[bed[,3]-bed[,2] == median(bed[,3]-bed[,2]),] ## Smaller possible b/c of indels.
#  bed <- bed[bed[,6] == "-",] ## For sanity check the two return identical info.
  return(bed)
}

doplot <- function(bw_file, bed_file_conserved, bed_file_gain, bed_file_loss, bed.op = center.bed.one.off, stp=100, halfWindow=10000, ...) {
 bw <- load.bigWig(bw_file)
 bedg <- read.table(bed_file_gain)
 bedl <- read.table(bed_file_loss)
 bed_conserved <- read.table(bed_file_conserved)

 meta_g <- avg.metaprofile.bigWig(bed.op(bedg, halfWindow), bw, bw, step=stp, abs.value = FALSE) # bedg[,6] == "+",1:3
 meta_l <- avg.metaprofile.bigWig(bed.op(bedl, halfWindow), bw, bw, step=stp, abs.value = FALSE) # bedl[,6] == "+",1:3
 meta_conserved <- avg.metaprofile.bigWig(bed.op(bed_conserved, halfWindow), bw, bw, step=stp, abs.value = FALSE) # bed_conserved[bed_conserved[,6] == "+",1:3

 signal_g <- meta_g$middle
 signal_l <- meta_l$middle
 signal_conserved <- meta_conserved$middle

 ylim=c(min(c(signal_g, signal_l, signal_conserved)), max(c(signal_g, signal_l, signal_conserved)))

 N = length(signal_conserved)
 x = ((1:N) - N/2)* stp #1:N*stp

 if(stp==1 & halfWindow<30) typ="b"
 else typ="l"

 plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal", ...) ## ylab= "Divergences/ 100 bp"
 lines(x, signal_g, col="dark red", type=typ, pch=19)
 lines(x, signal_l, col="dark blue", type=typ, pch=19)
 lines(x, signal_conserved, col="dark green", type=typ, pch=19)

}


tfbsdoplot <- function(bw, tf_name) {
  doplot(bw, paste(tf_name,".conserved.TFBS.bed.gz", sep=""), paste(tf_name,".Hgain.TFBS.bed.gz", sep=""), paste(tf_name,".Hloss.TFBS.bed.gz", sep=""), halfWindow=25, stp=1, main=tf_name, bed.op= exact.bed)
}

pdf("TFBS.pdf")
doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "conserved.TFBS.bed.gz", "Hgain.TFBS.bed.gz", "Hloss.TFBS.bed.gz", halfWindow= 100, stp=1, main="All TFs, 100-way")
doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "conserved.TFBS.bed.gz", "Hgain.TFBS.bed.gz", "Hloss.TFBS.bed.gz", halfWindow= 100, stp=1, main="All TFs, Primate PhyloP")

tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "ELF1")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "GABPA")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "ARNT")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "YY1")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "MYCN")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "STAT2")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "CREB1")

dev.off()

pdf("TFBS.motif.strandtest.pdf")
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "STAT2")

dev.off()


pdf("TFBS.byDist.pdf")

doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "tfbs.all.dist_0-10kb.bed.gz", "tfbs.all.dist_10-100kb.bed.gz", "tfbs.all.dist_100-1000kb.bed.gz", halfWindow= 100, stp=1, main="All TFs, 100-way")
doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "tfbs.all.dist_0-10kb.bed.gz", "tfbs.all.dist_10-100kb.bed.gz", "tfbs.all.dist_100-1000kb.bed.gz", halfWindow= 100, stp=1, main="All TFs, 100-way")

dev.off()


pdf("PhyloP_100way.pdf")

doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "conserved.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")
doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "conserved-distal.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")

#doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")
#doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "C-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz")
#doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "M-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz")

dev.off()

pdf("PrimatePhyloP.pdf")

doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "conserved.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")
doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "conserved-distal.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")

#doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")
#doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "C-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="PhyloP at Chimp LS dREG sites")
#doplot("/local/storage/data/hg19/all/phylopprimate/phyloP46way.primate.bw", "M-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="PhyloP at R.Macaque LS dREG sites")

dev.off()

pdf("DNASequence.pdf")

doplot("/home/cgd24/cbsudanko/data/hg19/all/neanderthal_ss_green/ntSssZScorePMVar.bw", "conserved.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="Neanderthal Green et. al. S scores at Human LS dREG sites", stp=500, halfWindow=25000)

doplot("hs.ls.diff.bigWig", "conserved.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="Human SNPs at Human LS dREG sites", stp=300, halfWindow=25000)
doplot("hs.ls.diff.bigWig", "conserved-distal.dREG_HD.bed.gz", "H-U.gain.dREG_HD.bed.gz", "H-U.loss.dREG_HD.bed.gz", main="Human SNPs at Human LS dREG sites", stp=300, halfWindow=25000)

doplot("hs.ls.diff.bigWig", "conserved.TFBS.bed.gz", "Hgain.TFBS.bed.gz", "Hloss.TFBS.bed.gz", halfWindow= 100, stp=1, main="All TFs, 100-way")
tfbsdoplot("hs.ls.diff.bigWig", "ELF1")
tfbsdoplot("hs.ls.diff.bigWig", "STAT2")

## Changes.
#doplot("hs.ls.diff.bigWig", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="Human divergences, Human gains/ losses.") 
#doplot("pt.ls.diff.bigWig", "C-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="Chimp divergences, Chimp gains/ losses.")
#doplot("rm.ls.diff.bigWig", "M-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="R. macaque divergences, R. macaque gains/ losses.")
#doplot("rm.ls.diff.bigWig", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="R. macaque divergences, Human gains/ losses.") ## Control.

dev.off()

