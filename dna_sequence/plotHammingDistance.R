require(bigWig)

doplot <- function(bw_file, bed_file_conserved, bed_file_gain, bed_file_loss, stp=100, halfWindow=10000, ...) {
 bw <- load.bigWig(bw_file)
 bedg <- read.table(bed_file_gain)
 bedl <- read.table(bed_file_loss)
 bed_conserved <- read.table(bed_file_conserved)

 meta_g <- metaprofile.bigWig(center.bed(bedg, halfWindow, halfWindow), bw, step=stp)
 meta_l <- metaprofile.bigWig(center.bed(bedl, halfWindow, halfWindow), bw, step=stp)
 meta_conserved <- metaprofile.bigWig(center.bed(bed_conserved, halfWindow, halfWindow), bw, step=stp)

 signal_g <- meta_g$middle
 signal_l <- meta_l$middle
 signal_conserved <- meta_conserved$middle

 ylim=c(min(c(signal_g, signal_l, signal_conserved)), max(c(signal_g, signal_l, signal_conserved)))

 N = length(signal_conserved)
 x = ((1:N) - N/2)* stp #1:N*stp

 plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal", ...) ## ylab= "Divergences/ 100 bp"
 lines(x, signal_g, col="dark red")
 lines(x, signal_l, col="dark blue")
 lines(x, signal_conserved, col="dark green")

}


tfbsdoplot <- function(bw, tf_name) {
  doplot(bw, paste(tf_name,".Hgain-loss.TFBS.bed.gz", sep=""), paste(tf_name,".conserved.TFBS.bed.gz", sep=""), halfWindow=50, stp=1, main=tf_name)
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
tfbsdoplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "CREB")

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



## Changes.
#doplot("hs.ls.diff.bigWig", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="Human divergences, Human gains/ losses.") 
#doplot("pt.ls.diff.bigWig", "C-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="Chimp divergences, Chimp gains/ losses.")
#doplot("rm.ls.diff.bigWig", "M-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="R. macaque divergences, R. macaque gains/ losses.")
#doplot("rm.ls.diff.bigWig", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="R. macaque divergences, Human gains/ losses.") ## Control.

dev.off()

