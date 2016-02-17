require(bigWig)

doplot <- function(bw_file, bed_file, bed_file_conserved, stp=100, halfWindow=10000, ...) {
 bw <- load.bigWig(bw_file)
 bed <- read.table(bed_file)
 bed_conserved <- read.table(bed_file_conserved)

 meta_change <- metaprofile.bigWig(center.bed(bed, halfWindow, halfWindow), bw, step=stp)
 meta_conserved <- metaprofile.bigWig(center.bed(bed_conserved, halfWindow, halfWindow), bw, step=stp)

# plot.metaprofile(meta_change, X0=halfWindow/stp, ...)

 signal_change <- meta_change$middle#/ NROW(meta_change)
 signal_conserved <- meta_conserved$middle#/ NROW(meta_conserved)

 ylim=c(min(c(signal_change, signal_conserved)), max(c(signal_change, signal_conserved)))

 N = length(signal_change)
 x = ((1:N) - N/2)* stp #1:N*stp

 plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal", ...) ## ylab= "Divergences/ 100 bp"
 lines(x, signal_change, col="dark red")
 lines(x, signal_conserved, col="dark green")

}

doplot("/local/storage/data/hg19/all/phyloP100way/hg19.100way.phyloP100way.bw", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="PhyloP at Human LS dREG sites")

pdf("DNASequence.pdf")

## Changes.
doplot("hs.ls.diff.bigWig", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="Human divergences, Human gains/ losses.") 
doplot("pt.ls.diff.bigWig", "C-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="Chimp divergences, Chimp gains/ losses.")
doplot("rm.ls.diff.bigWig", "M-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="R. macaque divergences, R. macaque gains/ losses.")

doplot("rm.ls.diff.bigWig", "H-U.gain-loss.dREG_HD.bed.gz", "conserved.dREG_HD.bed.gz", main="R. macaque divergences, Human gains/ losses.") ## Control.

dev.off()
