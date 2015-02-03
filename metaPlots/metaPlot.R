require(bigWig)
DATA_PATH="/home/cgd24/NHP/AllData/All_Merge/"

HP <- load.bigWig(paste(DATA_PATH, "H-U_plus.bw", sep=""))
HM <- load.bigWig(paste(DATA_PATH, "H-U_minus.bw", sep=""))
CP <- load.bigWig(paste(DATA_PATH, "C-U_plus.hg19.bw", sep=""))
CM <- load.bigWig(paste(DATA_PATH, "C-U_minus.hg19.bw", sep=""))
MP <- load.bigWig(paste(DATA_PATH, "M-U_plus.hg19.bw", sep=""))
MM <- load.bigWig(paste(DATA_PATH, "M-U_minus.hg19.bw", sep=""))

doit <- function(bed, stp=100, halfWindow= 2500, ...) {
	bed <- bed[grep("Un|random", bed$V1, invert=TRUE),]
	bed_rev <- bed; bed_rev[bed[,6] == "+",6] <- "-"; bed_rev[bed[,6] == "-",6] <- "+"

	H_meta_p <- metaprofile.bigWig(bed, HP, HM, step=stp)
	H_meta_m <- metaprofile.bigWig(bed_rev, HP, HM, step=stp)

        C_meta_p <- metaprofile.bigWig(bed, CP, CM, step=stp)
        C_meta_m <- metaprofile.bigWig(bed_rev, CP, CM, step=stp)

        M_meta_p <- metaprofile.bigWig(bed, MP, MM, step=stp)
        M_meta_m <- metaprofile.bigWig(bed_rev, MP, MM, step=stp)

        N = length(H_meta_p$middle)
        x = ((1:N) - N/2)* stp
        ylim=c(-1*max(c(H_meta_m$top, C_meta_m$top, M_meta_m$top)), max(c(H_meta_p$top, C_meta_p$top, M_meta_p$top)))
	
	par(mfrow=c(1,3))
	plot.metaprofile(H_meta_p, minus.profile=H_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(C_meta_p, minus.profile=C_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(M_meta_p, minus.profile=M_meta_m, X0=halfWindow/stp, ylim=ylim)

}

bed <- center.bed(read.table("myc.bed.gz"), upstreamWindow=halfWindow, downstreamWindow=halfWindow)
doit(bed)

bed <- center.bed(read.table("HMBOX.bed.gz"), upstreamWindow=halfWindow, downstreamWindow=halfWindow)
doit(bed)

gc <- fiveprime.bed(read.table("../annotations/gencode.v18.transcript.tsv"), upstreamWindow=halfWindow, downstreamWindow=halfWindow)
gc <- gc[sample(which(gc$V7 == "protein_coding"), 500),]
doit(gc)
