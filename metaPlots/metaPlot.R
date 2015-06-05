require(bigWig)
DATA_PATH="/home/cgd24/NHP/AllData/All_Merge/"

HP <- load.bigWig(paste(DATA_PATH, "H-U_plus.rpkm.bw", sep=""))
HM <- load.bigWig(paste(DATA_PATH, "H-U_minus.rpkm.bw", sep=""))
CP <- load.bigWig(paste(DATA_PATH, "C-U_plus.hg19.rpkm.bw", sep=""))
CM <- load.bigWig(paste(DATA_PATH, "C-U_minus.hg19.rpkm.bw", sep=""))
MP <- load.bigWig(paste(DATA_PATH, "M-U_plus.hg19.rpkm.bw", sep=""))
MM <- load.bigWig(paste(DATA_PATH, "M-U_minus.hg19.rpkm.bw", sep=""))

HPpi <- load.bigWig(paste(DATA_PATH, "H-PI_plus.rpkm.bw", sep=""))
HMpi <- load.bigWig(paste(DATA_PATH, "H-PI_minus.rpkm.bw", sep=""))
CPpi <- load.bigWig(paste(DATA_PATH, "C-PI_plus.hg19.rpkm.bw", sep=""))
CMpi <- load.bigWig(paste(DATA_PATH, "C-PI_minus.hg19.rpkm.bw", sep=""))
MPpi <- load.bigWig(paste(DATA_PATH, "M-PI_plus.hg19.rpkm.bw", sep=""))
MMpi <- load.bigWig(paste(DATA_PATH, "M-PI_minus.hg19.rpkm.bw", sep=""))

doit <- function(bed, HP, HM, CP, CM, MP, MM, stp, halfWindow, ...) {
	bed <- bed[grep("Un|random", bed$V1, invert=TRUE),]

	H_meta_p <- metaprofile.bigWig(bed, HP, HM, step=stp)
	H_meta_m <- metaprofile.bigWig(bed, HM, HP, step=stp)

        C_meta_p <- metaprofile.bigWig(bed, CP, CM, step=stp)
        C_meta_m <- metaprofile.bigWig(bed, CM, CP, step=stp)

        M_meta_p <- metaprofile.bigWig(bed, MP, MM, step=stp)
        M_meta_m <- metaprofile.bigWig(bed, MM, MP, step=stp)

        N = length(H_meta_p$middle)
        x = 1:N*stp ## ((1:N) - N/2)* stp
        ylim=c(-1*max(c(H_meta_m$top, C_meta_m$top, M_meta_m$top)), max(c(H_meta_p$top, C_meta_p$top, M_meta_p$top)))
	
	par(mfrow=c(1,4))
	plot.metaprofile(H_meta_p, minus.profile=H_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(C_meta_p, minus.profile=C_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(M_meta_p, minus.profile=M_meta_m, X0=halfWindow/stp, ylim=ylim)

        plot.metaprofile(H_meta_p, minus.profile=H_meta_m, X0=halfWindow/stp, ylim=ylim)
	lines(x, C_meta_p$middle, col="#17b92b")
	lines(x, -1*C_meta_m$middle, col="#17b92b")
        lines(x, M_meta_p$middle, col="#5b6c0c")
        lines(x, -1*M_meta_m$middle, col="#5b6c0c")
}

stp=100; halfWindow= 2500

bed <- center.bed(read.table("myc.bed.gz"), upstreamWindow=halfWindow, downstreamWindow=halfWindow)
doit(bed, HP, HM, CP, CM, MP, MM, stp=stp, halfWindow=halfWindow)
doit(bed, HPpi, HMpi, CPpi, CMpi, MPpi, MMpi, stp=stp, halfWindow=halfWindow)

bed <- center.bed(read.table("HMBOX.bed.gz"), upstreamWindow=halfWindow, downstreamWindow=halfWindow)
doit(bed, stp=stp, halfWindow=halfWindow)
doit(bed, HPpi, HMpi, CPpi, CMpi, MPpi, MMpi, stp=stp, halfWindow=halfWindow)

gc <- fiveprime.bed(read.table("../annotations/gencode.v18.transcript.tsv"), upstreamWindow=halfWindow, downstreamWindow=halfWindow)
gc <- gc[sample(which(gc$V7 == "protein_coding"), 500),]
doit(gc, HP, HM, CP, CM, MP, MM, stp=stp, halfWindow=halfWindow)
doit(gc, HPpi, HMpi, CPpi, CMpi, MPpi, MMpi, stp=stp, halfWindow=halfWindow)

