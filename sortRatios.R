## Process command arguments.
args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
out_file <- args[2]
conf_th  <- as.double(args[3])

data <- read.table(in_file)
o.v4 <- as.double(as.character(data$V4))
o.v4[o.v4 < 0] <- 0
o.v4[is.nan(o.v4)] <- 0

o.v5 <- as.double(as.character(data$V5))
o.v5[o.v5 < 0] <- 0
o.v5[is.nan(o.v5)] <- 0

o.v6 <- as.double(as.character(data$V6))
o.v6[o.v6 < 0] <- 0
o.v6[is.nan(o.v6)] <- 0

maxCM <- sapply(c(1:NROW(data)), function(x) {return(max(o.v5[x], o.v6[x]))})
maxHC <- sapply(c(1:NROW(data)), function(x) {return(max(o.v4[x], o.v5[x]))})
maxHM <- sapply(c(1:NROW(data)), function(x) {return(max(o.v4[x], o.v6[x]))})

minCM <- sapply(c(1:NROW(data)), function(x) {return(min(o.v5[x], o.v6[x]))})
minHC <- sapply(c(1:NROW(data)), function(x) {return(min(o.v4[x], o.v5[x]))})
minHM <- sapply(c(1:NROW(data)), function(x) {return(min(o.v4[x], o.v6[x]))})

pc <- 0.01
ratioHg <- log((o.v4+pc)/(maxCM+pc))
ratioCg <- log((o.v5+pc)/(maxHM+pc))
ratioMg <- log((o.v6+pc)/(maxHC+pc))
ratioHl <- log((o.v4+pc)/(minCM+pc))
ratioCl <- log((o.v5+pc)/(minHM+pc))
ratioMl <- log((o.v6+pc)/(minHC+pc))

data <- cbind(data, ratioHg, ratioCg, ratioMg, ratioHl, ratioCl, ratioMl)

ord <- rev(order(ratioHg))
data <- data[ord,]

colnames(data) <- c("chrom", "chromStart", "chromEnd", "H_DREG", "C_DREG", "M_DREG", "Hgain_log_ratio", "Cgain_log_ratio", "Mgain_log_ratio", "Hloss_log_ratio", "Closs_log_ratio", "Mloss_log_ratio")

write.table(data, paste(out_file,".tsv", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

write.table(data[which((data$Hgain_log_ratio >  1) & (data$H_DREG > conf_th)), c(1:3)], 
	paste(out_file,".H-gain.",conf_th,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

write.table(data[which((data$Hloss_log_ratio < -1) & (data$C_DREG > conf_th) & (data$M_DREG > conf_th)), c(1:3)], 
	paste(out_file,".H-loss.",conf_th,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

write.table(data[which((data$H_DREG > conf_th) & (data$C_DREG > conf_th) & (data$M_DREG > conf_th)), c(1:3)], 
	paste(out_file,".conserved.",conf_th,".bed", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



