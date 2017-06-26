source("../lib/circplot.R")

pdf("SGPP2.pdf")
#pdf("18s.pdf")

dat <- read.table("data.tsv", sep="\t", header=TRUE)
data <- as.matrix(dat[,c("R1", "R2", "R3")])#, "R4", "R5", "R6")])

dat18s <- read.table("18s.tsv", sep="\t", header=TRUE)
data18s <- as.matrix(dat18s[,c("R1", "R2", "R3")])#, "R4", "R5", "R6")])

## Sanity check ... no 18s.
#data <- data18s

## Confirm that no significant correlation exists.
#cor.test(c(data), c(data18s), method="spearman", na.rm=TRUE)
data <- data/ data18s

## Another sanity check ... scores are really driven by the one outlier: G4, R3.  
## I believe the 40m and 0m conditions were switched inprocessing ... 
#data[4,3] <- NA; data[12,3] <- NA
#data[5,3] <- NA; data[15,3] <- NA

tvc <- 12 #9

#############################################################
## Combine data from gRNAs for each TRE of interest ...
IE0_indx <- dat$Site == "IE" & dat$Time == 0
IE0 <- c(log(data[IE0_indx,]/rbind(data[1,], data[1,], data[1,]), 2)); IE0 <- IE0[!is.na(IE0)]
t.test(x= IE0) ## Paired t.test for IE.

UE0_indx <- dat$Site == "UE" & dat$Time == 0
UE0 <- c(log(data[UE0_indx,]/rbind(data[1,], data[1,], data[1,], data[1,]), 2)); UE0 <- UE0[!is.na(UE0)]
t.test(x= UE0) ## Paired t.test for UE.

PR0_indx <- dat$Site == "PR" & dat$Time == 0
PR0 <- c(log(data[PR0_indx,]/rbind(data[1,], data[1,]), 2)); PR0 <- PR0[!is.na(PR0)]
t.test(x= PR0) ## Paired t.test for PR.

cd.circplot( c(as.double(IE0), as.double(UE0), as.double(PR0)), 
			c(rep("IE", NROW(IE0)), rep("UE", NROW(UE0)), rep("PR", NROW(PR0))), 
			fill="dark gray", title= "enhancer", jidder=TRUE, lims=c(4, -4))


#############################################################
## Combine data from 40 min.
IE40_indx <- dat$Site == "IE" & dat$Time == 40
IE40 <- c(log(data[IE40_indx,]/rbind(data[tvc,], data[tvc,], data[tvc,]), 2))
t.test(x= IE40) ## Paired t.test for IE.

UE40_indx <- dat$Site == "UE" & dat$Time == 40
UE40 <- c(log(data[UE40_indx,]/rbind(data[tvc,], data[tvc,], data[tvc,], data[tvc,]), 2))
t.test(x= UE40) ## Paired t.test for UE.

PR40_indx <- dat$Site == "PR" & dat$Time == 40
PR40 <- c(log(data[PR40_indx,]/rbind(data[tvc,], data[tvc,]), 2))
t.test(x= PR40) ## Paired t.test for PR.

cd.circplot( c(as.double(IE40), as.double(UE40), as.double(PR40)), 
			c(rep("IE", NROW(IE40)), rep("UE", NROW(UE40)), rep("PR", NROW(PR40))), 
			fill="dark gray", title= "enhancer", jidder=TRUE, lims=c(4, -4))


#############################################################
## Compare 0 and 40 min.
CTRL040 <- c(log(data[tvc,]/data[1,], 2)); CTRL040 <- CTRL040[!is.na(CTRL040)]
t.test(x= CTRL040) ## Paired t.test for IE.

IE040 <- c(log(data[IE40_indx,]/data[IE0_indx,], 2)); IE040 <- IE040[!is.na(IE040)]
t.test(x= IE040) ## Paired t.test for IE.

UE040 <- c(log(data[UE40_indx,]/data[UE0_indx,], 2)); UE040 <- UE040[!is.na(UE040)]
t.test(x= UE040) ## Paired t.test for UE.

PR040 <- c(log(data[PR40_indx,]/data[PR0_indx,], 2)); PR040 <- PR040[!is.na(PR040)]
t.test(x= PR040) ## Paired t.test for PR.

cd.circplot( c(as.double(CTRL040), as.double(IE040), as.double(UE040), as.double(PR040)), 
			c(rep("VC", NROW(CTRL040)), rep("IE", NROW(IE040)), rep("UE", NROW(UE040)), rep("PR", NROW(PR040))), 
			fill="dark gray", title= "enhancer", jidder=TRUE, lims=c(5, -5))

			
			
#############################################################
## Include all in the same plot.

cd.circplot( c(as.double(IE0), as.double(UE0), as.double(PR0), as.double(CTRL040), as.double(IE040), as.double(UE040), as.double(PR040)), 
			c(rep("IE-0", NROW(IE0)), rep("UE-0", NROW(UE0)), rep("PR-0", NROW(PR0)), rep("VC-40", NROW(CTRL040)), rep("IE-40", NROW(IE040)), rep("UE-40", NROW(UE040)), rep("PR-40", NROW(PR040))), 
			fill="dark gray", title= "enhancer", jidder=TRUE, lims=c(5, -5))

dev.off()

pdf("barplots.pdf")

#############################################################
## Now use barplots.
require(ggplot2)

ds <- data.frame(

## Names: 
enhancer= factor(rep(c("Vc", "Int Enhancer", "Ups Enhancer", "Promoter"), 2),
			levels=c("Vc", "Ups Enhancer", "Promoter", "Int Enhancer")),

timepoint= factor(c(0, 0, 0, 0, 40, 40, 40, 40), levels=c("0", "40")),

## Compute median.
medians= c( 1,
	median(data[IE0_indx,]/data[1,], na.rm=TRUE), 
	median(data[UE0_indx,]/data[1,], na.rm=TRUE),
	median(data[PR0_indx,]/data[1,], na.rm=TRUE),

	median(data[tvc,]/data[1,], na.rm=TRUE),
	median(data[IE40_indx,]/data[1,], na.rm=TRUE),
	median(data[UE40_indx,]/data[1,], na.rm=TRUE),
	median(data[PR40_indx,]/data[1,], na.rm=TRUE)),

## Compute SEM.
sem= c( 0,
	sqrt(var(as.double(data[IE0_indx,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[IE0_indx,])),
	sqrt(var(as.double(data[UE0_indx,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[UE0_indx,])),
	sqrt(var(as.double(data[PR0_indx,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[PR0_indx,])),

        sqrt(var(as.double(data[tvc,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[tvc,])),
        sqrt(var(as.double(data[IE40_indx,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[IE40_indx,])),
        sqrt(var(as.double(data[UE40_indx,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[UE40_indx,])),
        sqrt(var(as.double(data[PR40_indx,]/data[1,]), na.rm=TRUE))/sum(!is.na(data[PR40_indx,]))
	)
)

ggplot(ds, aes(enhancer, medians, fill= timepoint))+scale_y_continuous(limits=c(0,10))+
	theme_bw()+theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))+scale_fill_manual(values=c("#CCCCCC","#FFFFFF"))+
	geom_bar(position=position_dodge(0.9), colour="black", stat = "identity")+
	geom_errorbar(position=position_dodge(0.9), aes(ymin= medians-sem, ymax= medians+sem, width=0.3))+
	geom_hline(yintercept=1)+ylab("Fold over Vc")+xlab("gRNA Target")


ggplot(ds[ds$timepoint == 0,], aes(enhancer, medians, fill= timepoint))+scale_y_continuous(limits=c(0,1.5))+
	theme_bw()+theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"))+scale_fill_manual(values=c("#CCCCCC","#FFFFFF"))+
        geom_bar(position=position_dodge(0.9), colour="black", stat = "identity")+
        geom_errorbar(position=position_dodge(0.9), aes(ymin= medians-sem, ymax= medians+sem, width=0.3))+
        geom_hline(yintercept=1)+ylab("Fold over Vc")+xlab("gRNA Target")

dev.off()

#################################################################
			
			

indx <- dat$Site == "IE" & dat$Time == 40
t.test(x= c(log(as.double(data[indx,])/rep(data[9,], 2))) ) ## Paired t.test for IE.

indx <- dat$Site == "UE" & dat$Time == 40
#c(log(as.double(data[indx,])/rep(data[1,], 3)))
t.test(x= c(log(as.double(data[indx,])/rep(data[9,], 3))) ) ## Paired t.test for UE.

indx <- dat$Site == "PR" & dat$Time == 40
t.test(x= c(log(as.double(data[indx,])/rep(data[9,], 2))) ) ## Paired t.test for PR.



data.frame(dat, p.value= -log(sapply(1:NROW(data), function(x) {t.test(x= data[1,], y= data[x,], paired=TRUE)$p.value}), 10))

data.frame(dat, p.value= c(1,sapply(2:NROW(data), function(x) {t.test(x= data[1,], y= data[x,], paired=FALSE)$p.value})))

data <- cbind(data, data)

data.frame(dat, p.value= sapply(1:NROW(data), function(x) {t.test(x= data[12,], y= data[x,], paired=TRUE)$p.value}))

