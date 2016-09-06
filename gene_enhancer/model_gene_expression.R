## Returns indices in BED2 that intersect BED1.
source("../lib/getOverlap.R")

##Prepare loops.
dist<- 5000 # 15000 # 1000
loops <- read.table("/local/storage/data/hg19/cd4/chiapet_h3k4me2/H3K4me2_interact_hg19.bed.gz")
loopdist <- function(i) { ## Get the actual distance in the detected loop interaction.
 loop1 <- sapply(strsplit(as.character(loops$V4), split=";"), function(x) {x[[i]]})
 sapply(strsplit(loop1, split="-"), function(x) {as.double(x[[2]])-as.double(x[[1]])})
}
loops1<- loops; loops1[,2] <- loops[,2]+as.integer(loopdist(1)/2)-dist; loops1[,3] <- loops[,2]+as.integer(loopdist(1)/2)+dist
loops2<- loops; loops2[,2] <- loops[,3]+as.integer(loopdist(2)/2)-dist; loops2[,3] <- loops[,3]+as.integer(loopdist(2)/2)+dist
#loops1<- loops; loops1[,3] <- loops[,2]+loopdist(1) #+dist
#loops2<- loops; loops2[,2] <- loops[,3]-loopdist(2) #-dist
hist(log(loops[,3]-loops[,2], 10), 20) ## Summary stats on loop distances.

##Prepare TADS.
TADs <- read.table("/local/storage/data/hg19/gm12878/hic/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.bed.gz")


#######################################################################
## Get loops, nearby, and AS transcription
getLoopNearby <- function(prefix="H", column=21, post_pro= ".change-U.tsv", post_enh= ".change-U.all.tsv", dist= 25000) {
 ## Prune out loops < distance so they are not counted twice.
 loops1 <- loops1[loops[,3]-loops[,2] > dist, ]
 loops2 <- loops2[loops[,3]-loops[,2] > dist, ]

 ## Get TSS of changed, annotated protein coding genes.
 HC <- read.table(paste("../annotations/chage_expr/",prefix,post_pro, sep=""))
 genes <- HC[HC$V7 == "protein_coding",]
 tss   <- genes
 tss[tss[,6] == "+",2] <- tss[tss[,6] == "+",2]-250; tss[tss[,6] == "+",3] <- tss[tss[,6] == "+",2]+1
 tss[tss[,6] == "-",3] <- tss[tss[,6] == "-",3]+251; tss[tss[,6] == "-",2] <- tss[tss[,6] == "-",3]-1

 ## nearby ... 
 nearby <- tss
 nearby[,2] <- tss[,2]-dist
 nearby[,3] <- tss[,3]+dist

 ## dREG sites & antisense TUs.
 tres <- read.table(paste("../annotations/chage_expr/",prefix, post_enh, sep=""))
 tres <- tres[grep("dREG", tres$V10),]

 uas <- read.table(paste("../annotations/chage_expr/",prefix,post_enh, sep="")) ## post_enh b/c usually want regardless of change.
 uas <- uas[uas$V7 == "ups_antisense", ]
 uas[,2] <- uas[,2] - 2000 ## Look nearby ... up to 2kb.
 uas[,3] <- uas[,3] + 2000

 as <- read.table(paste("../annotations/chage_expr/",prefix,post_enh, sep="")) ## post_enh b/c usually want regardless of change.
 as <- as[as$V7 == "antisense", ]

 ## Find out which loops intersect.
 enh_pro_change <- NULL
 for(i in c(1:NROW(tss))) {
  RE_changes <- double(0)
  loop_ <- NA #0
  near_ <- NA #0
  uas_  <- NA #0
  as_   <- NA #0
  tad_  <- NA #0

  ## Get loop pairs 
  indx1 <- getOverlap(tss[i,], loops1)
  indx2 <- getOverlap(tss[i,], loops2)
  indxTRE <- integer(0)
  if(NROW(indx1)>0) {
        indxTRE <- c(indxTRE, getOverlap(loops2[indx1,], tres))
  }
  if(NROW(indx2)>0) {
        indxTRE <- c(indxTRE, getOverlap(loops1[indx2,], tres))
  }
  if(NROW(indxTRE)>0) {
	loop_ <- mean(tres[indxTRE, column])
        RE_changes <- c(tres[indxTRE, column])
  }

  ## Get TREs in the same TAD, excluding the TSS.
  indxtss <- getOverlap(tss[i,], tres)
  indxtad <- getOverlap(tss[i,], TADs)
  
  if(NROW(indxtad)>0) {
    indx <- getOverlap(TADs[indxtad[1],], tres) ## Assume one TAD per tss.
    truth_ <- rep(FALSE, NROW(tres)); truth_[indx] <- TRUE; truth_[indxtss] <- FALSE
    indx <- which(truth_)

    if(NROW(indx)>0) {
	tad_ <- mean(tres[indx, column])
	RE_changes <- c(RE_changes, c(tres[indx, column]))
    }
  }

  ## Add nearby TREs.
  indxtss <- getOverlap(tss[i,], tres)
  indx <- getOverlap(nearby[i,], tres)
  truth_ <- rep(FALSE, NROW(tres)); truth_[indx] <- TRUE; truth_[indxtss] <- FALSE
  indx <- which(truth_)

  if(NROW(indx)>0) {
	near_ <- mean(tres[indx, column])
#        RE_changes <- c(RE_changes, c(tres[indx, column])) ## Replaced by TREs in the same TAD on 8-31-2016
  }

  ## Add antisense direction.
  indx_uas <- getOverlap(tss[i,], uas)
  indx_as  <- getOverlap(genes[i,], as)

  if(NROW(indx_uas)>0) {
	uas_ <- mean(uas[indx_uas, column])
	RE_changes <- c(RE_changes, c(uas[indx_uas, column]))
  }
  if(NROW(indx_as)>0) {
	as_ <- mean(as[indx_as, column])
	RE_changes <- c(RE_changes, c(as[indx_as, column]))
  }

  enh_pro_change <- rbind(enh_pro_change, data.frame(pro=tss[i,column], enh=mean(RE_changes), tad= tad_, near= near_, loop= loop_, uas= uas_, as= as_))
 }

 return(enh_pro_change)
}

enh_pro_change <- rbind(getLoopNearby("H", 21),
                        getLoopNearby("C", 22),
                        getLoopNearby("M", 23))

save.image("Enhancer-Promoter-Loops.RData")

## Treat missing values as 0s?!
for(i in 2:NCOL(enh_pro_change)) enh_pro_change[is.na(enh_pro_change[,i]),i] <- 0

cor.test(enh_pro_change$pro, enh_pro_change$enh)
plot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Mean Enhancers", pch=19)

# Check correlation within TADs.
indx <- which(enh_pro_change$tad != 0)
cor.test(enh_pro_change$tad[indx], enh_pro_change$pro[indx], method="pearson")
plot(enh_pro_change$pro[indx], enh_pro_change$tad[indx], xlab= "Gene Expression", ylab="ncRNAs in TAD", pch=19)

# Check correlation of loops.
indx <- which(enh_pro_change$loop != 0)
cor.test(enh_pro_change$loop[indx], enh_pro_change$pro[indx], method="pearson")
plot(enh_pro_change$pro[indx], enh_pro_change$loop[indx], xlab= "Gene Expression", ylab="Loop", pch=19)

require(vioplot)
vioplot(enh_pro_change$enh[enh_pro_change$pro < 0 & !is.nan(enh_pro_change$enh)], enh_pro_change$enh[enh_pro_change$pro > 0 & !is.nan(enh_pro_change$enh)]); abline(h=0)

## LM:
gl <- glm(pro~near+loop+uas+enh+as, data=enh_pro_change)
cor.test(predict(gl, enh_pro_change), enh_pro_change$pro)
cor.test(enh_pro_change$pro, enh_pro_change$enh)

## LM:
indx<- rep(TRUE, NROW(enh_pro_change)) 
indx <- enh_pro_change$uas!=0 #|| enh_pro_change$near!=0
cor.test(enh_pro_change$pro[indx], enh_pro_change$enh[indx])
plot(enh_pro_change$pro[indx], enh_pro_change$enh[indx])

sm_cng <- enh_pro_change[indx,]

train <- sample(c(1:NROW(sm_cng)), NROW(sm_cng)*0.2)
test  <- rep(TRUE, NROW(sm_cng)); test[train] <- FALSE; test <- which(test)

gl <- glm(pro~near+uas+loop+as+tad+near:uas+near:tad+loop:tad+loop:near+uas:as+uas:near, data=sm_cng[train,])
#cor.test(sm_cng$enh[test], sm_cng$pro[test])
#cor.test(sm_cng$uas[test], sm_cng$pro[test])
cor.test(predict(gl, sm_cng[test,]), sm_cng$pro[test])
cor.test(predict(gl, sm_cng[train,]), sm_cng$pro[train])
plot(predict(gl, sm_cng[test,]), sm_cng$pro[test], pch=19); abline(0,1)

cor.test(predict(gl, sm_cng[test,]), sm_cng$pro[test], method="spearman")
source("../lib/densScatterplot.R")
densScatterplot(sm_cng$pro[test], predict(gl, sm_cng[test,]), xlab= "Fold Change in Gene Transcription", ylab="Predicted Change in Gene Transcription")
abline(a=0, b=1)


pdf("Predicting_Changes_in_transcription.pdf")
  densScatterplot(sm_cng$pro[test], predict(gl, sm_cng[test,]), xlab= "Fold Change in Gene Transcription", ylab="Predicted Change in Gene Transcription")
  abline(a=0, b=1)
dev.off()

#require(e1071)
#sv <- svm(pro~near+uas+loop+as+tad+enh+uas:enh+enh:loop, data=sm_cng[train,])
#cor.test(as.double(predict(sv, sm_cng[test,])), sm_cng$pro[test])


