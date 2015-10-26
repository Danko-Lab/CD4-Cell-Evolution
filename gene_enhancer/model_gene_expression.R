## Returns indices in BED2 that intersect BED1.
source("../lib/getOverlap.R")

##Prepare loops.
dist<- 1000 # 15000
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

 as <- read.table(paste("../annotations/chage_expr/",prefix,post_enh, sep="")) ## post_enh b/c usually want regardless of change.
 as <- as[as$V7 == "ups_antisense", ]
 as[,2] <- as[,2] - 1000 ## Look nearby ... up to 1kb.
 as[,3] <- as[,3] + 1000

 ## Find out which loops intersect.
 enh_pro_change <- NULL
 for(i in c(1:NROW(tss))) {
  RE_changes <- double(0)
  loop_ <- 0
  near_ <- 0
  uas_  <- 0

  #print(i)
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

  ## Add nearby TREs.
  indxtss <- getOverlap(tss[i,], tres)
  indx <- getOverlap(nearby[i,], tres)
  truth_ <- rep(FALSE, NROW(tres)); truth_[indx] <- TRUE; truth_[indxtss] <- FALSE
  indx <- which(truth_)

  if(NROW(indx)>0) {
	near_ <- mean(tres[indx, column])
        RE_changes <- c(RE_changes, c(tres[indx, column]))
  }

  ## Add antisense direction.
  indx_uas <- getOverlap(tss[i,], as)

  if(NROW(indx_uas)>0) {
	uas_ <- mean(as[indx_uas, column])
	RE_changes <- c(RE_changes, c(as[indx_uas, column]))
  }

  enh_pro_change <- rbind(enh_pro_change, data.frame(pro=tss[i,column], enh=mean(RE_changes), near= near_, loop= loop_, uas= uas_))
 }

 return(enh_pro_change)
}

enh_pro_change <- rbind(getLoopNearby("H", 21),
                        getLoopNearby("C", 22),
                        getLoopNearby("M", 23))
cor.test(enh_pro_change$pro, enh_pro_change$enh)
plot(enh_pro_change$pro, enh_pro_change$enh, xlab= "Gene Expression", ylab="Mean Enhancers", pch=19)

require(vioplot)
vioplot(enh_pro_change$enh[enh_pro_change$pro < 0 & !is.nan(enh_pro_change$enh)], enh_pro_change$enh[enh_pro_change$pro > 0 & !is.nan(enh_pro_change$enh)]); abline(h=0)

## LM:
gl <- glm(pro~near+loop+uas, data=enh_pro_change)
cor.test(predict(gl), enh_pro_change$pro)

cor.test(enh_pro_change$pro, enh_pro_change$enh)

## Data for all...
indx <- (enh_pro_change$near!=0 & enh_pro_change$uas!=0)# & enh_pro_change$loop!=0)
indx2<- (enh_pro_change$near!=0 & enh_pro_change$uas!=0 & enh_pro_change$loop!=0)

indx2<- (enh_pro_change$near!=0 & enh_pro_change$uas!=0 & enh_pro_change$loop!=0)
cor.test(enh_pro_change$pro[indx2], enh_pro_change$enh[indx2])
plot(enh_pro_change$pro[indx2], enh_pro_change$enh[indx2])

## LM:
indx<- enh_pro_change$uas!=0
cor.test(enh_pro_change$pro[indx], enh_pro_change$enh[indx])
plot(enh_pro_change$pro[indx], enh_pro_change$enh[indx])

sm_cng <- enh_pro_change[indx,]

train <- sample(c(1:NROW(sm_cng)), NROW(sm_cng)*0.9)
test  <- rep(TRUE, NROW(sm_cng)); test[train] <- FALSE; test <- which(test)

gl <- glm(pro~near+loop+uas, data=sm_cng[train,])
cor.test(predict(gl, sm_cng[test,]), sm_cng$pro[test])




