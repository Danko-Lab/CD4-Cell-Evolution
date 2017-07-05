## Returns indices in BED2 that intersect BED1.
source("../lib/getOverlap.R")
load("Enhancer-Promoter-Loops.RData")
require(boot)

## Treat missing values as 0s?!
for(i in 2:NCOL(enh_pro_change)) enh_pro_change[is.na(enh_pro_change[,i]),i] <- 0

indx <- enh_pro_change$uas!=0 #|| enh_pro_change$near!=0
sm_cng <- enh_pro_change[indx,]

## FIrst get at R^2 values for all sites.
getDiff <- function(data, indices) {
  train <- sample(c(1:NROW(sm_cng)), NROW(sm_cng)*0.8)
  test  <- rep(TRUE, NROW(sm_cng)); test[train] <- FALSE; test <- which(test)
  #train <- rep(TRUE, NROW(sm_cng))
  #test  <- rep(TRUE, NROW(sm_cng))

  gl.tad <- glm(pro~near+uas+loop+as+tad+near:uas+near:tad+loop:tad+loop:near+uas:as+uas:near, data=sm_cng[train,])
  c.tad <- cor.test(predict(gl.tad, sm_cng[test,]), sm_cng$pro[test])$estimate*cor.test(predict(gl.tad, sm_cng[test,]), sm_cng$pro[test])$estimate

  gl.no_tad <- glm(pro~near+uas+loop+as+near:uas+near:near+uas:as+uas:near, data=sm_cng[train,])
  c.no_tad <- cor.test(predict(gl.no_tad, sm_cng[test,]), sm_cng$pro[test])$estimate*cor.test(predict(gl.no_tad, sm_cng[test,]), sm_cng$pro[test])$estimate

  c.tad - c.no_tad
}

results <- boot(data= sm_cng, statistic= getDiff, R=1000)
#boot.ci(results, type="basic")
median(results$t)
quantile(results$t)
quantile(results$t, probs= c(0.025, 0.975))

## Now estimate difference in model accuracy.  Use all sites for this...
## Use bootstrap to estimate difference.
getAIC <- function(data, indices) {
  gl.tad <- glm(pro~near+uas+loop+as+tad+near:uas+near:tad+loop:tad+loop:near+uas:as+uas:near, data=data[indices,])
  gl.no_tad <- glm(pro~near+uas+loop+as+near:uas+near:near+uas:as+uas:near, data=data[indices,])

  gl.tad$aic
  gl.no_tad$aic

  ## According to: https://en.wikipedia.org/wiki/Akaike_information_criterion
  ## interpretation is that this reflects the degree to which the gl.tad model is 
  ## more likely than gl.no_tad
 # exp((gl.tad$aic - gl.no_tad$aic)/2)

  gl.tad$aic - gl.no_tad$aic
}

results <- boot(data= sm_cng, statistic= getAIC, R=1000)

exp(results$t0/2)
exp(median(results$t)/2)
1/exp(quantile(results$t, probs= c(0.025, 0.975))/2)

#boot.ci(results, type="basic")


