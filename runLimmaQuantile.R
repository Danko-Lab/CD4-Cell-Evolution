##
## Obtained here: https://bitbucket.org/soccin/seqc/src/3c971e74c9a5df35f87880914e5168767870037c/src/run_limma.R?at=master
##
## Reference: 
##  Rapaport, F., Khanin, R., Liang, Y., Pirun, M., Krek, A., Zumbo, P., … Betel, D. (2013). 
##  Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data. 
##  Genome Biology, 14(9), R95. doi:10.1186/gb-2013-14-9-r95
## 

runLimmaQuantile <- function(count.dat, conditions, genes,
                             useVoom=FALSE,
                             q.cut=0.05, lfc=0.0, condA='condA', condB='condB'){
  require(limma)
  ## updated limma from 3.12.1 to 3.12.3 to 3.14.1 CGD: Now to 3.18.13
  if(! packageDescription("limma")$Version == "3.18.13"){
    stop("Wrong version of limma package. This script requires  limma-3.14.1")
  }

  groups=as.factor(conditions)
  design <- model.matrix(~0+groups)
  colnames(design) = levels(groups)
  
  if(useVoom){ ### THIS PART REQUIRES UPDATING!
    require(edgeR)
    nf <- calcNormFactors(count.dat)
    dat <- voom(count.dat, design, plot=FALSE, lib.size=colSums(count.dat) * nf)
  } else{
    counts.log.dat=as.matrix(log2(count.dat+1)) ## CGD: Cast appears to be required in Limma 3.18.13
    counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')

    dat=counts.log.norm.dat
#    plotDensities(dat)

  ## design <- model.matrix(~0+groups)
  ## colnames(design) = levels(groups)
  }
  
  fit=lmFit(dat,design)

  contrast.matrix <- makeContrasts("condB - condA", levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  res=decideTests(fit2,p.value=q.cut,lfc=lfc)
  tab<-topTable(fit2, adjust = "BH", number=nrow(fit2), sort.by='logFC')
#  rownames(tab)=tab[,1]
  if(useVoom){
    counts.limma=2^dat$E[match(rownames(tab), rownames(dat$E)),] ## CGD: Changed to re-order using rownames ...
  }else{
    counts.limma=2^dat[match(rownames(tab), rownames(dat$E)),]
  }
  val1=apply(counts.limma[,which(conditions==condA)],1,mean)
  val2=apply(counts.limma[,which(conditions==condB)],1,mean)

  tab=tab[,c("P.Value","adj.P.Val","logFC")]
  tab=as.matrix(cbind(tab,val1,val2))
  nam1=paste("mean_",condA,sep="")
  nam2=paste("mean_",condB,sep="")
  colnames(tab)[4]=nam1
  colnames(tab)[5]=nam2

  return(list(tab=tab,res=res,counts=counts.limma))
}
