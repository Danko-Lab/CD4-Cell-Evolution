##
## Obtained here: https://bitbucket.org/soccin/seqc/src/3c971e74c9a5df35f87880914e5168767870037c/src/run_limma.R?at=master
##
## Reference: 
##  Rapaport, F., Khanin, R., Liang, Y., Pirun, M., Krek, A., Zumbo, P., … Betel, D. (2013). 
##  Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data. 
##  Genome Biology, 14(9), R95. doi:10.1186/gb-2013-14-9-r95
## 

removePC <- function(dat, i, plot.pc=TRUE) {
  pca <- prcomp(dat, center=FALSE, scale=FALSE) ## UNT

  cols <- c(rep("red",3), rep("green",2), rep("blue", 3), rep("dark red", 3), rep("dark green", 2), rep("dark blue", 2), "black", "black")
  pch <- c(rep(19,8), rep(6,7), 9, 24)

  summary(pca) # Prints variance summary for all principal components.
  #plot(pca$rotation[,1], pca$rotation[,2], pch=19, col=cols)
  if(plot.pc) {
    pairs(pca$rotation[,1:5], col=cols, pch=pch)
  }

  n <- NROW(pca$rotation)
  adj_df <- ((pca$x[,c(1:(i-1),(i+1):n)] %*% t(pca$rotation[,c(1:(i-1),(i+1):n)])))

  return(adj_df)
}

runLimmaQuantile <- function(count.dat, conditions, genes, condA, condB, 
                             q.cut=0.01, lfc=0.0, lib.size=colSums(count.dat),
                             useVoom=FALSE, plotMA=FALSE, remove.pc=NULL){
  require(limma)
  ## updated limma from 3.12.1 to 3.12.3 to 3.14.1 CGD: Now to 3.18.13
  if(! packageDescription("limma")$Version >= "3.18.13"){
    stop("Wrong version of limma package. This script requires  limma-3.18.3")
  }

  groups=as.factor(conditions)
  design <- model.matrix(~0+groups)
  colnames(design) = levels(groups)

  if(useVoom){ ### THIS PART REQUIRES UPDATING!
    require(edgeR)
    nf <- calcNormFactors(count.dat)
    dat <- voom(count.dat, design, plot=FALSE, lib.size=lib.size * nf)
  } else{
    counts.log.dat=as.matrix(log2(count.dat+1)) ## CGD: Cast appears to be required in Limma 3.18.13
    counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')

    dat=counts.log.norm.dat
#    plotDensities(dat)
  }
  ## PCA.
#  if(!is.null(remove.pc)) {
#    dat <- removePC(dat, remove.pc)
#  }

  ## DGE.  
  fit=lmFit(dat,design)
  contrast.matrix <- makeContrasts(contrasts= paste(condA,"-",condB), levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  res=decideTests(fit2,p.value=q.cut,lfc=lfc)
  tab<-topTable(fit2, adjust = "BH", number=nrow(fit2), sort.by='none')

  if(plotMA) {
#   plotMA(fit2, array=2, status=as.factor(decideTests(fit2)[,2]), col=c("black", "red", "blue"))
   status <- as.character(p.adjust(fit2$p.value, method="fdr")<PVAL)
#   status[genes[,4] == "chr5_140005300_140013300"] <- "zCD14"
   plotMA(fit2, array=2, status=status, col=c("black", "red", "blue"))
   abline(h=0, col="blue")

#   plotMDS(dat,top=500,labels=species,gene.selection="common")
  }

#  rownames(tab)=tab[,1]
  if(useVoom){
    counts.limma=2^dat$E[match(rownames(tab), rownames(dat$E)),] ## CGD: Changed to re-order using rownames ...
  }else{
    counts.limma=2^dat[match(rownames(tab), rownames(dat)),]
  }
  genes= genes[match(rownames(tab), rownames(genes)),]
  val1=apply(counts.limma[,which(conditions==condA)],1,mean)
  val2=apply(counts.limma[,which(conditions==condB)],1,mean)

  tab=tab[,c("P.Value","adj.P.Val","logFC")]
  tab=as.matrix(cbind(tab,val1,val2))
  nam1=paste("mean_",condA,sep="")
  nam2=paste("mean_",condB,sep="")
  colnames(tab)[4]=nam1
  colnames(tab)[5]=nam2
  tab= cbind(genes, tab)

  return(list(tab=tab,res=res,counts=counts.limma))
}
