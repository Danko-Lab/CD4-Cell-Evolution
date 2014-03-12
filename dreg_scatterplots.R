##
## Compare dREG scores across species.
##
##

a <- read.table("HCM-U-PI.dREG-tss-clusters.tsv")
for(i in 7:12) a[is.na(a[,i]),i] <- 0 ## Assign 'NA' scores to 0.

## Use densCols() output to get density at each point.
## Thanks to StackOverflow for the quicky R lift: http://stackoverflow.com/questions/17093935/r-scatter-plot-symbol-color-represents-number-of-overlapping-points
densScatterplot <- function(x1, x2, ...) {
  df <- data.frame(x1, x2)

  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
 
  ## Map densities to colors
#  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <- colorRampPalette(c("light gray", "#000099", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(128)
  df$col <- cols[df$dens]

  ## Plot it, reordering rows so that densest points are plotted on top
  plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, ...) 
}

densScatterplot(a$V7, a$V8, xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")
densScatterplot(a$V7, a$V9, xlab="Human Unt.", ylab="Rhesus Macaque Unt.", main="dREG Scores") 

densScatterplot(a$V7, a$V10, xlab="Human Unt.", ylab="Human P/I", main="dREG Scores") 

