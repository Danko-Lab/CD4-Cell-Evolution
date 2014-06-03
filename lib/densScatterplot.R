
## Use densCols() output to get density at each point.
## Thanks to StackOverflow for the quicky R lift: http://stackoverflow.com/questions/17093935/r-scatter-plot-symbol-color-represents-number-of-overlapping-points
densScatterplot <- function(x1, x2, n=256, ...) {
  df <- data.frame(x1, x2)

  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L

  ## Map densities to colors
#  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <- colorRampPalette(c("light gray", "#000099", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(n)
  df$col <- cols[df$dens]

  ## Plot it, reordering rows so that densest points are plotted on top
  plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, ...)
}


