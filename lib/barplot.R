## Barplot

cd.barplot<- function(data, error, names, fill) {
  lab <- pretty(c(0,1), n=5)

  ## Reorder.
  ord <- order(data)
  data <- data[ord]
  names <- names[ord]
  error <- error[ord]
  
  ## Newpage.
  require(grid)
  grid.newpage()
  top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
  pushViewport(top.vp)

  ## Plot histone bars/ ... .
  width <- 1/NROW(data)/3
  centers <- seq(width, 1, 1/NROW(data))

  next.vp <- viewport(x= 0.1, y= 0.15, width = 0.9, height = 0.85, just=c("left","bottom"))
  pushViewport(next.vp)
  
  ## Gridlines.
  for(i in lab) {
    grid.lines(y= i, gp=gpar(col="light gray"))
  }
  
  for(i in 1:NROW(data)) {
	ym <- data[i]
    x <- c(centers[i]+width, centers[i]-width, centers[i]-width, centers[i]+width)
    y <- c(0, 0, ym, ym)
    grid.polygon(x, y, gp=gpar(fill=fill))

   ## Draw std. error.
    grid.lines(x=centers[i], y=c(ym, ym+error[i]))
    grid.lines(x=x[1:2], y=(ym+error[i]))
  }
  popViewport()
  
  ## Plot Y axis.
  next.vp <- viewport(x= 0.07, y= 0.15, width = 0.03, height = 0.85, just=c("left","bottom"))
  pushViewport(next.vp)
    grid.yaxis(at= seq(0, 1, length= length(lab)), label= lab)
  popViewport()
  ## Plot X axis labels.
  next.vp <- viewport(x= 0.1, y= 0, width = 0.9, height = 0.15, just=c("left","bottom"))
  pushViewport(next.vp)
  for(i in 1:NROW(data)) {
    grid.text(names[i], x= centers[i], y= 0.95, rot = 65, just=c("right", "top"))
  }
  popViewport()

}

