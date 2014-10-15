##
## Draw circles to plot data...
##


cd.circle <- function(x, y, r) {
	rads <- pi*c(32:0)/64
	qdx <- r*cos(rads)
	qdy <- r*sin(rads)
	dx <- c(qdx, rev(qdx), -1*qdx, rev(-1*qdx))
	dy <- c(qdy, rev(-1*qdy), -1*qdy, rev(qdy))
	return(list(x= x+dx, y= y+dy))
} ## TEST!


cd.circplot<- function(data, names, fill="black") {
	lab <- pretty(c(min(data),max(data)), n=5)
        data <- (data - min(lab)) / (max(lab) - min(lab))
	labels <- unique(names)

	## Newpage.
	require(grid)
	grid.newpage()
	top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
	pushViewport(top.vp)

	## Set up the plot.
	width <- 1/NROW(labels)/20
	centerseq <- seq(width, 1, 1/NROW(labels))
	centers <- centerseq[as.integer(as.factor(names))]
	next.vp <- viewport(x= 0.1, y= 0.15, width = 0.9, height = 0.85, just=c("left","bottom"))
	pushViewport(next.vp)

	## Gridlines.
	for(i in lab) {
		grid.lines(y= (i-min(lab)) / (max(lab) - min(lab)), gp=gpar(col="light gray"))
	}
	for(i in 1:NROW(data)) {
		ym <- data[i]
		cc1 <- cd.circle(centers[i], ym, width)
		grid.polygon(cc1$x, cc1$y, gp=gpar(fill=fill))
	}

	## Medians/ Means
	for(i in 1:NROW(labels)) {
		indx <- levels(as.factor(names))==labels[i]
		grid.lines(centerseq[i]+c(-1/12/NROW(labels),1/12/NROW(labels)), mean(data[names==labels[indx]])) ## Draw horizontal line at mean.
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
	for(i in 1:NROW(labels)) {
		grid.text(labels[levels(as.factor(names))==labels[i]], x= centerseq[i], y= 0.95, rot = 65, just=c("right", "top"))
	}
	popViewport()
}
## test:
## cd.circplot(c(1:10)/10, c(rep("A",4), rep("B", 2), rep("C", 4)))

