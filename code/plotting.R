#' Formatted time series plot
#'
#' @param x vector length n. x-axis values
#' @param y vector length n. y-axis values
#' @param X vector length m. x-axis values for additional lines plotted in gray
#' @param Y matrix mXk. y-axis values for `k` additional lines plotted in gray
#' @param mar 
#' @param ... 
plot_ts <- function(x, y, X=c(), Y=c(), mar=c(2, 0, 0, 0) + 0.1, bty="n", showAxis=T, ...) {
  par(mar=mar)
  plot(
    range(c(X, x), na.rm = T),
    range(c(Y, y), na.rm = T),
    type = "n", bty = bty, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
    ...
  )
  if (showAxis) axis.POSIXct(side=1,
                             at=as.POSIXct(c(seq(from=min(x),to=max(x),by=4*3600), max(x)), origin="1970-01-01"), 
                             format = "%H:%M")
  if (!is.null(X) & !is.null(Y)) {
    for(i in 1:dim(Y)[2]) lines(X, Y[,i],col="gray")
  }
  
  lines(x, y, lwd=3)
}

#' Formatted scatterplot
#'
#' @param x nX2 matrix or data frame
scatterplot <- function(x, xlab="", ylab="", bg="darkgray", pch=20, xaxt="n", yaxt="n", bty="n", cex = 0.3, 
                  mar=c(0.5,0.5,0,0)+0.1, asp=1, tick.num=c(4,4), showAxes=T, ...) {
  par(mar=mar)
  plot(x, xlab=xlab, ylab=ylab, bg=bg, pch=pch, xaxt=xaxt, yaxt=yaxt, cex=cex, asp=asp, bty=bty, ...)
  if (showAxes) {
    axis(side=1,
         at=seq(from=min(x[,1]), to=max(x[,1]), length.out=tick.num[1]),
         tck=-0.02
    )
    axis(side=2,
         at=seq(from=min(x[,2]), to=max(x[,2]), length.out=tick.num[2]),
         tck=-0.02
    )
  }
}

#' Formatted polygon
#'
#' @param pol 
plot_polygon <- function(pol, lwd=2, border=rgb(0.8, 0, 0, 0.5), ...) {
  polygon(pol, lwd=lwd, border=border, ...)
}

#' Formatted scatterplot with surrogates
#' Copied from kdd-poster-video-figures.Rmd on 2019-07-15.
#'
#' @param x 
#' @param surr 
#' @param ... 
#' @param mar 
#' @param cex 
#' @param pch 
#'
#' @return
scatterplot_surrogates <- function(x, surr=c(), ..., cex=0.3, pch=20) {
  scatterplot(x=x, ..., cex=cex, pch=pch)
  if (!is.null(surr)) {
    for (xy_surr in surr) points(xy_surr, cex=0.5*cex, pch=pch, col=rgb(0.6,0.6,0.6,0.5))
    points(x, cex=cex, pch=pch)
  }
}

#' Plot time series of particle concentration. Copied from plotting.R on 2019-07-25.
#'
#' @param data data.frame. First column is timestamp, other columns are particle concentrations (cm-3) for different particle sizes (nm)
#' If data is matrix without timestamp then workaround: banana.data[rownames(data),]
#' @param nlevels 
#' @param levels 
#' @param bare plot without axes
#' @param mar margins into par()
#' @param plot.title 
#' @param ... 
#'
#' @return contour plot x=time, y=particle sizes (nm), color=particle concentration (cm-3)
banana_plot <- function(data, nlevels=50, levels=seq(0, 10, 0.2),
                        bare=FALSE, mar=c(2,4,0.5,0.5) + 1, 
                        plot.title, ...) {
  x <- data$timestamp
  y <- name2num(colnames(data[-1])) # Particle sizes from colnames.
  z <- df2matrix_plus1_log(data[-1])
  
  clr_palette <- rev(RColorBrewer::brewer.pal(9, "YlGnBu"))
  clr <- colorRampPalette(clr_palette)
  
  if (is.null(levels) | (max(z, na.rm=T) > max(levels))) levels <- pretty(range(z, finite = TRUE), nlevels)
  
  if (bare) mar <- c(0,2,0,2)
  par(mar = mar)
  plot.new()
  plot.window(xlim=range(x, finite = TRUE), 
              ylim=range(log(y), finite = TRUE), 
              #log="y", 
              xaxs="i", yaxs="i")
  
  .filled.contour(x, log(y), z, as.double(levels), 
                  col = clr(length(levels)-1))
  if (!bare) {
    subtitle <- paste0("Concentrations range (1/cm3):", paste0(round(exp(range(levels))-1, 0), collapse=" "))
    if (missing(plot.title)) title(main="Particle concentration (cm-3)", ylab="Particle size (nm)") else plot.title
    title(sub=subtitle)
    xaxis_ticks <- c(min(x), pretty(x), max(x))
    xaxis_ticks <- pretty(x, n=6)
    axis(1, at=xaxis_ticks, label = format(as.POSIXct(xaxis_ticks), "%H:%M"))
    #axis(2, at=log(pretty(y)), label=pretty(y))
    axis(2, at=log(10^(0:5)), label=10^(0:5))
  }
  # key.axes = axis(4, at = log(10^(0:5)), labels=10^(0:5)) # legend
}

#' Plot lines2test on top of banana_plot().
#'
#' @param lines2test list of lines. See get_lines2test().
#' @param day data.frame. First column is timestamp.
#' @param plotLineNames
plot_lines2test <- function(lines2test, day, plotLineNames=F, col=rgb(1,1,1, 0.2), lwd=2, ...) {
  invisible(lapply(lines2test, plot_line, day=day, col=col, lwd=lwd, ...))
  if (plotLineNames) text(x=x[sapply(lines2test, function(L) L$x0)], 
                          y=y[sapply(lines2test, function(L) L$y0)], 
                          labels=names(lines2test), cex=0.5)
  plot_line(lines2test$user, day, lwd=2, col = rgb(0,0,0))
}

#' Draw a line L on banana_plot.
#'
#' @param L list. Line parameterized with x0, x1, y0, y1.
#' @param x,y numeric vectors. x-axis and y-axis values.
#' @param ... additional plot parameters
plot_line <- function(L, day, lwd=1, ...) {
  x <- day[, 1]
  y <- log(name2num(colnames(day[-1])))
  segments(x[L$x0], y[L$y0], x[L$x1], y[L$y1], lwd=lwd, ...)
  points(x[L$x0], y[L$y0], pch=19, ...)
}

#' Plot a rectangular grid on top of banana_plot().
#'
#' @param day data.frame. First column is timestamp.
#' @param grid list(x,y) See make_banana_grid().
banana_plot_grid <- function(day, grid=make_banana_grid(dim(day))) {
  x <- day[, 1]
  y <- log(name2num(colnames(day[-1])))
  abline(v=x[grid$x], h=y[grid$y])
}

#' Plot banana_plot and plot subset of lines2test. Used for plotting significant lines.
#'
#' @param day data frame. See banana_plot().
#' @param idx_lines2test logical vector. Subsets lines2test.
#' @param lines2test list of lines. See get.lines2test.
#' @param grid list(x,y) See create_banana_grid().
banana_plot_lines <- function(day, idx_lines2test=T, lines2test, plot.title={title(ylab="Particle size (nm)")}, grid=NULL) {
  banana_plot(day, plot.title=plot.title)
  if (!is.null(grid)) banana_plot_grid(day, grid)
  plot_lines2test(lines2test[idx_lines2test], day)
}

#' Wrapper for banana_plot_lines() for use with output of make_banana_test_param().
#'
#' @param D output from make_banana_test_param()
#' @param idx_lines2test see banana_plot_lines()
#' @param grid see banana_plot_lines()
#' @param signif_lines logical of length = length(lines2test). Subset of significant lines to plot in brown.
#'
#' @return plot
banana_plot_lines_wrapper <- function(test_params, idx_lines2test=T, grid=test_params$grid, idx_signif_lines=F) {
  banana_plot_lines(test_params$data, idx_lines2test=idx_lines2test, lines2test = test_params$lines2test, grid=grid)
  plot_lines2test(test_params$lines2test[idx_signif_lines], test_params$data, col="brown", lwd=1)
}