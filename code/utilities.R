check_if_packages_are_present <- function(required_packages) {
  missing_packages  <- setdiff(required_packages, installed.packages())
  
  if (length(missing_packages) > 0) {
    cat("Error! Cannot continue due to missing packages\n")
    cat("Please install the following packages:\n")
    for (i in missing_packages)
      cat(paste0("\t- ", i, "\n"))
    stop("Install packages and try again.")
  }
}

#' Read and preprocess air quality data
#'
#' @return data.frame
readAQ <- function() {
  AQ <- read.csv("data/AirQualityUCI.csv",
                 header=TRUE,sep=";",dec=",")
  AQ <- AQ[!is.na(AQ[,"CO.GT."]),]
  TI <- mapply(function(x,y) as.POSIXct(strptime(sprintf("%s %s",x,y),
                                                 format="%d/%m/%Y %H.%M.%S")),
               AQ[,"Date"],AQ[,"Time"])
  for(i in 3:15) AQ[AQ[,i]==-200,i] <- NA
  AQ <- AQ[!is.na(TI),]
  TI <- TI[!is.na(TI)]
  data.frame(time=mapply(function(x,y)
    as.POSIXct(strptime(sprintf("%s %s",x,y),format="%d/%m/%Y %H.%M.%S"),origin="1970-01-01"),
    AQ[,"Date"],AQ[,"Time"]),
    AQ[,3:15])
}

#' Save a figure as PDF
#'
#' @param fig function that outputs plot
#' @param filename 
#' @param width 
#' @param height 
savePDF <- function(fig, filename, width, height) {
  pdf(file = filename, width = width, height = height)
  fig
  invisible(dev.off())
}

#' Print formatted latex table
#'
#' @param df data frame
#' @param label 
#' @param caption 
#' @param include.colnames 
#'
#' @return latex table
print_xtable_formatted <- function(df, label=NULL, caption=NULL, 
                                   include.colnames=FALSE, 
                                   include.rownames=TRUE) {
  print(xtable::xtable(t(as.matrix(df)), 
                       label=label, 
                       caption=caption, align=c("l", rep("r", times=nrow(df)))), 
        include.colnames=include.colnames, 
        include.rownames=include.rownames,
        sanitize.text.function=function(x){x}, 
        booktabs=TRUE)
}

#' Returns rownames of data points inside a polygon selection. Utilizes sp library. 
#'
#' @param data      nX2 matrix with xy coordinates
#' @param selection kX2 matrix of xy coordinates for polygon vertices 
#'
#' @return
points_in_selection <- function(data, selection) {
  data <- data.frame(data)
  selection <- data.frame(selection)
  idx <- sp::point.in.polygon(point.x = data[ , 1], point.y = data[ , 2], 
                              pol.x = selection[ , 1], pol.y = selection[ , 2])
  if (length(idx) < 2) warning(length(idx), " points in selection.")
  as.integer(rownames(data)[as.logical(idx)])
}

#' Compute polygon of convex hull of points. 
#'
#' @param xy  matrix nX2
#' @param R   vector of m indices (for a subset of xy)
#'
#' @return    matrix kX2 with k points on convex hull of xy[R,]
convhull_subset <- function(xy, R) {
  xyR <- xy[R,]
  xyR[chull(xyR),]
}

#' Convert column names in SMEAR data into numeric, e.g. HYY_DMPS.d282e1 into 2.82
#'
#' @param x colnames(banana.data)
#'
#' @return numeric vector
name2num <- function(x) {as.numeric(sapply(x, function(x) strsplit(x, split = ".d", fixed=T)[[1]][2]))/1000}

char2posixct <- function(s) {as.POSIXct(s, tz = "UTC")}

fdate <- function(d) format(d, "%Y-%m-%d")

#' Find last TRUE in logical vector and change it into FALSE.
#'
#' @param i logical vector
change_last_true_into_false <- function(i) {
  idx <- which(i)
  i[idx[length(idx)]] <- FALSE
  i
}

#' Create a grid for banana_plot().
#'
#' @param dims numeric length 2. Dimensions of grid e.g. dim(z).
#' @param spacing numeric c(a,b). Keep every a:th point in x and every b:th point in y. e.g. c(3,1) means that the x-axis is sparser (every 3rd point is kept).
make_banana_grid <- function(dims, spacing=c(1,1)) {
  list(x=rep(c(T, rep(F, spacing[1] - 1)), length.out=dims[1]), 
       y=rep(c(T, rep(F, spacing[2] - 1)), length.out=dims[2]))
}

#' Convert data frame into matrix, add 1 and log.
#'
#' @param X data frame 
#'
#' @return matrix
df2matrix_plus1_log <- function(X) {
  X[X<0] <- NA # assume that negative particle concentrations are a mistake. NAs are plotted as white.
  as.matrix(log(X + 1))
}

#' Convert linenames into index inside lines2test 
#' e.g. "R38" becomes 39 because user.line is first.
#'
#' @param Lines line generated from get_lines2test. Has name attribute of the form "R12"
#'
#' @return int
linename2idx <- function(Lines) {
  tmp <- 1 + as.numeric(unlist(strsplit(names(Lines), split="R")))
  tmp[!is.na(tmp)]
}

#' Get points of line segment L, rounded to a discrete grid.
#' This is a simplified "vector-to-raster" conversion (geoinformatics).
#'
#' @param L list with x0, y0, x1, y1
get_line_points_between_start_end <- function(L) {
  x0<-L$x0
  x1<-L$x1
  y0<-L$y0
  y1<-L$y1
  xL <- x0:x1
  yL <- (xL - x0) * (y1 - y0) / (x1 - x0) + y0 
  if (x0==x1) yL <- y0:y1
  list(x=xL, y=round(yL))
}

#' Generate list of lines.
#'
#' @param user_line 
#' @param n_lines 
#' @param sample_lines function that returns list(x0,y0,x1,y1)
#'
#' @return list of lines
get_lines2test <- function(user_line, n_lines, sample_lines) {
  c(user = list(user_line), 
    R = replicate(n=n_lines-1, sample_lines(), simplify = F))
}
