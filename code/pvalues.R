#' Empirical p-value (one-sided)
#'
#' @param test: output from tester. test$t0: observed test statistic. test$t: vector of replicate test statistics.
#'
#' @return p-value
pvalue <- function(test, na.rm=T) {
  tt <- test$t
  t0 <- test$t0
  if (na.rm) tt <- tt[!is.na(tt)]
  n <- length(tt)
  r <- sum(tt >= t0)
  (r + 1) / (n + 1)
}

#' Empirical p-value (one sided) for multiple test statistics.
#'
#' @param test tester2 output. 
#' test$t0: vector of m test statistics on the observed data. 
#' test$t: nXm matrix of m test statistics on n replicate data.
#'
#' @return vector of m p-values (one for each test statistic)
pvalue2 <- function(test, na.rm=T, na2minus_inf=T) {
  t0 <- test$t0
  if (na2minus_inf) t0[is.na(t0)] <- -Inf # If NA then make it insignificant.
  tt <- t(test$t)
  n <- ncol(tt)  
  r <- rowSums(tt >= t0, na.rm = na.rm)
  (r + 1) / (n + 1)
} 

#' Bonferroni adjusted p-value
#'
#' @param test 
#' @param na.rm 
#'
#' @return
pvalue_bonf <- function(test, na.rm=F) {
  p.adjust(pvalue2(test, na.rm = na.rm), method = "bonf")
}

#' Iteration adjustment for p-values
#'
#' @param pval  p-value
#' @param t     iteration step
#'
pval_iter_adj <- function(pval, t) {
  min(1, 2^t * pval)
}

#' minP adjusted p-value
#'
#' @param test 
#' @param na2minus_inf 
#'
#' @return
pvalue_minP <- function(test, na2minus_inf=F) {
  if (na2minus_inf) {
    test$t[is.na(test$t)] <- -Inf
    test$t0[is.na(test$t0)] <- -Inf
  }
  minP(test$t0, t(test$t))
}

#' Compute maxT adjusted FWER-controlled p-values
#' Author: Kai Puolamäki
#' 
#' @param x vector of length n, test statistic in the original data
#' @param y matrix nXm, columns are the test statistic in m samples from null hypothesis
#' @return Vector of length n, maxT adjusted mid-P p-values for the test statistic
#' 
#' @export
maxT <- function(x,y) {
  n <- length(x)
  if(!is.matrix(y)) stop("maxT: !is.matrix(y)")
  if(n!=dim(y)[1]) stop("maxT: n!=dim(y)[1]")
  m <- dim(y)[2]
  
  i <- order(x)
  j <- rev(i)
  xm <- matrix(rep(x,m),n,m)
  ym <- if(n>1) apply(y[i,,drop=FALSE],2,cummax)[order(i),,drop=FALSE] else y
  p <- (1+rowSums(xm<=ym))/(1+m)
  names(p) <- names(x)
  cummax(p[j])[order(j)]
}

#' Compute minP adjusted FWER-controlled p-values
#' Author: Kai Puolamäki
#' 
#' @param x vector of length n, test statistic in the original data
#' @param y matrix nXm, columns are the test statistic in m samples from null hypothesis
#' @return Vector of length n, maxT adjusted mid-P p-values for the test statistic
#' 
#' @export
minP <- function(x,y) {
  n <- length(x)
  if(!is.matrix(y)) stop("maxT: !is.matrix(y)")
  if(n!=dim(y)[1]) stop("maxT: n!=dim(y)[1]")
  m <- dim(y)[2]
  
  z <- matrix(c(x,y),n,m+1)
  p <- t(apply(z,1,topvalues))
  rownames(p) <- names(x)
  maxT(-p[,1],-p[,-1,drop=FALSE])
}

#' Convert vector of test statistic to the vector of p-values
#' Author: Kai Puolamäki
#' 
#' @param x vector of length n, test statistic
#' @return Vector of length n, p-values.
#' 
#' @export
topvalues <- function(x) {
  n <- length(x)
  idx <- is.na(x)
  if(any(idx)) {
    p <- rep(NA,n)
    p[!idx] <- topvalues(x[!idx])
  } else {
    i <- order(x)
    p <- (n+1-order(i))/n
    if(n>1) {
      ## If there are equally sized test statistic...
      for(j in 2:n) {
        if(x[i[j]]==x[i[j-1]]) p[i[j]] <- p[i[j-1]]
      }
    }
  }
  p
}

#' Create data frame with test statistics and pvalues (raw, bonferroni, minP).
#'
#' @param test list with test$t, test$t0 (output from tester)
#' @param test_points vector length n. points at which test statistics are computed.
#' @param names_point_value character vector length 2. colnames for points and values.
#' @param seed integer. random seed (for minP) 
#' @param include_t0 logical. Include test statistic values.
#'
#' @return data frame with 3-5 columns: test statistic points and values, raw pvalue, bonferroni pvalue, minP pvalue
pvalue_table <- function(test_points=NA, 
                         test, 
                         seed=42, 
                         name_points = "time", 
                         name_t0 = "test_statistic", 
                         include_t0=FALSE) {
  set.seed(seed)
  df <- data.frame(test_points,
                   round(test$t0, 2),
                   pval_unadjusted = round(pvalue2(test, na.rm=T), 2), 
                   pval_bonferroni = round(pvalue_bonf(test, na.rm=T), 2), 
                   pval_minP = round(pvalue_minP(test, na2minus_inf = T), 2))
  colnames(df) <- c(name_points, name_t0, "$p_{raw}$", "$p_{bonf}$", "$p_{minP}$")
  if (!include_t0) df <- df[-2]
  if (all(is.na(test_points))) df <- df[-1]
  df
}