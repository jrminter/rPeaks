#' FitSingleLogNormal
#'
#' Given an isolated, background-subtracted peak with the independent
#' variable axis, \code{vecX}, and the dependent variable axis,
#' \code{vecY}, fit the data to a single mode lognormal distribution
#' and output the a list with the coefficients and a data frame with a
#' smooth curve with \code{n.pts.out}, suitable for plotting.
#'
#' @param vecX A vector containing the independent variable.
#' @param vecY A vector containing the dependent variable.
#' @param n.pts.out The number of points in the output curve (default 100),
#' @param b.debug Default FALSE, a flag to print debugging intormation.
#'
#' @return A list containing a vector of coefficients and a data frame with a smooth curve of vectors \code{xc} and for \code{yc} plotting.
#' @export
#'
#' @examples
#' # not run
FitSingleLogNormal <- function(vecX, vecY, n.pts.out=100, b.debug=TRUE){
  x <- sapply(vecX, log)
  amp.est <- max(vecX)
  y <- vecY
  n.pts <- length(x)
  mu.est <- PeakEstimateMu(x, y)
  # need to underestimate sigma
  sd.est <- 0.5*PeakEstimateSigma(x, y, mu.est)

  est <- c(mu.est, sd.est)
  names(est) <- c("mu.est", "sd.est")

  print(est)

  # fit a single mode lognormal
  fit <- nls( y~a*exp(-(x-c)^2/(2*b^2)),
              start=list( a=amp.est,
                          b=sd.est,
                          c=mu.est),
              control=nls.control(tol=1E-5,
                                  minFactor=1/1024),
              trace=b.debug)

  s <- summary(fit)
  print(s$coef)

  a <- s$coef[1]
  b <- s$coef[2]
  c <- s$coef[3]
  print (c(a,b,c))

  coef = c(a, exp(c), exp(b))
  names(coef) <- c("amplitude", "gmd", "gsd")
  myFunc <- function(x){
    val=a*exp(-(x-c)^2/(2*b^2))
    val
  }

  meanX <- mean(x)

  deltaX <- 1.5*max((max(x)- meanX), (meanX-min(x)))

  lxc <- seq(from=(meanX-deltaX), to=(meanX+deltaX), length.out=n.pts.out)
  xc  <- sapply(lxc, exp)
  yc  <- sapply(lxc,myFunc)

  smCurve <- data.frame(xc=xc,yc=yc)

  ret <- list(coef, smCurve)
  ret
}
