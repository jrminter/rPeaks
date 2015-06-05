#' PeakEstimateSigma
#'
#' Given a background-subtracted, isolated peak with frequency \code{y}
#' and scale \code{x}, and an estimate of the the mean, \code{mu},
#' calculate an estimate of the standard deviation. This is useful for
#' parameter estimates for fitting to a Gaussian.
#'
#' @param x The scale axis (vector).
#' @param y The frequency axis (vector).
#' @param mu An estimate of the mean (numeric).
#'
#' @return sdEst An estimate of the standard deviation.
#' @export
#'
#' @examples
#' # not run
PeakEstimateSigma <- function(x,y, mu){
  l <- length(x)
  sy  <- 0.
  sd2 <- 0.
  for (i in 1:l){
    sy <- sy + y[i]
    sd2 <- sd2 +  y[i] * (x[i]-mu)^2
  }
  varEst <- sd2/ (sy - 1.)
  sdEst = sqrt(varEst)
  sdEst
}
