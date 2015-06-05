#' PeakEstimateMu
#'
#' Given a background-subtracted, isolated peak with frequency \code{y}
#' and scale \code{x}, estimate the mean. This is useful for fitting to
#' a Gaussian.
#'
#' @param x The scale axis.
#' @param y The frequency axis.
#'
#' @return muEst An estimate of the mean.
#'
#' @export
#'
#' @examples
#' # not run
PeakEstimateMu <- function(x,y){
  l <- length(x)
  sy  <- 0.
  sxy <- 0.
  for (i in 1:l){
    sy <- sy + y[i]
    sxy <- sxy +  x[i]*y[i]
  }
  muEst <- sxy /sy
  muEst
}
