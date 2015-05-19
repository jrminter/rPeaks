#' FitpeakToGaussian
#'
#' Using the output from EstimateGaussianParameters, do a nonlinear
#' fit to a Gaussian model
#'
#' @param lData - the output (a list) from EstimateGaussianParameters
#' @param doPlot - A boolean (default FALSE) to plot the data, estimate, and the fit.
#' @param pTitle - An optional plot title
#' @param xTitle - An optional X-axis label
#' @param yTitle  - An optional Y-axis label
#'
#' @return sum - A summary of the fit
#' @export
#'
#' @examples
#' # Not run
FitPeakToGaussian <- function(lData, doPlot=FALSE, pTitle='Peak',xTitle='x', yTitle='y'){
  # extract the data we need
  dat      <- lData$peakData
  muEst    <- lData$muEst
  sigmaEst <- lData$sigmaEst
  htEst    <- lData$htEst
  x <- dat$x
  y <- dat$y

  (res <- nls(y ~ scale*exp(-0.5*(x-mu)^2/sigma^2),
              start=c(mu=muEst, sigma=sigmaEst, scale=htEst)))

  sum <- summary(res)
  ht <- sum$coefficients[3]
  mu <- sum$coefficients[1]
  sigma <- sum$coefficients[2]

  if(doPlot==TRUE){
    yc <- htEst*exp(-0.5*(x-muEst)^2/sigmaEst^2)
    yc2 <- ht*exp(-0.5*(x-mu)^2/sigma^2)
    xMin <- min(x)
    xMax <- max(x)
    yMax <- max(y)
    plot(c(xMin, xMax), c(0, yMax), type='n',
         xlab=xTitle, ylab=yTitle, main=pTitle)
         points(x,y,pch=19)
         lines(x,yc, col='red', lw=2)
         lines(x,yc2, col='blue', lw=2)
         legend("topright", bty="n",
         legend=c("points", "estimate", "fit"),
         col=c('black', 'red', 'blue'),
         lty=c(0, 1, 1),
         lwd=c(0, 1,1),
         pch=c(19, NA, NA))
  }
  return(sum)
  }
