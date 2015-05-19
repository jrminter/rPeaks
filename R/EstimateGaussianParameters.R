#' Estimate the Gaussian parameters for a peak
#'
#' Given the vectors, x and y, describing a peak and the centroid,
#' estimate the parameters for a Gaussian for input to a nonlinear
#' fit.
#'
#' @param x The vector containing the values for the independent axis
#' @param y The vector containing the values for the dependent axis.
#' @param centroid The centroid for the peak. Typically determined using \code{SpectrumSearch}
#' @param cutoff The fraction of intensity (default 0.25) from the maximum to consider.
#' @param bDebug Default False - print debugging messages
#'
#' @return A list with muEst, sigmaEst, htEst, lCO (lower cut off), uC (upper cutoff) and peakData - a data frame (x,y).
#'
#' @export
#'
#'
#' @examples
#' # Not run
EstimateGaussianParameters <- function(x, y, centroid, cutoff=0.25, bDebug=FALSE){
  l <- length(x)
  i <- which(x==centroid)
  if(bDebug){
    print(l)
    print(i)
  }
  bOK <- TRUE
  ymax <- y[i]
  if(bDebug){
    print(ymax)
  }

  j <- 0
  bOK <- TRUE
  while(bOK==TRUE){
    j <- j + 1
    it <- i - j
    yt <- y[it]
    if(yt >= 0.5*ymax){
      iL50 <- it
    }
    if(yt >= .25*ymax){
      iL25 <- it
    } else {
      bOK <- FALSE
    }
  }
  if(bDebug){
    print(c(iL25, iL50))
  }

  j <- 0
  bOK <- TRUE
  while(bOK==TRUE){
    j <- j + 1
    it <- i + j
    yt <- y[it]
    if(yt >= 0.5*ymax){
      iU50 <- it
    }
    if(yt >= .25*ymax){
      iU25 <- it
    } else {
      bOK <- FALSE
    }
  }
  if(bDebug){
    print(c(iU25, iU50))
  }

  fwhm <- x[iU50] - x[iL50]
  # from http://mathworld.wolfram.com/GaussianFunction.html
  sigmaEst <- fwhm/2.3548
  muEst <- centroid
  if(bDebug){
    print(round(c(muEst, sigmaEst), 3))
  }

  # make a data frame with a cropped peak
  peakData <- data.frame(x=x, y=y)
  peakData <- peakData[iL25:iU25,]

  # create a list to return
  ret <- list(muEst=muEst,
              sigmaEst=sigmaEst,
              htEst=y[i],
              lCO=as.integer(iL25),
              uCO=as.integer(iU25),
              peakData=peakData)
}
