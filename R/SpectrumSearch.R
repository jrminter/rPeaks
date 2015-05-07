#' Automatically detect peaks on a continuous background with noise.
#'
#' @param y Numeric vector of source spectrum
#' @param sigma Sigma of searched peaks
#' @param threshold Threshold value in \% for selected peaks, peaks with amplitude less than \code{threshold*highest_peak/100} are ignored
#' @param background Remove background. Logical variable, set to \code{TRUE} if the removal of background before deconvolution is desired.
#' @param iterations Number of iterations in deconvolution operation.
#' @param markov Logical variable, if it is \code{TRUE}, first the source spectrum is replaced by new spectrum calculated using Markov chains method.
#' @param window Averanging window of searched peaks, applies only for Markov smoothing
#'
#' Algorithm is straightforward. The function removes background and smooths (if requested) source vector \code{y}, then deconvolves it using Gaussian with \code{sigma} as responce vector and after that searches for peaks in deconvoluted vector which are above \code{threshold}.
#'
#' @return List with two vectors: \code{y} Deconvoluted source vector and \code{pos} Indexes of found peaks in spectrum
#'
#' @export
#'
#' @useDynLib rPeaks R_SpectrumSearchHighRes
#'
#' @examples
#' # Not run
SpectrumSearch <-  function(y,
                            sigma=3.0,
                            threshold=10.0,
                            background=FALSE,
                            iterations=13,
                            markov=FALSE,
                            window=3){
  p <- .Call("R_SpectrumSearchHighRes",
             as.vector(y),
             as.numeric(sigma),
             as.numeric(threshold),
             as.integer(background),
             as.integer(iterations),
             as.integer(markov),
             as.integer(window) )
  return(p)
}
