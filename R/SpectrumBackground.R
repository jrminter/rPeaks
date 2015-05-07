#' Compute the spectrum background
#'
#' The method is based on Sensitive Nonlinear Iterative Peak (SNIP)
#' clipping algorithm.
#'
#' References:
#'
#' C. G Ryan et al.: SNIP, a statistics-sensitive background treatment
#' for the quantitative analysis of PIXE spectra in geoscience
#' applications. NIM, B34 (1988), 396-402.
#'
#' M. Morhac, J. Kliman, V. Matoucek, M. Veselsky, I. Turzo.:
#' Background elimination methods for multidimensional gamma-ray
#' spectra. NIM, A401 (1997) 113-132.
#'
#' D. D. Burgess, R. J. Tervo: Background estimation for gamma-ray
#' spectroscopy. NIM 214 (1983), 431-434.
#'
#'
#' @param y The vector of source spectrum
#' @param iterations Maximal width of clipping window
#' @param decreasing The direction of change of clipping window.
#' If \code{TRUE} the window is decreasing, otherwise the window is
#' increasing.
#' @param order The order of clipping filter
#' @param smoothing Logical variable whether the smoothing operation
#' in the estimation of background will be included.
#' @param window Width of smoothing window
#' @param compton Logical variable whether the estimation of Compton
#' edge (step-like feature at the peaks positions) will be included.
#'
#' @return The background
#'
#' @export
#'
#' @useDynLib rPeaks R_SpectrumBackground
#'
#' @examples
#' # Not run
#'
SpectrumBackground <- function(y,
              iterations=100,
              decreasing=FALSE,
              order=c("2","4","6","8"),
              smoothing=FALSE,
              window=c("3","5","7","9","11","13","15"),
              compton=FALSE){

  p <- .Call("R_SpectrumBackground",
             as.vector(y),
             as.integer(iterations),
             as.integer(decreasing),
             as.integer(as.integer(match.arg(order))/2-1),
             as.integer(smoothing),
             as.integer(as.integer(match.arg(window))),
             as.integer(compton))
  return(p)
}
