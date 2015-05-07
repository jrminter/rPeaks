#' Supress noise with a discrete Markov chain.
#'
#' @param y Numeric vector of source spectrum
#' @param window Width of averaging smoothing window
#'
#' @return p The smoothed spectrum
#'
#' @export
#'
#' @useDynLib rPeaks R_SpectrumSmoothMarkov
#'
#' @examples
#' # Not run
SpectrumSmoothMarkov <- function(y,window=3){
  p <- .Call("R_SpectrumSmoothMarkov",
             as.vector(y),
             as.integer(window))
  return(p)
}
