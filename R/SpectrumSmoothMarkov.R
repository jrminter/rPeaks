#' Supress noise with a discrete Markov chain.
#'
#' This function calculates smoothed spectrum from source spectrum
#' based on Markov chain method.
#'
#' The algorithm is based on discrete Markov chain, which has very
#' simple invariant distribution:
#'
#' \deqn{U_2=\frac{p_{1,2}}{p_{2,1}}U_1}
#' \deqn{U_3=\frac{p_{2,3}}{p_{3,2}}U_2 U_1 }
#' \deqn{\ldots}
#' \deqn{U_n=\frac{p_{n-1,n}}{p_{n,n-1}}U_{n-1} \ldots U_2 U_1}
#'
#' and \eqn{U_1} being defined from the normalization condition:
#'
#' \deqn{\sum_{i=1}^{n}U_i=1}
#'
#' \eqn{n} is the length of the smoothed spectrum.
#'
#' The probability of the change of the peak position from channel
#' \eqn{i} to the channel \eqn{i+1} is :
#'
#' \deqn{p_{i,i \pm 1}=A_i \sum_{k=1}^{m}exp \left( \frac{y(i \pm k)
#' -y(i)}{y(i \pm k)+y(i)}\right)}
#'
#' where \eqn{A_i} is the normalization constant so that:
#'
#' \deqn{p_{i,i-1}+p_{i,i+1}=1}
#'
#' and \eqn{m} is a width of smoothing window.
#'
#' References:
#'
#' Z.K. Silagadze, A new algorithm for automatic photopeak searches.
#' NIM A 376 (1996), 451.
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
