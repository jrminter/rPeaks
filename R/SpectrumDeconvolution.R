#' Deconvolute spectrum
#'
#' This function is used to strip-off a known instrumental function from
#' a source spectrum. It is achieved by  deconvolution of source
#' spectrum based upon a response spectrum using the Gold or
#' Richardson-Lucy algorithms. Both methods provides less oscillating
#' solutions than Fourier or VanCittert algorithms.
#'
#' Both methods search iteratively for solution of deconvolution
#' problem
#'
#' \deqn{y(i)=\sum_{j=1}^{n}h(i-j)x(j)+e(i)}
#'
#'   in the form
#'
#' \deqn{x^{(k)}(i)=M^{(k)}(i)x^{(k-1)}(i)}
#'
#' For Gold method:
#'
#' \deqn{M^{(k)}(i)=\frac{x^{(k-1)}(i)}{\sum_{j=1}^{n}h(i-j)x^{(k-1)}(j)}}
#'
#' For Richardson-Lucy:
#'
#' \deqn{M^{(k)}(i)=\sum_{l=0}^{n}h(i-l)\frac{x^{(k-1)}(l)}{\sum_{j=1}^{n}h(l-j) x^{(k-1)}(j)}}
#'
#' Boosting is the exponentiation of iterated value with boosting
#' coefficient/exponent. It is generally improve stability.
#'
#' References:
#'
#' Abreu M.C. et al., A four-dimensional deconvolution method to correct
#' NA38 experimental data, NIM A 405 (1998) 139.
#'
#' Lucy L.B., A.J. 79 (1974) 745.
#'
#' Richardson W.H., J. Opt. Soc. Am. 62 (1972) 55.
#'
#' Gold R., ANL-6984, Argonne National Laboratories, Argonne Ill, 1964.
#'
#' Coote G.E., Iterative smoothing and deconvolution of one- and
#' two-dimensional elemental distribution data, NIM B 130 (1997) 118.
#'
#' M. Morhac, J. Kliman, V. Matousek, M. Veselsky, I. Turzo.:
#' Efficient one- and two-dimensional Gold deconvolution and its
#' application to gamma-ray spectra decomposition. NIM, A401 (1997)
#' 385-408.
#'
#' Morhac M., Matousek V., Kliman J., Efficient algorithm of
#' multidimensional deconvolution and its application to nuclear data
#' processing, Digital Signal Processing 13 (2003) 144.
#'
#'
#'
#' @param y Numeric vector of source spectrum
#' @param response Vector of response spectrum. Its length should be less or equal the length of \code{y}
#' @param iterations Number of iterations (parameter L in the Gold deconvolution algorithm) between boosting operations
#' @param repetitions Number of repetitions of boosting operations. It must be greater or equal to one. So the total number of iterations is \code{repetitions*iterations}
#' @param boost Boosting coefficient/exponent. Applies only if \code{repetitions} is greater than one. Recommended range [1..2].
#' @param method Method selected for deconvolution. Either Gold or Richardson-Lucy.
#'
#' @return p The deconvoluted spectrum
#'
#' @export
#'
#' @useDynLib rPeaks R_SpectrumDeconvolution R_SpectrumDeconvolutionRL
#'
#' @examples
#' # not run
SpectrumDeconvolution <- function(y,response,iterations=10,repetitions=1,boost=1.0,method=c("Gold","RL")){
  method <- match.arg(method)
  if (length(as.vector(response))<length(as.vector(y))){
    response <- c(response,rep(0,length(y)-length(response)))
  }
  if (length(as.vector(response))>length(as.vector(y))){
    stop("response length should be shorter or equal y length")
  }
  switch(method,
         Gold={
           p1 <- "R_SpectrumDeconvolution"
         },
         RL={
           p1 <- "R_SpectrumDeconvolutionRL"
         })
  p <- .Call(p1,
             as.vector(y),
             as.vector(response),
             as.integer(iterations),
             as.integer(repetitions),
             as.numeric(boost))

  return(p)
}
