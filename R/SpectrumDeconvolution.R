#' Deconvolute spectrum
#'
#' @param y Numeric vector of source spectrum
#' @param response Vector of response spectrum. Its length shold be less or equal the length of \code{y}
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
