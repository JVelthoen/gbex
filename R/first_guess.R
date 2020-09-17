#' Make the first guess of GPD parameters
#'
#' Estimate parameters of a GPD distribution on a set of observations
#'
#' @param y numeric vector of positive observations
#' @param gamma_positive boolean indicating whether gamma should be positive for boosting procedure
#' @details The maximum likelihood estimator is used for the to obtain the parameters.
#' When gamma < 0 and gamma_positive= True gamma is set to 0.01 and a warning is issued.
#' @return Returning a data.frame with length(y) identical rows with the estimated parameters for sigma and gamma.
#' @export
first_guess <- function(y,gamma_positive){
  n = length(y)
  theta = POT::fitgpd(y,threshold=0)$fitted.values
  if(theta[2] < 0 & gamma_positive){
    warning("Unconditional estimate of gamma is negative (initial estimate is set to 0.01")
    theta[2] = 0.01
  }
  theta_init = data.frame(s = rep(theta[1],n), g= rep(theta[2],n))
  return(theta_init)
}
