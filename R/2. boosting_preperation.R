#' Make the first guess of GPD parameters
#'
#'
#' @param y A vector of observations
#' @details Estimation is done based on maximum likelihood
#' @return Returning a data frame with in each of the length(y) rows the estimated parameters.
#' @export
first_guess <- function(y){
  n = length(y)
  theta = POT::fitgpd(y,threshold=0)$fitted.values # estimate GPD parameters
  theta_init = data.frame(s = rep(theta[1],n), g= rep(theta[2],n)) # Copy the parameters for n rows
  return(theta_init)
}

#' Boosting DF
#'
#' Create an data frame with the observations, current value of the parameters and the derivitives of the power divergence.
#' Used to perform one step the the boosting procedure.
#'
#' @param data A data frame with the response and covariates
#' @param theta The parameter values for each observation
#' @param alpha The power for the divergence if 0 then maximum likelihood is used (default 0)
#'
#' @return A data frame with the observations, parameter values and derivatives of the power divergence
#' @export
get_boosting_df <- function(data,theta,alpha=0){
  if(nrow(theta) != nrow(data)) stop("Number of rows of data and theta are not equal.")

  divergence = compute_divergence(data$y,theta,alpha)
  boosting_df = cbind(data,divergence)

  return(boosting_df)
}

#' Divergence
#'
#' Computes the power divergence for the generalized pareto distribution together with the first and second derivatives.
#'
#' @param y Vector of observations
#' @param theta Named data frame of parameters for the generelized pareto distribution
#' @param alpha The power for the divergence if 0 then maximum likelihood is used (default 0)
#'
#' @return A data frame GPD parameters calculated power divergence and derivatives of power divergence
#' @export
compute_divergence <- function(y,theta,alpha){
  divergence_input = cbind(theta,y)
  if(alpha == 0){ # Compute everything for maximum likelihood
    divergence_input_ML = divergence_input
    dev = apply(divergence_input_ML,1,GP_dev)
    r_s = apply(divergence_input_ML,1,GP_dev_diff_s)
    r_g = apply(divergence_input_ML,1,GP_dev_diff_g)
    r2_s = apply(divergence_input_ML,1,GP_dev_diff2_s)
    r2_g = apply(divergence_input_ML,1,GP_dev_diff2_g)
  } else{ # Compute everything with power divergence
    A = apply(divergence_input,1,A_func)
    B = apply(divergence_input,1,B_func)
    divergence_input_PD = cbind(divergence_input,A,B)

    dev = apply(divergence_input,1,function(x) PD_dev(x[1:3],x[4],x[5],alpha))
    r_s = apply(divergence_input,1,function(x) PD_dev_diff_s(x[1:3],x[4],x[5],alpha))
    r_g = apply(divergence_input,1,function(x) PD_dev_diff_g(x[1:3],x[4],x[5],alpha))
    r2_s = apply(divergence_input,1,function(x) PD_dev_diff2_s(x[1:3],x[4],x[5],alpha))
    r2_g = apply(divergence_input,1,function(x) PD_dev_diff2_g(x[1:3],x[4],x[5],alpha))
  }
  output = data.frame(s=theta$s,g=theta$g,dev=dev,r_s=r_s, r_g=r_g, r2_s=r2_s, r2_g=r2_g)
  return(output)
}
