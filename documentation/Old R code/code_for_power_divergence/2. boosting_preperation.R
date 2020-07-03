#' Make the first guess of GPD parameters
#'
#'
#' @param y A vector of observations
#' @param gamma_positive boolean indicating whether gamma should be positive
#' @details Estimation is done based on unconditional maximum likelihood.
#' Where the density is reparameterized by beta = log(sigma).
#' @return Returning a data frame with in each of the length(y) rows the estimated parameters.
#' @export
first_guess <- function(y,gamma_positive){
  n = length(y)
  theta = POT::fitgpd(y,threshold=0)$fitted.values # estimate GPD parameters
  if(theta[2] < 0 & gamma_positive){
    warning("Unconditional estimate of gamma is negative (initial estimate is set to 0.01")
    theta[2] = 0.01
  }
  theta_init = data.frame(s = rep(theta[1],n), g= rep(theta[2],n))
  return(theta_init)
}


#' Transform parameters to constrain them to be positive or not
#'
#' @param theta data frame with parameters
#' @param gamma_positive boolean indicating whether gamma should be positive
#' @param inverse_transform boolean indicating wheter the tranformation should be inverted
#' @return tranformed parameters
#' @export
transform_parameters <- function(theta,gamma_positive,inverse_transform=F){
  if(inverse_transform){
    if(gamma_positive){
      theta_transformed = data.frame(st = log(theta$s),gt = log(theta$g))
    } else{
      theta_transformed = data.frame(st = log(theta$s),gt = theta$g)
    }
  } else{
    if(gamma_positive){
      theta_transformed = data.frame(s = exp(theta$st), g = exp(theta$g))
    } else{
      theta_transformed = data.frame(s = exp(theta$st), g = theta$gt)
    }
  }
}

#' Boosting DF
#'
#' Create an data frame with the observations, current value of the parameters and the derivitives of the power divergence.
#' Used to perform one step the the boosting procedure.
#'
#' @param data A data frame with the response and covariates
#' @param theta The parameter values for each observation
#' @param alpha The power for the divergence if 0 then maximum likelihood is used (default 0)
#' @param gamma_positive boolean indicating whether gamma should be positive
#'
#' @return A data frame with the observations, parameter values and derivatives of the power divergence
#' @export
get_boosting_df <- function(data,theta,alpha,gamma_positive){
  if(nrow(theta) != nrow(data)) stop("Number of rows of data and theta are not equal.")

  divergence = compute_divergence(data$y,theta,alpha,gamma_positive)

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
#' @param gamma_positive boolean indicating whether gamma should be positive
#'
#' @return A data frame GPD parameters calculated power divergence and derivatives of power divergence
#' @details the sigma parameter is reparematerized as exp(beta) hence the output are the derivatives in terms of beta and gamma
#' @export
compute_divergence <- function(y,theta,alpha,gamma_positive){
  divergence_input = cbind(transform_parameters(theta,gamma_positive,inverse_transform = F),y)
  if(alpha == 0){ # Compute everything for maximum likelihood
    divergence_input_ML = divergence_input
    dev = apply(divergence_input_ML,1,GP_dev)
    r_st = exp(theta$st)*apply(divergence_input_ML,1,GP_dev_diff_s)
    r2_st = exp(2*theta$st)*apply(divergence_input_ML,1,GP_dev_diff2_s) + r_st
    if(gamma_positive){
      r_gt = exp(theta$gt)*apply(divergence_input_ML,1,GP_dev_diff_g)
      r2_gt = exp(2*theta$gt)*apply(divergence_input_ML,1,GP_dev_diff2_g) + r_gt
    } else{
      r_gt = apply(divergence_input_ML,1,GP_dev_diff_g)
      r2_gt = apply(divergence_input_ML,1,GP_dev_diff2_g)
    }
  } else{ # Compute everything with power divergence
    A = apply(divergence_input,1,A_func)
    B = apply(divergence_input,1,B_func)
    divergence_input_PD = cbind(divergence_input,A,B)

    dev = apply(divergence_input_PD,1,function(x) PD_dev(x[1:3],x[4],x[5],alpha))
    r_st = exp(theta$st)*apply(divergence_input_PD,1,function(x) PD_dev_diff_s(x[1:3],x[4],x[5],alpha))
    r2_st = exp(2*theta$st)*apply(divergence_input_PD,1,function(x) PD_dev_diff2_s(x[1:3],x[4],x[5],alpha)) + r_st
    if(gamma_positve){
      r_gt = exp(theta$gt)*apply(divergence_input_PD,1,function(x) PD_dev_diff_g(x[1:3],x[4],x[5],alpha))
      r2_gt = exp(2*theta$gt)*apply(divergence_input_PD,1,function(x) PD_dev_diff2_g(x[1:3],x[4],x[5],alpha)) + r_gt
    } else{
      r_gt = apply(divergence_input_PD,1,function(x) PD_dev_diff_g(x[1:3],x[4],x[5],alpha))
      r2_gt = apply(divergence_input_PD,1,function(x) PD_dev_diff2_g(x[1:3],x[4],x[5],alpha))
    }
  }
  output = data.frame(st=theta$st,gt=theta$gt,dev=dev,r_st=r_st, r_gt=r_gt, r2_st=r2_st, r2_gt=r2_gt)
  return(output)
}
