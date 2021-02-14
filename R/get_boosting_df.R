#' Boosting DF
#'
#' Given the data current parameter estimates update compute the likelihood and derivatives of the parameters
#' for the next tree to be fitted in the boosting procedure of gbex.
#'
#' @param data A data frame with the response and covariates
#' @param theta The parameter values for each observation
#' @param gamma_positive boolean indicating whether gamma should be positive
#'
#' @return A data frame with the data parameter values and derivatives of the GPD likelihood
#' @export
get_boosting_df <- function(data,theta,gamma_positive){
  deviance_input = cbind(transform_parameters(theta,gamma_positive,inverse_transform = F),data$y)

  dev = apply(deviance_input,1,GP_dev)
  r_st = exp(theta$st)*apply(deviance_input,1,GP_dev_diff_s)
  r2_st = exp(2*theta$st)*apply(deviance_input,1,GP_dev_diff2_s) + r_st

  if(gamma_positive){
    r_gt = exp(theta$gt)*apply(deviance_input,1,GP_dev_diff_g)
    r2_gt = exp(2*theta$gt)*apply(deviance_input,1,GP_dev_diff2_g) + r_gt
  } else{
    r_gt = apply(deviance_input,1,GP_dev_diff_g)
    r2_gt = apply(deviance_input,1,GP_dev_diff2_g)
  }

  update = data.frame(st=theta$st,gt=theta$gt,dev=dev,r_st=r_st, r_gt=r_gt, r2_st=r2_st, r2_gt=r2_gt)

  boosting_df = cbind(data,update)

  return(boosting_df)
}
