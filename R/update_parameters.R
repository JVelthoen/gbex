#' Update parameters
#'
#' Function is used in the boosting procedure of gbex to update the parameters after a set of gradient trees is fitted.
#'
#' @param tree_sigma A gradient tree object fitted for sigma parameter
#' @param tree_gamma A gradient tree object fitted for gamma parameter
#' @param boosting_df the boosting_df
#' @param lambda a vector with the learning rate of sigma and gamma
#' @return A data frame with the udpated parameters sigma and gamma
#' @export
update_parameters <- function(tree_sigma,tree_gamma,boosting_df,lambda){
  update_sigma = -lambda[1]*predict(tree_sigma,boosting_df)
  update_gamma = -lambda[2]*predict(tree_gamma,boosting_df)

  theta <- data.frame(st = boosting_df$st + update_sigma, gt = boosting_df$gt + update_gamma)
  return(theta)
}
