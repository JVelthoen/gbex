#' GPD boosting
#'
#' gbex is used to fit the parameters of a generalized pareto distribution as a function of covariates to data.
#' In a boosting procedure two sequences of trees are fitted, one for the sigma parameter (scale) and one for the gamma parameter (shape).
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param lambda numeric, learning rate for the sigma and gamma parameter
#' @param lambda_ratio numeric, the ratio of the lambda_sigma/lambda_gamma (Use this together with lambda_scale, see details)
#' @param lambda_scale numeric learning rate of lambda_sigma (Use this together with lambda_ratio, see details)
#' @param B numeric Number of gradient boosting steps
#' @param depth numeric, the maximum depth for the sigma and gamma trees
#' @param min_leaf_size numeric, of length two indicatin the minimum of observations in a leaf node for sigma and gamma trees.
#' @param sf numeric, sample fraction between 0 and 1
#' @param initial_values numeric, the initial values of the sigma and gamma parameter. If NULL the parameters are estimated by fitting the GPD to Y.
#' @param pred_sigma charcter vector with names of the predictors to use for sigma (default NULL all covariates will be used)
#' @param pred_gamma charcter vector with names of the predictors to use for gamma (default NULL all covariates will be used)
#' @param silent boolean if true progress is given during the fitting procedure.
#'
#' @return gbex returns an object of class "gbex" which contains the following components:
#' \item{theta}{Data frame with the estimated gamma and sigma parameter for each observation}
#' \item{dev}{Numeric with deviance of model}
#' \item{theta_init}{numeric with inital parameters for sigma and gamma}
#' \item{gamma_positive}{boolean indicating if the model is fit restricting the gamma parameter to be positive}
#' \item{trees_sigma}{List with gradient_tree objects for sigma}
#' \item{trees_gamma}{List with gradient_tree objects for gamma}
#' \item{lambda}{Numeric with the learining rate used for sigma and gamma}
#' \item{B}{Numeric with number of trees}
#' \item{depth}{Numeric with maximum tree depth for sigma and gamma}
#' \item{sf}{Numeric subsample fraction used in fitting}
#' \item{min_leaf_size}{numeric minimum number of observations per leaf used in fitting}
#' \item{data}{The data used in fitting}
#' @details
#' The learning rates of the two two sequence of trees are interconnected, therefore it is more natural to only define the learning rate for the scale, and
#' indirectly specifying the other lambda by specifying lambda_ratio = lambda_sigma/lambda_gamma.
#' The tuning parameters can be chosen by means of cross validation via CV_gbex.
#'
#' @export
gbex <- function(y,X,B=100,lambda=NULL,
                 lambda_ratio = 7, lambda_scale = 0.05,
                 depth=c(2,2),min_leaf_size=rep(max(10,length(y)/100),2),
                 sf=0.75, gamma_positive = T,initial_values = NULL,
                 silent=F, pred_sigma = NULL, pred_gamma = NULL){
  if(!is.data.frame(X)) X = data.frame(X=X)
  if(!silent){
    cat("Fit gbex\n")
    pb = txtProgressBar(style = 3,width = 76)
  }
  if(is.null(lambda)){
    lambda = lambda_scale*c(1,1/lambda_ratio)
  }
  if(B == 0){
    warning("Setting B=0 is the same as fitting an unconditional GPD on y")
  }
  if(is.null(pred_sigma)){
    pred_sigma = colnames(X)
  } else if(!(all(pred_sigma %in%  colnames(X)))){
    stop("pred_sigma should only contain predictors names which are columns in X")
  }

  if(is.null(pred_gamma)){
    pred_gamma = colnames(X)
  } else if(!(all(pred_gamma %in%  colnames(X)))){
    stop("pred_gamma should only contain predictors names which are columns in X")
  }

  n = length(y)
  covariates = colnames(X)
  data = cbind(y,X)
  colnames(data) = c("y",covariates)

  # First parameters are the unconditional tail parameters
  if(is.null(initial_values)){
    theta_init = first_guess(y,gamma_positive)
  } else{
    theta_init = data.frame(s= rep(initial_values[1],n), g= rep(initial_values[2],n))
  }
  # Transofrmation of the paramters in order to fit the model with the right constrains
  theta_init_t = transform_parameters(theta_init,gamma_positive,inverse_transform=T)

  # Create a data.frame used for the boosting procedure with data, parameters, first and second derivatives
  boosting_df = get_boosting_df(data,theta_init_t,gamma_positive)

  # Save the results of the boosting steps
  # TREES contains the boosting trees
  # dev contains the deviance for each iteration
  trees_sigma = trees_gamma = list()
  dev=rep(mean(boosting_df$dev),B+1)
  if(B>0){
    for(b in 1:B){
      subsample_index = sample(1:n,round(sf*n),replace=F)
      # Fit gradient trees for sigma and gamma parameter
      tree_sigma = gradient_tree(boosting_df[subsample_index,pred_sigma,drop=FALSE], ##gradient_tree
                                 boosting_df$r_st[subsample_index],
                                 boosting_df$r2_st[subsample_index], depth[1], min_leaf_size[1])
      tree_gamma = gradient_tree(boosting_df[subsample_index,pred_gamma,drop=FALSE],
                                 boosting_df$r_gt[subsample_index],
                                 boosting_df$r2_gt[subsample_index], depth[2], min_leaf_size[2])



      # Use the estimated trees to update the parameters
      theta_hat = update_parameters(tree_sigma,tree_gamma,boosting_df,lambda)

      # create a new data.frame with updated derivatives
      boosting_df = get_boosting_df(data,theta_hat,gamma_positive)

      # Save the estimated trees and the deviance
      trees_sigma[[b]] = tree_sigma
      trees_gamma[[b]] = tree_gamma
      dev[b+1] = mean(boosting_df$dev)

      if(!silent){
        setTxtProgressBar(pb,b/(B+1))
      }
    }
  } else{
    theta_hat = theta_init_t
  }

  ## Transform the parameters back to orginal values of parameters
  theta_hat = transform_parameters(theta_hat,gamma_positive,inverse_transform=F)

  func_call = match.call()

  output = list(theta = theta_hat, dev = dev,
                trees_sigma = trees_sigma, trees_gamma = trees_gamma,
                theta_init= theta_init[1,],
                gamma_positive = gamma_positive,
                lambda=lambda,B=B,depth=depth,
                sf =sf, min_leaf_size = min_leaf_size,
                data = data,
                pred_sigma = pred_sigma,
                pred_gamma = pred_gamma,
                call = func_call)
  class(output) = "gbex"

  if(!silent){
    setTxtProgressBar(pb,(B+1)/(B+1))
    close(pb)
  }
  return(output)
}

#' Print function for gbex
#'
#' @param object A fitted gbex object
#' @return Prints the gbex object
#' @export
print.gbex <- function(object){
  cat(paste0(deparse(object$call),"\n"))
  cat("A gradient boosting model for extremes fitted by likelihood.\n")
  cat(paste0(object$B," trees are fitted.\n"))
  cat(paste0("Training error was equal to ",object$dev[length(object$dev)],".\n"))
}

#' Plot function for gbex
#'
#' @param object A fitted gbex object
#' @return Plots the deviance for each iteration
#' @export
plot.gbex <- function(object){
  data = data.frame(dev = object$dev, B = 0:object$B)
  g = ggplot2::ggplot(data,ggplot2::aes(y=dev,x=B)) +
    ggplot2::geom_line(size=1) +
    ggplot2::labs(title="Training Deviance", x = "Iteration (number of trees)", y = "Deviance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  return(g)
}
