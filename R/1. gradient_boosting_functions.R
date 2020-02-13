#' GPD boosting
#'
#' Estimate the Generalized Pareto distribution conditional on covariates using a boosting procedure.
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param lambda learning rate for the scale and shape parameter
#' @param lambda_ratio the ratio of the lambda_sigma/lambda_gamma (Use this together with lambda_scale, see details)
#' @param lambda_scale equal to the suze of lambda_sigma  (Use this together with lambda_scale, see details)
#' @param B Number of gradient boosting steps
#' @param depth Maximum depth of the trees
#' @param sf  sample fraction used for fitting the trees
#' @param alpha the power for power divergence (default alpha = 0 meaning maximum likelihood is used)
#' @param silent boolean indicating whether progress during fitting procedure should be printed.
#' @param save_data boolean indicating whether data should be saved.
#' @return gbex returns an object of class "gbex" which contains the following components:
#' \item{theta}{Data frame with the estimated gamma and sigma parameter for each observation}
#' \item{dev}{Numeric with deviance of model}
#' \item{trees_sigma}{List with gradient_tree objects for sigma}
#' \item{trees_gamma}{List with gradient_tree objects for gamma}
#' \item{lambda}{Numeric with the learining rate of sigma and gamma}
#' \item{B}{Numeric with number of trees}
#' \item{depth}{Numeric with maximum tree depth for sigma and gamma}
#' \item{alpha}{Power divergence parameter used}
#' @details Instead of specifying the learning rates lambda for sigma and gamma seperately the ratio of the two can be specified together with the size of the first one.
#' This can be done using lambda_ratio and lambda_scale. The used lambda is then given by lambda_scale*(1,1/lambda_ratio).
#' @export
gbex <- function(y,X,B=180,lambda=NULL,
                 lambda_ratio = 15, lambda_scale = 0.01,
                 depth=c(2,2),min_leaf_size=c(30,30),sf=0.5,
                 alpha = 0,silent=F, save_data = T){
  if(!is.data.frame(X)) X = data.frame(X=X)
  if(!silent) cat("Fitting Boosting Trees for Model:\n")
  n = length(y)
  data = cbind(y,X)
  if(is.null(lambda)){
    lambda = lambda_scale*c(1,1/lambda_ratio)
  }
  # First parameters are the unconditional tail parameters
  theta_init = first_guess(y)

  # Create a data.frame used for the boosting procedure with data, parameters, first and second derivatives
  boosting_df = get_boosting_df(data,theta_init,alpha)

  # Save the results of the boosting steps
  # TREES contains the boosting trees
  # dev contains the deviance for each iteration
  trees_beta = trees_gamma = list()
  dev=rep(mean(boosting_df$dev),B+1)
  for (b in 1:B){
    # Take a subsample from the entire data.frame
    tree_df = boosting_df[sample(1:n,sf*n,replace=F),]

    # Fit gradient trees for sigma and gamma parameter
    tree_beta = gradient_tree_beta(tree_df,depth[1],min_leaf_size[1])
    tree_gamma = gradient_tree_gamma(tree_df,depth[2],min_leaf_size[2])

    # Use the estimated trees to update the parameters
    theta_hat = update_parameters(tree_beta,tree_gamma,boosting_df,lambda)

    # create a new data.frame with updated derivatives
    boosting_df = get_boosting_df(data,theta_hat,alpha)

    # Save the estimated trees and the deviance
    trees_beta[[b]] = tree_beta
    trees_gamma[[b]] = tree_gamma
    dev[b+1] = mean(boosting_df$dev)

    if(!silent & b %in% round(((1:10)*(B/10)))){
      cat(paste0(round(b/B,1)*100,"% of trees fitted\n"))
    }
  }

  func_call = match.call()

  output = list(theta = theta_hat, dev = dev,
                 trees_beta = trees_beta, trees_gamma = trees_gamma,
                 theta_init= theta_init[1,],
                 lambda=lambda,B=B,depth=depth,alpha=alpha,
                 call = func_call)

  if(save_data){
    output$y = y
    output$X = X
  } else{
    output$y = NULL
    output$X = NULL
  }

  class(output) = "gbex"
  return(output)
}


#' Predict function for gbex
#'
#' @param object A fitted gbex object
#' @param newdata A data frame with covariates for which to predict the sigma and gamma parameter
#' @param probs A vector of probabilities for which to estimate quantiles (needs par="quant")
#' @param what Character indicating what to predict, currently only "par"
#' @param Blim Integer number indicating how many trees to use of the fitted object only used for newdata
#' @return A data.frame object with the estimated sigma and gamma parameters
#' @export
predict.gbex <- function(object, newdata = NULL,probs = NULL,what="par",Blim=NULL){
  if(is.null(Blim)){
    Blim = object$B
  } else if(Blim > object$B){
    Blim = object$B
  } else{
    Blim = round(Blim)
  }

  if(Blim == object$B){
    theta = object$theta
    theta$s <- exp(theta$b)
    theta <- theta[c("s","g")]
  } else if(!is.null(newdata)){
    beta_updates = lapply(object$trees_beta,predict,newdata=newdata)[1:Blim]
    gamma_updates = lapply(object$trees_gamma,predict,newdata=newdata)[1:Blim]

    theta = data.frame(s = exp(object$theta_init$b - object$lambda[1]*Reduce("+",beta_updates)),
                       g= object$theta_init$g - object$lambda[2]*Reduce("+",gamma_updates))
  } else{
    stop("Blim can only be used when new data is provided")
  }

  if(what == "par"){
    return(theta)
  } else if(what == "quant"){
    if(is.null(probs)){
      warning("probs needs specification for quantile prediction, using probs = 0.99")
      probs = 0.99
    }
      quant <- sapply(probs,function(p,s,g) (s/g) * ( (1-p)^(-g) - 1),s = theta$s,g = theta$g)
      return(quant)
  } else {
    stop("This specification of par is not defined")
  }
}

#' Print function for gbex
#'
#'
#' @param object A fitted gbex object
#' @return Prints the gbex object
#' @export
print.gbex <- function(object){
  cat(paste0(deparse(object$call),"\n"))
  cat("A gradient boosting model for extremes fitted by ")
  if(object$alpha == 0){
    cat("likelihood.\n")
  } else{
    cat(paste0("power divergence with alpha = ",object$alpha,".\n"))
  }
  cat(paste0(object$B," trees are fitted.\n"))
  cat(paste0("Training error was equal to ",object$dev[length(object$dev)],".\n"))
}


#' Deviance per step
#'
#' Compute the deviance per boosting step
#'
#' @param object A fitted gbex object
#' @param X A dataframe with the right column names
#' @param y A vector of observations
#' @return A vector with the deviance at each step
#' @export
dev_per_step <- function(object,y=NULL,X=NULL){
  if(is.null(X)){
    dev = object$dev
  } else{
    beta_updates = cbind(0,do.call('cbind',lapply(object$trees_beta,predict,newdata=X)))
    gamma_updates = cbind(0,do.call('cbind',lapply(object$trees_gamma,predict,newdata=X)))

    beta_per_step = object$theta_init$b - object$lambda[1]*apply(beta_updates,1,cumsum)
    gamma_per_step = object$theta_init$g - object$lambda[2]*apply(gamma_updates,1,cumsum)

    divergence_input = cbind(exp(as.vector(t(beta_per_step))),as.vector(t(gamma_per_step)),rep(y,object$B+1))
    if(object$alpha == 0){
      dev_vec = apply(divergence_input,1,GP_dev)
      dev_matrix = matrix(dev_vec,ncol=nrow(beta_per_step))
      dev = apply(dev_matrix,2,mean)
    } else{
      stop("this is not implemented yet for power divergence")
    }
  }
  return(dev)
}

#' Trim trees from gbex object
#'
#' @param object A fitted gbex object
#' @param B_new A new total number of trees
#' @return A gbex object with B_new trees
#' @export
trim_trees <- function(object,B_new){
  if(B_new >= object$B){
    warning("B_new is larger than number of trees in object")
    return(object)
  } else{
    object$B = B_new
    object$dev = object$dev[1:(B_new + 1)]
    object$trees_beta = object$trees_beta[1:B_new]
    object$trees_gamma = object$trees_gamma[1:B_new]
    object$theta = predict(object,newdata=object$X)
    return(object)
  }
}
