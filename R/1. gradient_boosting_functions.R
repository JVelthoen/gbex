#' GPD boosting
#'
#' Estimate the Generalized Pareto distribution conditional on covariates using a boosting procedure.
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param lambda learning rate for the scale and shape parameter
#' @param B Number of gradient boosting steps
#' @param depth Maximum depth of the trees
#' @param sf  sample fraction used for fitting the trees
#' @param alpha the power for power divergence (default alpha = 0 meaning maximum likelihood is used)
#' @param silent boolean indicating whether progress during fitting procedure should be printed.
#' @return gbex returns an object of class "gbex" which contains the following components:
#' \item{theta}{Data frame with the estimated gamma and sigma parameter for each observation}
#' \item{dev}{Numeric with deviance of model}
#' \item{trees_sigma}{List with gradient_tree objects for sigma}
#' \item{trees_gamma}{List with gradient_tree objects for gamma}
#' \item{lambda}{Numeric with the learining rate of sigma and gamma}
#' \item{B}{Numeric with number of trees}
#' \item{depth}{Numeric with maximum tree depth for sigma and gamma}
#' \item{alpha}{Power divergence parameter used}
#' @export
gbex <- function(y,X,B=180,lambda=c(0.025,0.0025),
                 depth=c(2,2),min_leaf_size=c(30,30),sf=0.5,
                 alpha = 0,silent=F){
  if(!is.data.frame(X)) X = data.frame(X=X)
  n = length(y)
  data = cbind(y,X)
  # First parameters are the unconditional tail parameters
  theta_init = first_guess(y)

  # Create a data.frame used for the boosting procedure with data, parameters, first and second derivatives
  boosting_df = get_boosting_df(data,theta_init,alpha)

  # Save the results of the boosting steps
  # TREES contains the boosting trees
  # dev contains the deviance for each iteration
  trees_sigma = trees_gamma = list()
  dev=rep(mean(boosting_df$dev),B+1)
  for (b in 1:B){
    # Take a subsample from the entire data.frame
    tree_df = boosting_df[sample(1:n,sf*n,replace=F),]

    # Fit gradient trees for sigma and gamma parameter
    tree_sigma = gradient_tree_sigma(tree_df,depth[1],min_leaf_size[1])
    tree_gamma = gradient_tree_gamma(tree_df,depth[2],min_leaf_size[2])

    # Use the estimated trees to update the parameters
    theta_hat = update_parameters(tree_sigma,tree_gamma,boosting_df,lambda)

    # create a new data.frame with updated derivatives
    boosting_df = get_boosting_df(data,theta_hat,alpha)

    # Save the estimated trees and the deviance
    trees_sigma[[b]] = tree_sigma
    trees_gamma[[b]] = tree_gamma
    dev[b+1] = mean(boosting_df$dev)

    if(!silent & b %in% round(((1:10)*(B/10)))){
      cat(paste0(round(b/B,1)*100,"% of trees fitted\n"))
    }
  }

  output = list(theta = theta_hat, dev = dev,
                 trees_sigma = trees_sigma, trees_gamma = trees_gamma,
                 theta_init= theta_init[1,],
                 lambda=lambda,B=B,depth=depth,alpha=alpha)
  class(output) = "gbex"
  return(output)
}


#' Predict function for gbex
#'
#' @param object A fitted gbex object
#' @param newdata A data frame with covariates for which to predict the sigma and gamma parameter
#' @param what Character indicating what to predict, currently only "par"
#' @return A data.frame object with the estimated sigma and gamma parameters
#' @export
predict.gbex <- function(object, newdata = NULL,what="par"){
  if(is.null(newdata)){
    theta = object$theta
  } else{
    sigma_updates = lapply(object$trees_sigma,predict,newdata=newdata)
    gamma_updates = lapply(object$trees_gamma,predict,newdata=newdata)

    theta = data.frame(s= object$theta_init$s - object$lambda[1]*Reduce("+",sigma_updates),
                        g= object$theta_init$g - object$lambda[2]*Reduce("+",gamma_updates))
  }
  return(theta)
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
    sigma_updates = cbind(0,do.call('cbind',lapply(object$trees_sigma,predict,newdata=X)))
    gamma_updates = cbind(0,do.call('cbind',lapply(object$trees_gamma,predict,newdata=X)))

    sigma_per_step = object$theta_init$s - object$lambda[1]*apply(sigma_updates,1,cumsum)
    gamma_per_step = object$theta_init$g - object$lambda[2]*apply(gamma_updates,1,cumsum)

    divergence_input = cbind(as.vector(t(sigma_per_step)),as.vector(t(gamma_per_step)),rep(y,object$B+1))
    if(object$alpha == 0){
      dev_vec = apply(divergence_input,1,GP_dev)
      dev_matrix = matrix(dev_vec,ncol=nrow(sigma_per_step))
      dev = apply(dev_matrix,2,mean)
    } else{
      stop("this is not implemented yet for power divergence")
    }
  }
}

par_grid = list(depth = list(c(1,1),c(2,2),c(3,3)),
                min_leaf_size = list(c(25,25),c(50,50),c(75,75)),
                sf = list(0.25,0.5,0.75))

CV_tune_hyper_parameters <- function(y,X,num_folds = 8,lambda_high,lambda_low,Bmax,par_fixed = NULL,par_grid = NULL){
  if(is.null(par_fixed)) par_fixed = list(depth=c(1,1),sf=0.5,min_leaf_size=rep(0.1*length(y),2),silent=F,alpha=0)
  par_fixed[["B"]] <- Bmax

  if(!is.null(par_grid)){
    par_fixed[["lambda"]] <- lambda_high
    B_search = CV_grid_search(y,X,NULL,NULL,num_folds,par_fixed)
    par_fixed[["B"]] = B_search$par

    if("depth" %in% names(par_grid)){
      par_fixed[["depth"]] = NULL
      depth_search = CV_grid_search(y,X,"depth",par_grid[["depth"]],num_folds,par_fixed)
      par_fixed[["depth"]] = depth_search$par
    }

    if("min_leaf_size" %in% names(par_grid)){
      par_fixed[["min_leaf_size"]] = NULL
      min_leaf_size_search = CV_grid_search(y,X,"min_leaf_size",par_grid[["min_leaf_size"]],num_folds,par_fixed)
      par_fixed[["min_leaf_size"]] = min_leaf_size_search$par
    }

    if("sf" %in% names(par_grid)){
      par_fixed[["sf"]] = NULL
      sf_search = CV_grid_search(y,X,"sf",par_grid[["sf"]],num_folds,par_fixed)
      par_fixed[["sf"]] = sf_search$par
    }
  }

  par_fixed[["lambda"]] <- lambda_low
  par_fixed[["B"]] <- Bmax
  B_search = CV_grid_search(y,X,NULL,NULL,num_folds,par_fixed)
  par_fixed[["B"]] = B_search$par

  return(par_fixed)
}
