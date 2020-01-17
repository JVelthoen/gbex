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
#' @param grid_search Boolean indicating if tuning parameters need to be optimized
#' @param par_grid Named list with a grid for each parameter to tune (Only when grid_search = True)
#' @param lambda_grid An extra lambda parameter used for grid_search (Only when grid_search = True)
#' @param num_folds Number of folds used for grid search Default = 5 (Only when grid_search = True)
#' @return gbex returns an object of class "gbex" which contains the following components:
#' \item{theta}{Data frame with the estimated gamma and sigma parameter for each observation}
#' \item{dev}{Numeric with deviance of model}
#' \item{trees_sigma}{List with gradient_tree objects for sigma}
#' \item{trees_gamma}{List with gradient_tree objects for gamma}
#' \item{lambda}{Numeric with the learining rate of sigma and gamma}
#' \item{B}{Numeric with number of trees}
#' \item{depth}{Numeric with maximum tree depth for sigma and gamma}
#' \item{alpha}{Power divergence parameter used}
#' @details The cross validation procedure is currently implemented for depth, min_leaf_size, sf and B.
#' The initial starting parameters are the ones that are supplied to the function. Then the procedure is as follows:
#' \itemize{
#'  \item{Find the optimal B with initial parameters and lambda_grid}
#'  \item{Find optimal depth over the grid of depth using the chosen B and lambda_grid}
#'  \item{Find optimal min_leaf_size over the grid of min_leaf_size using the chosen B and lambda_grid}
#'  \item{Find optimal sf over the grid of sf using the chosen B and lambda_grid}
#'  \item{Find optimal B with new chosen parameters and lambda}
#' }
#' Note that if a grid for a parameter is not supplied the step for this parameter is skipped.
#' For B no grid needs to be specified only a maximum initial value.
#' @export
gbex <- function(y,X,B=180,lambda=c(0.025,0.0025),
                 depth=c(2,2),min_leaf_size=c(30,30),sf=0.5,
                 alpha = 0,silent=F,
                 grid_search = F, par_grid = NULL, lambda_grid=NULL, num_folds = 5){
  if(!is.data.frame(X)) X = data.frame(X=X)
  if(grid_search){
    par_fixed = list(depth=depth,min_leaf_size=min_leaf_size,sf=sf,alpha=alpha,silent=silent)
    if(is.null(lambda_grid)) lambda_grid = lambda
    Bmax = B
    par_opt = CV_tuning_parameters(y,X,num_folds,lambda_grid,lambda,Bmax,par_fixed,par_grid)
    B = par_opt$B
    depth = par_opt$depth
    min_leaf_size = par_opt$min_leaf_size
    sf = par_opt$sf
  }

  if(!silent) cat("Fitting Boosting Trees for Model:\n")
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


#' Tuning parameters
#'
#' Obtain optimal tuning parameters for gbex by performing a grid search.
#'
#' @param y Numeric vector of observations
#' @param X Data frame with the right column names
#' @param num_folds Integer for the number of folds in the cross validation
#' @param lambda_grid Numeric with learning rate for sigma and gamma used for doing the grid search
#' @param lambda Numeric with learning rate for sigma and gamma used for final model
#' @param Bmax Numeric indicating the maximum number of trees
#' @param par_fixed Named list with for each parameter an initial value
#' @param par_grid Named list with for each parameter a grid to optimize over
#' @return A named list with the optimal parameters
#' @details The function uses two lambda parameters.
#' The parameters lambda_grid should be larger than lambda and is used to perform the gridsearch as less trees need to be fitted.
#' @export
CV_tuning_parameters <- function(y,X,num_folds = 8,lambda_grid,lambda,Bmax,par_fixed,par_grid = NULL){
  if(par_fixed$silent == F){
    par_fixed$silent = T
    silent = F
  } else{
    silent =T
  }

  if(!silent) cat("Cross Validation procedure:\n")

  par_fixed[["B"]] <- Bmax

  if(!is.null(par_grid)){
    par_fixed[["lambda"]] <- lambda_grid
    B_search = CV_grid_search(y,X,NULL,NULL,num_folds,par_fixed)
    par_fixed[["B"]] = B_search$par
    if(!silent) cat(paste0("B for lambda_grid set to: ",B_search$par,".\n"))

    if("depth" %in% names(par_grid)){
      par_fixed[["depth"]] = NULL
      depth_search = CV_grid_search(y,X,"depth",par_grid[["depth"]],num_folds,par_fixed)
      par_fixed[["depth"]] = depth_search$par
      if(!silent) cat(paste0("depth set to: c(",paste0(depth_search$par,collapse=", "),").\n"))
    }

    if("min_leaf_size" %in% names(par_grid)){
      par_fixed[["min_leaf_size"]] = NULL
      min_leaf_size_search = CV_grid_search(y,X,"min_leaf_size",par_grid[["min_leaf_size"]],num_folds,par_fixed)
      par_fixed[["min_leaf_size"]] = min_leaf_size_search$par
      if(!silent) cat(paste0("min_leaf_size set to: c(",paste0(min_leaf_size_search$par,collapse=", "),").\n"))
    }

    if("sf" %in% names(par_grid)){
      par_fixed[["sf"]] = NULL
      sf_search = CV_grid_search(y,X,"sf",par_grid[["sf"]],num_folds,par_fixed)
      par_fixed[["sf"]] = sf_search$par
      if(!silent) cat(paste0("sf set to: ",sf_search$par,".\n"))
    }
  }

  par_fixed[["lambda"]] <- lambda
  par_fixed[["B"]] <- Bmax
  B_search = CV_grid_search(y,X,NULL,NULL,num_folds,par_fixed)
  par_fixed[["B"]] = B_search$par
  if(!silent) cat(paste0("B for lambda set to: ",B_search$par,".\n"))

  return(par_fixed)
}
