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
#' @param gamma_positive boolean indicating whether gamma should be positive
#' @param alpha the power for power divergence (default alpha = 0 meaning maximum likelihood is used)
#' @param silent boolean indicating whether progress during fitting procedure should be printed.
#' @param VI_type character vector indicating which variable importance should be calculated (when VI_type = NULL no variable importance is calculated)
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
gbex <- function(y,X,B=100,lambda=NULL,
                 lambda_ratio = 10, lambda_scale = 0.01,
                 depth=c(1,1),min_leaf_size=c(30,30),sf=0.75,
                 gamma_positive = T,
                 alpha = 0,silent=F, VI_type = NULL){
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

  n = length(y)
  covariates = colnames(X)
  data = cbind(y,X)
  colnames(data) = c("y",covariates)

  # First parameters are the unconditional tail parameters
  theta_init = first_guess(y,gamma_positive)

  # Transofrmation of the paramters in order to fit the model with the right constrains
  theta_init_t = transform_parameters(theta_init,gamma_positive,inverse_transform=T)

  # Create a data.frame used for the boosting procedure with data, parameters, first and second derivatives
  boosting_df = get_boosting_df(data,theta_init_t,alpha,gamma_positive)

  # Save the results of the boosting steps
  # TREES contains the boosting trees
  # dev contains the deviance for each iteration
  trees_sigma = trees_gamma = list()
  subsample_indices = matrix(nrow=round(sf*n),ncol = B)
  dev=rep(mean(boosting_df$dev),B+1)
  if(B>0){
    for(b in 1:B){
      # Take a subsample from the entire data.frame
      subsample_index = sample(1:n,round(sf*n),replace=F)
      subsample_indices[,b] = subsample_index

      # Fit gradient trees for sigma and gamma parameter
      tree_sigma = gradient_tree(boosting_df[subsample_index,covariates,drop=FALSE], ##gradient_tree
                                              boosting_df$r_st[subsample_index],
                                              boosting_df$r2_st[subsample_index], depth[1], min_leaf_size[1])
      tree_gamma = gradient_tree(boosting_df[subsample_index,covariates,drop=FALSE],
                                              boosting_df$r_gt[subsample_index],
                                              boosting_df$r2_gt[subsample_index], depth[2], min_leaf_size[2])



      # Use the estimated trees to update the parameters
      theta_hat = update_parameters(tree_sigma,tree_gamma,boosting_df,lambda)

      # create a new data.frame with updated derivatives
      boosting_df = get_boosting_df(data,theta_hat,alpha,gamma_positive)

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
                lambda=lambda,B=B,depth=depth,alpha=alpha,
                subsamples = subsample_indices,
                data = data,
                call = func_call)
  class(output) = "gbex"

  if(!is.null(VI_v)){
    output$variable_importance = calc_VI(output,type=VI_type)
  }

  if(!silent){
    setTxtProgressBar(pb,(B+1)/(B+1))
    close(pb)
  }
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
    warning("The gbex object does not contain so many trees")
  } else{
    Blim = round(Blim)
  }

  if(is.null(newdata) & Blim == object$B){
    theta = object$theta
  } else if(Blim>0){
      theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)
      sigma_updates = lapply(object$trees_sigma,predict,newdata=newdata)[1:Blim]
      gamma_updates = lapply(object$trees_gamma,predict,newdata=newdata)[1:Blim]

      theta_transformed = data.frame(st = theta_init$st - object$lambda[1]*Reduce("+",sigma_updates),
                                     gt = theta_init$gt - object$lambda[2]*Reduce("+",gamma_updates))

      theta = transform_parameters(theta_transformed,object$gamma_positive,inverse_transform=F)
  } else{
    theta_init = object$theta_init
    theta = data.frame(s = rep(theta_init$s,nrow(newdata)),g=rep(theta_init$g,nrow(newdata)))
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

#' Plot function for gbex
#'
#' @param object A fitted gbex object
#' @return Plots the deviance for each iteration
#' @export
plot.gbex <- function(object){
  data = data.frame(dev = object$dev, B = 0:object$B)
  g = ggplot2::ggplot(data,ggplot2::aes(y=dev,x=B)) +
    ggplot2::geom_line(size=1) +
    ggplot2::labs(title="Training deviance", x = "Iteration", y = "Deviance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  return(g)
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
    if(object$B==0){
      theta_input = data.frame(s = rep(object$theta_init$s,length(y)),g=rep(object$theta_init$g,length(y)))
    } else{
      theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)

      sigma_updates = cbind(0,do.call('cbind',lapply(object$trees_sigma,predict,newdata=X)))
      gamma_updates = cbind(0,do.call('cbind',lapply(object$trees_gamma,predict,newdata=X)))

      sigma_per_step = theta_init$st - object$lambda[1]*apply(sigma_updates,1,cumsum)
      gamma_per_step = theta_init$gt - object$lambda[2]*apply(gamma_updates,1,cumsum)

      theta_input = transform_parameters(data.frame(st = as.vector(sigma_per_step),
                                                    gt = as.vector(gamma_per_step)),
                                         object$gamma_positive, inverse_transform=F)
    }

    divergence_input = cbind(theta_input,as.vector(sapply(y,rep,times=object$B+1)))
    if(object$alpha == 0){
      dev_vec = apply(divergence_input,1,GP_dev)
      dev_matrix = matrix(dev_vec,nrow=object$B+1)
      dev = apply(dev_matrix,1,mean)
    } else{
      alpha = object$alpha
      A = apply(divergence_input,1,A_func)
      B = apply(divergence_input,1,B_func)
      divergence_input_PD = cbind(divergence_input,A,B)
      dev_vec = apply(divergence_input_PD,1,function(x) PD_dev(x[1:3],x[4],x[5],alpha))
      dev_matrix = matrix(dev_vec,nrow=nrow(sigma_per_step))
      dev = apply(dev_matrix,1,mean)
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
