#' Deviance per step
#'
#' Compute the deviance per boosting iteration
#'
#' @param object A fitted gbex object
#' @param X data.frame with a testset of covariates
#' @param y vector with a test set of observations
#' @return A vector with the deviance after adding a set of trees for sigma and gamma parameters
#' @details Without specifying y and X the data in the gbex object are used
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

    deviance_input = cbind(theta_input,as.vector(sapply(y,rep,times=object$B+1)))
    dev_vec = apply(deviance_input,1,GP_dev)
    dev_matrix = matrix(dev_vec,nrow=object$B+1)
    dev = apply(dev_matrix,1,mean)
  }
  return(dev)
}
