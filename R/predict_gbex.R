#' Predict function for gbex
#'
#' Fpredicted values for a fitted gbe model.
#'
#' @param object A fitted gbex object
#' @param newdata A data.frame with covariates for which to predict if NULL it is done for insample data
#' @param what character with "par" for parameters, "quant" for quantiles
#' @param probs numeric with probabilities for which to estimate quantiles  (needs what == "quant")
#' @param Blim Numeric number of trees to use for the prediction
#' @return A data.frame object with the estimated sigma and gamma parameters or a data.frame with the estimated quantiles
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
