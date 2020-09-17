#' Gradient boosting for extremes with threshold estimation
#'
#' Estimate extreme quantiles of y|X by first estimating the p quantile and then estimating the tail model with gbex function
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param lambda learning rate for the scale and shape parameter
#' @param lambda_ratio the ratio of the lambda_sigma/lambda_gamma (Use this together with lambda_scale, see details)
#' @param lambda_scale equal to the suze of lambda_sigma  (Use this together with lambda_scale, see details)
#' @param B Number of gradient boosting steps
#' @param depth Maximum depth of the trees
#' @param sf  sample fraction used for fitting the trees
#' @param p probability level to estimate the threshold
#' @param gamma_positive boolean indicating whether gamma should be positive
#' @param silent boolean indicating whether progress during fitting procedure should be printed.
#' @param VI_type character vector indicating which variable importance should be calculated (when VI_type = NULL no variable importance is calculated)
#' @param ... additional arguments for the quantile forest for the threshold model
#' @return gbex returns an object of class "gbex" which contains the following components:
#' \item{theta}{Data frame with the estimated gamma and sigma parameter for each observation}
#' \item{dev}{Numeric with deviance of model}
#' \item{trees_sigma}{List with gradient_tree objects for sigma}
#' \item{trees_gamma}{List with gradient_tree objects for gamma}
#' \item{lambda}{Numeric with the learining rate of sigma and gamma}
#' \item{B}{Numeric with number of trees}
#' \item{depth}{Numeric with maximum tree depth for sigma and gamma}
#' @details Instead of specifying the learning rates lambda for sigma and gamma seperately the ratio of the two can be specified together with the size of the first one.
#' This can be done using lambda_ratio and lambda_scale. The used lambda is then given by lambda_scale*(1,1/lambda_ratio).
#' @export
gbex_threshold <- function(y, X, B=100, lambda=NULL,
                           lambda_ratio = 10, lambda_scale = 0.01,
                           depth=c(1,1), min_leaf_size=c(30,30), sf=0.75,
                           p=0.8, gamma_positive = T,
                           silent=F, VI_type = NULL,...){

  threshold_model = forest_threshold(y,X,p,...)

  z = y - threshold_model$threshold
  X_exceed = X[z>0,]
  z = z[z>0]

  gbex_model = gbex(z,X_exceed,B,lambda,lambda_ratio,lambda_scale,depth,
                    min_leaf_size,sf,gamma_positive,silent,VI_type)

  output = list(threshold_model = threshold_model,gbex_model = gbex_model, p = p)
  class(output) = "gbex_threshold"
  return(output)
}




#' Estimate a threshold model for gbex extrapolation
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param p probability level to estimate the threshold
#' @param ... additional arguments for the quantile forest
#' @return a list with the estimated quantile forest and the estimated threshold
#' @export
forest_threshold <- function(y,X,p=0.8,...){
  arguments_extra = list(...)
  arguments_grf = list(X=X,Y=y)
  forest = do.call(grf::quantile_forest,c(arguments_grf,arguments_extra))
  threshold = predict(forest,quantiles = p)[,1]

  output = list(forest = forest, threshold= threshold)
  return(output)
}


#' Predict function for gbex_threshold
#'
#' @param object A fitted gbex_threshold object
#' @param newdata A data frame with covariates for which to predict the sigma and gamma parameter
#' @param probs A vector of probabilities for which to estimate quantiles (needs par="quant")
#' @param what Character indicating what to predict, currently only "par"
#' @param Blim Integer number indicating how many trees to use of the fitted object only used for newdata for gbex data
#' @return A data.frame object with the estimated sigma and gamma parameters
#' @export
predict.gbex_threshold <- function(object, newdata = NULL,probs = NULL,what="par",Blim=NULL){
  if(what == "par"){
    output = predict(object$gbex_model,newdata=newdata,what = "par",Blim=Blim)
  } else if(what == "quant"){
    if(any(probs < object$p)){
      stop("Quantiles can only be computed for probabilities larger than p")
    }
    rescaled_probs = 1-(1-probs)/(1-object$p)
    gbex_quantile = predict(object$gbex_model,newdata=newdata,what = "quant",
                            probs = rescaled_probs,Blim=Blim)
    threshold_quantile = predict(output$threshold_model$forest,newdata=newdata,quantiles = object$p)
    output = threshold_quantile + gbex_quantile
  }
  return(output)
}



