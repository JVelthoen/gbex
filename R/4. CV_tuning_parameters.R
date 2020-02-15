#' Cross validation for gbex
#'
#' A cross validation procedure for gbex.
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param par_name parameter name for which to do the cross validation (see details)
#' @param par_grid the grid of parameter values for par to do a gridsearch, either a list or vector (see details)
#' @param Bmax the maximum number of trees
#' @param stratified whether the cross validation procedure should use stratified sampling
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return A CV_gbex object
#' @details Performs a cross validation grid search over the provided grid,
#' for one variable at the time. Because the optimal number of trees depends heavily on the chosen parameter values,
#' for each parameter in grid Bmax trees are estimated and the optimal number over all cross validation folds is chosen to compare with other parameter values.
#'
#' For opimizing only the number of trees set par_name = NULL.
#' Other possible values of "par" are:
#' \itemize{
#' \item{depth}{The depth of the trees grid is specified to be a list where each element contains the depth of beta and gamma}
#' \item{sf}{Sample fraction of the trees, grid is specified to be a vector.}
#' \item{min_leaf_size}{The minimum leafsize of a tree where grid is specified to be a list where each element contains the minimum leafsize of beta and gamma}
#' \item{lambda_ratio}{The ratio between lambda_sigma/lambda_gamma}
#' \item{lambda_scale}{The size of lambda_sigma when lambda ratio is specified.}
#' \item{lambda}{A vector with lambda for sigma and gamma}
#' \item{alpha}{The power of the power divergence}
#'
#' }
#'
#' In the ... the values of other tuning parameters are specified. If these are not supplied the standard values of gbex will be used
#' @export
CV_gbex <- function(y,X,num_folds,Bmax,par_name=NULL,par_grid=NULL,stratified=F,...){
  if(is.null(par_name)){
    CV_result = CV_normal(y,X,num_folds,Bmax,stratified,...)
  } else if(par_name %in% c("lambda","lambda_ratio","lambda_scale","depth","min_leaf_size","sf","alpha")){
    CV_result = CV_par(y,X,num_folds,par_name,par_grid,Bmax,stratified,...)
  } else if(par_name == "B"){
    stop("For optimizing B leave par_name = NULL")
  } else{
    stop("This is not a parameter from gbex that can be optimized by cross validation.")
  }
  CV_result$call = match.call()
  return(CV_result)
}

#' Divide the data up in folds for cross validation
#'
#' @param y Vector of observations
#' @param num_folds Integer number of folds used for cross validation
#' @param stratified Boolean indicating wheter stratified sampling should be used
#' @return A vector with integers corresponding to the fold
#' @details The stratified sampling creates stratified clusters of size num_folds from the ordered observations y.
#' For each fold one observations is sampled.
#' @export
divide_in_folds <- function(y,num_folds,stratified = F){
  n = length(y)
  if(stratified){
    folds_matrix <- sapply(1:ceiling(n/num_folds),function(i){sample(1:num_folds)})
    folds_vector <- folds_matrix[1:n]
    folds <- folds_vector[rank(-y)]
  } else{
    index_shuffled = sample(1:n)
    folds = cut(seq(1,length(index_shuffled)),breaks=num_folds,labels=F)[order(index_shuffled)]
  }
  return(folds)
}

#' Cross validation to inspect performance of the  model
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param Bmax maximum number of trees used for finding the optimal
#' @param stratified indicate whether stratified sampling should be used
#' @param ... Additional arguments supplied to gbex function
#' @return A CV_gbex object
#' @export
CV_normal <- function(y,X,num_folds,Bmax,stratified,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  dev_matrix = sapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax))
    fit = do.call(gbex,arguments_gbex)
    dev = dev_per_step(fit,y=ytest,X=Xtest)
    return(dev)
  })

  dev = apply(dev_matrix,1,mean)
  B_opt = which(dev == min(dev))-1

  output = list(par_CV = B_opt, par = "B", grid = 1:Bmax,
                grid_B = 0:Bmax, dev_all = dev, dev_folds = dev_matrix,
                num_folds = num_folds,folds=folds, y=y,X=X,
                stratified = stratified, call = match.call())
  class(output) = "CV_gbex"
  return(output)
}

#' Cross validation for a parameter in gbex function
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param par_name the name of the parameter
#' @param par_grid a grid to perform parameter optimization either a vector or a list
#' @param Bmax maximum number of trees used for finding the optimal
#' @param stratified indicate whether stratified sampling should be used
#' @param ... Additional arguments supplied to gbex function
#' @return A CV_gbex object
#' @export
CV_par <- function(y,X,num_folds,par_name,par_grid,Bmax,stratified,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  dev_matrix_list = lapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    dev_matrix <- sapply(par_grid,function(par){
      arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax))
      arguments_gbex[[par_name]] = par
      fit = do.call(gbex,arguments_gbex)
      dev = dev_per_step(fit,y=ytest,X=Xtest)
      return(dev)
    })
    return(dev_matrix)
  })

  dev = Reduce("+",dev_matrix_list)/num_folds
  index_opt = which(apply(dev,2,min) == min(dev))
  B_opt = which(dev[,index_opt] == min(dev[,index_opt]))
  par_opt = unlist(par_grid[index_opt])

  output = list(par_CV = par_opt, par = par_name, grid = par_grid,
                grid_B = 0:Bmax, B_opt = B_opt,
                dev_all = dev, dev_folds = dev_matrix_list,
                num_folds = num_folds, folds=folds,
                y=y, X=X, stratified = stratified,
                call = match.call())
  class(output) = "CV_gbex"
  return(output)
}
