#' Cross validation for gbex
#'
#' A cross validation procedure for gbex.
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param par parameter name for which to do the cross validation (see details)
#' @param grid the grid of parameter values for par to do a gridsearch, either a list or vector (see details)
#' @param Bmax the maximum number of trees
#' @param stratified whether the cross validation procedure should use stratified sampling
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return A CV_gbex object
#' @details Performs a cross validation grid search over the provided grid,
#' for one variable at the time. Because the optimal number of trees depends heavily on the chosen parameter values,
#' for each parameter in grid Bmax trees are estimated and the optimal number over all cross validation folds is chosen to compare with other parameter values.
#'
#' The possible values of "par" are:
#' \itemize{
#' \item{B}{The number of trees for this no grid needs to be specified only the maximum number of trees Bmax }
#' \item{depth}{The depth of the trees grid is specified to be a list where each element contains the depth of beta and gamma}
#' \item{sf}{Sample fraction of the trees, grid is specified to be a vector.}
#' \item{min_leaf_size}{The minimum leafsize of a tree where grid is specified to be a list where each element contains the minimum leafsize of beta and gamma}
#' }
#' In the ... the values of other tuning parameters are specified. If these are not supplied the standard values of gbex will be used
#' @export
CV_gbex <- function(y,X,num_folds,par,Bmax,grid=NULL,stratified=F,...){
  if(par == "B"){
    CV_result = CV_B(y,X,num_folds,Bmax,stratified,...)
  } else if(par == "depth"){
    CV_result = CV_depth(y,X,num_folds,grid,Bmax,stratified,...)
  } else if(par == "sf"){
    CV_result = CV_sf(y,X,num_folds,grid,Bmax,stratified,...)
  } else if(par == "min_leaf_size"){
    CV_result = CV_min_leaf_size(y,X,num_folds,grid,Bmax,stratified,...)
  } else{
    stop("No cross validation gridserach defined for this parameter.")
  }
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

#' Cross validation for number of trees B
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param Bmax maximum number of trees used for finding the optimal
#' @param stratified indicate whether stratified sampling should be used
#' @param ... Additional arguments supplied to gbex function
#' @return A CV_gbex object
#' @export
CV_B <- function(y,X,num_folds,Bmax,stratified,...){
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
                num_folds = num_folds,folds=folds, y=y)
  class(output) = "CV_gbex"
  return(output)
}

#' Cross validation for the deth of the trees
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param depth_list A list with parameter values for depth
#' @param Bmax maximum number of trees used for finding the optimal
#' @param stratified indicate whether stratified sampling should be used
#' @param ... Additional arguments supplied to gbex function
#' @return A CV_gbex object
#' @export
CV_depth <- function(y,X,num_folds,depth_list,Bmax,stratified,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  dev_matrix_list = lapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    dev_matrix <- sapply(depth_list,function(depth){
      arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax,depth=depth))
      fit = do.call(gbex,arguments_gbex)
      dev = dev_per_step(fit,y=ytest,X=Xtest)
      return(dev)
    })
    return(dev_matrix)
  })

  dev = Reduce("+",dev_matrix_list)/num_folds
  Bopt = apply(dev,2,function(dev) which(dev == min(dev))-1)
  depth_opt = depth_list[[which(apply(dev,2,min) == min(dev))]]

  output = list(par_CV = depth_opt, par = "depth", grid = depth_list,
                grid_B = 0:Bmax,
                dev_all = dev, dev_folds = dev_matrix_list,
                num_folds = num_folds, folds=folds, y=y)
  class(output) = "CV_gbex"
  return(output)
}


#' Cross validation for the sample fraction for the trees
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param sf_vec A list with parameter values for depth
#' @param Bmax maximum number of trees used for finding the optimal
#' @param stratified indicate whether stratified sampling should be used
#' @param ... Additional arguments supplied to gbex function
#' @return A CV_gbex object
#' @export
CV_sf <- function(y,X,num_folds,sf_vec,Bmax,stratified,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  dev_matrix_list = lapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    dev_matrix <- sapply(sf_vec,function(sf){
      arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax,sf=sf))
      fit = do.call(gbex,arguments_gbex)
      dev = dev_per_step(fit,y=ytest,X=Xtest)
      return(dev)
    })
    return(dev_matrix)
  })

  dev = Reduce("+",dev_matrix_list)/num_folds
  index_opt = which(apply(dev,2,min) == min(dev))
  B_opt = which(dev[,index_opt] == min(dev[,index_opt]))
  sf_opt = depth[index_opt,Bopt]

  output = list(par_CV = sf_opt, par = "sf", grid = sf_vec,
                grid_B = 0:Bmax,
                dev_all = dev, dev_folds = dev_matrix_list,
                num_folds = num_folds, folds=folds, y=y)
  class(output) = "CV_gbex"
  return(output)
}

#' Cross validation for the minimum leaf size of the trees
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param sf_vec A list with parameter values for depth
#' @param Bmax maximum number of trees used for finding the optimal
#' @param stratified indicate whether stratified sampling should be used
#' @param ... Additional arguments supplied to gbex function
#' @return A CV_gbex object
#' @export
CV_min_leaf_size <- function(y,X,num_folds,min_leaf_size_list,Bmax,stratified,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  dev_matrix_list = lapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    dev_matrix <- sapply(min_leaf_size_list,function(min_leaf_size){
      arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax,depth=depth))
      fit = do.call(gbex,arguments_gbex)
      dev = dev_per_step(fit,y=ytest,X=Xtest)
      return(dev)
    })
    return(dev_matrix)
  })

  dev = Reduce("+",dev_matrix_list)/num_folds
  index_opt = which(apply(dev,2,min) == min(dev))
  B_opt = which(dev[,index_opt] == min(dev[,index_opt]))
  min_leaf_size_opt = depth[index_opt,Bopt]

  output = list(par_CV = min_leaf_size_opt, par = "min_leaf_size", grid = min_leaf_size_list,
                grid_B = 0:Bmax,
                dev_all = dev, dev_folds = dev_matrix_list,
                num_folds = num_folds, folds=folds, y=y)
  class(output) = "CV_gbex"
  return(output)
}
