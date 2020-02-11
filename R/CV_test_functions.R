#' Function to find optimal number of trees with other parmaeters given
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param Bmax maximum number of trees used for finding the optimal
#' @param lambda the learning rate
#' @param depth the depth of the trees
#' @param min_leaf_size the minimum leaf node size
#' @param sf subample fraction used for fitting the trees
#' @return The optimal number of trees used
#' @export
find_optimal_B <- function(y,X,num_folds = 4,Bmax,lambda,depth,min_leaf_size,sf){
  n = length(y)
  folds = divide_in_folds(n,num_folds)

  dev_matrix = sapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    fit = gbex(ytrain,Xtrain,Bmax,lambda,depth,min_leaf_size,sf,alpha=0)
    dev = dev_per_step(fit,y=ytest,X=Xtest)
    return(dev)
  })

  dev = apply(dev_matrix,1,mean)
  Bopt = which(dev == min(dev))-1
  return(Bopt)
}


#' Function to find optimal number of trees with other parmaeters given
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param Bmax maximum number of trees used for finding the optimal
#' @param lambda the learning rate
#' @param depth the depth of the trees
#' @param min_leaf_size the minimum leaf node size
#' @param sf subample fraction used for fitting the trees
#' @return The optimal number of trees used
#' @export
find_optimal_depth <- function(y,X,num_folds = 4,depth_list,Bmax,lambda,min_leaf_size,sf){
  n = length(y)
  folds = divide_in_folds(n,num_folds)

  dev_matrix_list = lapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    dev_matrix <- sapply(depth_list,function(depth){
      fit = gbex(ytrain,Xtrain,Bmax,lambda,depth,min_leaf_size,sf,alpha=0)
      dev = dev_per_step(fit,y=ytest,X=Xtest)
      return(dev)
    })
    return(dev_matrix)
  })

  dev = Reduce("+",dev_matrix_list)/num_folds
  Bopt = apply(dev,2,function(dev) which(dev == min(dev))-1)
  dev_opt = apply(dev,2,min)
  return(list(depth= depth_list[[which(dev_opt == min(dev_opt))]],dev_per_depth = dev_opt))
}


#' Function to find optimal number of trees with other parmaeters given
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param num_folds the number of folds used to determine the optimal number of trees
#' @param Bmax maximum number of trees used for finding the optimal
#' @param lambda the learning rate
#' @param depth the depth of the trees
#' @param min_leaf_size the minimum leaf node size
#' @param sf subample fraction used for fitting the trees
#' @return The optimal number of trees used
#' @export
find_optimal_sf <- function(y,X,num_folds = 8,sf_vec,Bmax,lambda,depth,min_leaf_size){
  n = length(y)
  folds = divide_in_folds(n,num_folds)

  dev_matrix_list = lapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    dev_matrix <- sapply(sf_vec,function(sf){
      fit = gbex(ytrain,Xtrain,Bmax,lambda,depth,min_leaf_size,sf,alpha=0)
      dev = dev_per_step(fit,y=ytest,X=Xtest)
      return(dev)
    })
    return(dev_matrix)
  })

  dev = Reduce("+",dev_matrix_list)/num_folds
  Bopt = apply(dev,2,function(dev) which(dev == min(dev))-1)
  dev_opt = apply(dev,2,min)
  return(list(sf= sf_vec[which(dev_opt == min(dev_opt))],dev_per_sf = dev_opt))
}
