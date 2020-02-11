#' Cross Validation Grid Search
#'
#' For a given grid for a parameter find the optimal value such that the cross validation deviance is minimized
#'
#' @param X Data frame of covariates
#' @param y Numeric vector of response
#' @param par_name Character name of the parameter to adjust
#' @param grid Either a numeric vector or list of possible parameter values
#' @param num_folds Integer indicating the number of folds
#' @param par_fixed Named list of values of parameters that are kept fixed
#' @return The function returns a named list with the following elements:
#' \item{par}{Optimal parameter value for all parameter values in grid}
#' \item{par_name}{The paramter name}
#' \item{grid}{A vector or list with tried parameter values}
#' \item{dev}{Numeric vector with deviance for each parameter value in grid}
#' \item{B}{Integers with optimal B value for each value in grid}
#' @export
CV_grid_search <- function(y,X,par_name=NULL,grid=NULL,num_folds = 4,par_fixed = list()){
  if(is.null(grid)){
    grid = par_fixed$B
    par_name = "B"
    par_fixed[[par_name]] <- NULL
  }

  n = length(y)
  folds = divide_in_folds(n,num_folds)

  dev = numeric(length(grid))
  B_opt = numeric(length(grid))
  for(grid_cnt in 1:length(grid)){
   args = par_fixed
   args[[par_name]] = unlist(grid[grid_cnt])
   dev_per_fold = lapply(1:num_folds,CV_single_fold,folds=folds,data=list(y=y,X=X),args=args)
   dev_all_folds = Reduce("+",dev_per_fold)/num_folds
   B_opt[grid_cnt] = which(dev_all_folds == min(dev_all_folds)) - 1
   dev[grid_cnt] = dev_all_folds[B_opt[grid_cnt]]
  }

  if(par_name!="B"){
    index_opt <- which(dev == min(dev))
    par_opt <- unlist(grid[index_opt])
  } else{
    par_opt=B_opt
  }
  return(list(par=par_opt,par_name=par_name,grid=grid,dev=dev,B_opt=B_opt))
}

#' Estimate a single fold for cross validation
#'
#' @param fold Integer indicating which fold should be computed
#' @param folds Integer vector of length n with the division in folds
#' @param data List where response in named y and covariates are name X
#' @param args Named list with arguments for the gbex model
#' @return numeric vector with deviance of each boosting iteration
#' @export
CV_single_fold <- function(fold,folds,data,args){
  bool_test = (folds == fold)
  y_test = data$y[bool_test]
  y_train = data$y[!bool_test]
  X_test = data$X[bool_test,]
  X_train = data$X[!bool_test,]
  args = c(args,list(X=X_train,y=y_train))
  fit = do.call(gbex,args)
  dev = dev_per_step(fit,y=y_test,X=X_test)
  return(dev)
}

#' Divide the data up in folds for cross validation
#'
#' @param n Integer number of observations
#' @param num_folds Integer number of folds used for cross validation
#' @return A vector with integers corresponding to the fold
#' @export
divide_in_folds <- function(n,num_folds){
  index_shuffled = sample(1:n)
  folds = cut(seq(1,length(index_shuffled)),breaks=num_folds,labels=F)[order(index_shuffled)]
  return(folds)
}
