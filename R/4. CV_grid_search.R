#' Optimize the model by doing a gridsearch over a parameter using cross validation
#'
#' @param X data.frame of covariates
#' @param y vector of response
#' @param par character name of the parameter to adjust, if it is NULL cross validation is performed using parameters in par_fixed
#' @param grid_values numeric vector of values to perform gridsearch ove, if it is NULL cross validation is performed using parameters in par_fixed
#' @param num_folds integer number of folds used for cross validation
#' @param par_fixed named list of values of parameters that should remain fixed
#' @return the optimal parameter value and a dataframe with cross validation performance for each of the parameter values together with the chosen B
#' @export
CV_grid_search <- function(X,y,par=NULL,grid_values=NULL,num_folds = 4,par_fixed){
  folds <- divide_in_folds(y,num_folds)

  CV_results <- 1:num_folds %>%
    purrr::map(get_CV_data,X=X,y=y,folds=folds) %>%
    purrr::map(get_argument_list,par=par,grid_values=grid_values,par_fixed=par_fixed) %>%
    purrr::map(performance_per_fold) %>%
    add_folds_togehter()

  if(all(!is.null(par),!is.null(grid_values))){
    CV_df <- get_CV_df(CV_results,grid_values,par)
    par_opt <- get_optimal_parameters(CV_df)
    return(list(par_opt=par_opt,CV_df=CV_df))
  } else{
    CV_df <- get_CV_df(CV_results,1,"par")
    par_opt <- CV_df$B
    return(list(par_opt = par_opt))
  }
}

#' Divide the data up in folds for cross validation
#'
#' @param y vector of response
#' @param num_folds integer number of folds used for cross validation
#' @return a vector with integers corresponding to the fold
#' @export
divide_in_folds <- function(y,num_folds){
  n <- length(y)
  index_shuffled <- sample(1:n)
  folds <- cut(seq(1,length(index_shuffled)),breaks=num_folds,labels=F)[order(index_shuffled)]
  return(folds)
}

#' Run the gradient boosting algorithm for a specific set of parameters
#'
#' @param Xtrain data.frame of covariates for training
#' @param ytrain vector of response for training
#' @param Xtest data.frame of covariates for testing
#' @param ytest vector of response for testing
#' @param B maximum number of trees
#' @param lambda vector with learning rates for sigma and gamma
#' @param depth vector with depth of the trees for sigma and gamma
#' @param min_leaf_size vector with the minimum leaf size of trees for sigma and gamma
#' @param sf the sample fraction for each tree
#' @return the deviance for each step of the boosting algorithm
#' @export
CV_model_run <- function(Xtrain,ytrain,Xtest,ytest,B,lambda,depth,min_leaf_size,sf){
  fit <- gbex(y=ytrain,X=Xtrain,
              B=B,lambda=lambda,depth=depth,
              min_leaf_size=min_leaf_size,sf=sf,
              silent=T,alpha=0)
  dev_test <- dev_per_step(fit,Xtest,ytest)
  return(dev_test)
}

#' Get a list of training and testing data for a single fold
#'
#' @param fold a numeric indicating which of the folds in the cross validation
#' @param folds a vector with for each data point the fold to which they belong
#' @param X data.frame of covariates
#' @param y vector of response
#' @return a list with training and testing data
#' @export
get_CV_data <- function(fold,folds,X,y){
  data <- list(Xtrain = X[folds!=fold,],
         ytrain = y[folds!=fold],
         Xtest = X[folds==fold,],
         ytest = y[folds==fold])
  return(data)
}

#' Get a list of function arguments for gbex for each of the grid_values
#'
#' @param data a list with training and testing data
#' @param par a character with the name of the paramter over which to do a grid search
#' @param grid_values a numeric vector with values over which to do a grid_search
#' @param par_fixed a named list with fixed values for the other parameters
#' @return a list with all arguments for the gbex function
#' @export
get_argument_list <- function(data,par,grid_values,par_fixed){
  if(all(!is.null(par),!is.null(grid_values))){
    argument_list <- grid_values %>%
      purrr::array_branch(m=1) %>%
      purrr::map(function(value){
        data[[par]] <- value
        return(c(data,par_fixed))
      })
  } else{
    argument_list = list(c(data,par_fixed))
  }
  return(argument_list)
}

#' Calculate deviance per fold in the cross validation
#'
#' @param arg_lists a list with argument lists to call the gbex function
#' @return a matrix where each row matches the the elements of arg_lists and the columns the number of trees
#' @export
performance_per_fold <- function(arg_lists){
  performance_matrix <- purrr::invoke_map(CV_model_run,arg_lists) %>%
    purrr::invoke('rbind',.)
  return(performance_matrix)
}

#' Add the cross validation of different folds together
#'
#' @param fold_results list of results of each fold
#' @return a pointwise mean over the different folds
#' @export
add_folds_togehter <- function(fold_results){
  num_folds <- length(fold_results)
  aggregated_folds <- purrr::reduce(fold_results,`+`)/num_folds
  return(aggregated_folds)
}

#' Calculate the optimal parameter
#'
#' @param CV_df a dataframe obtained by get_CV_df
#' @return a numeric vector with the optimal parameters
#' @export
get_optimal_parameters <- function(CV_df){
  par_opt <- CV_df %>%
    dplyr::filter(dev == min(dev)) %>%
    dplyr::select(-c(dev,B)) %>%
    unlist(use.names=F)
  return(par_opt)
}

#' Calculate the CV dataframe with for each parameter value the deviance and the best chosen number of trees
#'
#' @param performance a matrix as returned by performance_by_fold
#' @param grid_values the values of the parameter to try
#' @param par character vecotr with the name of a single parameter
#' @return A dataframe with for each parameter value the deviance and the best chosen number of trees
#' @export
get_CV_df <- function(performance,grid_values,par){
  CV_df <- performance %>%
    apply(1,function(x){
      data.frame(B=which(x == min(x)) -1, dev=min(x))
      }) %>%
    purrr::invoke("rbind",.) %>%
    cbind(par = grid_values)
  return(CV_df)
}
