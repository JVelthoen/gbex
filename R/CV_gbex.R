#' Cross validation for gbex
#'
#' A cross validation procedure for gbex that can be used for any tuning parameter by specifying a grid of parameter values over which to perform the cross validation.
#'
#' @param X data.frame of covariates
#' @param y numeric vector of response
#' @param num_folds numeric number of folds to use
#' @param Bmax numeric maximum number of trees to fit for each fold
#' @param par_name character name of the parameter for which to perform a cross validation (see details)
#' @param par_grid list or vector a grid for the parameter in par_name (see details)
#' @param stratified boolean indicating to use stratified sampling for the folds
#' @param ncores number of cores to used for parallelization of the cross validation, if not specified cores are chosen by detectCores function from the parallel package (Does not work on windows)
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return CV_gbex returns an object of class "CV_gbex" which contains the following components:
#' \item{par_CV}{Optimal parameter based on cross validation}
#' \item{par_name}{character parameter name for which CV is performed}
#' \item{par_grid}{parameter grid used in CV}
#' \item{Bmax}{Maximum number of trees estimated per fold}
#' \item{dev}{matrix or vector (if par_name = NULL) of deviance averaged over the folds for each iteration 0:Bmax}
#' \item{dev_folds}{list of matrices or matrix (if par_name = NULL) of deviance for each fold for each iteration 0:Bmax}
#' \item{B_opt}{Vector of the optimal B value}
#' \item{num_folds}{The number of folds used in the CV}
#' \item{folds}{vector of integers of length length(y) assigning each observation to a fold}
#' \item{y}{The response}
#' \item{X}{The covariates}
#' \item{stratified}{boolean inidcating if folds are obtained by stratified sampling}
#' @details
#' If there is no par_name is set to NULL the cross validation is over the values of B.
#' If par_name is specified togehter with par_grid the cross validation is two dimensional over the specified parameter and B.
#' This is necesary the optimal B changes as a function of the other tuning parameters.
#'
#' Sometimes not for each fold Bmax trees are fitted because the test data falls outisde the support of the GPD for the specific parameters.
#' This generally happens due to overfitting, it could also be an indication to take a smaller learning rate.
#'
#' Other possible values of "par_name" are:
#' \itemize{
#' \item{depth}{The depth of the trees grid is specified to be a list where each element contains the depth of beta and gamma}
#' \item{sf}{Sample fraction of the trees, grid is specified to be a vector.}
#' \item{min_leaf_size}{The minimum leafsize of a tree where grid is specified to be a list where each element contains the minimum leafsize of beta and gamma}
#' \item{lambda_ratio}{The ratio between lambda_sigma/lambda_gamma}
#' \item{lambda_scale}{The size of lambda_sigma when lambda ratio is specified.}
#' \item{lambda}{A vector with lambda for sigma and gamma}
#' }
#'
#' @export
CV_gbex <- function(y,X,num_folds,Bmax,
                    par_name=NULL,par_grid=NULL,stratified=F, ncores = parallel::detectCores(),...){
  if(is.null(par_name)){
    CV_result = CV_gbex_B(y,X,num_folds,Bmax,stratified,ncores,...)
  } else if(par_name %in% c("lambda","lambda_ratio","lambda_scale","depth","min_leaf_size","sf","alpha")){
    CV_result = CV_gbex_par(y,X,num_folds,par_name,par_grid,Bmax,stratified,ncores,...)
  } else if(par_name == "B"){
    stop("For optimizing B leave par_name = NULL")
  } else{
    stop("This is not a parameter from gbex that can be optimized by cross validation.")
  }
  CV_result$call = match.call()
  return(CV_result)
}
