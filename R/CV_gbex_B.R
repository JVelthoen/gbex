#' Cross validation for the number of trees B
#'
#' @param X data.frame of covariates
#' @param y numeric vector of response
#' @param num_folds numeric number of folds to use
#' @param repeat_cv numeric indicating how many times to repeat the cross validation (repeating makes the cross validation less sensitive to how the folds are split)
#' @param Bmax numeric maximum number of trees to fit for each fold
#' @param stratified boolean indicating to use stratified sampling for the folds
#' @param ncores number of cores to used for parallelization of the cross validation, if not specified cores are chosen by detectCores function from the parallel package (Does not work on windows)
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return A CV_gbex object see description of CV_gbex
#' @export
CV_gbex_B <- function(y,X,num_folds,repeat_cv,Bmax,stratified,ncores,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,repeat_cv,stratified)
  func = function(test_index){
    n = length(y)
    train_index = setdiff(1:n,test_index)
    ytrain = y[train_index]
    ytest = y[test_index]
    Xtrain = X[train_index,]
    Xtest = X[test_index,]
  
    arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax))
    fit = do.call(gbex,arguments_gbex)
    dev = dev_per_step(fit,y=ytest,X=Xtest)
    return(dev)
  }
  # Compute the cross validation in parallel
  cl <- parallel::makePSOCKcluster(ncores)
  parallel::clusterEvalQ(cl, c(library(gbex)))
  parallel::clusterExport(cl,varlist=list('y','X','Bmax','arguments'),envir=environment())
  dev_list = parallel::parLapply(cl, folds, func)
  parallel::stopCluster(cl)
  
  dev_matrix = do.call("cbind",dev_list)

  dev = apply(dev_matrix,1,mean)
  if(any(is.na(dev))){
    dev = dev[1:(which(is.na(dev))[1]-1)]
    warning("In CV validation observations out of range of tail fit. Optimal choice of B might not be correct.")
  }
  B_opt = which(dev == min(dev))-1

  output = list(par_CV = B_opt, par_name = "B", par_grid = 0:Bmax,
                Bmax =  Bmax, dev_all = dev, dev_folds = dev_matrix,
                num_folds = num_folds, folds=folds,  repeat_cv = repeat_cv,
                y=y,X=X,stratified = stratified, call = match.call())
  class(output) = "CV_gbex"
  return(output)
}
