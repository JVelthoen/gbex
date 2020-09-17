#' Cross validation for the number of trees B
#'
#' @param X data.frame of covariates
#' @param y numeric vector of response
#' @param num_folds numeric number of folds to use
#' @param Bmax numeric maximum number of trees to fit for each fold
#' @param stratified boolean indicating to use stratified sampling for the folds
#' @param ncores number of cores to used for parallelization of the cross validation, if not specified cores are chosen by detectCores function from the parallel package (Does not work on windows)
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return A CV_gbex object see description of CV_gbex
#' @export
CV_gbex_B <- function(y,X,num_folds,Bmax,stratified,ncores,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  dev_list = parallel::mclapply(1:num_folds,function(fold){
    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax))
    fit = do.call(gbex,arguments_gbex)
    dev = dev_per_step(fit,y=ytest,X=Xtest)
    return(dev)
  },mc.cores=ncores)
  dev_matrix = do.call("cbind",dev_list)

  dev = apply(dev_matrix,1,mean)
  if(any(is.na(dev))){
    dev = dev[1:(which(is.na(dev))[1]-1)]
    warning("In CV validation observations out of range of tail fit. Optimal choice of B might not be correct.")
  }
  B_opt = which(dev == min(dev))-1

  output = list(par_CV = B_opt, par_name = "B", par_grid = 0:Bmax,
                Bmax =  Bmax, dev_all = dev, dev_folds = dev_matrix,
                num_folds = num_folds, folds=folds, y=y,X=X,
                stratified = stratified, call = match.call())
  class(output) = "CV_gbex"
  return(output)
}
