#' Cross validation for a parameter in gbex function
#'
#' @param X data.frame of covariates
#' @param y numeric vector of response
#' @param num_folds numeric number of folds to use
#' @param repeat_CV numeric indicating how many times to repeat the cross validation (repeating makes the cross validation less sensitive to how the folds are split)
#' @param Bmax numeric maximum number of trees to fit for each fold
#' @param par_name character name of the parameter for which to perform a cross validation (see details of CV_gbex)
#' @param par_grid list or vector a grid for the parameter in par_name (see details of CV_gbex)
#' @param stratified boolean indicating to use stratified sampling for the folds
#' @param ncores number of cores to used for parallelization of the cross validation, if not specified cores are chosen by detectCores function from the parallel package (Does not work on windows)
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return A CV_gbex object see CV_gbex
#' @export
CV_gbex_par <- function(y,X,num_folds,repeat_cv,par_name,par_grid,Bmax,stratified,ncores,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,repeat_cv,stratified)

  parallelization_list = unlist(lapply(1:length(folds),function(fold_nr){
    lapply(par_grid,function(par_value){
      list(par_value = par_value,fold_nr = fold_nr, test_index = folds[[fold_nr]])
    })
  }),recursive = F)
  func = function(job){
    n = length(y)
    test_index = job$test_index
    train_index = setdiff(1:n,test_index)
    par_value = job$par_value
    
    ytrain = y[train_index]
    ytest = y[test_index]
    Xtrain = X[train_index,]
    Xtest = X[test_index,]
    
    arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax))
    arguments_gbex[[par_name]] = par_value
    fit = do.call(gbex,arguments_gbex)
    dev = dev_per_step(fit,y=ytest,X=Xtest)
    return(dev)
  }
  
  # Compute in parallel
  cl <- parallel::makePSOCKcluster(ncores)
  parallel::clusterEvalQ(cl, c(library(gbex)))
  parallel::clusterExport(cl,varlist=list('y','X','Bmax','par_name','arguments','folds'),envir=environment())
  dev_list = parallel::parLapply(cl, folds, func)
  parallel::stopCluster(cl)

  job_fold_nr = sapply(parallelization_list,function(job){job$fold_nr})
  dev_matrix_list = tapply(dev_list,job_fold_nr,function(dev){do.call("cbind",dev)})

  dev = Reduce("+",dev_matrix_list)/len(folds)
  index_opt = which(apply(dev,2,min,na.rm = T) == min(dev,na.rm=T))
  B_opt = which(dev[,index_opt] == min(dev[,index_opt],na.rm=T))
  par_opt = unlist(par_grid[index_opt])

  if(any(is.na(dev))) warning("In some folds NA values were generated be carefull interpreting the results!")
  output = list(par_CV = par_opt, par_name = par_name, par_grid = par_grid,
                Bmax = Bmax, B_opt = B_opt,
                dev_all = dev, dev_folds = dev_matrix_list,
                num_folds = num_folds, folds=folds,
                repeat_cv = repeat_cv,
                y=y, X=X, stratified = stratified,
                call = match.call())
  class(output) = "CV_gbex"
  return(output)
}
