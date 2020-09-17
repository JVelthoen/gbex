#' Cross validation for a parameter in gbex function
#'
#' @param X data.frame of covariates
#' @param y numeric vector of response
#' @param num_folds numeric number of folds to use
#' @param Bmax numeric maximum number of trees to fit for each fold
#' @param par_name character name of the parameter for which to perform a cross validation (see details of CV_gbex)
#' @param par_grid list or vector a grid for the parameter in par_name (see details of CV_gbex)
#' @param stratified boolean indicating to use stratified sampling for the folds
#' @param ncores number of cores to used for parallelization of the cross validation, if not specified cores are chosen by detectCores function from the parallel package (Does not work on windows)
#' @param ... Other arguments to be passed to the gbex function (details)
#' @return A CV_gbex object see CV_gbex
#' @export
CV_gbex_par <- function(y,X,num_folds,par_name,par_grid,Bmax,stratified,ncores,...){
  arguments = list(...)
  folds = divide_in_folds(y,num_folds,stratified)

  parallelization_list = unlist(lapply(1:num_folds,function(fold){
    lapply(par_grid,function(par_value){
      list(par_value = par_value,fold = fold)
    })
  }),recursive = F)


  dev_list = parallel::mclapply(parallelization_list,function(job){
    fold = job$fold
    par_value = job$par_value

    ytrain = y[folds!=fold]
    ytest = y[folds==fold]
    Xtrain = X[folds!=fold,]
    Xtest = X[folds==fold,]

    arguments_gbex = c(arguments,list(y=ytrain,X=Xtrain,B=Bmax))
    arguments_gbex[[par_name]] = par_value
    fit = do.call(gbex,arguments_gbex)
    dev = dev_per_step(fit,y=ytest,X=Xtest)
    return(dev)
  },mc.cores = ncores)

  job_fold = sapply(parallelization_list,function(job){job$fold})
  dev_matrix_list = tapply(dev_list,job_fold,function(dev){do.call("cbind",dev)})

  dev = Reduce("+",dev_matrix_list)/num_folds
  index_opt = which(apply(dev,2,min,na.rm = T) == min(dev,na.rm=T))
  B_opt = which(dev[,index_opt] == min(dev[,index_opt],na.rm=T))
  par_opt = unlist(par_grid[index_opt])

  if(any(is.na(dev))) warning("In some folds NA values were generated be carefull interpreting the results!")
  output = list(par_CV = par_opt, par_name = par_name, par_grid = par_grid,
                Bmax = Bmax, B_opt = B_opt,
                dev_all = dev, dev_folds = dev_matrix_list,
                num_folds = num_folds, folds=folds,
                y=y, X=X, stratified = stratified,
                call = match.call())
  class(output) = "CV_gbex"
  return(output)
}
