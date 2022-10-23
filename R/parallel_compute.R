#' Compute list in parallel
#'
#' @param l elements to parallize
#' @param func function to apply on elements of l
#' @param vars list of variable names to be distributed over the nodes
#' @param envir environment in which the varibales live
#' @param ncores number of cores to parallize
#' @return A CV_gbex object see description of CV_gbex
#' @export
parallel_compute <- function(l,func, vars, envir, ncores){
  cl <- parallel::makePSOCKcluster(ncores)
  parallel::clusterEvalQ(cl, c(library(gbex)))
  parallel::clusterExport(cl,varlist=vars,envir=envir)
  result = parallel::parLapply(cl, l, func)
  parallel::stopCluster(cl)
  return(result)
}
