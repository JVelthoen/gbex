#' Print function for CV_gbex
#'
#' @param object CV_gbex object
#' @return prints a basic description of the CV_gbex object
#' @export
print.CV_gbex <- function(object){
  cat(paste0(deparse(object$call),"\n"))
  cat(paste0("A cross validation object for gbex for parameter ",object$par_name,".\n"))
  cat(paste0("The number of folds is ",object$num_folds," which are"))
  if(object$stratified){
    cat(" obtained by stratified sampling.\n")
  } else{
    cat(" obtained by random sampling.\n")
  }
  if(object$par_name %in% c("B","sf","lambda_ratio")){
    cat(paste0("The optimal parameter value is ",object$par_name," = ",object$par_CV,".\n"))
  } else{
    cat(paste0("The optimal parameter value is ",object$par_name," = (",paste(object$par_CV,collapse=","),").\n"))
  }
}
