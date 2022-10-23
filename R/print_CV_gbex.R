#' Print function for CV_gbex
#'
#' @param object CV_gbex object
#' @return prints a basic description of the CV_gbex object
#' @export
print.CV_gbex <- function(object){
  cat(paste0(deparse(object$call),"\n"))
  cat(paste0("A cross validation object for gbex for parameter ",object$par_name,".\n"))
  cat(paste0("Results obtained by repeating ",object$repeat_cv, " times a cross validation with ",object$num_folds,"  folds.\n"))
  cat(paste0("Folds are sampled "))
  if(object$stratified){
    cat(" by stratified sampling.\n")
  } else{
    cat(" by random sampling.\n")
  }
  if(object$par_name %in% c("B","sf","lambda_ratio")){
    cat(paste0("The optimal parameter value is ",object$par_name," = ",object$par_CV,".\n"))
  } else{
    cat(paste0("The optimal parameter value is ",object$par_name," = (",paste(object$par_CV,collapse=","),").\n"))
  }
}
