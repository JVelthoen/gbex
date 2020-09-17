#' Permutation variable importance
#'
#' @param var_name Variable name for which to calculate variable importance
#' @param object A gbex object
#' @return A named vector with importance scores in decrease of deviance
#' @export
VI_permutation <- function(var_name, object){
  dev_model = object$dev[length(object$dev)]
  permuted_data = object$data[-1]
  permuted_data[[var_name]] = permuted_data[[var_name]][sample(1:nrow(permuted_data))]
  dev_perm = dev_per_step(object,y=object$data$y,X=permuted_data)[length(object$dev)]
  VI = dev_perm-dev_model
  return(VI)
}
