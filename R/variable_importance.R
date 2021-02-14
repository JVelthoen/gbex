#' Variable importance for a gbex object
#'
#' @param object A gbex object
#' @param type either "relative" for relative importance or "permutation" for permutation importancee
#' @param scaled boolean if True the variable importance is scaled such that the most important variable is scaled to 100.
#' @return A VI_gbex object with
#' \item{type}{the type of variable importance}
#' \item{VI_df}{data.frame with the variable importance measures for each variable}
#' \item{scaled}{If the variable importance is scaled after computing}
#' @details The relative importance score is calculated sepearately for sigma and gamma parameter and the permutation importance measure is calculated for deviance.
#' @export
variable_importance <- function(object,type=c("relative","permutation"),scaled = T){
  var_names = colnames(object$data[-1])
  VI_df = data.frame(variables = var_names)

  if(type == "relative"){
    if(object$depth[1] > 0){
      VI_per_tree = sapply(object$trees_sigma,VI_relative,var_names=var_names)
      VI_sigma = apply(VI_per_tree,1,sum)
      if(scaled){
        VI_sigma = VI_sigma/max(VI_sigma) * 100
      }
    } else {
      VI_sigma = NULL
    }

    if(object$depth[2] > 0){
      VI_per_tree = sapply(object$trees_gamma,VI_relative,var_names=var_names)
      VI_gamma = apply(VI_per_tree,1,sum)
      if(scaled){
        VI_gamma = VI_gamma/max(VI_gamma) * 100
      }
    } else{
      VI_gamma = NULL
    }
    if(!is.null(VI_sigma)) VI_df$sigma = VI_sigma
    if(!is.null(VI_gamma)) VI_df$gamma = VI_gamma
  } else if(type == "permutation"){
    VI_dev = sapply(var_names,VI_permutation,object= object)
    if(scaled){
      VI_dev = VI_dev/max(VI_dev)*100
    }
    VI_df$permutation = VI_dev
  } else{
    stop("This type of variable importance is not implemented")
  }
  output = list(type=type, VI_df =VI_df, scaled=scaled)
  class(output) = "VI_gbex"
  return(output)
}
