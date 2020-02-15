#' Calculate variable importance for a gbex object
#'
#' @param object A gbex object
#' @param type a character vector either relative for "relative" importance or "permutation" for permutation importance
#' @return A plot with an ordered histogram of the scaled variable importance measures
#' @details In the case of relative importance two figures are supplied with the importance for sigma and for gamma.
#' For permutation importance only one figure is supplied as importance for deviance.
#' @export
variable_importance <- function(object,type){
  VI = calc_VI(object,type)
  if(type=="permutation"){
    data = data.frame(dev=VI$dev,var=factor(names(VI$dev),levels=names(VI$dev[order(dev,decreasing=T)])))
    g = ggplot2::ggplot(data,ggplot2::aes(x=var,y=dev)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::labs(title="Permutation variable importance", x = "Variable", y = "Importance") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.border = ggplot2::element_rect(fill=NA))
  } else if(type == "relative"){
    data = data.frame(VI=c(VI$sigma,VI$gamma),
                            par = c(rep("sigma",length(VI$sigma)),rep("gamma",length(VI$gamma))),
                            var = c(names(VI$sigma),names(VI$gamma)))
    g = ggplot2::ggplot(data,ggplot2::aes(x=var,y=VI)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::labs(title="Relative importance", x = "Variable", y = "Importance") +
      ggplot2::facet_wrap(vars(par)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.border = ggplot2::element_rect(fill=NA))
  }
  return(g)
}

#' Calculate variable importance for a gbex object
#'
#' @param object A gbex object
#' @param type a character vector either relative for "relative" importance or "permutation" for permutation importance
#' @return A named vector with scaled importance scores for each variable
#' @details The relative importance score is calculated sepearately for sigma and gamma parameter.
#' the permutation importance measure is calculated for deviance.
#' @export
calc_VI <- function(object,type=c("relative","permutation")){
  if(is.null(object$X)){
    stop("Data is not saved into gbex object")
  }

  var_names = colnames(object$X)
  if(type == "relative"){
    VI_per_tree = sapply(object$trees_beta,VI_relative_tree,var_names=var_names)
    VI_sigma = apply(VI_per_tree,1,sum)
    VI_sigma = VI_sigma/max(VI_sigma) * 100

    VI_per_tree = sapply(object$trees_gamma,VI_relative_tree,var_names=var_names)
    VI_gamma = apply(VI_per_tree,1,sum)
    VI_gamma = VI_gamma/max(VI_gamma) * 100
    VI = list(sigma = VI_sigma, gamma = VI_gamma, type=type)
  } else if(type == "permutation"){
    VI_dev = sapply(var_names,VI_permutation,object= object)
    VI_dev = VI_dev/max(VI_dev)*100
    VI = list(dev = VI_dev,type=type)
  } else{
    stop("This type of variable importance is not implemented")
  }
  return(VI)
}

#' Calculate permutation variable importance
#'
#' @param var_name The variable name
#' @param object A gbex object
#' @return A named vector with importance scores relative to deviance
#' @export
VI_permutation <- function(var_name, object){
  dev_model = object$dev[length(object$dev)]
  permuted_data = object$X
  permuted_data[[var_name]] = permuted_data[[var_name]][sample(1:nrow(permuted_data))]
  dev_perm = dev_per_step(object,y=object$y,X=permuted_data)[length(object$dev)]
  VI = dev_perm-dev_model
  return(VI)
}

#' Calculate relative variable importance for a gradient tree
#'
#' @param tree a gradient tree object
#' @param var_names The set of variable names used for fitting
#' @return A named vector with importance scores of the variables for the tree
#' @export
VI_relative_tree <- function(tree,var_names){
  tree_frame = tree$tree$frame
  tree_frame = tree_frame[tree_frame$var != "<leaf>",]
  MSE_per_split = tree_frame$dev/tree_frame$wt
  MSE_per_var = aggregate(MSE_per_split,by=list(tree_frame$var),FUN=sum)
  VI = numeric(length(var_names))
  VI[match(MSE_per_var$Group.1,var_names)] = MSE_per_var$x
  names(VI) = var_names
  return(VI)
}
