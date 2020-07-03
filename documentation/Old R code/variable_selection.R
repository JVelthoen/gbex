variable_importance <- function(object,type){
  VI = calc_VI(object,type)
  if(type=="permutation"){
    data = data.frame(dev=VI$dev,var=factor(names(VI$dev),levels=names(VI$dev[order(VI$dev,decreasing=T)])))
    g = ggplot2::ggplot(data,ggplot2::aes(x=var,y=dev)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::labs(title="Permutation variable importance", x = "Variable", y = "Importance") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.border = ggplot2::element_rect(fill=NA))
  } else if(type == "relative"){
    if(!(is.null(VI$gamma) & is.null(VI$sigma))){
      data = data.frame(VI=c(VI$sigma,VI$gamma),
                        par = c(rep("sigma",length(VI$sigma)),rep("gamma",length(VI$gamma))),
                        var = c(names(VI$sigma),names(VI$gamma)))
    } else if(!is.null(VI$sigma)){
      data = data.frame(VI= VI$sigma,
                        par = rep("sigma",length(VI$sigma)),
                        var = names(VI$sigma))
    } else if(!is.null(VI$sigma)){
      data = data.frame(VI= VI$gamma,
                        par = rep("gamma",length(VI$gamma)),
                        var = names(VI$gamma))
    } else{
      stop("Importance is not defined for a model with (0,0) depth")
    }

    g = ggplot2::ggplot(data,ggplot2::aes(x=var,y=VI)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::labs(title="Relative importance", x = "Variable", y = "Importance") +
      ggplot2::facet_wrap(ggplot2::vars(par)) +
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
    if(object$depth[1] > 0){
      VI_per_tree = sapply(object$trees_sigma,VI_relative_tree,var_names=var_names)
      VI_sigma = apply(VI_per_tree,1,sum)
      VI_sigma = VI_sigma/max(VI_sigma) * 100
    } else {
      VI_sigma = NULL
    }

    if(object$depth[2] > 0){
      VI_per_tree = sapply(object$trees_gamma,VI_relative_tree,var_names=var_names)
      VI_gamma = apply(VI_per_tree,1,sum)
      VI_gamma = VI_gamma/max(VI_gamma) * 100
    } else{
      VI_gamma = NULL
    }
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


#' Variable importance for gbex fitted object
#'
#' @param object A gbex object
#' @return A plot with an ordered histogram of the scaled variable importance measures
#' @details The permutation importance measure is calculated for each tree by
#' the difference of the deviance of the out-of-bag data and the deviance for
#' the out of bag data with a single variable permuted
#' @export
variable_importance <- function(object){
  data = object$data
  data_perm = cbind(data.frame(y=data$y),data[sample(1:nrow(data)),-1])

  ## First calulate the parameter values after fitting each tree
  theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)

  sigma_updates = cbind(0,do.call('cbind',lapply(object$trees_sigma,predict,newdata=X)))
  gamma_updates = cbind(0,do.call('cbind',lapply(object$trees_gamma,predict,newdata=X)))

  s_t_per_step = theta_init$st - object$lambda[1]*apply(sigma_updates,1,cumsum)
  g_t_per_step = theta_init$gt - object$lambda[2]*apply(gamma_updates,1,cumsum)

  theta_per_step = transform_parameters(data.frame(st = as.vector(s_t_per_step),
                                                   gt = as.vector(g_t_per_step)),
                                        object$gamma_positive, inverse_transform=F)
  s_per_step = matrix(theta_per_step$s,nrow=object$B + 1)
  g_per_step = matrix(theta_per_step$g,nrow=object$B + 1)

  VI = matrix(0,ncol=ncol(data)-1,nrow=object$B)
  for(b in 1:object$B){
    ib_index = object$subsamples[,b]
    oob_index = setdiff(1:nrow(data),ib_index)

    ## Calculate normal deviance
    dev_input = cbind(s_per_step[b+1,oob_index],g_per_step[b+1,oob_index],data$y[oob_index])
    dev_oob = mean(apply(dev_input,1,GP_dev))

    ## Calculate permutation deviance
    vars_used = unique(c(names(object$trees_sigma[[b]]$tree$variable.importance),
                         names(object$trees_gamma[[b]]$tree$variable.importance)))
    vars_used = vars_used[order(match(vars_used,colnames(data)))]
    for(var in vars_used){
      permuted_data = data[oob_index,]
      permuted_data[var] = data_perm[oob_index,var]

      update_sigma = -object$lambda[1]*predict(object$trees_sigma[[b]],permuted_data)
      update_gamma = -object$lambda[2]*predict(object$trees_gamma[[b]],permuted_data)
      theta_new_t = data.frame(st = s_t_per_step[b,oob_index] + update_sigma,
                               gt = g_t_per_step[b,oob_index] + update_gamma)
      theta_new = transform_parameters(theta_new_t,object$gamma_positive,inverse_transform=F)
      dev_permuted = mean(apply(cbind(theta_new,permuted_data$y),1,GP_dev))
      VI[b,match(var,colnames(data))-1] = dev_permuted - dev_oob
    }
  }
  VI = apply(VI,2,cumsum)
  names(VI) = colnames(data)[-1]
  return(VI)
}






