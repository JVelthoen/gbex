#' Calculate partial dependence for gbex object
#'
#' @param object A gbex object
#' @param var_name a vector with one variable or two variable names to calculate partial depnedence
#' @return A list with for each variable the partial dependence of sigma and gamma with a
#' numeric value corresponding to the variable.
#' @export
calc_PD <- function(object,var_name){
  if(is.null(object$data)){
    stop("Data is not saved into gbex object")
  }

  if(length(var_name) == 1){
    values = sort(unique(as.vector(object$data[[var_name]])))

    theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)
    PD_per_tree_sigma = sapply(object$trees_sigma,PD_tree,var_name=var_name,values=values)
    PD_per_tree_gamma = sapply(object$trees_gamma,PD_tree,var_name=var_name,values=values)
    PD_transformed = data.frame(st = theta_init$st - apply(PD_per_tree_sigma,1,sum)*object$lambda[1],
                                gt = theta_init$gt -apply(PD_per_tree_gamma,1,sum)*object$lambda[2])
    PD = transform_parameters(PD_transformed,object$gamma_positive,inverse_transform=F)
    PD$values = values
  } else if(length(var_name) == 2){
    values = unique(object$data[var_name])

    theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)
    PD_per_tree_sigma = sapply(object$trees_sigma,PD_tree2,var_names=var_name,values=values)
    PD_per_tree_gamma = sapply(object$trees_gamma,PD_tree2,var_names=var_name,values=values)
    PD_transformed = data.frame(st = theta_init$st - apply(PD_per_tree_sigma,1,sum)*object$lambda[1],
                                gt = theta_init$gt -apply(PD_per_tree_gamma,1,sum)*object$lambda[2])
    PD = transform_parameters(PD_transformed,object$gamma_positive,inverse_transform=F)
    PD = cbind(PD,values)
  } else{
    stop("partial dependence can only be calculated for one or two variables")
  }
  return(PD)
}

#' Make partial dependence plot for given variable
#'
#' @param object A gbex object
#' @param variable variable name or index for which to make the dependence plot
#' @return Two figures side by side with on the left the partial dependence plot for sigma and on the right the partial dependence plot for gamma.
#' @export
partial_dependence <- function(object,variable){
  if(is.null(variable)) stop("variable is not specified")
  if(length(variable) > 1) stop("plotting is only possible for one variable")
  if(is.numeric(variable)){
    if(ncol(object$X) < round(variable)){
      stop("Column index out of range")
    } else{
      variable = colnames(object$data)[round(variable)]
    }
  } else if(is.null(object$data[[variable]])){
    stop("Variable name is not recognized.")
  }
  PD = calc_PD(object,variable)
  data1 = data.frame(values = PD$values,PD = PD$s)
  data2 = data.frame(values = PD$values,PD = PD$g)

  g1 = ggplot2::ggplot(data1,ggplot2::aes(x=values,y=PD)) +
    ggplot2::geom_line(lwd=1.5) +
    ggplot2::labs(title="sigma", x = variable, y = "Partial Dependence") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.border = ggplot2::element_rect(fill=NA))
  g2 = ggplot2::ggplot(data2,ggplot2::aes(x=values,y=PD)) +
    ggplot2::geom_line(lwd=2) +
    ggplot2::labs(title="gamma", x = variable, y = "Partial Dependence") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.border = ggplot2::element_rect(fill=NA))

  g = patchwork::wrap_plots(g1,g2)
  return(g)
}

#' Calculate partial dependence for gradient tree
#'
#' @param tree A gradient_tree object
#' @param var_name variable name for which to calulate the partial dependence
#' @param values Numeric vector with the values for var_name to calculate partial dependence
#' @return A vector with partial dependence from the gradient tree for vector values
#' @details The implementation makes use of the structure of the trees.
#' By keeping the observations in each node the tree can be followed down one time for the computation.
#' @export
PD_tree <- function(tree,var_name,values){

  tree_frame = tree$tree$frame
  tree_splits = tree$tree$splits
  update_table = tree$update_table

  PD_df = data.frame(tree_frame[c("var","n")],
                     node_nr = as.numeric(rownames(tree_frame)),
                     node_value = rep(NA,nrow(tree_frame)),
                     min = min(values),
                     max = max(values)+0.1)
  split_counter = 0

  for(row in 1:nrow(PD_df)){
    if(PD_df$var[row] == "<leaf>"){
      PD_df$node_value[row] = update_table$update[match(row,update_table$leaf)]
    } else if(PD_df$var[row] == var_name){
      split_counter = split_counter + 1
      child_index = match(PD_df$node_nr[row]*2 + (0:1),PD_df$node_nr)
      if(tree_splits[split_counter,2] == 1) child_index = rev(child_index)
      PD_df$n[child_index] = PD_df$n[row]
      PD_df$min[child_index] = c(PD_df$min[row],tree_splits[split_counter,4])
      PD_df$max[child_index] = c(tree_splits[split_counter,4],PD_df$max[row])
    } else{
      split_counter = split_counter + 1
      child_index = match(PD_df$node_nr[row]*2 + (0:1),PD_df$node_nr)
      if(tree_splits[split_counter,2] == 1) child_index = rev(child_index)
      PD_df$n[child_index] = PD_df$n[row]*(tree_frame$n[child_index]/tree_frame$n[row])
      PD_df$min[child_index] = rep(PD_df$min[row],2)
      PD_df$max[child_index] = rep(PD_df$max[row],2)
    }
  }

  PD_tree = numeric(length(values))
  for(row in which(PD_df$var == "<leaf>")){
      index = which(values >= PD_df$min[row] & values < PD_df$max[row])
      PD_tree[index] <- PD_tree[index] + PD_df$n[row]*PD_df$node_value[row]
  }
  PD_tree = PD_tree/tree_frame$n[1]
  return(PD_tree)
}


#' Calculate partial dependence for gradient tree for two variables
#'
#' @param tree A gradient_tree object
#' @param var_names two variable names for which to calulate the partial dependence
#' @param values dataframe with for each variable in var_names the values to calculate partial dependence
#' @return A vector with partial dependence from the gradient tree for vector values
#' @details The implementation makes use of the structure of the trees.
#' By keeping the observations in each node the tree can be followed down one time for the computation.
#' @export
PD_tree2 <- function(tree,var_names,values){

  tree_frame = tree$tree$frame
  tree_splits = tree$tree$splits
  update_table = tree$update_table

  PD_df = data.frame(tree_frame[c("var","n")],
                     node_nr = as.numeric(rownames(tree_frame)),
                     node_value = rep(NA,nrow(tree_frame)),
                     min.var1 = min(values[,1]),
                     max.var1 = max(values[,1])+0.1,
                     min.var2 = min(values[,2]),
                     max.var2 = max(values[,2])+0.1)

  split_counter = 0

  for(row in 1:nrow(PD_df)){
    if(PD_df$var[row] == "<leaf>"){
      PD_df$node_value[row] = update_table$update[match(row,update_table$leaf)]
    } else if(PD_df$var[row] == var_names[1]){
      split_counter = split_counter + 1
      child_index = match(PD_df$node_nr[row]*2 + (0:1),PD_df$node_nr)
      if(tree_splits[split_counter,2] == 1) child_index = rev(child_index)
      PD_df$n[child_index] = PD_df$n[row]
      PD_df$min.var1[child_index] = c(PD_df$min.var1[row],tree_splits[split_counter,4])
      PD_df$max.var1[child_index] = c(tree_splits[split_counter,4],PD_df$max.var1[row])
      PD_df$min.var2[child_index] = rep(PD_df$min.var2[row],2)
      PD_df$max.var2[child_index] = rep(PD_df$max.var2[row],2)
    } else if(PD_df$var[row] == var_names[2]){
      split_counter = split_counter + 1
      child_index = match(PD_df$node_nr[row]*2 + (0:1),PD_df$node_nr)
      if(tree_splits[split_counter,2] == 1) child_index = rev(child_index)
      PD_df$n[child_index] = PD_df$n[row]
      PD_df$min.var2[child_index] = c(PD_df$min.var2[row],tree_splits[split_counter,4])
      PD_df$max.var2[child_index] = c(tree_splits[split_counter,4],PD_df$max.var2[row])
      PD_df$min.var1[child_index] = rep(PD_df$min.var1[row],2)
      PD_df$max.var1[child_index] = rep(PD_df$max.var1[row],2)
    } else{
      split_counter = split_counter + 1
      child_index = match(PD_df$node_nr[row]*2 + (0:1),PD_df$node_nr)
      if(tree_splits[split_counter,2] == 1) child_index = rev(child_index)
      PD_df$n[child_index] = PD_df$n[row]*(tree_frame$n[child_index]/tree_frame$n[row])
      PD_df$min.var1[child_index] = rep(PD_df$min.var1[row],2)
      PD_df$max.var1[child_index] = rep(PD_df$max.var1[row],2)
      PD_df$min.var2[child_index] = rep(PD_df$min.var2[row],2)
      PD_df$max.var2[child_index] = rep(PD_df$max.var2[row],2)
    }
  }


  PD_tree = numeric(nrow(values))
  for(row in which(PD_df$var == "<leaf>")){
    index = which(values[,1] >= PD_df$min.var1[row] & values[,1] < PD_df$max.var1[row] &
                    values[,2] >= PD_df$min.var2[row] & values[,2] < PD_df$max.var2[row])
    PD_tree[index] <- PD_tree[index] + PD_df$n[row]*PD_df$node_value[row]
  }
  PD_tree = PD_tree/tree_frame$n[1]
  return(PD_tree)
}
