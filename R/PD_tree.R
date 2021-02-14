#' Calculate partial dependence for gradient tree
#'
#' @param tree A gradient_tree object
#' @param var_name variable name for which to calulate the partial dependence
#' @param var_names variable names for which to calculate the partial dependence
#' @param values Numeric vector with the values for var_name to calculate partial dependence
#' @return A vector with partial dependence from the gradient tree for vector values
#' @details The implementation makes use of the structure of the trees.
#' By keeping the observations in each node the tree can be followed down one time for the computation.
#' The function PD_tree is used for calculating one dimensional partial dependence functions and PD_tree2 is used for 2 dimensional partial dependence functions.
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

#' @rdname PD_tree
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
