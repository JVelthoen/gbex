#' Relative variable importance for a gradient tree
#'
#' @param tree a gradient tree object
#' @param var_names The set of variable names used for fitting
#' @return A named vector with importance scores of the variables for the tree
#' @export
VI_relative <- function(tree,var_names){
  tree_frame = tree$tree$frame
  tree_frame = tree_frame[tree_frame$var != "<leaf>",]
  MSE_per_split = tree_frame$dev/tree_frame$wt
  MSE_per_var = aggregate(MSE_per_split,by=list(tree_frame$var),FUN=sum)
  VI = numeric(length(var_names))
  VI[match(MSE_per_var$Group.1,var_names)] = MSE_per_var$x
  names(VI) = var_names
  return(VI)
}
