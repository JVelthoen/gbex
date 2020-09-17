#' Fit gradient tree
#'
#' Fit a gradient tree and fill the leaf nodes with updates based on a newton raphson optimization step
#'
#' @param X the covariates a data.frame
#' @param rr the derivative of likelihood or powerdivergence
#' @param rr2 the second derivative of likelihood or powerdivergence
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A gradient_tree object containing
#' \item{tree}{The fitted tree object of class rpart}
#' \item{update_table}{data.frame with for each leafnode a newton raphson optimization step}
#' @details The procedure relies on the rpart function of the rpart package.
#' Only here the leafnode predictions are not average of rr of observations falling into the leaf, but instead a newton raphson optimization
#' step.
#' @export
gradient_tree <- function(X,rr,rr2,depth,min_leaf_size){
  if(depth == 0){
    ctrl = rpart::rpart.control(maxdepth = 1, minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = nrow(X))
  } else {
    ctrl = rpart::rpart.control(maxdepth = depth, minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = min_leaf_size)
  }
  data = cbind(data.frame(rr=rr),X)
  tree = rpart::rpart(rr~.,data=data, method='anova',control=ctrl)

  leaf_node_derivatives = split(data.frame(rr=rr,rr2=rr2),tree$where)
  new_update = unlist(lapply(leaf_node_derivatives,function(x) return(sum(x$rr)/sum(x$rr2))))
  update_table = data.frame(leaf = as.numeric(names(new_update)),update = unname(new_update))
  update_table$update = ifelse(abs(update_table$update) > 1,sign(update_table$update),update_table$update) ## NOTE HERE WE BOUND THE UPDATE
  gradient_tree = list(tree=tree,update_table=update_table)
  class(gradient_tree) = "gradient_tree"
  return(gradient_tree)
}
