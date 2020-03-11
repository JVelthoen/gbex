#' Fit gradient Tree
#'
#' Estimate a gradient tree then calculate for the leafnodes the optimal gradient based on a single newton raphson step
#'
#' @param X the covariates a data.frame
#' @param rr the derivative of likelihood or powerdivergence
#' @param rr2 the second derivative of likelihood or powerdivergence
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A list with tree object and a named vector with for each leafnode the update value
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

#' Get the gradient update for each leafnode
#'
#' @param object A gradient tree object
#' @param newdata A New data frame with columns used for the gradient tree; If NULL the updates for the model data is computed

#' @return A vector with gradient steps
#' @export
predict.gradient_tree <- function(object,newdata=NULL){
  update_table = object$update_table
  tree = object$tree
  if(is.null(newdata)){
    leaf = treeClust::rpart.predict.leaves(tree)
    update = update_table$update[match(leaf,update_table$leaf)]
  } else{
    leaf = treeClust::rpart.predict.leaves(tree,newdata=newdata)
    update = update_table$update[match(leaf,update_table$leaf)]
  }
  return(update)
}


#' Update parameters
#'
#' Take a single gradient step for both gamma and sigma parameters
#'
#' @param tree_sigma A gradient tree object fitted for sigma parameter
#' @param tree_gamma A gradient tree object fitted for gamma parameter
#' @param boosting_df the boosting_df
#' @param lambda a vector with the learning rate of sigma and gamma
#' @return A data frame with the udpated parameters sigma and gamma
#' @export
update_parameters <- function(tree_sigma,tree_gamma,boosting_df,lambda){
  update_sigma = -lambda[1]*predict(tree_sigma,boosting_df)
  update_gamma = -lambda[2]*predict(tree_gamma,boosting_df)

  theta <- data.frame(st = boosting_df$st + update_sigma, gt = boosting_df$gt + update_gamma)
  return(theta)
}



