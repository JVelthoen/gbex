#' Gradient Tree beta
#'
#' Estimate a tree for the beta variable then calculate for the leafnodes how to update.
#' Beta is a reparametrization of sigma in the sense of beta = log(sigma)
#'
#' @param tree_df A dataframe to estimate the tree of the same format as the boosting_df
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A list with tree object and a named vector with for each leafnode the update value
#' @export
gradient_tree_beta <- function(tree_df,depth,min_leaf_size){
  ctrl = rpart::rpart.control(maxdepth = depth, minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = min_leaf_size)
  variable_names = setdiff(names(tree_df),c("y","dev","b","g","r_b","r_g","r2_b","r2_g"))
  formula = as.formula( paste( "r_b",paste( variable_names, collapse = " + "), sep = " ~ ") )
  tree = rpart::rpart(formula,data=tree_df, method='anova',control=ctrl)

  leaf_node_derivatives = split(tree_df[c('r_b','r2_b')],tree$where)
  update = unlist(lapply(leaf_node_derivatives,function(x) return(sum(x$r_b)/sum(x$r2_b))))
  update_table = data.frame(leaf = as.numeric(names(update)),update = unname(update))
  update_table$update = ifelse(abs(update_table$update) > 1,sign(update_table$update),update_table$update) ## NOTE HERE WE CAP THE UPDATE
  gradient_tree = list(tree=tree,update_table=update_table)
  class(gradient_tree) = "gradient_tree"
  return(gradient_tree)
}

#' Gradient Tree gamma
#'
#' Estimate a tree for the gamma variable then calculate for the leafnodes how to update
#'
#' @param tree_df A dataframe to estimate the tree of the same format as the boosting_df
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A list with tree object and a named vector with for each leafnode the update value
#' @export
gradient_tree_gamma <- function(tree_df,depth,min_leaf_size){
  ctrl = rpart::rpart.control(maxdepth = depth, minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = min_leaf_size)
  variable_names = setdiff(names(tree_df),c("y","dev","b","g","r_b","r_g","r2_b","r2_g"))
  formula = as.formula( paste( "r_g",paste( variable_names, collapse = " + "), sep = " ~ ") )
  tree = rpart::rpart(formula,data=tree_df, method='anova',control=ctrl)

  leaf_node_derivatives = split(tree_df[c('r_g','r2_g')],tree$where)
  update = unlist(lapply(leaf_node_derivatives,function(x) return(sum(x$r_g)/sum(x$r2_g))))
  update_table = data.frame(leaf = as.numeric(names(update)),update = unname(update))
  update_table$update = ifelse(abs(update_table$update) > 1,sign(update_table$update),update_table$update) ## NOTE HERE WE CAP THE UPDATE
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
#' @param tree_beta A gradient tree object fitted for sigma parameter
#' @param tree_gamma A gradient tree object fitted for gamma parameter
#' @param boosting_df the boosting_df
#' @param lambda a vector with the learning rate of sigma and gamma
#' @return A data frame with the udpated parameters sigma and gamma
#' @export
update_parameters <- function(tree_beta,tree_gamma,boosting_df,lambda){
  update_beta = -lambda[1]*predict(tree_beta,boosting_df)
  update_gamma = -lambda[2]*predict(tree_gamma,boosting_df)

  theta <- data.frame(b = boosting_df$b + update_beta ,g = boosting_df$g + update_gamma)
  return(theta)
}


