#' Predict the gradient tree
#'
#' predicting the gradient updates for the parameter
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
