#' Estimate the unconditional tail as a first guess
#'
#' @param y A set of observations
#' @param newdata A data frame with covariates for which to predict the sigma and gamma parameter
#' @return A data.frame object with the estimated sigma and gamma parameters
#' @export
first_guess <- function(y){
  GP_fit=POT::fitgpd(y,threshold=0)
  theta_init=GP_fit$fitted.values
  theta_init=matrix(rep(theta_init,length(y)),ncol=2,byrow=T)
  colnames(theta_init)=c('s','g')
  return(theta_init)
}

#' For a gradient tree fill the leafnodes with the gradient update gamma (CURRENTLY NOT USED)
#'
#' @param tree a gradient tree
#' @return A data.frame object with for each leaf the gradient update
#' @export
get_leaf_values <- function(tree){
  gamma <- predict(tree)
  leaf <- treeClust::rpart.predict.leaves(tree)
  update_df <- unique(data.frame(gamma = gamma,leaf = leaf))
  return(update_df)
}

#' Fill the leafnodes of the sigma tree using linesearch
#'
#' @param tree a gradient tree
#' @param DF_boost the boosting dataframe used for obtaining the tree
#' @return A data.frame object with for each leaf the gradient update
#' @export
fill_leafs_s <- function(tree,DF_tree){
  update_df <- data.frame(leaf = tree$where, r = DF_tree$r_s, r2 = DF_tree$r2_s) %>%
      dplyr::group_by(leaf) %>% dplyr::summarise(r=sum(r),r2=sum(r2)) %>%
      dplyr::mutate(gamma = r/r2) %>%
      dplyr::select(gamma,leaf) %>%
      dplyr::mutate(gamma = ifelse(abs(gamma)<1,gamma,gamma/abs(gamma)))
  return(update_df)
}

#' Fill the leafnodes of the gamma tree using linesearch
#'
#' @param tree a gradient tree
#' @param DF_tree the boosting dataframe used for obtaining the tree
#' @return A data.frame object with for each leaf the gradient update
#' @export
fill_leafs_g <- function(tree,DF_tree){
  update_df <- data.frame(leaf = tree$where, r = DF_tree$r_g, r2 = DF_tree$r2_g) %>%
    dplyr::group_by(leaf) %>% dplyr::summarise(r=mean(r),r2=mean(r2)) %>%
    dplyr::mutate(gamma = r/r2) %>%
    dplyr::select(gamma,leaf) %>%
    dplyr::mutate(gamma = ifelse(abs(gamma)<1,gamma,gamma/abs(gamma)))
  return(update_df)
}

#' Update the current parameters by one gradient step
#'
#' @param TREE_s a list of gradient tree for sigma together with a df with gradient step per leafnode
#' @param TREE_s a list of gradient tree for sigma together with a df with gradient step per leafnode
#' @param DF_boost the current boosting data frame
#' @param lambda a vector of length two with learning rate for sigma and gamma
#' @return A data.frame object with for each leaf the gradient update
#' @export
update_parameter <- function(TREE_s,TREE_g,DF_boost,lambda){
  gamma_s <- data.frame(leaf = treeClust::rpart.predict.leaves(TREE_s$tree,newdata=DF_boost)) %>%
    dplyr::left_join(TREE_s$values,by="leaf") %>% .$gamma

  gamma_g <- data.frame(leaf = treeClust::rpart.predict.leaves(TREE_g$tree,newdata=DF_boost)) %>%
    dplyr::left_join(TREE_g$values,by="leaf") %>% .$gamma

  theta <- data.frame(s=DF_boost$s_hat -  lambda[1]*gamma_s ,g= DF_boost$g_hat - lambda[2]*gamma_g)
  return(theta)
}

#' Create or update the DF_boost data.frame
#'
#' @param DF A data frame with the response and covariates or when updating the current DF_boost data frame
#' @param theta the new parameter updates for sigma and gamma
#'
#' @return A data.frame object with the data the parameter estimates the negative likelihood and the first and second derivatives of the negative likelihood towards sigma and gamma
#' @export
get_DF_boost <- function(DF,theta,alpha=0){
  data <- DF[ colnames(DF)[!(colnames(DF) %in% c("s_hat","g_hat","dev","r_s","r_g","r2_s","r2_g"))]]
  s_hat=theta[,'s']
  g_hat=theta[,'g']
  if(alpha==0){
    dev=apply(cbind(s_hat,g_hat,DF$y),1,GP_dev)
    r_s=apply(cbind(s_hat,g_hat,DF$y),1,GP_dev_diff_s)
    r_g=apply(cbind(s_hat,g_hat,DF$y),1,GP_dev_diff_g)
    r2_s=apply(cbind(s_hat,g_hat,DF$y),1,GP_dev_diff2_s)
    r2_g=apply(cbind(s_hat,g_hat,DF$y),1,GP_dev_diff2_g)
    DF_dev = data.frame(s_hat=s_hat, g_hat=g_hat, dev=dev, r_s=r_s, r_g=r_g, r2_s=r2_s, r2_g=r2_g)
  } else{
    A <- apply(cbind(s_hat,g_hat,DF$y),1,A_func)
    B <- apply(cbind(s_hat,g_hat,DF$y),1,B_func)
    dev=apply(cbind(s_hat,g_hat,DF$y,A,B),1,function(x) PD_dev(x[1:3],x[4],x[5],alpha))
    r_s=apply(cbind(s_hat,g_hat,DF$y,A,B),1,function(x) PD_dev_diff_s(x[1:3],x[4],x[5],alpha))
    r_g=apply(cbind(s_hat,g_hat,DF$y,A,B),1,function(x) PD_dev_diff_g(x[1:3],x[4],x[5],alpha))
    r2_s=apply(cbind(s_hat,g_hat,DF$y,A,B),1,function(x) PD_dev_diff2_s(x[1:3],x[4],x[5],alpha))
    r2_g=apply(cbind(s_hat,g_hat,DF$y,A,B),1,function(x) PD_dev_diff2_g(x[1:3],x[4],x[5],alpha))
  }
  DF_boost=data.frame(data, s_hat=s_hat, g_hat=g_hat, dev=dev, r_s=r_s, r_g=r_g, r2_s=r2_s, r2_g=r2_g)
  return(DF_boost)
}
