#' Estimate scale and shape for the GPD using gradient boosting
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param lambda learning rate for the scale and shape parameter
#' @param B Number of gradient boosting steps
#' @param depth Maximum depth of the trees
#' @param sf  sample fraction used for fitting the trees
#' @param alpha the power for power divergence (default alpha = 0 meaning maximum likelihood is used)
#' @param silent boolean indicating whether progress during fitting procedure should be printed.
#' @return A list with the estimates of scale and shape the deviance for each optimaization step a list of trees and a data frame with all parameters of the final optimization step.
#' @export
gbex <- function(y,X,B=180,lambda=c(0.025,0.0025),depth=c(2,2),min_leaf_size=c(30,30),sf=0.5,alpha = 0,silent=F){
  n <- length(y)
  if(!is.data.frame(X)) X = data.frame(X=X)
  cov_names <- colnames(X)
  # Estimate the unconditional tail first and set the estimates as the first guess
  theta_init <- first_guess(y)

  # Create a data.frame used for the boosting procedure with data, parameters, first and second derivatives
  DF_boost <- get_DF_boost(cbind(y,X),theta_init,alpha)

  # Save the results of the boosting steps
  # TREES contains the boosting trees
  # dev contains the deviance for each iteration
  TREES <- list()
  dev=rep(mean(DF_boost$dev),B+1)
  for (b in 1:B){
    # Take a subsample from the entire data.frame
    DF_tree <- DF_boost[sample(1:n,sf*n,replace=F),]

    # fit a tree for sigma the final TREE object is the tree itself and the new values
    ctrl_s=rpart::rpart.control(maxdepth = depth[1], minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = min_leaf_size[1])
    formula_s <- as.formula(
      paste("r_s",
            paste(cov_names, collapse = " + "),
            sep = " ~ "))
    tree_s= rpart::rpart(formula_s,data=DF_tree, method='anova',control=ctrl_s)
    leaf_values_s <- fill_leafs_s(tree_s,DF_tree)
    TREE_s <- list(tree=tree_s,values = leaf_values_s)

    # fit a tree for gamma the final TREE object is the tree itself and the new values
    ctrl_g=rpart::rpart.control(maxdepth = depth[2], minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = min_leaf_size[2])
    formula_g <- as.formula(
      paste("r_g",
            paste(cov_names, collapse = " + "),
            sep = " ~ "))
    tree_g= rpart::rpart(formula_g,data=DF_tree, method='anova',control=ctrl_g)
    leaf_values_g <- fill_leafs_g(tree_g,DF_tree)
    TREE_g <- list(tree=tree_g,values = leaf_values_g)

    # Update the parameter vectors by one gradient step new gradient step
    theta_hat <- update_parameter(TREE_s,TREE_g,DF_boost,lambda)
    # create a new data.frame with updated derivatives
    DF_boost <- get_DF_boost(DF_boost,theta_hat,alpha)

    # Save the estimated trees and the deviance
    TREES[[b]] <- list(tree_s = TREE_s, tree_g=TREE_g)
    dev[b+1] <- mean(DF_boost$dev)

    if(!silent & b %in% round(((1:10)*(B/10)))){
      cat(paste0(round(b/B,1)*100,"% of trees fitted\n"))
    }
  }

  output <- list(s_hat= DF_boost$s_hat, g_hat = DF_boost$g_hat,
                 dev=dev,
                 TREES=TREES, theta_init= theta_init[1,],
                 DF_boost = DF_boost,
                 lambda=lambda,B=B,depth=depth,alpha=alpha)
  class(output) <- "gbex"
  return(output)
}


#' Predict function for gbex object
#'
#' @param object A fitted gbex object
#' @param newdata A data frame with covariates for which to predict the sigma and gamma parameter
#' @return A data.frame object with the estimated sigma and gamma parameters
#' @export
predict.gbex <- function(object, newdata = NULL){
  if(is.null(newdata)){
    pred = data.frame(s = object$s_hat,g=object$g_hat)
  } else{
    TREES_s <- lapply(object$TREES,function(trees) trees$tree_s)
    TREES_g <- lapply(object$TREES,function(trees) trees$tree_g)

    updates_s <- lapply(TREES_s,function(tree){
      gamma = data.frame(leaf = treeClust::rpart.predict.leaves(tree$tree,newdata)) %>%
        dplyr::left_join(tree$values,by="leaf")  %>% .$gamma
      return(gamma*object$lambda[1])
    })

    s_hat <- object$theta_init[1] - Reduce('+',updates_s)

    updates_g <- lapply(TREES_g,function(tree){
      gamma = data.frame(leaf = treeClust::rpart.predict.leaves(tree$tree,newdata)) %>%
        dplyr::left_join(tree$values,by="leaf")  %>% .$gamma
      return(gamma*object$lambda[2])
    })

    g_hat <- object$theta_init[2] - Reduce('+',updates_g)

    pred = data.frame(s=s_hat,g=g_hat)
  }
  return(pred)
}

#' Compute the deviance at each iteration for a given dataset (CURRENTLY only working for likelihood)
#'
#' @param object A fitted gbex object
#' @param X A dataframe with the right column names
#' @param y A vector of observations
#' @return A vector with the deviance at iteration
#' @export
dev_per_step <- function(object,X=NULL,y=NULL){
  if(is.null(X)){
    dev = object$dev
  } else{
    TREES_s <- lapply(object$TREES,function(trees) trees$tree_s)
    TREES_g <- lapply(object$TREES,function(trees) trees$tree_g)

    updates_s <- lapply(TREES_s,function(tree){
      gamma = data.frame(leaf = treeClust::rpart.predict.leaves(tree$tree,newdata=X)) %>%
        dplyr::left_join(tree$values,by="leaf")  %>% .$gamma
      return(gamma*object$lambda[1])
    }) %>% do.call('cbind',.) %>% cbind(0,.) %>% apply(1,cumsum) %>% t()

    s_hat_matrix <- object$theta_init[1] - updates_s

    updates_g <- lapply(TREES_g,function(tree){
      gamma = data.frame(leaf = treeClust::rpart.predict.leaves(tree$tree,X)) %>%
        dplyr::left_join(tree$values,by="leaf")  %>% .$gamma
      return(gamma*object$lambda[2])
    }) %>% do.call('cbind',.) %>% cbind(0,.) %>% apply(1,cumsum) %>% t()

    g_hat_matrix <- object$theta_init[2] - updates_g

    dev <- c(sapply(1:ncol(s_hat_matrix),function(i){
      mean(apply(cbind(s_hat_matrix[,i],g_hat_matrix[,i],y),1,GP_dev))
    }))
  }
  return(dev)
}
