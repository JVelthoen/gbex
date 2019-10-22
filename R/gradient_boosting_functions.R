#' Estimate scale and shape for the GPD using gradient boosting
#'
#' @param y Response variables (vector of length n)
#' @param X Covariate matrix (matrix of dimension (n x d))
#' @param lambda learning rate for the scale and shape parameter
#' @param B Number of gradient boosting steps
#' @param depth Maximum depth of the trees
#' @param sf  sample fraction used for fitting the trees
#' @param constant_g indicate whether gamma should be taken constant
#' @return A list with the estimates of scale and shape the deviance for each optimaization step a list of trees and a data frame with all parameters of the final optimization step.
#' @export
gbex <- function(y,X,B=250,lambda=c(0.01,0.005),depth=c(2,2),min_leaf_size=c(10,10),sf=0.2){
  if(!is.data.frame(X)) X = data.frame(X=X)
  # Estimate the unconditional tail first and set the estimates as the first guess
  theta_init <- first_guess(y)

  # Create a data.frame used for the boosting procedure with data, parameters, first and second derivatives
  DF_boost <- get_DF_boost(cbind(y,X),theta_init)

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
    tree_s= rpart::rpart(r_s~X.1+X.2,data=DF_tree, method='anova',control=ctrl_s)
    leaf_values_s <- get_leaf_values(tree_s)
    TREE_s <- list(tree=tree_s,values = leaf_values_s)

    # fit a tree for gamma the final TREE object is the tree itself and the new values
    ctrl_g=rpart::rpart.control(maxdepth = depth[2], minsplit=2, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = min_leaf_size[2])
    tree_g= rpart::rpart(r_g~X.1+X.2,data=DF_tree, method='anova',control=ctrl_g)
    tree_values_g <- get_leaf_values(tree_g)
    TREE_g <- list(tree=tree_g,values = tree_values_g)

    # Update the parameter vectors by one gradient step new gradient step
    theta_hat <- update_parameter(TREE_s,TREE_g,DF_boost,lambda)
    #par(mfrow=c(1,3))
    #hist(theta_hat[,1])
    #hist(theta_hat[,2])

    # create a new data.frame with updated derivatives
    DF_boost <- get_DF_boost(DF_boost,theta_hat)
    #plot(dev)

    # Save the estimated trees and the deviance
    TREES[[b]] <- list(tree_s = TREE_s, tree_g=TREE_g)
    dev[b+1] <- mean(DF_boost$dev)
  }

  output <- list(s_hat= DF_boost$s_hat, g_hat = DF_boost$g_hat,
                 dev=dev,
                 TREES=TREES, theta_init= theta_init[1,],
                 DF_boost = DF_boost,
                 lambda=lambda,B=B,depth=depth)
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

#' Compute the deviance at each iteration for a given dataset
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
      gamma = data.frame(leaf = treeClust::rpart.predict.leaves(tree$tree,X)) %>%
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
