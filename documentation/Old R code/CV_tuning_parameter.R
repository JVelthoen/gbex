set.seed(25051979)

### Data generating process ###
n <- 1000
d <- 2
data <- get_data(n,d)
s <- data$s # the sigma parameter
g <- data$g # the gamma parameter
y <- data$y # the response vector
X <- data.frame(X=data$X) # the covariates

CV_parameters <- function(y,X){
  n <- length(y)
  index_shuffled <- sample(1:n)
  folds <- cut(seq(1,length(index_shuffled)),breaks=4,labels=F)[order(index_shuffled)]

  lambda=c(0.01,0.0025) # The learning rate for the sigma and gamma parameter
  min_leaf_size = c(30,30) # The maximum depth for a single sigma and gamma tree
  depth=c(2,2)
  sf = 0.75 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)
  B <- CV_find_B(y,X,folds,B_max=500,lambda=lambda,min_leaf_size=min_leaf_size,sf=sf,depth=depth)
  depth <- CV_find_depth(y,X,folds,depth_s=c(1,2,3,4),lambda=lambda,min_leaf_size=min_leaf_size,sf=sf,B=B)
  min_leaf_size <- CV
}



CV_one <- function(y,X,folds,B_max = 500,depth_s= c(1,2,3,4),depth_g = c(1,2,3),...){
  k <- length(unique(folds))
  depth_matrix <- lapply(depth_g,function(d_g) cbind(depth_s,d_g)) %>% do.call('rbind',.)
  p <- nrow(depth_matrix)

  B_per_depth <- dev_per_depth <- numeric(p)
  for(j in 1:p){
    dev_per_fold <- matrix(nrow=k,ncol=B_max+1)
    cat(paste("Start",j,"of",p,"depth configurations\n"))
    for(i in 1:k){
      test_index <- which(folds==i)
      train_index <- which(folds!=i)
      y_train <- y[train_index]
      X_train <- X[train_index,]
      y_test <- y[test_index]
      X_test <- X[test_index,]

      fit <- gbex(y_train,X_train,B=B_max,depth=depth_matrix[j,],lambda=lambda,min_leaf_size=min_leaf_size,sf=sf)
      dev_per_fold[i,] <- dev_per_step(fit,X=X_test,y=y_test)
      cat(paste("Fold",i,"of",k,"\n"))
    }
    dev <- dev_per_fold %>% apply(2,mean)
    dev_per_depth[j] <- min(dev)
    B_per_depth[j] <- (0:B_max)[which(dev == min(dev))]
  }
  index_optimal <- which(dev_per_depth == min(dev_per_depth))
  depth_opt <- depth_matrix[index_optimal,]
  B_optimal <- B_per_depth[index_optimal]

  return(list(B=B_opt,depth=depth_opt)
}




CV_find_B <- function(y,X,folds,B_max = 500,...){
  dev_per_fold <- matrix(nrow=B_max + 1,ncol=4)
  for(i in 1:4){
    test_index <- index_shuffled[folds==i]
    train_index <- index_shuffled[folds!=i]
    y_train <- y[train_index]
    X_train <- X[train_index,]
    y_test <- y[test_index]
    X_test <- X[test_index,]

    fit <- gbex(y_train,X_train,B=B_max,...)
    dev_per_fold[,i] <- dev_per_step(fit,X=X_test,y=y_test)
    cat(paste0("Optimizing B (",i,"/4)\n"))
  }
  dev <- dev_per_fold %>% apply(1,mean)
  B_opt <- (0:B_max)[which(dev==min(dev))[1]]
  return(B_opt)
}

CV_find_depth <- function(y,X,folds,depth_vec = c(1,2,3,4),...){
  depth_matrix <- cbind(depth_vec,depth_vec)
  dev_per_fold <-matrix(nrow=nrow(depth_matrix),ncol=4)

  for(i in 1:4){
    test_index <- index_shuffled[folds==i]
    train_index <- index_shuffled[folds!=i]
    y_train <- y[train_index]
    X_train <- X[train_index,]
    y_test <- y[test_index]
    X_test <- X[test_index,]
    for(j in 1:nrow(depth_matrix)){
      fit <- gbex(y_train,X_train,depth=depth_matrix[j,],...)
      parameters <- predict(fit,newdata=X_test)
      dev_per_fold[j,i] <- apply(cbind(parameters,y_test),1,GP_dev) %>% mean()
    }
    cat(paste0("Optimizing depth (",i,"/4)\n"))
  }
  dev <- dev_per_fold %>% apply(1,mean)
  depth_opt <- depth_matrix[which(dev==min(dev))[1],]
  return(depth_opt)
}


CV_find_min_leaf_size <- function(y,X,folds,min_leaf_vec = c(20,30,40,50,60,70,80),...){
  min_leaf_matrix <- cbind(min_leaf_vec,min_leaf_vec)
  dev_per_fold <-matrix(nrow=nrow(min_leaf_matrix),ncol=4)

  for(i in 1:4){
    test_index <- index_shuffled[folds==i]
    train_index <- index_shuffled[folds!=i]
    y_train <- y[train_index]
    X_train <- X[train_index,]
    y_test <- y[test_index]
    X_test <- X[test_index,]
    for(j in 1:nrow(min_leaf_matrix)){
      fit <- gbex(y_train,X_train,min_leaf_size = min_leaf_matrix[j,],depth=c(2,2),B=128,sf=sf,lambda=lambda)
      parameters <- predict(fit,newdata=X_test)
      dev_per_fold[j,i] <- apply(cbind(parameters,y_test),1,GP_dev) %>% mean()
    }
    cat(paste0("Optimizing minimum leaf size (",i,"/4)\n"))
  }
  dev <- dev_per_fold %>% apply(1,mean)
  plot(dev)
  min_leaf_opt <- min_leaf_matrix[which(dev==min(dev))[1],]
  return(min_leaf_opt)
}
