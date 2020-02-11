rm(list=ls())
library(gbex);library(parallel)

sim_nr <- 1

### Data generating process ###
n <- 1000
d <- 2

results <- mclapply(1:4,function(i){
  data <- get_data(n,d)
  s <- data$s # the sigma parameter
  g <- data$g # the gamma parameter
  y <- data$y # the response vector
  X <- data.frame(X=data$X) # the covariates

  result <- data.frame()
  for(num_folds in c(3,6,9)){
    #### WITHOUT CV ON TUNING PARAMETERS
    Bmax <-  150 # The number of gradient boosting steps
    lambda_ratio <- 15# The learning rate for the sigma and gamma parameter
    lambda_size <- 0.01#
    lambda <- lambda_size*c(1,1/lambda_ratio)
    depth=c(2,2) # The maximum depth for a single sigma and gamma tree
    min_leaf_size=c(5,5) # The minimum number of observations contained in the leafnodes
    sf=0.75 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)

    B <- sapply(c(1,2), function(i){find_optimal_B(y=y,X=X,num_folds = num_folds,
           Bmax = 150,lambda=lambda,depth=depth,
           min_leaf_size=min_leaf_size,sf=sf)})

    result <- rbind(result,data.frame(Bmax= max(B),Bmin=min(B),folds=num_folds))
  }

  return(result)
},mc.cores = 2)

results <- do.call('rbind',results)

saveRDS(results,paste0("sim_res_",sim_nr,".rds"))
