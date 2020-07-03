## example script
library(gbex)
set.seed(25051979)

### Data generating process ###
n <- 1000
d <- 5
data <- get_data(n,2)
s <- data$s # the sigma parameter
g <- data$g # the gamma parameter
y <- data$y # the response vector
X <- data.frame(X=data$X,N=matrix(runif(n*(d-2),-1,1),ncol=d-2)) # the covariates
colnames(X) = paste0("X",1:ncol(X))

B <-  100 # The number of gradient boosting steps
lambda=c(0.01,0.0025) # The learning rate for the sigma and gamma parameter
depth=c(2,2) # The maximum depth for a single sigma and gamma tree
min_leaf_size=c(30,30) # The minimum number of observations contained in the leafnodes
sf=1 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)


fit = gbex(y,X,B=B,lambda=lambda,depth=c(1,1),min_leaf_size = min_leaf_size,sf =sf,gamma_positive=T)




variable_importance(fit)

par(mfrow=c(1,2))
sigma_variables= fit$trees_sigma %>% map("tree") %>% map("splits") %>%
  map(rownames) %>% unlist() %>% table()
barplot(sigma_variables)
gamma_variables= fit$trees_gamma %>% map("tree") %>% map("splits") %>%
  map(rownames) %>% unlist() %>% table()
barplot(gamma_variables)

sigma_variables= fit$trees_sigma %>% map("tree") %>% map("splits") %>%
  map(rownames) %>% unlist()
paths = sapply(paste0("X",1:ncol(X)),function(par){cumsum(sigma_variables == par)})
data.frame(paths,index=1:B) %>% pivot_longer(1:5,names_to="var",values_to = "chosen") %>%
  ggplot(aes(x=index,y=chosen,colour=var)) +
  geom_line()

gamma_variables= fit$trees_gamma %>% map("tree") %>% map("splits") %>%
  map(rownames) %>% unlist()
paths = sapply(paste0("X",1:ncol(X)),function(par){cumsum(gamma_variables == par)})
data.frame(paths,index=1:B) %>% pivot_longer(1:5,names_to="var",values_to = "chosen") %>%
  ggplot(aes(x=index,y=chosen,colour=var)) +
  geom_line()

barplot(sigma_variables)
gamma_variables= fit$trees_gamma %>% map("tree") %>% map("splits") %>%
  map(rownames) %>% unlist() %>% table()
barplot(gamma_variables)

plot(sigma_p)
arguments = list(lambda=lambda,
                 depth=depth,
                 min_leaf_size = min_leaf_size,
                 sf=sf)

CV_fit = CV_gbex(y,X,num_folds = 3,Bmax = 200,loss_func = "likelihood", tau_loss = 0.5,ncores = 4,
                 lambda = lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf)

plot(CV_fit)

CV_fit_quant = CV_gbex(y,X,num_folds = 4,Bmax = 200,loss_func = "quant_loss", tau_loss = 0.9,ncores = 4,
                       lambda = lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf)

plot(CV_fit_quant)

CV_fit_crps = CV_gbex(y,X,num_folds = 4,Bmax = 200,loss_func = "crps", tau_loss = 0.5,ncores = 4,
                       lambda = lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf)

plot(CV_fit_crps)
plot(CV_fit,what="per_fold")

