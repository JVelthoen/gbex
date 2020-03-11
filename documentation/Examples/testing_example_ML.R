## example script
rm(list=ls())
library(gbex)
set.seed(30071993)

### Data generating process ###
n <- 1000
d <- 2
data <- get_data(n,d)
s <- data$s # the sigma parameter
g <- data$g # the gamma parameter
y <- data$y # the response vector
X <- cbind(data.frame(X=data$X)) # the covariates
colnames(X) = paste0("X",1:ncol(X))
#### WITHOUT CV ON TUNING PARAMETERS
B <- 250 # The number of gradient boosting steps
lambda_ratio <- 15# The learning rate for the sigma and gamma parameter
lambda_size <- 0.01#
lambda <- lambda_size*c(1,1/lambda_ratio)
depth=c(2,2) # The maximum depth for a single sigma and gamma tree
min_leaf_size=c(5,5) # The minimum number of observations contained in the leafnodes
sf=0.75 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)


## ESTIMATING THE MODEL AND PRINTING OUR SIMPLE FIGURES
fit <- gbex(y,X,B=B,lambda=lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf,gamma_positive = T)

print(fit)
plot(fit)
partial_dependence(fit,variable = 1)
partial_dependence(fit,variable = 2)
variable_importance(fit,type="permutation")
variable_importance(fit,type="relative")


theta <- predict(fit,newdata=X)

x_bar=apply(X^2, 1,mean)
par(mfrow=c(1,3))
plot(x_bar,s,pch='.', col='red',ylim=c(min(c(s,theta$s)),max(c(s,theta$s))))
points(x_bar,theta$s,pch='.',col='blue')
plot(x_bar,g,pch='.', col='red',ylim=c(min(c(g,theta$g)),max(c(g,theta$g))))
points(x_bar,theta$g,pch='.',col='blue')
plot(fit$dev,xlab="iteration",ylab="deviance",type="l")


## TESTING THE PREDICT FUNCTION
y_test <- s*(runif(n)^(-g)-1)/g


## TESTING THE DEVIANCE PER STEP FUNCTION
dev <- dev_per_step(fit)
dev_test <- dev_per_step(fit,X=X,y=y_test)

minimum_iteration <- (0:B)[which(dev_test==min(dev_test))]
par(mfrow=c(1,1))
plot(0:B,dev,type="l",xlab="iteration 0,..,B",col="blue",ylim=c(min(dev,dev_test),max(dev,dev_test)))
lines(0:B,dev_test,col="red")
abline(v=minimum_iteration,lty=2,lwd=2)

## Optimal number of trees to minimize test sample deviance
Blim = minimum_iteration-1
theta <- predict(fit,newdata=X,Blim=Blim)

x_bar=apply(X^2, 1,mean)
par(mfrow=c(1,2))
plot(x_bar,s,pch='.', col='red',ylim=c(min(c(s,theta$s)),max(c(s,theta$s))))
points(x_bar,theta$s,pch='.',col='blue')
plot(x_bar,g,pch='.', col='red',ylim=c(min(c(g,theta$g)),max(c(g,theta$g))))
points(x_bar,theta$g,pch='.',col='blue')

#####################################################################
################### CROSS VALIDATION TESTS ##########################
#####################################################################

##### CV for B ####
n <- 1000
d <- 2
data <- get_data(n,d)
s <- data$s # the sigma parameter
g <- data$g # the gamma parameter
y <- data$y # the response vector
X <- data.frame(X=data$X) # the covariates

lambda = c(0.01,0.0006666666)
min_leaf_size = c(10,10)
sf = 0.75
depth = c(2,2)

num_folds = 4

CV_fit = CV_gbex(y,X,num_folds,300,stratified=T,ncores = 8,lambda=lambda,min_leaf_size=min_leaf_size,sf=sf,depth=depth)

print(CV_fit)

plot(CV_fit,what="general")
plot(CV_fit,what="per_fold")
plot(CV_fit,what="hist_data")


##### CV for lambda_ratio ####
par_grid = c(2,5,10)
par_name = "lambda_ratio"
lambda_scale = 0.01
num_folds = 3

CV_fit = CV_gbex(y,X,num_folds,300,
                 par_name=par_name, par_grid=par_grid, stratified=T,
                 min_leaf_size= min_leaf_size,sf=sf,depth=depth)

print(CV_fit)

plot(CV_fit,what="general")
plot(CV_fit,what="hist_data")




