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
X <- data.frame(X=data$X) # the covariates

#### WITHOUT CV ON TUNING PARAMETERS
B <-  250 # The number of gradient boosting steps
lambda_ratio <- 15# The learning rate for the sigma and gamma parameter
lambda_size <- 0.01#
lambda <- lambda_size*c(1,1/lambda_ratio)
depth=c(2,2) # The maximum depth for a single sigma and gamma tree
min_leaf_size=c(5,5) # The minimum number of observations contained in the leafnodes
sf=0.75 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)


## ESTIMATING THE MODEL AND PRINTING OUR SIMPLE FIGURES
fit <- gbex(y,X,B=B,lambda=lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf)

print(fit)

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


## Estimation of B with cross validation
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

depth_list <- list(c(1,1),c(2,2))
sf_vec <- c()
num_folds = 12
CV_fit = CV_gbex(y,X,num_folds,"B",250,stratified=T,depth=depth,lambda=lambda,sf=sf,min_leaf_size=min_leaf_size)

plot(CV_fit,what="all")
plot(CV_fit,what="folds")
plot(CV_fit,what="data")


initial_dev = CV_fit$dev_all[1]
dev = CV_fit$dev_all - initial_dev
dev_fold = apply(CV_fit$dev_folds,2,function(dev_temp){dev_temp - dev_temp[1] })

grid_B = CV_fit$grid_B
layout(c(1),widths = c(5,1))
par(mai=rep(0.5, 4))
plot(grid_B,dev,type="l",lwd=2,
     xlab="B",ylab="dev",ylim=range(dev_fold))
abline(v=grid_B[dev==min(dev)],lty = "dashed",lwd=3)
for(fold in 1:num_folds){
  lines(grid_B,dev_fold[,fold],type="l",lwd=2,lty=3)

}

object <- CV_fit
folds = object$folds
y = object$y

num_folds = object$num_folds
nrow = floor(sqrt(num_folds))
ncol = ceiling(num_folds/nrow)
dev = object$dev_folds

layout(matrix(c(1:num_folds,rep(0,nrow*ncol - num_folds)), ncol=ncol), widths=c(rep(1,ncol)))
par(mai=rep(0.5, 4))
for(fold in 1:num_folds){
  plot(X[folds==fold,1],X[folds==fold,2],main=paste("Fold",fold), cex=0.7,pch=16,
       xlab="y",ylab="counts")
  points(X[(folds==fold)&(y>10),1],X[(folds==fold)&(y>10),2],main=paste("Fold",fold), cex=1.2,pch=16,col="red",
       xlab="y",ylab="counts")
  abline(h=c(0.5,-0.5),v=c(0.5,-0.5))
}

