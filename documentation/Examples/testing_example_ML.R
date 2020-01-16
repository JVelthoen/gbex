## example script
rm(list=ls())
library(gbex)
set.seed(25051979)

### Data generating process ###
n <- 1000
d <- 2
data <- get_data(n,d)
s <- data$s # the sigma parameter
g <- data$g # the gamma parameter
y <- data$y # the response vector
X <- data.frame(X=data$X) # the covariates


B <-  250 # The number of gradient boosting steps
lambda=c(0.01,0.0025) # The learning rate for the sigma and gamma parameter
depth=c(2,2) # The maximum depth for a single sigma and gamma tree
min_leaf_size=c(30,30) # The minimum number of observations contained in the leafnodes
sf=0.5 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)


## ESTIMATING THE MODEL AND PRINTING OUR SIMPLE FIGURES
fit <- gbex(y,X,B=B,lambda=lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf)
theta <- fit$theta

x_bar=apply(X^2, 1,mean)
par(mfrow=c(1,3))
plot(x_bar,s,pch='.', col='red',ylim=c(min(c(s,theta$s)),max(c(s,theta$s))))
points(x_bar,theta$s,pch='.',col='blue')
plot(x_bar,g,pch='.', col='red',ylim=c(min(c(g,theta$g)),max(c(g,theta$g))))
points(x_bar,theta$g,pch='.',col='blue')
plot(fit$dev,xlab="iteration",ylab="deviance")


## TESTING THE PREDICT FUNCTION
y_test <- s*(runif(n)^(-g)-1)/g
theta_test <- predict(fit,newdata=X)

par(mfrow=c(1,2))
plot(x_bar,s,pch='.', col='red',ylim=c(min(c(s,theta_test$s)),max(c(s,theta_test$s))))
points(x_bar,theta_test$s,pch='.',col='blue')
plot(x_bar,g,pch='.', col='red',ylim=c(min(c(g,theta_test$g)),max(c(g,theta_test$g))))
points(x_bar,theta_test$g,pch='.',col='blue')


## TESTING THE DEVIANCE PER STEP FUNCTION
dev <- dev_per_step(fit)
dev_test <- dev_per_step(fit,X=X,y=y_test)

minimum_iteration <- (0:B)[which(dev_test==min(dev_test))]
par(mfrow=c(1,1))
plot(0:B,dev,type="l",xlab="iteration 0,..,B",col="blue",ylim=c(min(dev,dev_test),max(dev,dev_test)))
lines(0:B,dev_test,col="red")
abline(v=minimum_iteration,lty=2,lwd=2)

