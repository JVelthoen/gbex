## example script
library(gbex)
set.seed(25051979)

### Data generating process ###
n <- 1000
d <- 2
X = data.frame(X1 = runif(n),X2=runif(n))
s = X$X1*0.25 + 0.5
g = rep(0.25,n)
y <- sim_gpd_data(n,s,g)

B <-  250 # The number of gradient boosting steps
lambda=c(0.01,0.0025) # The learning rate for the sigma and gamma parameter
depth=c(2,0) # The maximum depth for a single sigma and gamma tree
min_leaf_size=c(30,30) # The minimum number of observations contained in the leafnodes
sf=0.75 # The subsample fraction used for estimating each tree (in a single step only one subsample is drawn and used for both the sigma and gamma tree)

fit <- gbex(y,X,B=B,lambda=lambda,depth=depth,min_leaf_size=min_leaf_size,sf=sf)

par(mfrow=c(1,3))
plot(X$X1,s,pch='.', col='red',ylim=c(min(c(s,fit$theta$s)),max(c(s,fit$theta$s))))
points(X$X1,fit$theta$s,pch='.',col='blue')
plot(X$X1,g,pch='.', col='red',ylim=c(min(c(g,fit$theta$g)),max(c(g,fit$theta$g))))
points(X$X1,fit$theta$g,pch='.',col='blue')
plot(fit$dev,xlab="iteration",ylab="deviance")

partial_dependence(fit,"X1")


y_test <- s*(runif(n)^(-g)-1)/g
theta_hat <- predict(fit,newdata=X)
s_hat_test <- theta_hat$s
g_hat_test <- theta_hat$g

dev <- dev_per_step(fit)
dev_test <- dev_per_step(fit,X=X,y=y_test)
minimum_iteration <- (0:B)[which(dev_test==min(dev_test))]
par(mfrow=c(1,1))
plot(0:B,dev,type="l",xlab="iteration 0,..,B",col="blue",ylim=c(min(dev,dev_test),max(dev,dev_test)))
lines(0:B,dev_test,col="red")
abline(v=minimum_iteration,lty=2,lwd=2)



data = data.frame(x1 = c(-1,-1,1,1),
                  x2 = c(-1,1,-1,1),
                  y = c(-3,-1,1,3))

## TESTS FOR SINGLE SPLIT
ctrl = rpart::rpart.control(maxdepth = 1, minsplit=1, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = 1)
tree = rpart::rpart(y~.,data=data, method='anova',control=ctrl)
update_table = data.frame(leaf = c(2,3),
                          update = c(-2,2))
tree = list(tree = tree, update_table = update_table)
class(tree) = "gradient_tree"

PD_tree(tree,"x1",c(-1,1)) ## SHOULD BE c(-2,2)
PD_tree(tree,"x2",c(-1,1)) ## SHOULD BE c(0,0)

## TESTS For Double split
ctrl = rpart::rpart.control(maxdepth = 2, minsplit=1, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = 1)
tree = rpart::rpart(y~.,data=data, method='anova',control=ctrl)
update_table = data.frame(leaf = c(3,4,6,7),
                          update = c(-3,-1,1,3))
tree = list(tree = tree, update_table = update_table)
class(tree) = "gradient_tree"

PD_tree(tree,"x1",c(-1,1)) ## SHOULD BE c(-2,2)
PD_tree(tree,"x2",c(-1,1)) ## SHOULD BE c(-1,1)
