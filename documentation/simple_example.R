## example script
rm(list=ls())
library(grf)
library(gbex)

set.seed(13)

### Data generating process ###
n <- 2000
d <- 3

X = data.frame(matrix(runif(n*d,-1,1),ncol=d))
Xbar = apply(X[,1:2]^2,1,mean)
sigma = exp(Xbar)
df = 4 - Xbar
y = sigma*rt(n,df=4)

## Parameters
tau_c = 0.8

## Fit a quantile forest for the threshold
fit_grf = quantile_forest(X,y)
x_threshold = predict(fit_grf,quantiles = tau_c)

## fit gbex model
Z = y-x_threshold
X_z = X[Z>0,]
Z = Z[Z>0]
fit_gbex = gbex(Z,X_z,B=150,lambda_ratio = 10, lambda_scale = 0.01,
                depth = c(2,2), min_leaf_size = c(25,25), sf = 0.75,
                gamma_positive=T)

## plot insample deviance
print(fit_gbex)
plot(fit_gbex)

## variable importance figure
variable_importance(fit_gbex,type = "relative",scaled = T)
variable_importance(fit_gbex,type = "permutation",scaled = T)

## partial dependence plot
partial_dependence(fit_gbex,"X1")
partial_dependence(fit_gbex,"X2")


## use cross validation for optimal B
CV_gbex = CV_gbex(Z,X_z,num_folds = 5, Bmax = 500,stratified = T,ncores = 2)
B_CV = CV_gbex$par_CV

## plot CV deviance
print(CV_gbex)
plot(CV_gbex)

## plot CV deviance for each fold
plot(CV_gbex, what="per_fold")

