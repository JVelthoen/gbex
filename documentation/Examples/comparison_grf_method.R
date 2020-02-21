library(gbm)
library(grf)

n = 2000
d = 40
taus = c(0.7,0.8,0.9)

X = data.frame(X = matrix(runif(n*d,-1,1),ncol=d))
y = rnorm(n)*(1+(X[,1]>0))

fit = quantile_forest(X,y,mtry=d)
fit70 = quantile_forest(X,y,mtry=d,quantiles = 0.7)
fit80 = quantile_forest(X,y,mtry=d,quantiles = 0.8)
fit90 = quantile_forest(X,y,mtry=d,quantiles = 0.9)

tau = 0.7
MSE_all = mean((as.vector(predict(fit,quantiles=tau)) - qnorm(tau,0,(1+(X[,1]>0))))^2)
MSE_spe = mean((as.vector(predict(fit70,quantiles=tau)) - qnorm(tau,0,(1+(X[,1]>0))))^2)
print(c(MSE_all,MSE_spe))

tau = 0.8
quant80_all = mean((as.vector(predict(fit,quantiles=tau)) - qnorm(tau,0,(1+(X[,1]>0))))^2)
quant80_spe = mean((as.vector(predict(fit80,quantiles=tau)) - qnorm(tau,0,(1+(X[,1]>0))))^2)
print(c(quant80_all,quant80_spe))

tau = 0.9
quant90_all = mean((as.vector(predict(fit,quantiles=tau)) - qnorm(tau,0,(1+(X[,1]>0))))^2)
quant90_spe = mean((as.vector(predict(fit90,quantiles=tau)) - qnorm(tau,0,(1+(X[,1]>0))))^2)
print(c(quant90_all,quant90_spe))


