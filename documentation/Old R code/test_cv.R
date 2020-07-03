library(gbex)
set.seed(1234)

n = 2000
d = 2
df = 4
reps = 5

## data parameters
df0 = 4
a= 1
b=1

## gbex parameters
tau_c = 0.8
sf = 1
min_leaf_size = rep(0.1*n*(1-tau_c),2)
lambda_scale = 0.01
lambda_ratio = 7
B = 25
depth = c(2,2)

### Data generating process ###
n <- 2000
d <- 10
X = data.frame(matrix(runif(n*d,-1,1),ncol=d))
Xbar = apply(X[,1:2]^2,1,mean)
sigma = exp(a*Xbar)
df = df0 - b*Xbar
y = sigma*rt(n,df=df)
colnames(X) = paste0("X",1:ncol(X))


fit = gbex(y,X,B=B,lambda=lambda,depth=c(2,2),min_leaf_size = min_leaf_size,sf =sf,gamma_positive=T)

plot(fit)



