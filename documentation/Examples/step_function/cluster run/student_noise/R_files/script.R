rm(list=ls())
library(gbex)
library(grf)
library(tidyverse)
library(parallel)

sim_nr = 1

### Simulation parameters
n = 1000
tau_thresh = 0.75
df = 4

## Fixed parameters
m = 100
tau_extreme = c(0.99,0.995,0.999)
d = 40
depth = c(1,0)
B = 150
min_leaf_size = c(10,10)
sf = 0.75


output = list()
for(i in 1:m){
  ## Data simulation
  X = data.frame(X = matrix(runif(n*d,-1,1),ncol=d))
  y = rt(n,df=df)*(1+(X[,1]>0))

  ### FIT THE GRF
  fit_grf = grf::quantile_forest(X,y,mtry=d,num.threads = 1)

  ### Compute the OOB threshold important
  threshold = predict(fit_grf,quantiles=tau_thresh)

  ## Compute data for gbex and fit the gbex model
  Z = y-threshold
  y_gbex = Z[Z>0]
  X_gbex = X[Z>0,]

  fit_gbex = gbex(y_gbex,X_gbex,B = B,depth=depth,sf=sf,min_leaf_size = min_leaf_size,silent=T)

  ## Define testset
  Xtest = lapply(c(-0.5,0.5),function(x1){
    res = data.frame(X = matrix(runif(250*d,-1,1),ncol=d))
    res[,1] = x1
    return(res)
  }) %>% bind_rows()
  colnames(Xtest) = colnames(X)

  ## Predict quantile with gbex
  threshold_test = predict(fit_grf,newdata=Xtest,quantiles=tau_thresh)
  quantiles_gbex = apply(predict(fit_gbex,newdata=Xtest,probs =1-(1-tau_extreme)/(1-tau_thresh),what="quant"),2,
                         function(q){
                           q + threshold_test
                         })

  ## Predict quantile with grf
  quantiles_grf = predict(fit_grf,newdata = Xtest ,quantiles=tau_extreme)

  ## True quantiles
  quantiles_true = sapply(tau_extreme,function(tau){ qt(tau,df=df)*ifelse(Xtest[,1]>0,2,1)})


  sim_data = data.frame(X=Xtest[,1],
                        gbex_MSE=(quantiles_gbex-quantiles_true)^2,
                        grf_MSE = (quantiles_grf-quantiles_true)^2,
                        gbex_bias = (quantiles_gbex - quantiles_true),
                        grf_bias = (quantiles_grf-quantiles_true)) %>%
    group_by(X) %>%
    summarize_all(mean)
  print(i)
  output[[length(output) + 1]] <- sim_data
}

result = output %>%
  do.call("rbind",.) %>%
  group_by(X) %>%
  summarize_all(mean)


saveRDS(result,paste0("sim_res_",sim_nr,".rds"))
