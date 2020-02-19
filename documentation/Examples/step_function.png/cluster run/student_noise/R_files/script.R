library(gbex)
library(grf)
library(tidyverse)
library(parallel)

sim_nr = 1

### Simulation parameters
n = 1000
tau_thresh = 0.80
df = 4

## Fixed parameters
m = 25
tau_extreme = c(0.99,0.999,0.9999)
d = 40
depth = c(1,0)
B = 150
min_leaf_size = c(10,10)
sf = 0.75


output = mclapply(1:m,function(i){
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

  fit_gbex = gbex(y_gbex,X_gbex,B = B,depth=c(1,0),sf=0.75,min_leaf_size = c(10,10))

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


  sim_data = data.frame(X=Xtest[,1],gbex=quantiles_gbex,grf = quantiles_grf,true = quantiles_true)  %>%
    pivot_longer(2:(length(tau_extreme)*2+1),names_to="method",values_to="quantile") %>%
    mutate(probs = tau_extreme[as.numeric(substr(method,nchar(method),nchar(method)))],
           method = substr(method,1,nchar(method)-2),
           truth = ifelse(probs == tau_extreme[1],true.1,ifelse(probs == tau_extreme[2],true.2,true.3))) %>%
    group_by(X,probs,method) %>%
    summarize(MSE = mean((quantile-truth)^2),bias=mean(quantile-truth))

  return(sim_data)
},mc.cores = 2)

output = do.call('rbind',output)

saveRDS(output,paste0("sim_res_",sim_nr,".rds"))
