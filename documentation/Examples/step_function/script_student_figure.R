rm(list=ls())
library(gbex)
library(grf)
library(tidyverse)

sim_nr = 1

### Simulation parameters
n = 1000
tau_thresh = 0.75
df = 4

## Fixed parameters
tau_extreme = c(0.99,0.995,0.999)
d = 40
depth = c(1,0)
B = 150
min_leaf_size = c(10,10)
sf = 0.8
alpha = 0

## Data simulation
X = data.frame(X = matrix(runif(n*d,-1,1),ncol=d))
y = rt(n,df=df)*(1+(X[,1]>0))

### FIT THE GRF
fit_grf = grf::quantile_forest(X,y,mtry=d)

### Compute the OOB threshold important
threshold = predict(fit_grf,quantiles=tau_thresh)
#threshold = qt(tau_thresh,df=df)*(1+(X[,1]>0))

## Compute data for gbex and fit the gbex model
Z = y-threshold
y_gbex = Z[Z>0]
X_gbex = X[Z>0,]


fit_gbex = gbex(y_gbex,X_gbex,B = B,depth=depth,sf=sf,min_leaf_size = min_leaf_size,alpha=alpha)

## Define testset
Xtest = lapply(seq(-1,1,length.out=100)[-c(1,100)],function(x1){
  res = data.frame(X = matrix(runif(25*d,-1,1),ncol=d))
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


plot_data = data.frame(X=Xtest[,1],gbex=quantiles_gbex,grf = quantiles_grf,true = quantiles_true)  %>%
  pivot_longer(2:(length(tau_extreme)*3+1),names_to="method",values_to="quantile") %>%
  mutate(probs = tau_extreme[as.numeric(substr(method,nchar(method),nchar(method)))],
         method = substr(method,1,nchar(method)-2)) %>%
  group_by(X,probs,method) %>%
  summarize(quant= mean(quantile),
            quant_min= quantile(quantile,probs=0.25),
            quant_max= quantile(quantile,probs=0.75))

## Make a figure
ggplot(plot_data,aes(x=X,y=quant,col=method)) +
  geom_line() +
  geom_ribbon(data=plot_data,aes(x=X,ymin=quant_min,ymax=quant_max,fill=method),alpha=0.25)+
  facet_grid(vars(probs)) +
  theme_minimal()


#g + ggsave("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/step_function.png/step_fun.png")

