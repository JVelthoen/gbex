library(gbex)
library(grf)
library(tidyverse)

### Simulation parameters
n = 500

## Fixed parameters
tau_extreme = 1 - c(5/n,1/n)
d = 5
depth = c(1,0)
B = 250
min_leaf_size = c(10,10)
sf = 0.75

## Data simulation
X = data.frame(X = matrix(runif(n*d,-1,1),ncol=d))
s = (1+(X[,1]>0))
g = 0.5
y = sim_gpd_data(n,s,g)

### FIT THE GRF
fit_grf = quantile_forest(X,y,mtry=d)

fit_gbex = gbex(y,X,B = B,depth=depth,sf=sf,min_leaf_size = min_leaf_size)

## Define testset
Xtest = lapply(seq(-1,1,length.out=100)[-c(1,100)],function(x1){
  res = data.frame(X = matrix(runif(25*d,-1,1),ncol=d))
  res[,1] = x1
  return(res)
}) %>% bind_rows()
colnames(Xtest) = colnames(X)

## Predict quantile with gbex
#gbex_quant = qnorm(tau_thresh,mean=0,sd=(1+Xtest[,1]>0)) + predict(fit,newdata=Xtest,probs =(1-tau_extreme)/(1-tau_thresh),what="quant")
quantiles_gbex = predict(fit_gbex,newdata=Xtest,probs =tau_extreme,what="quant")

## Predict quantile with grf
quantiles_grf = predict(fit_grf,newdata = Xtest ,quantiles=tau_extreme)

## True quantiles
quantiles_true = sapply(tau_extreme,function(tau){
  g = 0.5
  s = (1+(Xtest[,1]>0))
  s*((1-tau)^(-g)-1)/g
})


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

