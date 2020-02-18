library(gbex)
library(grf)
library(quantregForest)
library(tidyverse)


n = 2000
d = 40
X = data.frame(X = matrix(runif(n*d,-1,1),ncol=d))
y = rnorm(n)*(1+(X[,1]>0))

fit_grf = grf::quantile_forest(X,y,mtry=d)

Xtest = lapply(seq(-1,1,length.out=100)[-c(1,100)],function(x1){
  res = data.frame(X = matrix(runif(100*d,-1,1),ncol=d))
  res[,1] = x1
  return(res)
}) %>% bind_rows()
colnames(Xtest) = colnames(X)

tau_thresh = 0.9
threshold = qnorm(tau_thresh,mean=0,sd=(1+X[,1]>0))
Z = y-threshold
y_gbex = Z[Z>0]
X_gbex = X[Z>0,]
fit = CV_gbex(y_gbex,X_gbex,num_folds = 8,Bmax = 300,stratified=T,depth=c(1,0),min_leaf_size = c(10,10),sf=0.75)




probabilities = c(0.9,0.999)
names(probabilities) = paste0("q",1:length(probabilities))
quantiles_grf = predict(fit_grf,newdata = Xtest ,quantiles=probabilities)
colnames(quantiles_grf) = names(probabilities)

data_grf = cbind(data.frame(X=Xtest[,1]),quantiles_grf)  %>%
  pivot_longer(2:(length(probabilities)+1),names_to = "probs",values_to="qu") %>%
  mutate(probs = recode(probs,!!!probabilities)) %>%
  group_by(X,probs) %>%
  summarize(quant= mean(qu),quant_sd= sd(qu)) %>%
  mutate(quant_min = quant - quant_sd,quant_max = quant + quant_sd)

data_true = data_qrf %>%
  mutate(quant = qnorm(probabilities,mean=0,sd=ifelse(X>0,2,1)))

g = ggplot(data_grf,aes(x=X,y=quant,group=probs)) +
  geom_ribbon(aes(ymin=quant_min,ymax=quant_max),alpha=0.25,fill="red") +
  geom_line(col="red",lwd=2) +
  geom_line(data=data_true,aes(x=X,y=quant),col="black") +
  theme_minimal()

print(g)
#g + ggsave("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/step_function.png/step_fun.png")

