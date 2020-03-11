setwd("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/lambda_ratio/cluster run/results")
library(tidyverse)

files <- list.files("third simulation",full.names=T,recursive=T)


results <- lapply(files,readRDS)

lambda_grid = results[[1]]$model1$CV_fit$grid

lambda_choice_model1 = data.frame(
  lambda = sapply(results,function(res){ res$model1$CV_fit$par_CV }),
  B = sapply(results,function(res){ res$model1$CV_fit$B_opt }),
  model = "constant")
lambda_choice_model2 = data.frame(
  lambda = sapply(results,function(res){ res$model2$CV_fit$par_CV }),
  B = sapply(results,function(res){ res$model2$CV_fit$B_opt }),
  model="slow")
lambda_choice_model3 = data.frame(
  lambda = sapply(results,function(res){ res$model3$CV_fit$par_CV }),
  B = sapply(results,function(res){ res$model3$CV_fit$B_opt }),
  model="quick")

lambda_choice = rbind(lambda_choice_model1,lambda_choice_model2,lambda_choice_model3) %>%
  mutate(index = rep(1:(n()/3),3))

optimal_dev_model1 = lapply(results,function(res){
  apply(res$model1$CV_fit$dev_all,2,min)
}) %>% do.call("rbind",.) %>%
  as.data.frame() %>%
  setNames(paste("lambda",lambda_grid)) %>%
  mutate(index = 1:n(),
         model = "constant")

optimal_dev_model2 = lapply(results,function(res){
  apply(res$model2$CV_fit$dev_all,2,min)
}) %>% do.call("rbind",.) %>%
  as.data.frame() %>%
  setNames(paste("lambda",lambda_grid)) %>%
  mutate(index = 1:n(),
         model= "slow")

optimal_dev_model3 = lapply(results,function(res){
  apply(res$model3$CV_fit$dev_all,2,min)
}) %>% do.call("rbind",.) %>%
  as.data.frame() %>%
  setNames(paste("lambda",lambda_grid)) %>%
  mutate(index = 1:n(),
         model="quick")

optimal_dev = rbind(optimal_dev_model1,optimal_dev_model2,optimal_dev_model3) %>%
  pivot_longer(1:8,names_to = "lambda_ratio",values_to = "dev") %>%
  mutate(lambda_ratio = as.numeric(substr(lambda_ratio,8,nchar(lambda_ratio))))


relative_dev_model1 = lapply(results,function(res){
  x= res$model1$CV_fit$dev_all[1,] -apply(res$model1$CV_fit$dev_all,2,min)
  return(x/max(x)*100)
}) %>% do.call("rbind",.) %>%
  as.data.frame() %>%
  setNames(paste("lambda",lambda_grid)) %>%
  mutate(index = 1:n(),
         model = "constant")

relative_dev_model2 = lapply(results,function(res){
  x= res$model2$CV_fit$dev_all[1,] -apply(res$model2$CV_fit$dev_all,2,min)
  return(x/max(x)*100)
}) %>% do.call("rbind",.) %>%
  as.data.frame() %>%
  setNames(paste("lambda",lambda_grid)) %>%
  mutate(index = 1:n(),
         model = "slow")

relative_dev_model3 = lapply(results,function(res){
  x= res$model3$CV_fit$dev_all[1,] -apply(res$model3$CV_fit$dev_all,2,min)
  return(x/max(x)*100)
}) %>% do.call("rbind",.) %>%
  as.data.frame() %>%
  setNames(paste("lambda",lambda_grid)) %>%
  mutate(index = 1:n(),
         model = "quick")


relative_dev = rbind(relative_dev_model1,relative_dev_model2,relative_dev_model3) %>%
  pivot_longer(1:8,names_to = "lambda_ratio",values_to = "dev") %>%
  mutate(lambda_ratio = as.numeric(substr(lambda_ratio,8,nchar(lambda_ratio))))



simulation_result = list(lambda_choice = lambda_choice,
                         optimal_dev = optimal_dev,
                         relative_dev = relative_dev,
                         lambda_grid = lambda_grid)

saveRDS(simulation_result,"simulation_result.rds")


  