rm(list=ls())
library(tidyverse)
setwd("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/lambda_ratio/cluster run/")


sim_res = readRDS("results/simulation_result.rds")


### chosen lambda
plot_data = sim_res$lambda_choice
ggplot(plot_data,aes(x=lambda)) +
  geom_histogram(bins=8) +
  facet_wrap(vars(model))

### optimal deviance per lambda
plot_data = sim_res$optimal_dev %>%
  group_by(lambda_ratio,model) %>%
  summarise(dev_min = quantile(dev,probs=0.1),
            dev_max = quantile(dev,probs=0.9),
            dev=mean(dev))
ggplot(plot_data,aes(x=lambda_ratio,y=dev)) +
  geom_line() +
  facet_wrap(vars(model))

### optimal deviance per lambda
plot_data = sim_res$relative_dev %>%
  group_by(lambda_ratio,model) %>%
  summarise(dev_min = quantile(dev,probs=0.1),
            dev_max = quantile(dev,probs=0.9),
            dev=mean(dev)) %>%
  mutate(model = factor(model,levels=c("constant","slow","quick")))
ggplot(plot_data,aes(x=lambda_ratio,y=dev)) +
  geom_line() +
  geom_ribbon(aes(ymin=dev_min,ymax=dev_max),alpha=0.3) +
  facet_wrap(vars(model))

### consequence of misspecifying the model
plot_data = sim_res$relative_dev %>%
  group_by(index,model) %>%
  mutate(lambda_chosen = sim_res$lambda_grid[which(dev==100)]) %>%
  ungroup() %>%
  group_by(lambda_chosen,lambda_ratio,model) %>%
  summarise(dev_min = quantile(dev,probs=0.1),
            dev_max = quantile(dev,probs=0.9),
            dev=mean(dev))

ggplot(plot_data,aes(x = lambda_ratio, y = dev)) +
  geom_line(aes(col=model)) +
  geom_ribbon(aes(ymin=dev_min,ymax=dev_max,fill=model),alpha=0.3) +
  facet_wrap(vars(lambda_chosen))







