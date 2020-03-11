rm(list=ls())
setwd("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/stratified sampling/cluster run/results/folds_8_d_2")
library(tidyverse)


results = readRDS("all_results.rds")

results = list()
files = list.files(full.names=T)
results_new = lapply(files,readRDS)
results = list(results,results_new) %>% unlist(recursive=F)
saveRDS(results,"all_results.rds")

mean_strat = results %>%
  map("strat") %>%
  map(mean) %>%
  unlist()

mean_rand = results %>%
  map("rand") %>%
  map(mean) %>%
  unlist()

sd_strat = results %>%
  map("strat") %>%
  map(sd) %>%
  unlist()

sd_rand = results %>%
  map("rand") %>%
  map(sd) %>%
  unlist()

data.frame(strat = mean_strat,rand=mean_rand) %>%
  pivot_longer(1:2,names_to = "sampling",values_to = "B") %>%
  ggplot(aes(x=B,fill=sampling)) +
  geom_histogram(bins=10,position="identity",alpha = 0.3)

data.frame(strat = sd_strat,rand=sd_rand) %>%
  pivot_longer(1:2,names_to = "sampling",values_to = "B") %>%
  ggplot(aes(x=B,fill=sampling)) +
  geom_histogram(bins=10,position="identity",alpha = 0.3)



