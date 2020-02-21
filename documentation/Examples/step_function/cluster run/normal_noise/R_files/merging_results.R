rm(list=ls())
library(tidyverse)
setwd("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/step_function/cluster run/normal_noise/results/")


file_names = list.files(full.names = T)
sim_res = lapply(file_names,readRDS)
n = c(rep(1000,10),rep(2000,10))
thresh = rep(c(rep(0.8,5),rep(0.9,5)),2)

index_1000_80 = which(n==1000 & thresh == 0.8)
index_1000_90 = which(n==1000 & thresh == 0.9)
index_2000_80 = which(n==2000 & thresh == 0.8)
index_2000_90 = which(n==2000 & thresh == 0.9)

table_1000_80 = Reduce("+",sim_res[index_1000_80])/5
table_1000_90 = Reduce("+",sim_res[index_1000_90])/5
table_2000_80 = Reduce("+",sim_res[index_2000_80])/5
table_2000_90 = Reduce("+",sim_res[index_2000_90])/5

## making latex output
temp = table_2000_90[-1] %>% t()

latex_matrix = as.matrix(cbind(cbind(temp[1:3,1],temp[7:9,1],temp[4:6,1],temp[10:12,1]),cbind(temp[1:3,2],temp[7:9,2],temp[4:6,2],temp[10:12,2])))
rownames(latex_matrix) = NULL
latex_matrix = round(latex_matrix,3)
print_statement = apply(latex_matrix,1,paste0,collapse = " & ")
for(i in print_statement) cat(paste0(i,"\n"))
