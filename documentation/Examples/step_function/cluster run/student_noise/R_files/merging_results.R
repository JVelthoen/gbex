rm(list=ls())
library(tidyverse)
setwd("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/step_function/cluster run/student_noise/results/")


file_names = list.files(full.names = T)
sim_res = lapply(file_names,readRDS)
n = c(rep(1000,15),rep(2000,15))
thresh = rep(c(rep(0.75,5),rep(0.8,5),rep(0.85,5)),2)

index_1000_75 = which(n==1000 & thresh == 0.75)
index_1000_85 = which(n==1000 & thresh == 0.85)
index_2000_75 = which(n==2000 & thresh == 0.75)
index_2000_85 = which(n==2000 & thresh == 0.85)

table_1000_75 = Reduce("+",sim_res[index_1000_75])/5
table_1000_85 = Reduce("+",sim_res[index_1000_85])/5
table_2000_75 = Reduce("+",sim_res[index_2000_75])/5
table_2000_85 = Reduce("+",sim_res[index_2000_85])/5

## making latex output
temp = table_2000_85[-1] %>% t()

latex_matrix = as.matrix(cbind(cbind(temp[1:3,1],temp[7:9,1],temp[4:6,1],temp[10:12,1]),cbind(temp[1:3,2],temp[7:9,2],temp[4:6,2],temp[10:12,2])))
rownames(latex_matrix) = NULL
latex_matrix = round(latex_matrix,3)
print_statement = apply(latex_matrix,1,paste0,collapse = " & ")
for(i in print_statement) cat(paste0(i,"\n"))
