## Simulation for the right learning rate size and ratio
library(gbex)

sim_nr = 1
num_folds = 8
reps = 25

## Simulation parameters
n = 1000
d = 2

## Setting the fixed parameters
lambda_scale = 0.025
lambda_ratio = 8
depth = c(1,1)
min_leaf_size = c(10,10)
sf = 0.75
Bmax = 250
ncores = 4

## Simulating data
X = as.data.frame(matrix(runif(n*d,-1,1),ncol=d))
s = 1 + 0.5*X[,1]
g = 0.4 + 0.1*X[,2]
y = sim_gpd_data(n,s,g)

B_CV_strat = B_CV_rand = numeric(reps)
for(i in 1:reps){
  B_CV_strat[i] = CV_gbex(y,X,num_folds=num_folds,Bmax = Bmax,
                       stratified = T,ncores = ncores,
                       lambda_scale = lambda_scale, lambda_ratio=lambda_ratio,
                       depth=depth, min_leaf_size = min_leaf_size,sf=sf,
                       gamma_positive=F,silent=T)$par_CV
  B_CV_rand[i] = CV_gbex(y,X,num_folds=num_folds,Bmax = Bmax,
                         stratified = F,ncores = ncores,
                         lambda_scale = lambda_scale, lambda_ratio=lambda_ratio,
                         depth=depth, min_leaf_size = min_leaf_size,sf=sf,
                         gamma_positive=F,silent=T)$par_CV
  print(i)
}


simulation_result = list(strat= B_CV_strat,rand = B_CV_rand, num_folds=num_folds)

file_save = paste0("sim_8_res_",sim_nr,".rds")
saveRDS(simulation_result,file_save)
