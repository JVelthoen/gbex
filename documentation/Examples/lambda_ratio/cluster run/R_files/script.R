## Simulation for the right learning rate size and ratio
library(gbex)

sim_nr = 1

## Simulation parameters
n = 1000
d = 10

## Setting the fixed parameters
lambda_scale = 0.025
depth = c(1,1)
min_leaf_size = c(10,10)
sf = 0.75
num_folds = 4
Bmax = 500
ncores = 2


## Setting grid for lambda ratio
grid_lambda_ratio = c(1,2,3,5,8,12,16,20)

## Simulating data
X = as.data.frame(matrix(runif(n*d,-1,1),ncol=d))

## Definition model 1
s1 = 1 + 0.5*X[,1]
g1 = 0.4
y_model1 = sim_gpd_data(n,s1,g1)

CV_fit_model1 = CV_gbex(y_model1,X,num_folds=num_folds,Bmax = Bmax,
                        par_name = "lambda_ratio",par_grid = grid_lambda_ratio,
                        stratified = T,ncores = ncores,
                        lambda_scale = lambda_scale,depth=depth,
                        min_leaf_size = min_leaf_size,sf=sf,gamma_positive=T,silent=T)

CV_fit_model1$X = CV_fit_model1$y = NULL

## Definition model 2
s2 = 1 + 0.5*X[,1]
g2 = 0.4 + 0.1*X[,1]
y_model2 = sim_gpd_data(n,s2,g2)

CV_fit_model2 = CV_gbex(y_model2,X,num_folds=num_folds,Bmax = Bmax,
        par_name = "lambda_ratio",par_grid = grid_lambda_ratio,
        stratified = T,ncores = ncores,
        lambda_scale = lambda_scale,depth=depth,
        min_leaf_size = min_leaf_size,sf=sf,gamma_positive=T,silent=T)

CV_fit_model2$X = CV_fit_model2$y = NULL

## Definition model 3
s3 = 1 + 0.5*X[,1]
g3 = 0.4 + 0.2*X[,1]
y_model3 = sim_gpd_data(n,s3,g3)

CV_fit_model3 = CV_gbex(y_model3,X,num_folds=num_folds,Bmax = Bmax,
                        par_name = "lambda_ratio",par_grid = grid_lambda_ratio,
                        stratified = T,ncores = ncores,
                        lambda_scale = lambda_scale,depth=depth,
                        min_leaf_size = min_leaf_size,sf=sf,gamma_positive=T,silent=T)

CV_fit_model3$X = CV_fit_model3$y = NULL

simulation_result = list(
  X = X,
  model1 = list(y=y_model1,CV_fit = CV_fit_model1,s=s1,g=g1),
  model2 = list(y=y_model2,CV_fit = CV_fit_model2,s=s2,g=g2),
  model3 = list(y=y_model3,CV_fit = CV_fit_model3,s=s3,g=g3),
  lambda_scale = lambda_scale,
  depth = depth,
  min_leaf_size = min_leaf_size,
  sf = sf,
  num_folds = num_folds,
  Bmax = Bmax,
  ncores = ncores
)

file_save = paste0("sim_res_",sim_nr,".rds")
saveRDS(simulation_result,file_save)





