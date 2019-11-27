optimize_tuning_parameters <- function(X,y,grid_start,grid_search,lambda_high,lambda_low,num_folds){
  parameter_names <- names(grid_start)
  if(!dplyr::setequal(parameter_names,c("sf","min_leaf_size","depth"))){
    warning("Not all tuning parameters are specified. Beware that default values are used instead")
  }

  par_fixed <- grid_start
  if(length(grid_search) > 0){
    par_fixed$lambda <- lambda_high
  } else{
    par_fixed$lambda <- lambda_low
  }

  par_fixed <- optimize_parameter(X,y,par_fixed,num_folds = num_folds)
  cat(paste("Optimal B for high learning rate by cross validation:",par_fixed$B,"\n"))

  if(!is.null(grid_search$depth)){
    par_fixed <- optimize_parameter(X,y,par_fixed,num_folds = num_folds,par="depth",grid_values=grid_search$depth)
    cat(paste("Optimal depth by cross validation:","(",paste(par_fixed$depth,collapse=","),")","\n"))
  }

  if(!is.null(grid_search$min_leaf_size)){
    par_fixed <- optimize_parameter(X,y,par_fixed,num_folds = num_folds,par="min_leaf_size",grid_values=grid_search$min_leaf_size)
    cat(paste("Optimal min_leaf_size by cross validation:","(",paste(par_fixed$min_leaf_size,collapse=","),")","\n"))
  }

  if(!is.null(grid_search$sf)){
    par_fixed <- optimize_parameter(X,y,par_fixed,num_folds = num_folds,par="sf",grid_values=grid_search$sf)
    cat(paste("Optimal sf by cross validation:",par_fixed$sf,")","\n"))
  }

  if(length(grid_search) > 0){
    par_fixed$lambda <- lambda_low
    par_fixed$B <- grid_start$B
    par_fixed <- optimize_parameter(X,y,par_fixed,num_folds = num_folds)
    cat(paste("Optimal B for low learning rate by cross validation:",par_fixed$B,"\n"))
  }

  return(par_fixed)
}

optimize_parameter <- function(X,y,par_fixed,num_folds = 4,par=NULL,grid_values=NULL){
  if(!is.null(par)) par_fixed[par] <- NULL

  opt_par <- CV_grid_search(X,y,par,grid_values,num_folds,par_fixed) %>%
    extract2("par_opt")

  if(all(!is.null(par),!is.null(grid_values))){
    opt_par_set <- par_fixed
    opt_par_set[[par]] <- opt_par
  } else{
    opt_par_set <- par_fixed
    opt_par_set$B <- opt_par
  }

  return(opt_par_set)
}
