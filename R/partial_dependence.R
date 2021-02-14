#' Partial dependence function
#'
#' @param object A gbex object
#' @param var_name character variable names for which to calculate the partial dependence function (maximum number of variables is 2)
#' @param grid_points numeric number of equally spaced grid points for a single dimension
#' @param grid named dataframe with grid on which to compute partial dependence (default NULL)
#' @return A list with for each variable the partial dependence of sigma and gamma with a
#' numeric value corresponding to the variable.
#' @details For a two dimensional partial dependence function the number of grid points used is grid_points^2.
#' @export
partial_dependence <- function(object,var_name,grid_points = 50,grid = NULL){
  if(is.null(object$data)){
    stop("Data is not saved into gbex object")
  }
  if(!is.null(grid)){
    if(!all(var_name %in% names(grid))){
      stop("grid supplied is not complete for specified var_names")
    }
  }

  if(length(var_name) == 1){
    if(is.null(grid)){
      values = seq(min(object$data[[var_name]]),max(object$data[[var_name]]),length.out=grid_points)
    } else{
      values = grid[[var_name]]
    }

    theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)
    PD_per_tree_sigma = sapply(object$trees_sigma,PD_tree,var_name=var_name,values=values)
    PD_per_tree_gamma = sapply(object$trees_gamma,PD_tree,var_name=var_name,values=values)
    PD_transformed = data.frame(st = theta_init$st - apply(PD_per_tree_sigma,1,sum)*object$lambda[1],
                                gt = theta_init$gt -apply(PD_per_tree_gamma,1,sum)*object$lambda[2])
    PD = transform_parameters(PD_transformed,object$gamma_positive,inverse_transform=F)
    PD[[var_name]] = values
  } else if(length(var_name) == 2){
    if(is.null(grid)){
      values = expand.grid(seq(min(object$data[[var_name[1]]]),max(object$data[[var_name[[1]]]]),length.out = grid_points),
                           seq(min(object$data[[var_name[2]]]),max(object$data[[var_name[[2]]]]),length.out = grid_points))
      colnames(values) = var_name
    } else{
      values = grid[var_name]
    }

    theta_init = transform_parameters(object$theta_init,object$gamma_positive,inverse_transform=T)
    PD_per_tree_sigma = sapply(object$trees_sigma,PD_tree2,var_names=var_name,values=values)
    PD_per_tree_gamma = sapply(object$trees_gamma,PD_tree2,var_names=var_name,values=values)
    PD_transformed = data.frame(st = theta_init$st - apply(PD_per_tree_sigma,1,sum)*object$lambda[1],
                                gt = theta_init$gt -apply(PD_per_tree_gamma,1,sum)*object$lambda[2])
    PD = transform_parameters(PD_transformed,object$gamma_positive,inverse_transform=F)
    PD = cbind(PD,values)
  } else{
    stop("partial dependence can only be calculated for one or two variables")
  }
  output = list(var_name = var_name,
                PD_df = PD,
                dimensions = length(var_name))
  class(output) = "PD_gbex"
  return(output)
}


#' Partial dependence function for quantiles
#'
#' @param object A gbex object
#' @param var_name character variable names for which to calculate the partial dependence function (maximum number of variables is 2)
#' @param grid_points numeric number of equally spaced grid points for a single dimension
#' @param grid named dataframe with grid on which to compute partial dependence (default NULL)
#' @return A list with for each variable the partial dependence of sigma and gamma with a
#' numeric value corresponding to the variable.
#' @details For a two dimensional partial dependence function the number of grid points used is grid_points^2.
#' @export
partial_dependence_quantile <- function(object,tau,var_name,grid_points = 50,grid = NULL){
  if(is.null(object$data)){
    stop("Data is not saved into gbex object")
  }
  if(!is.null(grid)){
    if(!all(var_name %in% names(grid))){
      stop("grid supplied is not complete for specified var_names")
    }
  }

  if(length(var_name) == 1){
    if(is.null(grid)){
      values = seq(min(object$data[[var_name]]),max(object$data[[var_name]]),length.out=grid_points)
    } else{
      values = grid[[var_name]]
    }

  quant_PD = numeric(length(values))
  for( icnt   in 1:length(values)){
    value = values[icnt]
    df_value = object$data[-1]
    df_value[var_name] = value
    quant = predict(object,df_value,what="quant",probs=tau)
    quant_PD[icnt] = mean(quant)
  }

  quant_PD = data.frame(Q = quant_PD)
  quant_PD[var_name] = values
  } else if(length(var_name) == 2){
    if(is.null(grid)){
      values = expand.grid(seq(min(object$data[[var_name[1]]]),max(object$data[[var_name[[1]]]]),length.out = grid_points),
                           seq(min(object$data[[var_name[2]]]),max(object$data[[var_name[[2]]]]),length.out = grid_points))
      colnames(values) = var_name
    } else{
      values = grid[var_name]
    }

    quant_PD = numeric(nrow(values))
    for( icnt in 1:nrow(values)){
      df_value = object$data[-1]
      df_value[var_name] = values[icnt,]
      quant = predict(object,df_value,what="quant",probs=tau)
      quant_PD[icnt] = mean(quant)
    }
    quant_PD = data.frame(Q = quant_PD)
    quant_PD = cbind(quant_PD,values)

  } else{
    stop("partial dependence can only be calculated for one or two variables")
  }
  output = list(var_name = var_name,
              PD_df = quant_PD,
              tau = tau)
  class(output) = "PD_gbex_quant"
  return(output)
}

