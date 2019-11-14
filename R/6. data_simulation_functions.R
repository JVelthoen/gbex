#' Simulate a dataset
#'
#' @param n number of observations
#' @param d dimension of the covariate space
#' @return the second derivative of the negative log likelihood to gamma for a single observation
#' @export
get_data <- function(n,d){
  # parameter functions
  fun_s=function(x,alpha=1){exp(alpha*mean(x^2))}
  fun_g=function(x,gamma_0=1/3,beta=1/10){gamma_0+beta*mean(x^2)}

  X=array(runif(n*d,-1,1),dim=c(n,d)) %>% as.data.frame() %>% setNames(paste0("X",1:d))

  # parameter vectors
  s=apply(X,1,fun_s)
  g=apply(X,1,fun_g)

  # generating GPD random observations
  y= s*(runif(n)^(-g)-1)/g

  return(list(y=y,X=X,s=s,g=g))
}

#' Simulate a dataset
#'
#' @param n number of observations
#' @param d dimension of the covariate space
#' @return the second derivative of the negative log likelihood to gamma for a single observation
#' @export
get_data2 <- function(n,d){
  # parameter functions
  fun_s=function(x,alpha=1){alpha*(1+ifelse(abs(max(x))<0.5,max(x),sign(max(x))*0.5))}
  fun_g=function(x,gamma_0=1/3,beta=1/10){gamma_0+beta*(1+ifelse(abs(max(x))<0.5,max(x),sign(max(x))*0.5))}

  X=array(runif(n*d,-1,1),dim=c(n,d)) %>% as.data.frame() %>% setNames(paste0("X",1:d))

  # parameter vectors
  s=apply(X,1,fun_s)
  g=apply(X,1,fun_g)

  # generating GPD random observations
  y= s*(runif(n)^(-g)-1)/g

  return(list(y=y,X=X,s=s,g=g))
}
