#' Likelihood functions for the generalized pareto density
#'
#' @param v vector with first entry the scale second entry shape and the third entry the observation
#' @return the negative log likelihood for a single observation
#' @export
GP_dev=function(v){
  (1+1/v[2])*log(1+v[2]*v[3]/v[1])+log(v[1])
}

#' @rdname GP_dev
GP_dev_diff_s=function(v){
  (1-((1+v[2])*v[3])/(v[1]+v[2]*v[3]))/v[1]
}

#' @rdname GP_dev
GP_dev_diff_g=function(v){
  -log(1+v[2]*v[3]/v[1])*v[2]^(-2)+(1+1/v[2])*v[3]/(v[1]+v[2]*v[3])
}

#' @rdname GP_dev
GP_dev_diff2_s=function(v){
  ((v[1]+v[2]*v[3])*v[1])^(-1)*(v[3]/v[1]+(v[3]-v[1])/(v[1]+v[2]*v[3]))
}

#' @rdname GP_dev
GP_dev_diff2_g=function(v){
  2*v[2]^(-3)*log(1+v[2]*v[3]/v[1])-2*v[3]*v[2]^(-2)/(v[1]+v[2]*v[3])-(1+1/v[2])*v[3]^2/(v[1]+v[2]*v[3])^2
}
