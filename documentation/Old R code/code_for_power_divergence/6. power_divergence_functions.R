#' The density power divergence
#'
#' @param v vector with first entry the scale second entry shape and the third entry the observation
#' @param alpha the power for the power divergence
#' @return Power divergence evaluated in the point v
#' @export
PD_dev = function(v,A,B,alpha=0.5){
  Int(c(v[1],v[2],alpha)) - (1 + 1/alpha)*GP_density(v,A,B)^alpha
}

#' The density power divergence differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape and the third entry the observation
#' @param alpha the power for the power divergence
#' @return the derivative evaluated in point v
#' @export
PD_dev_diff_s = function(v,A,B,alpha=0.5){
  Int_diff_s(c(v[1],v[2],alpha)) - (1+alpha)*GP_density(v,A,B)^(alpha -1)*GP_density_diff_s(v,A,B)
}

#' The density power divergence differentiated to gamma
#'
#' @param v vector with first entry the scale second entry shape and the third entry the observation
#' @param alpha the power for the power divergence
#' @return the derivative in point v
#' @export
PD_dev_diff_g = function(v,A,B,alpha=0.5){
  Int_diff_g(c(v[1],v[2],alpha)) - (1+alpha) * GP_density(v,A,B)^(alpha - 1) * GP_density_diff_g(v,A,B)
}

#' The density power divergence differentiated two times to sigma
#'
#' @param v vector with first entry the scale second entry shape and the third entry the observation
#' @param alpha the power for the power divergence
#' @return the derivative in point v
#' @export
PD_dev_diff2_s = function(v,A,B,alpha=0.5){
  Int_diff2_s(c(v[1],v[2],alpha)) -
    (alpha^2 -1)*GP_density(v,A,B)^(alpha - 2)*GP_density_diff_s(v,A,B)^2 -
    (1+alpha)*GP_density(v,A,B)^(alpha - 1)*GP_density_diff2_s(v,A,B)
}

#' The density power divergence differentiated two times to gamma
#'
#' @param v vector with first entry the scale second entry shape and the third entry the observation
#' @param alpha the power for the power divergence
#' @return the derivative in the point v
#' @export
PD_dev_diff2_g = function(v,A,B,alpha=0.5){
  Int_diff2_g(c(v[1],v[2],alpha)) -
    (alpha^2 -1)*GP_density(v,A,B)^(alpha -2)*(GP_density_diff_g(v,A,B))^2 -
    (1+alpha)*GP_density(v,A,B)^(alpha -1)*GP_density_diff2_g(v,A,B)
}

#' The integral in the first part of empirical power divergence
#'
#' @param v vector with first entry the scale second entry shape and the third parameter the alpha
#' @return the integral for parameter vector v
#' @export
Int = function(v){
  1/(v[1]^v[3] * (1 + v[3]*(1 + v[2])))
}

#' The integral in the first part of empirical power divergence differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape and the third parameter the alpha
#' @return the derivative in point v
#' @export
Int_diff_s = function(v){
  -v[3]/(v[1]^(v[3] + 1) * (1 + v[3]*(1 + v[2])))
}

#' The integral in the first part of empirical power divergence two times differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape and the third parameter the alpha
#' @return the derivative in point v
#' @export
Int_diff2_s = function(v){
  v[3]*(v[3] + 1)/(v[1]^(v[3] + 2) * (1 + v[3]*(1 + v[2])))
}

#' The integral in the first part of empirical power divergence differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape and the third parameter the alpha
#' @return the derivative in point v
#' @export
Int_diff_g = function(v){
  -v[3]/(v[1]^v[3] * (1 + v[3]*(1 + v[2]))^2)
}

#' The integral in the first part of empirical power divergence two times differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape and the third parameter the alpha
#' @return the derivative in point v
#' @export
Int_diff2_g = function(v){
  2*v[3]^2 / (v[1]^v[3] * (1 + v[3]*(1 + v[2]))^3)
}

#' The A function in the split of the GP density
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the function evaluated at v
#' @export
A_func = function(v){
  1/v[1] * (1 + v[2]*v[3]/v[1])^(-1)
}

#' The B function in the split of the GP density
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the function evaluated at v
#' @export
B_func = function(v){
  (1 + v[2]*v[3]/v[1])^(-1/v[2])
}

#' The A function in the split of the GP density differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
A_diff_s = function(v,A){
  A = A_func(v)
  -A/v[1] + v[2]*v[3]/v[1] * A^2
}

#' The A function in the split of the GP density differentiated two times to sigma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
A_diff2_s = function(v,A){
  2*A/v[1]^2 * (1-v[2]*v[3]*A)^2
}

#' The A function in the split of the GP density differentiated to gamma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
A_diff_g = function(v,A){
  -v[3]*A^2
}

#' The A function in the split of the GP density differentiated two times to gamma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
A_diff2_g = function(v,A){
  2 * v[3]^2 * A^3
}

#' The B function in the split of the GP density differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
B_diff_s = function(v,A,B){
 v[3]/v[1] * A * B
}

#' The B function in the split of the GP density differentiated two times to sigma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
B_diff2_s = function(v,A,B){
  v[3]*A*B/v[1]^2 * (v[3]*A*(v[2] + 1) - 2)
}

#' The B function in the split of the GP density differentiated to gamma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
B_diff_g = function(v,A,B){
  -1/v[2] * B * log(B) - v[3]/v[2] * A * B
}

#' The B function in the split of the GP density differentiated two times to gamma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the derivative in the given point v
#' @export
B_diff2_g = function(v,A,B){
  B/v[2]^2 *( (log(B) + v[3]*A)^2 + 2* (log(B) + v[3]*A) ) + v[3]^2/v[2] * A^2 * B
}

#' The density of the generalized pareto distribution
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the density in point v
#' @export
GP_density = function(v,A,B){
  A*B
}

#' The density of the generalized pareto distribution differentiated to sigma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the density in point v
#' @export
GP_density_diff_s = function(v,A,B){
  A_diff_s(v,A)*B + A*B_diff_s(v,A,B)
}

#' The density of the generalized pareto distribution differentiated two times to sigma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the density in point v
#' @export
GP_density_diff2_s = function(v,A,B){
  A_diff2_s(v,A)*B +2*A_diff_s(v,A)*B_diff_s(v,A,B) + B_diff2_s(v,A,B)*A
}

#' The density of the generalized pareto distribution differentiated to gamma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the density in point v
#' @export
GP_density_diff_g = function(v,A,B){
  A_diff_g(v,A)*B + A*B_diff_g(v,A,B)
}

#' The density of the generalized pareto distribution differentiated two times to gamma
#'
#' @param v vector with first entry the scale second entry shape third parameter the observations
#' @return the density in point v
#' @export
GP_density_diff2_g = function(v,A,B){
  A_diff2_g(v,A)*B +2*A_diff_g(v,A)*B_diff_g(v,A,B) + B_diff2_g(v,A,B)*A
}



