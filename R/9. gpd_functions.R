#' Quantile function GPD
#'
#' @param p A probability level
#' @param sigma a vector of sigma values
#' @param gamma a vecotr of gamma values =
#' @return quantile of the GPD
#' @export
qgpd=function(p,sigma,gamma){
  quant = sigma/gamma * ((1-p)^(-gamma) -1)
  return(quant)
}

#' Distribution function GPD
#'
#' @param x A probability level
#' @param sigma a vector of sigma values
#' @param gamma a vecotr of gamma values =
#' @return distribution of the GPD
#' @export
pgpd=function(x,sigma,gamma){
  prob = 1- (1+gamma*x/sigma)^(-1/gamma)
  return(prob)
}

#' density function GPD
#'
#' @param x A probability level
#' @param sigma a vector of sigma values
#' @param gamma a vecotr of gamma values =
#' @return density of the GPD
#' @export
dgpd=function(x,sigma,gamma){
  density = 1/sigma*(1+gamma*x/sigma)^(-1/gamma +1)
  return(density)
}

#' quantile loss for gpd estimates
#'
#' @param sigma a vector of sigma values
#' @param gamma a vecotr of gamma values
#' @param y a vector of observations
#' @param tau a probability level for the quantile loss
#' @return quantile loss for the tau probability level
#' @export
gpd_quant_loss = function(sigma,gamma,y,tau){
  quant = qgpd(tau,sigma,gamma)
  u = y-quant
  quant_loss = u*(tau - (u<0))
  return(quant_loss)
}

#' crps loss for gpd estimates
#'
#' @param sigma a vector of sigma values
#' @param gamma a vecotr of gamma values
#' @param y a vector of observations
#' @return crps loss for the gpd
#' @export
gpd_crps_loss = function(sigma,gamma,y){
  crps = scoringRules::crps_gpd(y,gamma,location=0,scale=sigma)
  return(crps)
}
