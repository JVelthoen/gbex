#' Divide observations randomly over folds
#'
#' Function assigns randomly assigns each observation to a fold for cross validation
#'
#' @param n numeric number of observations
#' @param num_folds Integer number of folds used for cross validation
#' @param stratified Boolean indicating wheter sampling should be done in a stratified way
#' @return A vector with integers corresponding to the fold
#' @details The stratified sampling creates stratified clusters of size num_folds from the ordered observations y.
#' For each fold one observations is sampled.
#' @export
divide_in_folds <- function(y,num_folds,stratified = F){
  n = length(y)
  if(stratified){
    folds_matrix <- sapply(1:ceiling(n/num_folds),function(i){sample(1:num_folds)})
    folds_vector <- folds_matrix[1:n]
    folds <- folds_vector[rank(-y)]
  } else{
    index_shuffled = sample(1:n)
    folds = cut(seq(1,length(index_shuffled)),breaks=num_folds,labels=F)[order(index_shuffled)]
  }
  return(folds)
}
