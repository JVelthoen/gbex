#' Divide observations randomly over folds
#'
#' Function assigns randomly assigns each observation to a fold for cross validation
#'
#' @param n numeric number of observations
#' @param num_folds Integer number of folds used for cross validation
#' @param repeat_cv Integer number of repetitions of cross validation
#' @param stratified Boolean indicating wheter sampling should be done in a stratified way
#' @return A vector with integers corresponding to the fold
#' @details The stratified sampling creates stratified clusters of size num_folds from the ordered observations y.
#' For each fold one observations is sampled.
#' @export
divide_in_folds <- function(y,num_folds,repeat_cv, stratified = F){
  n = length(y)
  ordered_y = order(y)
  all_folds = list()
  for(rep in 1:repeat_cv){
    if(stratified){
      folds_vector <- sapply(1:ceiling(n/num_folds),function(i){sample(1:num_folds)})[1:n]
      folds <- lapply(1:num_folds,function(fold){ ordered_y[folds_vector == fold]})
    } else{
      index_shuffled = sample(1:n)
      folds = unname(split(sample(1:n),cut(1:n,breaks=num_folds,labels=F)))
    }
    all_folds[[length(all_folds)+1]] <- folds
  }
  all_folds = unlist(all_folds,recursive=F)
  return(all_folds)
}
