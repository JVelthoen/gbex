#' Plotting function for CV_par object
#'
#' @param object par_CV object
#' @param what string taking "general" average over all folds, "per_fold" for a figure for each fold or "hist_data" for histograms of the different folds.
#' @return plots the cross validation results for the CV_gbex object
#' @export
plot.CV_gbex <- function(object,what = "general"){
  grid_B = object$grid_B
  if(what == "general"){
    if(object$par == "B"){
      data = data.frame(loss=object$loss_all, B = object$grid_B)
      g = ggplot2::ggplot(data,ggplot2::aes(y=loss,x=B)) +
        ggplot2::geom_line(size=1) +
        ggplot2::geom_vline(xintercept = which(object$loss_all == min(object$loss_all)) + 1, lty=2,size=0.8) +
        ggplot2::labs(title="CV loss", x = "Iteration", y = "Loss") +
        ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5))
    } else if(object$par %in% c("lambda","lambda_ratio","lambda_scale","depth","min_leaf_size","sf","alpha")){
      data = data.frame(B=grid_B,loss=object$loss_all)
      data = reshape(data, direction = "long", varying = colnames(data)[-1],
                     v.names = "loss", timevar = "par_name")
      grid_names = unlist(lapply(object$grid,
                                 function(x){
                                   if(length(x) == 1) return(as.character(x))
                                   else return(paste0("(",paste0(x,collapse=","),")"))
                                               }
                                 ))
      data$par_name = factor(grid_names[data$par_name],levels=grid_names)
      data_vline = data.frame(par_name = factor(grid_names[which(sapply(object$grid,
                                                           function(grid_value){
                                                             all(grid_value == object$par_CV)
                                                           }))],levels=grid_names),
                              int = object$B_opt)
      g = ggplot2::ggplot(data,ggplot2::aes(y=loss,x=B,lty=par_name)) +
        ggplot2::geom_line(size=1) +
        ggplot2::geom_vline(data = data_vline, ggplot2::aes(xintercept = int), linetype = 2,size=0.8) +
        ggplot2::labs(title=paste("CV loss per",object$par), x = "Iteration", y = "Loss",lty=object$par) +
        ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5))
    }
  } else if(what == "per_fold"){
    if(object$par == "B"){
      loss_temp = apply(object$loss_folds,2,function(loss){ loss-loss[1]})
      data = data.frame(B = grid_B,loss_temp)
      data = reshape(data, direction = "long", varying = colnames(data)[-1],
              v.names = "loss", timevar = "fold")
      data$fold = paste("Fold",data$fold)
      data_all = data.frame(loss=object$loss_all-object$loss_all[1], B = object$grid_B,fold="all")
      data_abline = data.frame(fold = unique(data$fold),
                               int = apply(object$loss_folds,2,function(loss){which(loss== min(loss)) + 1}))
      g = ggplot2::ggplot(data,ggplot2::aes(y=loss,x=id,group = fold)) +
        ggplot2::geom_line(size=1,lty=2) +
        ggplot2::geom_line(data=data_all,ggplot2::aes(y=loss,x=B),size=1,lty=1) +
        ggplot2::labs(title="CV loss every fold", x = "Iteration", y = "Loss") +
        ggplot2::geom_vline(xintercept = which(object$loss_all == min(object$loss_all))+1,size=0.8) +
        ggplot2::theme_minimal() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       panel.border = ggplot2::element_rect(fill=NA))
    } else{
      stop("this type of plot can not be made for this variable")
    }
  } else if(what == "hist_data"){
    data = data.frame(fold = paste("Fold",object$folds),y=object$y)
    g = ggplot2::ggplot(data,ggplot2::aes(x=y)) +
      ggplot2::geom_histogram(binwidth = diff(range(data$y))/20) +
      ggplot2::labs(title="Observation histogram", x = "Observations", y = "Frequency") +
      ggplot2::facet_wrap(ggplot2::vars(fold)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.border = ggplot2::element_rect(fill=NA))
  } else{
    stop("This plot is not specified for a CV_gbex object.")
  }
  return(g)
}


#' printing function for CV_gbex object
#'
#' @param object par_CV object
#' @return prints a basic description of the CV_gbex object
#' @export
print.CV_gbex <- function(object){
  cat(paste0(deparse(object$call),"\n"))
  cat(paste0("A cross validation object for gbex for parameter ",object$par,".\n"))
  cat(paste0("The number of folds is ",object$num_folds," which are"))
  if(object$stratified){
    cat(" obtained by stratified sampling.\n")
  } else{
    cat(" obtained by random sampling.\n")
  }
  if(object$par %in% c("B","sf","lambda_ratio")){
    cat(paste0("The optimal parameter value is ",object$par," = ",object$par_CV,".\n"))
  } else{
    cat(paste0("The optimal parameter value is ",object$par," = (",paste(object$par_CV,collapse=","),").\n"))
  }
}
