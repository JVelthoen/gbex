#' Plotting function for CV_par object
#'
#' @param object par_CV object
#' @param what string taking "all" average over all folds, "folds" for a figure for each fold or "data" for histograms of the different folds.
#' @return plots the cross validation results for the CV_gbex object
#' @export
plot.CV_gbex <- function(object,what = "all"){
  par = object$par
  dev = object$dev_all
  grid_B = object$grid_B
  if(what == "all"){
    if(par=="B"){
      layout(c(1),widths = c(5,1))
      par(mai=rep(0.5, 4))
      plot(grid_B,dev,type="l",lwd=2,
           xlab="B",ylab="dev",main=paste("Grid Search",par))
      abline(v=grid_B[dev==min(dev)],lty = "dashed",lwd=3)
    } else{
      layout(matrix(c(1,2), ncol=2, byrow=TRUE),widths = c(5,1))
      par(mai=rep(0.5, 4))
      plot(grid_B,grid_B,type="n",lwd=2,
           xlab="B",ylab="dev",main=paste("Grid Search",par),
           ylim=range(dev))
      for(par_value in 1:ncol(dev)){
        lines(grid_B,dev[,par_value],lty=par_value,lwd=2)
        abline(v=grid_B[dev[,par_value] == min(dev[,par_value])],lty = par_value,lwd=3)
      }
      par(mai=c(0,0,0,0))
      plot.new()
      grid = object$grid
      if(par %in% c("depth","min_leaf_size")){
        legend_text = sapply(grid,function(par){paste0("(",par[1],",",par[2],")")})
      } else{
        legend_text = as.character(grid)
      }
      legend(x="center", ncol=1,legend=legend_text,
             lty=1:length(legend_text), title=par,bty="n",lwd=2)
    }
  } else if(what=="folds"){
    num_folds = object$num_folds
    nrow = floor(sqrt(num_folds))
    ncol = ceiling(num_folds/nrow)
    dev = object$dev_folds

    if(par=="B"){
      dev = object$dev_all - object$dev_all[1]
      dev_fold = apply(object$dev_folds,2,function(dev_temp){dev_temp - dev_temp[1] })
      grid_B = CV_fit$grid_B

      layout(c(1),widths = c(5,1))
      par(mai=rep(0.5, 4))
      plot(grid_B,dev,type="l",lwd=2,
           xlab="B",ylab="dev",ylim=range(dev_fold))
      abline(v=grid_B[dev==min(dev)],lty = "dashed",lwd=3)
      for(fold in 1:num_folds){
        lines(grid_B,dev_fold[,fold],type="l",lwd=2,lty=3)
      }
    } else{
      layout(matrix(c(1:num_folds,rep(0,nrow*ncol - num_folds)), ncol=ncol), widths=c(rep(1,ncol)))
      for(fold in 1:num_folds){
        par(mai=rep(0.5, 4))
        plot(grid_B,dev[[fold]][,1],type="n",
             xlab="B",ylab="dev",main=paste("fold",fold),
             ylim = c(min(dev[[fold]]),max(dev[[fold]])))
        for(par_value in 1:ncol(dev[[fold]])){
          lines(grid_B,dev[[fold]][,par_value],lty=par_value,lwd=2)
          abline(v=grid_B[dev[[fold]][,par_value] == min(dev[[fold]][,par_value])],lty = par_value,lwd=3)
        }
      }
      par(mai=c(0,0,0,0))
      plot.new()
      grid = object$grid
      if(par %in% c("depth","min_leaf_size")){
        legend_text = sapply(grid,function(par){paste0("(",par[1],",",par[2],")")})
      } else{
        legend_text = as.character(grid)
      }
      legend(x="right", ncol=1,legend=legend_text,
             lty=1:length(legend_text), title=par,bty="n",lwd=2,cex=1.5)
      }
    } else if(what == "data"){
    folds = object$folds
    y = object$y

    num_folds = object$num_folds
    nrow = floor(sqrt(num_folds))
    ncol = ceiling(num_folds/nrow)
    dev = object$dev_folds

    layout(matrix(c(1:num_folds,rep(0,nrow*ncol - num_folds)), ncol=ncol), widths=c(rep(1,ncol)))
    par(mai=rep(0.5, 4))
    for(fold in 1:num_folds){
      hist(y[folds==fold],seq(min(y)-0.1,max(y) + 0.1,length.out=20),main=paste("Fold",fold),
           xlab="y",ylab="counts",xlim=range(y),10)
    }
  } else{
    stop("This plot is not specified for a par_CV object.")
  }
}
