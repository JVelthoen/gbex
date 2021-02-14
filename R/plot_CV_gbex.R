#' Plot function for CV_gbex
#'
#' @param object par_CV object
#' @return linegraph of the the cross validation deviance for the different parameter values in the par_grid
#' @export
plot.CV_gbex <- function(object){
  grid_B = 0:object$Bmax
  if(object$par_name == "B"){
    data = data.frame(dev=object$dev_all, B = grid_B)
    data = data[stats::complete.cases(data),]
    g = ggplot2::ggplot(data,ggplot2::aes(y=dev,x=B)) +
        ggplot2::geom_line(size=1) +
        ggplot2::geom_vline(xintercept = which(object$dev_all == min(object$dev_all)) - 1, lty=2,size=0.8) +
        ggplot2::labs(title="CV deviance", x = "Iteration", y = "Deviance") +
        ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5))
    } else if(object$par_name %in% c("lambda","lambda_ratio","lambda_scale","depth","min_leaf_size","sf")){
      data = data.frame(B=grid_B,dev=object$dev_all)
      data = stats::reshape(data, direction = "long", varying = colnames(data)[-1],
                     v.names = "dev", timevar = "par_name")
      data = data[stats::complete.cases(data),]
      grid_names = unlist(lapply(object$par_grid,
                                 function(x){
                                   if(length(x) == 1) return(as.character(x))
                                   else return(paste0("(",paste0(x,collapse=","),")"))
                                 }
      ))
      data$par_name = factor(grid_names[data$par_name],levels=grid_names)
      data_vline = data.frame(par_name = factor(grid_names[which(sapply(object$par_grid,
                                                                        function(grid_value){
                                                                          all(grid_value == object$par_CV)
                                                                        }))],levels=grid_names),
                              int = object$B_opt)
      g = ggplot2::ggplot(data,ggplot2::aes(y=dev,x=B,lty=par_name)) +
        ggplot2::geom_line(size=1) +
        ggplot2::geom_vline(data = data_vline, ggplot2::aes(xintercept = int), linetype = 2,size=0.8) +
        ggplot2::labs(title=paste("CV deviance per",object$par_name), x = "Iteration", y = "Deviance",lty=object$par_name) +
        ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5))
    }
  return(g)
}
