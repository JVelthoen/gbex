#' Plot function for PD_gbex
#'
#' @param object A PD_gbex object
#' @return Two figures side by side with on the left the partial dependence plot for sigma and on the right the partial dependence plot for gamma.
#' @export
plot.PD_gbex <- function(object){
  dimension = object$dimensions
  if("Q" %in% colnames(object$PD_df)){
    plot_quant  = ggplot2::ggplot(object$PD_df,ggplot2::aes_string(x=object$var_name,y="Q")) +
      ggplot2::geom_line(lwd=1.5) +
      ggplot2::labs(title=paste0("Partial dependence Quantile",tau), x = object$var_name, y = "Partial dependence") +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.border = ggplot2::element_rect(fill=NA))
  } else{
    if(dimension == 1){
      plot_sigma  = ggplot2::ggplot(object$PD_df,ggplot2::aes_string(x=object$var_name,y="s")) +
        ggplot2::geom_line(lwd=1.5) +
        ggplot2::labs(title="Sigma", x = object$var_name, y = "Partial dependence") +
        ggplot2::theme_minimal() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       panel.border = ggplot2::element_rect(fill=NA))

      plot_gamma  = ggplot2::ggplot(object$PD_df,ggplot2::aes_string(x=object$var_name,y="g")) +
        ggplot2::geom_line(lwd=2) +
        ggplot2::labs(title="Gamma", x = object$var_name, y = "Partial dependence") +
        ggplot2::theme_minimal() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       panel.border = ggplot2::element_rect(fill=NA))
    } else if(dimension == 2){
      plot_sigma  = ggplot2::ggplot(object$PD_df,ggplot2::aes_string(x=object$var_name[1],y=object$var_name[2],fill="s")) +
        ggplot2::geom_raster() +
        ggplot2::labs(title="Sigma", x = object$var_name[1], y = object$var_name[2],fill="Partial dependence") +
        ggplot2::scale_fill_gradient(low="blue",high="red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       panel.border = ggplot2::element_rect(fill=NA))

      plot_gamma  = ggplot2::ggplot(object$PD_df,ggplot2::aes_string(x=object$var_name[1],y=object$var_name[2],fill="g")) +
        ggplot2::geom_raster() +
        ggplot2::labs(title="Gamma", x = object$var_name[1], y = object$var_name[2],fill="Partial dependence") +
        ggplot2::scale_fill_gradient(low="blue",high="red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       panel.border = ggplot2::element_rect(fill=NA))
    }
    g = patchwork::wrap_plots(plot_sigma,plot_gamma)
  }
  return(g)
}
