#'  # @import tidyr readr ggplot2
#' NULL
#'
#'
#' #' Plot observed p in the most recent run
#' #'
#' #' @param true_path path to the true p.
#' #' @param slice slice outputs
#' #' @param burn_in discard specified number of samplings
#' #'
#' #' @export
#' plot_observed_p <- function(true_path=0, slice=5, burn_in=100){
#' 	options(warn=-1)
#'
#' 	if(true_path != 0){
#' 		p_k <- read.csv("output_p_k.txt")
#' 		colnames(p_k) <- gsub("Col", "Topic ", colnames(p_k))
#'
#' 		p_k_true <- read_delim(true_path, delim=" ", col_names=F)
#' 		p_k_true <- data.frame(value = as.vector(t(as.matrix(p_k_true))))
#' 		p_k_true$Topic <- colnames(p_k)
#'
#' 		p_k$iter  <-  1:nrow(p_k)
#' 		p_k <- p_k[seq(burn_in, nrow(p_k), slice),]
#' 		p_k <- gather(p_k, key=Topic, value=value, -iter)
#'
#' 		p_k_fig <- ggplot(data=p_k, aes(x=iter, y=value, group=Topic)) +
#' 				 geom_line() +
#' 				 geom_point(size=0.3) +
#' 				 facet_wrap(~Topic, ncol=2, scales = "free") +
#' 				 geom_hline(data = p_k_true, aes(yintercept = value), size=0.6, color="red") +
#' 				 theme_bw()
#' 	} else{
#' 			p_k <- read.csv("output_p_k.txt")
#' 			colnames(p_k) <- gsub("Col", "Topic ", colnames(p_k))
#'
#' 			p_k$iter  <-  1:nrow(p_k)
#' 			if(nrow(p_k) >= 150){
#' 				p_k <- p_k[seq(burn_in, nrow(p_k), slice),]
#' 			}
#' 			p_k <- gather(p_k, key=Topic, value=value, -iter)
#'
#' 			p_k_fig <- ggplot(data=p_k, aes(x=iter, y=value, group=Topic)) +
#' 					 geom_line() +
#' 					 geom_point(size=0.3) +
#' 					 facet_wrap(~Topic, ncol=2, scales = "free") +
#' 					 theme_bw()
#' 	}
#'
#' 	options(warn=0)
#' 	return (p_k_fig)
#'
#' }
