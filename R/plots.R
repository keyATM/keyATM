#' Show a diagnosis plot of alpha
#'
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param start Slice iteration
#' @param show_topic a vector to specify topic indexes to show
#' @param thinning a integer for thinning
#' @param scale a parameter to control the scale of y-axis: 'free' adjusts y-axis for parameters
#'
#' @return ggplot2 object
#' @importFrom stats as.formula
#' @import ggplot2
#' @import magrittr
#' @export
plot_alpha <- function(x, start = 0, show_topic = NULL,
                       thinning = 5,
                       scale = "fixed")
{

  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)

  if (!"alpha_iter" %in% names(x$values_iter)) {
    stop("`alpha` is not stored. Please check the settings of the model.")  
  }

  thinning <- as.integer(thinning)
  enq_thinning <- enquo(thinning)

  if (is.null(show_topic)) {
    show_topic <- 1:ncol(x$theta)  
  }
  check_arg_type(show_topic, "numeric")
  enq_show_topic <- enquo(show_topic)

  x$values_iter$alpha_iter %>%
    dplyr::filter(Iteration %% (!!enq_thinning) == 0) %>%
    dplyr::filter(Iteration >= start) %>%
    dplyr::filter(Topic %in% (!!show_topic)) %>%
    dplyr::mutate(Topic = paste0("Topic", Topic)) -> res_alpha

  if (nrow(res_alpha) == 0) {
    stop("Nothing left to plot. Please check arguments.")  
  }

  if (modelname %in% c("base", "lda", "label")) {
    p <- ggplot(res_alpha, aes(x = Iteration, y = alpha, group = Topic)) +
          geom_line() +
          geom_point(size = 0.3) +
          facet_wrap(~ Topic, ncol = 2, scales = scales) +
          ylab("Value") +
          ggtitle("Estimated alpha") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
  } else if (modelname %in% c("hmm", "ldahmm")) {
    res_alpha %>% mutate(State = as.character(State)) -> res_alpha
    p <- ggplot(res_alpha, aes(x = Iteration, y = alpha, group = State, colour = State)) +
          geom_line() +
          geom_point(size = 0.3) +
          facet_wrap(~ Topic, ncol = 2, scales = scales) +
          ylab("Value") +
          ggtitle("Estimated alpha") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))  
  }
  return(p)
}


#' Show a diagnosis plot of log-likelihood and perplexity
#'
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param start starting value of the plot
#'
#' @return ggplot2 object
#' @import ggplot2
#' @importFrom stats as.formula
#' @export
plot_modelfit <- function(x, start = 1)
{

  check_arg_type(x, "keyATM_output")

  modelfit <- x$model_fit

  if (!is.numeric(start) | length(start) != 1) {
    message("`start` argument is invalid. Using the default (=1)")  
    start <- 1
  }

  if (!is.null(start)) {
    modelfit <- modelfit[ modelfit$Iteration >= start, ]
  }

  modelfit <- tidyr::gather(modelfit, key = Measures, value = value, -Iteration)

  p <- ggplot(data = modelfit, aes_string(x='Iteration', y='value',
                                          group='Measures', color='Measures')) +
     geom_line(show.legend = F) +
     geom_point(size = 0.3, show.legend = F) +
     facet_wrap(as.formula(paste("~", "Measures")), ncol = 2, scales = "free") +
     ylab("Value")

  p <- p + ggtitle("Model Fit") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


#' Show a diagnosis plot of p
#'
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param show_topic A vector to indicate topics to visualize
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @export
plot_p <- function(x, show_topic = NULL)
{
  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)

  if (modelname %in% c("lda", "ldacov", "ldahmm")) {
    stop(paste0("`", x$model, "` is not a model with keywords.")) 
  }

  num <- length(unique(x$p$Topic))
  if (is.null(show_topic)) {
    shoe_topic <- 1:num
  }

  check_arg_type(show_topic, "numeric")
  enq_show_topic <- enquo(show_topic)

  x$p %>%
    dplyr::filter(Topic %in% (!!show_topic)) %>%
    dplyr::mutate(Topic = paste0("Topic", Topic)) -> temp

  g  <- ggplot(temp, aes_string(x='Topic', y='Proportion')) +
      geom_bar(stat="identity") +
      theme_bw() +
      scale_x_discrete(limits = paste0("Topic", get("show_topic"))) +
      ylab("Proportion (%)") +
      xlab("Topic") +
      ggtitle("Proportion of words drawn from topic-word distribution") +
      theme(plot.title = element_text(hjust = 0.5))

  return(g)
}


#' Plot document-topic distribution by strata 
#'
#' @param x a strata_doctopic object (see \code{by_strata_DocTopic()})
#' @param topics topics to show
#' @param quantile_vec quantiles to show
#' @param ... additional arguments not used
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import magrittr
#' @export
plot.strata_doctopic <- function(x, topics = NULL, quantile_vec = c(0.05, 0.5, 0.95), ...)
{
  tables <- summary.strata_doctopic(x, quantile_vec = quantile_vec)
  by_name <- x$by_name
  by_values <- x$by_values

  if (is.null(topics)) {
    topics <- 1:nrow(tables[[1]]) 
  }

  tables <- dplyr::bind_rows(tables) %>%
              dplyr::filter(TopicId %in% topics)

  variables <- unique(tables$by)


  p <- ggplot(tables) +
        geom_linerange(aes(x = by,
                           ymin = Lower, ymax = Upper,
                           group = Topic, colour = Topic),
                       position = position_dodge(width = -1/2)) +
        coord_flip() +
        scale_x_discrete(limits = rev(variables)) +
        xlab(paste0("Value of ", by_name)) +
        ylab(expression(paste("Mean of ", theta))) +
        theme_bw()

  return(p)
}
