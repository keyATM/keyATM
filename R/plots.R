#' Show a diagnosis plot of alpha
#'
#' 
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param start integer. The start of slice iteration. Default is 0.
#' @param show_topic a vector to specify topic indexes to show. Default is \code{NULL}.
#' @param scale character. Control the scale of y-axis (the parameter in \code{facet_wrap()}): 'free' adjusts y-axis for parameters. Default is "fixed". 
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
#' @export
plot_alpha <- function(x, start = 0, show_topic = NULL, scale = "fixed")
{

  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)

  if (modelname %in% c("lda", "ldacov", "ldahmm")) {
    stop(paste0("This is not a model with keywords.")) 
  }

  if (!"alpha_iter" %in% names(x$values_iter)) {
    stop("`alpha` is not stored. Please check the settings of the model.")  
  }

  if (is.null(show_topic)) {
    show_topic <- 1:x$keyword_k
  }

  if (!is.numeric(start) | length(start) != 1) {
    stop("`start` argument is invalid.")  
  }

  tnames <- c(names(x$keywords_raw))[show_topic]
  x$values_iter$alpha_iter %>%
    dplyr::filter(.data$Iteration >= start) %>%
    dplyr::filter(.data$Topic %in% (!!show_topic)) %>%
    tidyr::pivot_wider(names_from = "Topic", values_from = "alpha") -> temp

  if (modelname %in% c("base", "lda", "label")) {
    temp %>% dplyr::rename_at(vars(-"Iteration"), ~tnames) %>%
      tidyr::pivot_longer(-"Iteration", names_to = "Topic", values_to = "alpha") -> res_alpha

    p <- ggplot(res_alpha, aes(x = .data$Iteration, y = .data$alpha, group = .data$Topic)) +
          geom_line() +
          geom_point(size = 0.3) +
          facet_wrap(~ .data$Topic, ncol = 2, scales = scale) +
          xlab("Iteration") + ylab("Value") +
          ggtitle("Estimated alpha") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
  } else if (modelname %in% c("hmm", "ldahmm")) {
    temp %>% dplyr::rename_at(vars(-.data$Iteration, -.data$State), ~tnames) %>%
      tidyr::pivot_longer(-c("Iteration", "State"), names_to = "Topic", values_to = "alpha") -> res_alpha
    res_alpha$State <- factor(res_alpha$State, levels = 1:max(res_alpha$State))

    p <- ggplot(res_alpha, aes(x = .data$Iteration, y = .data$alpha, group = .data$State, colour = .data$State)) +
          geom_line() +
          geom_point(size = 0.3) +
          facet_wrap(~ .data$Topic, ncol = 2, scales = scale) +
          xlab("Iteration") + ylab("Value") +
          ggtitle("Estimated alpha") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))  
  }

  return(p)
}


#' Show a diagnosis plot of log-likelihood and perplexity
#'
#' 
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param start integer. The starting value of iteration to use in plot. Default is 1.
#'
#' @return ggplot2 object
#' @import ggplot2
#' @importFrom stats as.formula
#' @importFrom rlang .data
#' @export
plot_modelfit <- function(x, start = 1)
{

  check_arg_type(x, "keyATM_output")

  modelfit <- x$model_fit

  if (!is.numeric(start) | length(start) != 1) {
    stop("`start` argument is invalid.")  
  }

  if (!is.null(start)) {
    modelfit <- modelfit[ modelfit$Iteration >= start, ]
  }

  modelfit <- tidyr::gather(modelfit, key = "Measures", value = "value", -"Iteration")

  p <- ggplot(data = modelfit, aes(x = .data$Iteration, y = .data$value, group = .data$Measures, color = .data$Measures)) +
     geom_line(show.legend = FALSE) +
     geom_point(size = 0.3, show.legend = FALSE) +
     facet_wrap(~ .data$Measures, ncol = 2, scales = "free") +
     xlab("Iteration") + ylab("Value")

  p <- p + ggtitle("Model Fit") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


#' Show a diagnosis plot of pi
#'
#' 
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param show_topic an integer or a vector. Indicate topics to visualize. Default is \code{NULL}.
#' @param start integer. The starting value of iteration to use in the plot. Default is 0.
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
#' @export
plot_pi <- function(x, show_topic = NULL, start = 0)
{
  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)

  if (modelname %in% c("lda", "ldacov", "ldahmm")) {
    stop(paste0("This is not a model with keywords.")) 
  }

  if (is.null(show_topic)) {
    show_topic <- 1:x$keyword_k
  } else if (sum(!show_topic %in% 1:x$keyword_k) != 0) {
    stop("`plot_pi` only visualize keyword topics.") 
  }

  if (!is.numeric(start) | length(start) != 1) {
    stop("`start` argument is invalid.")  
  }

  tnames <- c(names(x$keywords_raw))[show_topic]

  if (!is.null(x$values_iter$pi_iter)) {
    pi_mat <- t(sapply(x$values_iter$pi_iter, unlist, use.names = FALSE))[, show_topic]
    pi_mat %>%
      tibble::as_tibble(.name_repair = ~ tnames) %>%
      dplyr::mutate(Iteration = x$values_iter$used_iter) %>%
      dplyr::filter(.data$Iteration >= start) %>% 
      dplyr::select(-.data$Iteration) -> pi_mat

    if (nrow(pi_mat) == 0) {
      stop("Nothing left to plot. Please check arguments.")
    }

    pi_mat %>%
      tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Topic") %>%
      dplyr::group_by(.data$Topic) %>%
      dplyr::summarise(mean = mean(.data$value), uq = stats::quantile(.data$value, .975), 
                       lq = stats::quantile(.data$value, 0.025)) -> temp
    
    g <- ggplot(temp, aes(y = .data$mean, x = .data$Topic, color = .data$Topic)) + 
         theme_bw() +
         geom_errorbar(aes(ymin = .data$lq, ymax = .data$uq), data = temp, width = 0.01, size = 1) + 
         xlab("Topic") + ylab("Probability") +
         ggtitle("Probability of words drawn from keyword topic-word distribution") +
         theme(plot.title = element_text(hjust = 0.5))
  } else {
    x$pi %>%
      dplyr::mutate(Probability = .data$Proportion / 100) %>%
      dplyr::filter(.data$Topic %in% (!!show_topic)) %>%
      dplyr::mutate(Topic = tnames) -> temp

    g  <- ggplot(temp, aes(x = .data$Topic, y = .data$Probability)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        xlab("Topic") + ylab("Probability") +
        ggtitle("Probability of words drawn from keyword topic-word distribution") +
        theme(plot.title = element_text(hjust = 0.5))    
  }
  return(g)
}


#' Plot document-topic distribution by strata (for covariate models)
#' 
#'
#' @param x a strata_doctopic object (see \code{by_strata_DocTopic()})
#' @param topics a vector or an integer. Indicate topics to visualize.
#' @param quantile_vec a numeric. Quantiles to visualize
#' @param ... additional arguments not used
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
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
              dplyr::filter(.data$TopicId %in% topics)

  variables <- unique(tables$by)

  p <- ggplot(tables) +
        geom_linerange(aes(x = .data$by,
                           ymin = .data$Lower, ymax = .data$Upper,
                           group = .data$Topic, colour = .data$Topic),
                       position = position_dodge(width = -1/2)) +
        coord_flip() +
        scale_x_discrete(limits = rev(variables)) +
        xlab(paste0("Value of ", by_name)) +
        ylab(expression(paste("Mean of ", theta))) +
        guides(color = guide_legend(title = "Topic")) +
        theme_bw()

  return(p)
}
