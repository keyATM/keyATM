#' Save a figure
#' 
#' @param x the keyATM_fig object.
#' @param filename file name to create on disk.
#' @param ... other arguments passed on to the [ggplot2::ggsave()][ggplot2::ggsave] function.
#' @seealso [visualize_keywords()], [plot_alpha()], [plot_modelfit()], [plot_pi()], [plot_timetrend()], [by_strata_DocTopic()], [values_fig()]
#' @export
save_fig <- function(x, filename, ...) {
  UseMethod("save_fig")
}


#' Get values used to create a figure
#' 
#' @param x the keyATM_fig object.
#' @seealso [save_fig()], [visualize_keywords()], [plot_alpha()], [plot_modelfit()], [plot_pi()], [plot_timetrend()], [by_strata_DocTopic()]
#' @export
values_fig <- function(x) {
  UseMethod("values_fig")
}

#' @noRd
#' @export
values_fig.keyATM_fig <- function(x)
{
  return(x$values) 
}


#' @noRd
#' @export
save_fig.keyATM_fig <- function(x, filename, ...)
{
  ggplot2::ggsave(filename = filename, plot = x$figure, ...)
}


#' @noRd
#' @export
print.keyATM_fig <- function(x, ...)
{
  print(x$figure) 
}


#' Show a diagnosis plot of alpha
#'
#' @param x the output from a keyATM model (see [keyATM()]).
#' @param start integer. The start of slice iteration. Default is \code{0}.
#' @param show_topic a vector to specify topic indexes to show. Default is \code{NULL}.
#' @param scales character. Control the scale of y-axis (the parameter in [ggplot2::facet_wrap()][ggplot2::facet_wrap]): \code{free} adjusts y-axis for parameters. Default is \code{fixed}. 
#' @return keyATM_fig object
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
#' @seealso [save_fig()]
#' @export
plot_alpha <- function(x, start = 0, show_topic = NULL, scales = "fixed")
{
  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)

  if (modelname %in% c("lda", "ldacov", "ldahmm")) {
    stop(paste0("This is not a model with keywords."))  # only plot keywords later
  }
  if (!"alpha_iter" %in% names(x$values_iter)) {
    stop("`alpha` is not stored. Please check the options.\nNote that the covariate model does not have `alpha`.\nPlease check our paper for details.")
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
          facet_wrap(~ .data$Topic, ncol = 2, scales = scales) +
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
          facet_wrap(~ .data$Topic, ncol = 2, scales = scales) +
          xlab("Iteration") + ylab("Value") +
          ggtitle("Estimated alpha") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))  
  }

  p <- list(figure = p, values = res_alpha)
  class(p) <- c("keyATM_fig", class(p))
  return(p)
}


#' Show a diagnosis plot of log-likelihood and perplexity
#' 
#' @param x the output from a keyATM model (see [keyATM()]).
#' @param start integer. The starting value of iteration to use in plot. Default is \code{1}.
#' @return keyATM_fig object.
#' @import ggplot2
#' @importFrom stats as.formula
#' @importFrom rlang .data
#' @seealso [save_fig()]
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

  p <- list(figure = p, values = modelfit)
  class(p) <- c("keyATM_fig", class(p))
  return(p)
}


#' Show a diagnosis plot of pi
#'
#' @param x the output from a keyATM model (see [keyATM()]).
#' @param show_topic an integer or a vector. Indicate topics to visualize. Default is \code{NULL}.
#' @param start integer. The starting value of iteration to use in the plot. Default is \code{0}.
#' @param ci value of the credible interval (between 0 and 1) to be estimated. Default is \code{0.9} (90%). This is an option when calculating credible intervals (you need to set \code{store_pi = TRUE} in [keyATM()]).
#' @param method method for computing the credible interval. The Highest Density Interval (\code{hdi}, default) or Equal-tailed Interval (\code{eti}). This is an option when calculating credible intervals (you need to set \code{store_pi = TRUE} in [keyATM()]).
#' @param point method for computing the point estimate. \code{mean} (default) or \code{median}. This is an option when calculating credible intervals (you need to set \code{store_pi = TRUE} in [keyATM()]).
#' @return keyATM_fig object.
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
#' @seealso [save_fig()]
#' @export
plot_pi <- function(x, show_topic = NULL, start = 0, ci = 0.9, method = c("hdi", "eti"), point = c("mean", "median"))
{
  method <- match.arg(method)
  point <- match.arg(point)
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
    pi_mat <- t(sapply(x$values_iter$pi_iter, unlist, use.names = FALSE))[, show_topic, drop = FALSE]
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
      dplyr::summarise(x = list(tibble::enframe(calc_ci(.data$value, ci, method, point), "q", "value")), .groups = "drop_last") %>%
      tidyr::unnest(x) %>% tidyr::pivot_wider(names_from = .data$q, values_from = .data$value) -> temp
    
    p <- ggplot(temp, aes(y = .data$Point, x = .data$Topic)) + 
         theme_bw() + geom_point() +
         geom_errorbar(aes(ymin = .data$Lower, ymax = .data$Upper), data = temp, width = 0.01, size = 1) + 
         xlab("Topic") + ylab("Probability") +
         ggtitle("Probability of words drawn from keyword topic-word distribution") +
         theme(plot.title = element_text(hjust = 0.5))
  } else {
    x$pi %>%
      dplyr::mutate(Probability = .data$Proportion / 100) %>%
      dplyr::filter(.data$Topic %in% (!!show_topic)) %>%
      dplyr::mutate(Topic = tnames) -> temp

    p <- ggplot(temp, aes(x = .data$Topic, y = .data$Probability)) +
         geom_bar(stat = "identity") +
         theme_bw() +
         xlab("Topic") + ylab("Probability") +
         ggtitle("Probability of words drawn from keyword topic-word distribution") +
         theme(plot.title = element_text(hjust = 0.5))    
  }
  p <- list(figure = p, values = temp)
  class(p) <- c("keyATM_fig", class(p))
  return(p)
}


#' Plot document-topic distribution by strata (for covariate models)
#' 
#' @param x a strata_doctopic object (see [by_strata_DocTopic()]).
#' @param show_topic a vector or an integer. Indicate topics to visualize.
#' @param var_name the name of the variable in the plot.
#' @param by `topic` or `covariate`. Default is by `topic`.
#' @param ci value of the credible interval (between 0 and 1) to be estimated. Default is \code{0.9} (90%).
#' @param method method for computing the credible interval. The Highest Density Interval (\code{hdi}, default) or Equal-tailed Interval (\code{eti}).
#' @param point method for computing the point estimate. \code{mean} (default) or \code{median}.
#' @param width numeric. Width of the error bars.
#' @param show_point logical. Show point estimates. The default is \code{TRUE}.
#' @param ... additional arguments not used.
#' @return keyATM_fig object.
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
#' @seealso [save_fig()], [by_strata_DocTopic()]
#' @export
plot.strata_doctopic <- function(x, show_topic = NULL, var_name = NULL, by = c("topic", "covariate"),
                                 ci = 0.9, method = c("hdi", "eti"), point = c("mean", "median"),
                                 width = 0.1, show_point = TRUE, ...)
{
  by <- match.arg(by)
  method <- match.arg(method)
  point <- match.arg(point)

  tables <- summary.strata_doctopic(x, ci, method, point)
  by_var <- x$by_var
  by_values <- x$by_values
  if (!is.null(var_name)) {
    by_var <- var_name
  }

  if (is.null(show_topic)) {
    show_topic <- 1:nrow(tables[[1]]) 
  }

  tables <- dplyr::bind_rows(tables) %>%
              dplyr::filter(.data$TopicId %in% show_topic)
  variables <- unique(tables$label)

  p <- ggplot(tables) +
        coord_flip() +
        scale_x_discrete(limits = rev(variables)) +
        xlab(paste0(by_var)) +
        ylab(expression(paste("Mean of ", theta))) +
        guides(color = guide_legend(title = "Topic")) +
        theme_bw()

  if (by == "topic") {
    p <- p + geom_errorbar(width = width, aes(x = .data$label, ymin = .data$Lower, ymax = .data$Upper,
                group = .data$Topic), position = position_dodge(width = -1/2)) + facet_wrap(~Topic) 
    if (show_point)
      p <- p + geom_point(aes(x = .data$label, y = .data$Point))
  } else {
    p <- p + geom_errorbar(width = width, aes(x = .data$label, ymin = .data$Lower, ymax = .data$Upper,
                group = .data$Topic, colour = .data$Topic), position = position_dodge(width = -1/2))   
    if (show_point)
      p <- p + geom_point(aes(x = .data$label, y = .data$Point, colour = .data$Topic), position = position_dodge(width = -1/2))
  }

  p <- list(figure = p, values = tables)
  class(p) <- c("keyATM_fig", class(p))
  return(p)
}


#' Plot time trend
#' 
#' @param x the output from the dynamic keyATM model (see [keyATM()]).
#' @param show_topic an integer or a vector. Indicate topics to visualize. Default is \code{NULL}.
#' @param time_index_label a vector. The label for time index. The length should be equal to the number of documents (time index provided to [keyATM()]). 
#' @param ci value of the credible interval (between 0 and 1) to be estimated. Default is \code{0.9} (90%). This is an option when calculating credible intervals (you need to set \code{store_theta = TRUE} in [keyATM()]).
#' @param method method for computing the credible interval. The Highest Density Interval (\code{hdi}, default) or Equal-tailed Interval (\code{eti}). This is an option when calculating credible intervals (you need to set \code{store_theta = TRUE} in [keyATM()]).
#' @param point method for computing the point estimate. \code{mean} (default) or \code{median}. This is an option when calculating credible intervals (you need to set \code{store_theta = TRUE} in [keyATM()]).
#' @param xlab a character.
#' @param scales character. Control the scale of y-axis (the parameter in [ggplot2::facet_wrap()][ggplot2::facet_wrap]): \code{free} adjusts y-axis for parameters. Default is \code{fixed}. 
#' @param width numeric. Width of the error bars.
#' @param show_point logical. The default is \code{TRUE}. This is an option when calculating credible intervals.
#' @param ... additional arguments not used.
#' @return keyATM_fig object.
#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
#' @seealso [save_fig()]
#' @export
plot_timetrend <- function(x, show_topic = NULL, time_index_label = NULL, 
                           ci = 0.9, method = c("hdi", "eti"), point = c("mean", "median"),
                           xlab = "Time", scales = "fixed", width = 0.5, show_point = TRUE, ...)
{
  method <- match.arg(method)
  point <- match.arg(point)
  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)
  if (!modelname %in% c("hmm", "ldahmm")) {
    stop(paste0("This is not a model with keywords.")) 
  }

  if (!is.null(time_index_label)) {
    if (length(x$values_iter$time_index) != length(time_index_label)) {
      stop("The length of `time_index_label` does not match with the number of documents.")
    }
    time_index <- time_index_label
  } else {
    time_index <- x$values_iter$time_index
  }

  if (is.null(show_topic)) {
    show_topic <- 1:x$keyword_k
  }

  format_theta <- function(theta, time_index, tnames) {
    theta[ , show_topic, drop = FALSE] %>%
      tibble::as_tibble(.name_repair = ~tnames) %>%
      dplyr::mutate(time_index = time_index) %>%
      tidyr::pivot_longer(-.data$time_index, names_to = "Topic", values_to = "Proportion") %>%
      dplyr::group_by(.data$time_index, .data$Topic) %>%
      dplyr::summarize(Proportion = base::mean(.data$Proportion), .groups = "drop_last") -> res
    return(res)
  }

  tnames <- colnames(x$theta)

  if (is.null(x$values_iter$theta_iter)) {
    dat <- format_theta(x$theta, time_index, tnames[show_topic])
    p <- ggplot(dat, aes(x = .data$time_index, y = .data$Proportion, group = .data$Topic)) + 
          geom_line(size = 0.8, color = "blue") + geom_point(size = 0.9)
  } else {
    dplyr::bind_rows(lapply(x$values_iter$theta_iter, format_theta, time_index, tnames[show_topic])) %>%
      dplyr::group_by(.data$time_index, .data$Topic) %>%
      dplyr::summarise(x = list(tibble::enframe(calc_ci(.data$Proportion, ci, method, point), "q", "value"))) %>% 
      tidyr::unnest(.data$x) %>% dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = .data$q, values_from = .data$value) %>%
      stats::setNames(c("time_index", "Topic", "Lower", "Point", "Upper")) -> dat
    p <- ggplot(dat, aes(x = .data$time_index, y = .data$Point, group = .data$Topic)) +
          geom_ribbon(aes(ymin = .data$Lower, ymax = .data$Upper), fill = "gray75") +
          geom_line(size = 0.8, color = "blue")

    if (show_point)
      p <- p + geom_point(size = 0.9)
  }
  p <- p + xlab(xlab) + ylab(expression(paste("Mean of ", theta))) +
        facet_wrap(~.data$Topic, scales = scales) + theme_bw() + theme(panel.grid.minor = element_blank())
  p <- list(figure = p, values = dat)
  class(p) <- c("keyATM_fig", class(p))
  return(p)
}
