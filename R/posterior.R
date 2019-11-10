#' Get posterior quantities from model output
#'
#' \code{keyATM_output()} makes various quantities that help interpret the model.
#'
#' @param model a fitted keyATM model (an output of \code{keyATM_fit()})
#'
#' @return A keyATM_output containing:
#'   \describe{
#'     \item{keyword_k}{Number of keyword topics}
#'     \item{no_keyword_topics}{Number of regular unseeded topics}
#'     \item{V}{Number of word types}
#'     \item{N}{Number of documents}
#'     \item{theta}{Normalized topic proportions for each document}
#'     \item{phi}{Normalized topic specific word generation probabilities}
#'     \item{topic_counts}{Number of tokens assigned to each topic}
#'     \item{word_counts}{Number of times each word type appears}
#'     \item{doc_lens}{Length of each document in tokens}
#'     \item{vocab}{Words in the vocabulary}
#'     \item{model_fit}{Perplexity and log-likelihood}
#'     \item{p}{Estimated p}
#'     \item{values_iter}{Organized values stored during iterations}
#'   }
#'
#' @import magrittr
#'
keyATM_output <- function(model)
{
  message("Creating an output object. It may take time...")

  check_arg_type(model, "keyATM_fitted")
  values_iter <- list()  # store values by iteration
  model$model <- extract_full_model_name(model)

  # Make info
  info <- list()
  info$allK <- model$no_keyword_topics + length(model$keywords)
  info$V <- length(model$vocab)
  info$N <- length(model$Z)
  info$doc_lens <- sapply(model$Z, length)

  if (model$no_keyword_topics > 0 & length(model$keywords) != 0) {
    info$tnames <- c(names(model$keywords_raw), paste0("Other_", 1:model$no_keyword_topics))
  } else if (model$no_keyword_topics > 0 & length(model$keywords) == 0) {
    # No keywords (= lda models)
    info$tnames <- paste0("Topic_", 1:model$no_keyword_topics)
  } else {
    # Keywords only
    info$tnames <- c(paste0("", 1:length(model$keywords)))
  }


  # theta (document-topic distribution)
  theta <- keyATM_output_theta(model, info)

  # theta iter
  if (model$options$store_theta) {
    values_iter$theta_iter <- keyATM_output_theta_iter(model, info)  
    total_iter <- 1:(model$options$iterations)
    thinning <- model$options$thinning
    values_iter$used_iter <- total_iter[(total_iter %% thinning == 0) | (total_iter == 1) | total_iter == max(total_iter)]
  }

  # Phi (topic-word distribution)
  res <- keyATM_output_phi(model, info)
  phi <- res$phi
  topic_counts <- res$topic_counts
  word_counts <- res$word_counts
  

  # alpha_iter
  if (model$model %in% c("hmm", "ldahmm")) {
    values_iter$alpha_iter <- keyATM_output_alpha_iter_hmm(model, info)
  }

  if ((model$model %in% c("base", "lda"))) {
    if (model$options$estimate_alpha)
      values_iter$alpha_iter <- keyATM_output_alpha_iter_base(model, info)  
  }

  # model fit
  modelfit <- NULL
  if (length(model$model_fit) > 0) {
    model$model_fit %>%
      purrr::set_names(1:length(.)) %>%
      dplyr::bind_rows() %>%
      t() %>%
      tibble::as_tibble(., .name_repair = ~c("Iteration", "Log Likelihood", "Perplexity")) -> modelfit
  }

  # p
  if (model$model %in% c("base", "cov", "hmm")){
    p_estimated <- keyATM_output_p(model) 
  } else {
    p_estimated <- NULL 
  }

  # Make an object to return
  ll <- list(keyword_k = length(model$keywords), no_keyword_topics = model$no_keyword_topics,
             V = length(model$vocab), N = length(model$Z),
             model = abb_model_name(model$model),
             theta = theta, phi = phi,
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = info$doc_lens, vocab = model$vocab,
             keywords_raw = model$keywords_raw,
             model_fit = modelfit, p = p_estimated,
             values_iter = values_iter)
  class(ll) <- c("keyATM_output", model$model, class(ll))
  return(ll)
}

#' @noRd
#' @import magrittr
keyATM_output_p <- function(model)
{
  data <- tibble::tibble(Z = unlist(model$Z, use.names = F),
                         X = unlist(model$X, use.names = F))
  data %>%
    dplyr::mutate(Topic = Z+1L) %>%
    dplyr::select(-starts_with("Z")) %>%
    dplyr::group_by(Topic) %>%
    dplyr::summarize(count = (dplyr::n()), sumx = sum(X)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Proportion = sumx/count*100) %>%
    dplyr::select(-sumx) -> p_estimated

  return(p_estimated)
}


#' @noRd
#' @import magrittr
keyATM_output_theta <- function(model, info)
{

  # Theta
  if (model$model %in% c("cov", "ldacov")) {
    Alpha <- exp(model$model_settings$covariates_data %*% t(model$stored_values$Lambda_iter[[length(model$stored_values$Lambda_iter)]]))

    posterior_z <- function(docid){
      zvec <- model$Z[[docid]]
      alpha <- Alpha[docid, ]
      tt <- table(factor(zvec, levels = 1:(info$allK) - 1L))
      (tt + alpha) / (sum(tt) + sum(alpha)) # posterior mean
    }

    theta <- do.call(dplyr::bind_rows, lapply(1:length(model$Z), posterior_z))

  } else if (model$model %in% c("base", "lda")) {
    if (model$options$estimate_alpha) {
      alpha <- model$stored_values$alpha_iter[[length(model$stored_values$alpha_iter)]]  
    } else {
      alpha <- model$priors$alpha  
    }

    posterior_z <- function(zvec){
      tt <- table(factor(zvec, levels = 1:(info$allK) - 1L))
      (tt + alpha) / (sum(tt) + sum(alpha)) # posterior mean
    }  

    theta <- do.call(dplyr::bind_rows, lapply(model$Z, posterior_z))

  } else if (model$model %in% c("hmm", "ldahmm")) {
    S <- model$stored_values$S_iter[[length(model$stored_values$S_iter)]] + 1L  # adjust index for R
    S <- S[model$model_settings$time_index]  # retrieve doc level state info
    alphas <- matrix(model$stored_values$alpha_iter[[length(model$stored_values$alpha_iter)]][S],
                     nrow = length(model$Z), ncol = info$allK)

    Z_table <- do.call(dplyr::bind_rows, 
                       lapply(model$Z, 
                        function(zvec){table(factor(zvec, levels = 1:(info$allK) - 1L))}))

    tt <- Z_table + alphas
    theta <- tt / Matrix::rowSums(tt)
  }

  theta <- as.matrix(theta)
  colnames(theta) <- info$tnames # label seeded topics

  return(theta)
}


#' @noRd
#' @import magrittr
keyATM_output_phi <- function(model, info)
{
  all_words <- model$vocab[as.integer(unlist(model$W, use.names = F)) + 1L]
  all_topics <- as.integer(unlist(model$Z, use.names = F))
  
  res_tibble <- data.frame(
                        Word = all_words,
                        Topic = all_topics
                       ) %>%
                dplyr::group_by(Topic, Word) %>%
                dplyr::summarize(Count = dplyr::n())
  
  res_tibble %>%
    tidyr::spread(key = Word, value = Count)  -> beta
  beta <- apply(beta, 2, function(x){ifelse(is.na(x), 0, x)})
  beta <- beta[, 2:ncol(beta)] + model$priors$beta
  beta <- beta[, model$vocab]

  topic_counts <- Matrix::rowSums(beta)
  word_counts <- Matrix::colSums(beta)

  phi <- beta / topic_counts
  rownames(phi) <- info$tnames

  return(list(phi = phi, topic_counts = topic_counts, word_counts = word_counts))
}


keyATM_output_theta_iter <- function(model, info)
{
  if (model$model %in% c("cov", "ldacov")) {
    posterior_theta <- function(x){
      Z_table <- model$stored_values$Z_tables[[x]]
      lambda <- model$stored_values$Lambda_iter[[x]]
      Alpha <- exp(model$model_settings$covariates_data %*% t(lambda))

      tt <- Z_table + Alpha
      row.names(tt) <- NULL

      return(tt / Matrix::rowSums(tt))
    }
  } else if (model$model %in% c("hmm", "ldahmm")) {
    posterior_theta <- function(x){
      Z_table <- model$stored_values$Z_tables[[x]]
      S <- model$stored_values$S_iter[[x]] + 1L  # adjust index for R
      S <- S[model$model_settings$time_index]  # retrieve doc level state info

      alphas <- matrix(model$stored_values$alpha_iter[[x]][S],
                       nrow = length(model$Z), ncol = info$allK)
    
      tt <- Z_table + alphas
      theta <- tt / Matrix::rowSums(tt)
      return(theta)
    }
  } else {
    posterior_theta <- function(x){
      Z_table <- model$stored_values$Z_tables[[x]]
      alpha <- model$stored_values$alpha_iter[[x]]

      return((sweep(Z_table, 2, alpha, "+")) / 
              (Matrix::rowSums(Z_table) + sum(alpha)))
    }
  }  

  theta_iter <- lapply(1:length(model$stored_values$Z_tables),
                        posterior_theta)
  return(theta_iter)
}


#' @noRd
#' @import magrittr
keyATM_output_alpha_iter_base <- function(model, info)
{
  topics <- paste0(1:(info$allK))
  model$stored_values$alpha_iter %>%
    purrr::set_names(1:length(.))   %>%
    dplyr::bind_rows() %>%
    t() %>%
    tibble::as_tibble(., .name_repair = ~topics) %>%
    dplyr::mutate(Iteration = 1:(dplyr::n())) %>%
    tidyr::gather(key = Topic, value = alpha, -Iteration) %>%
    dplyr::mutate(Topic = as.integer(Topic)) -> alpha_iter
  return(alpha_iter)
}


keyATM_output_alpha_iter_hmm <- function(model, info)
{
  topics <- paste0(1:(info$allK))
  model$stored_values$alpha_iter %>%
    purrr::imap_dfr(., function(x, i){
                          x %>%
                            tibble::as_tibble(.,
                                              .name_repair = ~topics) %>%
                            dplyr::mutate(State = 1:(dplyr::n()),
                                          Iteration = i) %>%
                            tidyr::gather(key = Topic, value = alpha, -State, -Iteration)
                        }) %>%
     dplyr::mutate(Topic = as.integer(Topic)) -> alpha_iter
  return(alpha_iter)
}


#' @noRd
#' @export
print.keyATM_output <- function(x)
{
  cat(
      paste0(
             "keyATM_output object for the ",
             x$model,
             " model. ",
             "\n"
            )
     )
}


#' @noRd
#' @export
summary.keyATM_output <- function(x)
{
  cat(
      paste0(
             "keyATM_output object for the ",
             x$model,
             " model. ",
             "\n"
            )
     )
}


#' @noRd
#' @export
save.keyATM_output <- function(x, file = stop("'file' must be specified"))
{
  saveRDS(x, file = file)
}


#' @noRd
#' @export
plot.keyATM_output <- function(x)
{
  print(plot_modelfit(x))
}



#' Show the top words for each topic
#'
#' If \code{show_keyword} is true then words in their seeded categories
#' are suffixed with a check mark. Words from another seeded category
#' are labeled with the name of that category.
#'
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param n How many terms to show. Default: NULL, which shows all
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#' @param show_keyword Mark keywords. (default: TRUE)
#'
#' @return An n x k table of the top n words in each topic
#' @export
#'
top_words <- function(x, n = 10, measure = c("probability", "lift"),
                      show_keyword = TRUE)
{
  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)

  if (modelname %in% c("lda", "ldacov", "ldahmm"))
     show_keyword <- FALSE

  if (is.null(n))
    n <- nrow(x$theta)
  measure <- match.arg(measure)
  if (measure == "probability") {
     measuref <- function(xrow){
       colnames(x$phi)[order(xrow, decreasing = TRUE)[1:n]]
     }
  } else if (measure == "lift") {
     wfreq <- x$word_counts / sum(x$word_counts)
     measuref <- function(xrow){
       colnames(x$phi)[order(xrow / wfreq, decreasing = TRUE)[1:n]]
     }
  }
  res <- apply(x$phi, 1, measuref)

  if (show_keyword) {
    for (i in 1:ncol(res)) {
      for (j in 1:length(x$keywords_raw)) {
         inds <- which(res[,i] %in% x$keywords_raw[[j]])
         label <- ifelse(i == j,
                         paste0("[", "\U2713" ,"]"),
                         paste0("[", as.character(j), "]"))
         res[inds, i] <- paste(res[inds, i], label)
      }
    }
  }
  res <- as.data.frame(res)

  return(res)
}



#' Show the top topics for each document
#'
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param n How many topics to show. Default: 2
#'
#' @return An n x k table of the top n topics in each document
#' @import magrittr
#' @export
#'
top_topics <- function(x, n = 2)
{
  check_arg_type(x, "keyATM_output")
  check_arg_type(n, "numeric")

  if (n > ncol(x$theta))
    n <- ncol(x$theta)

  measuref <- function(xrow){
    colnames(x$theta)[order(xrow, decreasing = TRUE)[1:n]]
  }

  res <- t(apply(x$theta, 1, measuref)) %>%
          tibble::as_tibble(., .name_repair = ~paste0("Rank", 1:n))
  return(res)
}



#' Show the top documents for each topic
#'
#' @param x the output from a keyATM model (see \code{keyATM_output()})
#' @param n How many documents to show. Default: 10
#'
#' @return An n x k table of the top n documents for each topic, each number is a document index
#' @import magrittr
#' @export
top_docs <- function(x, n = 10)
{
  check_arg_type(x, "keyATM_output")
  if (is.null(n))
    n <- nrow(x$theta)

  measuref <- function(xcol){
    order(xcol, decreasing = TRUE)[1:n]
  }
  
  res <- apply(x$theta, 2, measuref) %>%
          tibble::as_tibble(.)
  return(res) 
}



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
                       scales = "fixed")
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

  if (modelname %in% c("base", "lda")) {
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
#' @param start 
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


#' Regress theta on covariates
#'
#' @param x the output from a keyATM model (see \code{keyATM()})
#' @param covariates_data matrix, data.frame, or tibble
#' @param covariates_formula
#'
#' @return keyATM_coefficients object
#' @import dplyr
#' @import magrittr
#' @export
estimate_coefficients <- function(x, covariates_data, covariates_formula = NULL,
                             thinning = 5, burn_in = NULL,
                             parallel = TRUE, mc.cores = NULL)
{
  # Check inputs
  if (nrow(covariates_data) != nrow(x$theta)) {
    stop("The row of `covariates_data` should be the same as the number of documents.")  
  }


  if (is.null(covariates_formula)) {
    formula <- as.formula("~ .")
  }

  if (is.null(burn_in)) {
    burn_in <- floor(max(x$model_fit$Iteration) / 2) 
  }

  # Check if it works as a valid regression 
  temp <- as.data.frame(covariates_data)
  temp$y <- rnorm(nrow(covariates_data))
  fit <- lm(y ~ 0 + ., data = temp)

  if (NA %in% fit$coefficients) {
    stop("Covariates are invalid.")    
  }


  # Run regression
  message("Fitting regression...") 
  outcome_name <- paste0("Y_",
                         as.character(format(Sys.time(), "%x%X")))  # unique name
  outcome_name <- gsub("/", "_", outcome_name)
  outcome_name <- gsub(":", "_", outcome_name)
  tname <- colnames(x$theta)
  covariates_data <- as.data.frame(covariates_data)

  if (is.null(x$values_iter$theta_iter)) {
    warning("`options$store_theta` in `keyATM()` was FALSE. keyATM cannot calculate credible intervals.") 
    res <- fit_regression(x$theta, covariates_data, formula, outcome_name, tname)
    obj <- list(res = res, topic_names = tname)
    class(obj) <- c("keyATM_coefficients_point", class(obj))
    return(obj)
  } else {
    # Run for stored theta
    used_iter <- x$values_iter$used_iter
    used_iter <- used_iter[used_iter > burn_in]

    if (parallel) {
      if (is.null(mc.cores)){
        num_core <- parallel::detectCores(all.tests = FALSE, logical = T) - 2L
      } else {
        num_core <- mc.cores 
      }
    } else {
      num_core <- 1L
    }

    res <- do.call(dplyr::bind_rows,
                   parallel::mclapply(1:length(used_iter), 
                                      function(i){
                                        return(fit_regression(x$values_iter$theta_iter[[i]], 
                                                              covariates_data, formula,
                                                              outcome_name, tname,
                                                              used_iter[i])) 
                                      },
                                      mc.cores = num_core
                                      ))
  }

  obj <- list(res = res, topic_names = tname)
  
  class(obj) <- c("keyATM_coefficients", class(obj))
  return(obj)

}

#' @noRd
#' @import magrittr
fit_regression <- function(theta, cov, formula, outcome_name, tname, iter = NULL)
{
  theta <- as.data.frame(theta)
  num_topics <- ncol(theta)

  res <- sapply(1:num_topics,
                function(k){
                   cov[, outcome_name] <- theta[, k]
                   fit <- lm(as.formula(paste(c(outcome_name, as.character(formula)), collapse = " ")),
                             data = cov)
                   return(fit$coefficients)
                })
  colnames(res) <- tname
  res %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Variable = row.names(res)) %>%
    tidyr::gather(key = Topic, value = Coefficient, -Variable) -> res

  if (!is.null(iter)) {
    res$Iteration <- iter 
  }

  return(res)
}


#' @noRd
#' @import magrittr
#' @import ggplot2
#' @export
plot.keyATM_coefficients <- function(obj, topics = NULL, prob_vec = c(0.05, 0.5, 0.95), variables = NULL)
{

  res <- obj$res
  tnames <- obj$topic_names

  coeff <- summary.keyATM_coefficients(obj, topics, prob_vec)$coeff

  if (is.null(variables)) {
    variables <- unique(coeff$Variable) 
  }

  p <- ggplot() +
            coord_flip() +
            geom_linerange(data = coeff,
                            aes(x = Variable, ymin = Lower, ymax = Upper,
                                group = Topic, colour = Topic), 
                            size=1.0, position = position_dodge(width = -1/2)) +
            scale_x_discrete(limits = rev(variables)) +
            geom_hline(aes(yintercept=0), linetype="dashed") +
            xlab("Variable") + ylab("Value") + ggtitle("Coefficient") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5),
                  text = element_text(size=12))

  return(p)
}


#' @noRd
#' @import magrittr
#' @import ggplot2
#' @export
plot.keyATM_coefficients_point <- function(obj, topics = NULL, variables = NULL)
{
  warning("`options$store_theta` in `keyATM()` was FALSE. keyATM cannot calculate credible intervals.")
  res <- obj$res
  tnames <- obj$topic_names

  if (is.null(variables)) {
    variables <- unique(res$Variable) 
  }

  p <- ggplot() +
            coord_flip() +
            geom_point(data = res,
                       aes(x = Variable, y = Coefficient,
                           group = Topic, colour = Topic), 
                       size=2.5, position = position_dodge(width = -1/2)) +
            scale_x_discrete(limits = rev(variables)) +
            geom_hline(aes(yintercept=0), linetype="dashed") +
            xlab("Variable") + ylab("Value") + ggtitle("Coefficient") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5),
                  text = element_text(size=12))

  return(p)
}



#' @noRd
#' @import magrittr
#' @export
summary.keyATM_coefficients <- function(obj, topics = NULL, prob_vec = c(0.05, 0.5, 0.95))
{
  res <- obj$res
  if (is.null(topics)) {
    topics <- 1:length(unique(res$Topic)) 
  }

  if (min(prob_vec) < 0) {
    stop("Check `prob_vec`. The minimum should not be smaller than 0.") 
  }

  if (max(prob_vec) > 1) {
    stop("Check `prob_vec`. The maximum shouhld not be greater than 1.") 
  }

  if (length(prob_vec) != 3) {
    stop("Check `prob_vec`. The length should be 3.") 
  }

  tnames <- obj$topic_names

  coeff <- do.call(dplyr::bind_rows,
                   lapply(topics,
                          function(k){
                             res %>%
                               dplyr::filter(Topic == tnames[k]) %>%
                               dplyr::select(-Topic) %>%
                               tidyr::spread(key = Variable, value = Coefficient) -> temp
                         
                             if ("(Intercept)" %in% colnames(temp)) {
                               temp %>% dplyr::rename(Intercept = `(Intercept)`) -> temp
                             }

                             temp <- data.frame(apply(temp[, 2:ncol(temp)], 2, quantile, prob = prob_vec))
                             temp$Percentile <- c("Lower", "Point", "Upper")
                             temp %>% 
                               tidyr::gather(key = Variable, value = Value, -Percentile) %>%
                               tidyr::spread(key = Percentile, value = Value) %>%
                               dplyr::mutate(Topic = tnames[k]) %>%
                               tibble::as_tibble() -> temp
                          })
                   )

  obj <- list(coeff = coeff, topics = topics)

  class(obj) <- c("summary.keyATM_coefficients", class(obj))
  return(obj)
}


#' @noRd
#' @import magrittr
#' @export
print.summary.keyATM_coefficients <- function(obj)
{
  print(data.frame(obj$coeff))
}

