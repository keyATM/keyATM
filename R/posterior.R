#' Create an output object
#'
#' @keywords internal
#' @import magrittr
keyATM_output <- function(model, keep)
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
  info$model <- model$model

  if (model$no_keyword_topics > 0 & length(model$keywords) != 0) {
    info$tnames <- c(names(model$keywords_raw), paste0("Other_", 1:model$no_keyword_topics))
  } else if (model$no_keyword_topics > 0 & length(model$keywords) == 0) {
    # No keywords (= lda models)
    info$tnames <- paste0("Topic_", 1:model$no_keyword_topics)
  } else {
    # Keywords only
    info$tnames <- c(names(model$keywords_raw))
  }

  # theta (document-topic distribution)
  theta <- keyATM_output_theta(model, info)

  # theta iter
  if (model$options$store_theta) {
    values_iter$theta_iter <- keyATM_output_theta_iter(model, info)
  }
  # used_iter is useful quantity
  total_iter <- 1:(model$options$iterations)
  thinning <- model$options$thinning
  values_iter$used_iter <- total_iter[(total_iter %% thinning == 0) | (total_iter == 1) | total_iter == max(total_iter)]
  info$used_iter <- values_iter$used_iter

  # Phi (topic-word distribution)
  res <- keyATM_output_phi(model, info)
  phi <- res$phi
  topic_counts <- res$topic_counts
  word_counts <- res$word_counts

  # alpha_iter
  if (model$model %in% c("hmm", "ldahmm")) {
    values_iter$alpha_iter <- keyATM_output_alpha_iter_hmm(model, info)
    values_iter$time_index <- model$model_settings$time_index
  }

  if ((model$model %in% c("base", "lda", "label"))) {
    if (model$options$estimate_alpha)
      values_iter$alpha_iter <- keyATM_output_alpha_iter_base(model, info)
  }

  # model fit
  modelfit <- NULL
  if (length(model$model_fit) > 0) {
    names(model$model_fit) <- 1:length(model$model_fit)
    model$model_fit %>%
      dplyr::bind_rows() %>%
      t() %>%
      tibble::as_tibble(.name_repair = ~c("Iteration", "Log Likelihood", "Perplexity")) -> modelfit
  }

  # pi
  if (model$model %in% c("base", "cov", "hmm", "label")) {
    pi_estimated <- keyATM_output_pi(model$Z, model$S, model$priors$gamma)
  } else {
    pi_estimated <- NULL
  }

  if (length(model$stored_values$pi_vectors) > 0) {
    values_iter$pi_iter <- model$stored_values$pi_vectors
  }

  # Rescale lambda
  if (model$model %in% c("cov", "ldacov")) {
    values_iter$Lambda_iter <- model$stored_values$Lambda_iter
  }

  # Add information
  pd <- utils::packageDescription("keyATM")
  information <- list(date_output_made = Sys.time(), version_keyATM = pd$Version)

  # Check what to keep
  if (model$model %in% c("cov", "ldacov")) {
    if (!"stored_values" %in% keep)
      keep <- c("stored_values", keep)

    if (!"model_settings" %in% keep)
      keep <- c("model_settings", keep)
  }

  if (model$model %in% c("hmm", "ldahmm")) {
    keep <- c("model_settings", keep)
  }

  kept_values <- list(doc_index_used = model$stored_values$doc_index)
  if (length(keep) != 0) {
    use_elements <- keep[keep %in% names(model)]
    for (i in 1:length(use_elements)) {
      kept_values[use_elements[i]]  <- model[use_elements[i]]
    }
  }

  if (model$options$store_theta & "stored_values" %in% keep) {
    kept_values$stored_values$Z_tables <- NULL  # duplicate information
  }

  # Make an object to return
  ll <- list(keyword_k = length(model$keywords), no_keyword_topics = model$no_keyword_topics,
             V = length(model$vocab), N = length(model$Z),
             model = abb_model_name(model$model),
             theta = theta, phi = phi,
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = info$doc_lens, vocab = model$vocab,
             priors = model$priors, options = model$options,
             keywords_raw = model$keywords_raw,
             model_fit = modelfit, pi = pi_estimated,
             values_iter = values_iter, information = information, kept_values = kept_values)
  class(ll) <- c("keyATM_output", model$model, class(ll))
  return(ll)
}


#' @noRd
#' @import magrittr
#' @importFrom rlang .data
keyATM_output_pi <- function(model_Z, model_S, prior)
{
  # p(p | S=s, n, a, b) \propto Be(a+s, b+(n-s))
  #   p(S=s | n, p) p(p | a, b)
  # Expectation is (a+s) / (a+b+n)

  data <- tibble::tibble(Z = unlist(model_Z, use.names = FALSE),
                         S = unlist(model_S, use.names = FALSE))
  data %>%
    dplyr::mutate(Topic = .data$Z+1L) %>%
    dplyr::select(-dplyr::starts_with("Z")) %>%
    dplyr::group_by(.data$Topic) %>%
    dplyr::summarize(count = (dplyr::n()), sums = sum(.data$S)) %>%
    dplyr::ungroup() -> temp

  # Check used topics
  if (nrow(temp) != nrow(prior)) {
    warning("Some of the topics are not used.")
    missing <- setdiff(1:nrow(prior), temp$Topic)
    temp %>%
      tibble::add_row(Topic = missing, count = 0, sums = 0) %>%
      dplyr::arrange(.data$Topic) -> temp
  }

  # Get p
  n <- temp$count
  s <- temp$sums
  a <- prior[, 1]
  b <- prior[, 2]
  p <- (a + s) / (a + b + n)
  temp %>%
    dplyr::mutate(Proportion = p * 100) %>%
    dplyr::select(-.data$sums) -> pi_estimated
  return(pi_estimated)
}


#' @noRd
#' @keywords internal
table_to_vector <- function(x) {
  purrr::set_names(as.numeric(x), names(x))
}


#' @noRd
#' @import magrittr
#' @keywords internal
bind_tables <- function(x) {
  dplyr::bind_rows(!!!lapply(x, table_to_vector)) %>%
    tibble::as_tibble()
}


#' @noRd
#' @import magrittr
#' @importFrom rlang .data
keyATM_output_theta <- function(model, info)
{

  # Theta
  if (model$model %in% c("cov", "ldacov")) {
    Alpha <- exp(model$model_settings$covariates_data_use %*% t(model$stored_values$Lambda_iter[[length(model$stored_values$Lambda_iter)]]))

    posterior_z_cov <- function(docid) {
      zvec <- model$Z[[docid]]
      alpha <- Alpha[docid, ]
      tt <- table(factor(zvec, levels = 1:(info$allK) - 1L))
      (tt + alpha) / (sum(tt) + sum(alpha)) # posterior mean
    }

    theta <- bind_tables(lapply(1:length(model$Z), posterior_z_cov))

  } else if (model$model %in% c("base", "lda", "label")) {
    if (model$options$estimate_alpha) {
      alpha <- model$stored_values$alpha_iter[[length(model$stored_values$alpha_iter)]]
    } else {
      alpha <- model$priors$alpha
    }

    posterior_z_base <- function(zvec) {
      tt <- table(factor(zvec, levels = 1:(info$allK) - 1L))
      (tt + alpha) / (sum(tt) + sum(alpha)) # posterior mean
    }

    theta <- bind_tables(lapply(model$Z, posterior_z_base))

  } else if (model$model %in% c("hmm", "ldahmm")) {
    R <- model$stored_values$R_iter[[length(model$stored_values$R_iter)]] + 1L  # adjust index for R
    R <- R[model$model_settings$time_index]  # retrieve doc level state info
    alphas <- matrix(model$stored_values$alpha_iter[[length(model$stored_values$alpha_iter)]][R],
                     nrow = length(model$Z), ncol = info$allK)

    Z_table <- bind_tables(lapply(model$Z, function(zvec) {
      table_to_vector(table(factor(zvec, levels = 1:(info$allK) - 1L)))
    }))

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
  all_words <- model$vocab[as.integer(unlist(model$W, use.names = FALSE)) + 1L]
  all_topics <- as.integer(unlist(model$Z, use.names = FALSE))

  if (model$model %in% c("base", "cov", "hmm", "label")) {
    pi_estimated <- keyATM_output_pi(model$Z, model$S, model$priors$gamma)
    all_s <- as.integer(unlist(model$S, use.names = FALSE))

    obj <- keyATM_output_phi_calc_key(all_words, all_topics, all_s, pi_estimated,
                                      keywords_raw = model$keywords_raw,
                                      vocab = model$vocab,
                                      priors = model$priors,
                                      tnames = info$tnames,
                                      model = model)
  } else if (model$model %in% c("lda", "ldacov", "ldahmm")) {
    obj <- keyATM_output_phi_calc_lda(all_words, all_topics,
                                      model$vocab, model$priors$beta, info$tnames)
  }
  return(obj)
}


#' @noRd
#' @import magrittr
#' @importFrom rlang .data
keyATM_output_phi_calc_key <- function(all_words, all_topics, all_s, pi_estimated,
                                       keywords_raw, vocab, priors, tnames, model)
{
  res_tibble <- tibble::tibble(
                        Word = all_words,
                        Topic = all_topics,
                        Switch = all_s
                       )

  prob1 <- pi_estimated %>% dplyr::pull(.data$Proportion) / 100
  prob0 <- 1 - prob1
  vocab_sorted <- sort(vocab)

  if ("beta_s0" %in% names(priors)) {
    beta_s0 <- priors$beta_s0
    colnames(beta_s0) <- vocab
    beta_s0 <- beta_s0[, vocab_sorted]
  } else {
    beta_s0 <- priors$beta
  }

  all_keywords <- unique(unlist(model$keywords_raw, use.names = FALSE))
  beta_s1 <- matrix(priors$beta_s, nrow = length(model$keywords), ncol = length(all_keywords))
  colnames(beta_s1) <- sort(all_keywords)
  if ("beta_s1" %in% names(priors)) {
    for (k in 1:length(model$keywords_raw)) {
      keywords_k <- model$keywords_raw[[k]]

      for (keyword in keywords_k) {
        index <- match(keyword, vocab)
        beta_s1[k, keyword] <- priors$beta_s1[k, index]
      }
    }
  }

  get_phi <- function(res_tibble, switch_val, model)
  {
    if (switch_val == 0) {
      # Use no-keyword topic-word dist
      prior <- beta_s0
    } else if (switch_val == 1) {
      # Use keyword topic-word dist
      prior <- beta_s1
    }

    temp <- res_tibble %>%
              dplyr::filter(.data$Switch == switch_val) %>%
              dplyr::group_by(.data$Topic, .data$Word) %>%
              dplyr::summarize(Count = dplyr::n())

    temp %>% tidyr::spread(key = "Word", value = "Count") -> phi

    # Check unused topic
    if (nrow(phi) != length(tnames)) {
      missing <- setdiff(0:(length(tnames)-1L), phi$Topic)
      phi %>%
        dplyr::ungroup() %>%
        dplyr::add_row(Topic = missing) %>%
        dplyr::arrange(.data$Topic) -> phi
    }

    # Deal with NAs
    phi <- apply(phi, 2, function(x) {ifelse(is.na(x), 0, x)})

    if (!is.matrix(phi)) {
      phi <- t(phi)
    }

    phi <- phi[, 2:ncol(phi), drop = FALSE]
    topic_counts <- Matrix::rowSums(phi)

    rownames(phi) <- tnames[1:nrow(phi)]

    if (switch_val == 1) {
      # keyword topic-word dist
      phi_ <- phi
      all_keywords <- unique(unlist(model$keywords_raw, use.names = FALSE))
      phi <- matrix(0.0, nrow = length(model$keywords), ncol = length(all_keywords))
      colnames(phi) <- sort(all_keywords)

      for (k in 1:length(model$keywords_raw)) {
        phi[k, which(colnames(phi) %in% colnames(phi_))] <- phi_[k, ]
      }

      phi <- phi[, sort(colnames(phi)), drop = FALSE]
      for (k in 1:length(keywords_raw)) {
        phi[k, ] <- phi[k, ] + prior[k, ]
      }
      phi <- phi / Matrix::rowSums(phi)
      phi <- apply(phi, c(1,2), function(x) {ifelse(is.na(x), 0, x)})

      # keyword topic-word dist should have the same dimension as no-keyword dist
      # for marginilization, but no-keyword elements are 0
      phi_ <- matrix(0, nrow = length(tnames),
                     ncol = length(vocab))
      colnames(phi_) <- vocab_sorted
      phi_[1:nrow(phi), which(colnames(phi_) %in% colnames(phi))] <-
          phi[, which(colnames(phi) %in% colnames(phi_))]
      phi <- phi_
    } else {
      # no-keyword topic-word dist
      # Should have the same dimension as vocab
      phi_ <- matrix(0, nrow = length(tnames),
                     ncol = length(vocab))
      colnames(phi_) <- vocab_sorted
      phi <- phi[, sort(colnames(phi)), drop = FALSE]
      rownames(phi_) <- tnames

      # phi with all words
      phi_[, which(colnames(phi_) %in% colnames(phi))] <-
            phi[, which(colnames(phi) %in% colnames(phi_))]


      phi <- phi_ + prior  # add scalar (don't use label or matrix)
      phi <- phi / Matrix::rowSums(phi)
    }
    return(phi)
  }

  # Regular
  phi0 <- get_phi(res_tibble, switch_val = 0, model)

  # Keyword
  phi1 <- get_phi(res_tibble, switch_val = 1, model)

  # Marginal out switch
  blank_vec <- rep(0, length(vocab))
  names(blank_vec) <- vocab_sorted

  phi <- sapply(1:length(tnames),
                function(k) {
                  regular <- blank_vec
                  regular[colnames(phi0)] <- phi0[k, ]

                  key <- blank_vec
                  key[colnames(phi1)] <- phi1[k, ]

                  res <- regular * prob0[k] + key * prob1[k]
                  return(res)
                }) %>% t()
  colnames(phi) <- vocab_sorted  # same as colnames(phi0), colnames(phi1)
  rownames(phi) <- tnames

  topic_counts <- res_tibble %>%
                    dplyr::group_by(.data$Topic) %>%
                    dplyr::summarize(Count = dplyr::n()) %>%
                    dplyr::pull(.data$Count)

  word_counts <- res_tibble %>%
                    dplyr::group_by(.data$Word) %>%
                    dplyr::summarize(Count = dplyr::n()) %>%
                    dplyr::arrange(match(.data$Word, vocab)) %>%  # same order as vocab
                    dplyr::pull(.data$Count)

  if (ncol(phi) == length(vocab)) {
    phi <- phi[, vocab]
  } else {
    # This can happen in `by_strata_TopicWord`, does nothing
  }
  return(list(phi = phi, topic_counts = topic_counts, word_counts = word_counts))
}


#' @noRd
#' @import magrittr
#' @importFrom rlang .data
keyATM_output_phi_calc_lda <- function(all_words, all_topics, vocab, priors, tnames)
{
  res_tibble <- data.frame(
                        Word = all_words,
                        Topic = all_topics
                       ) %>%
                dplyr::group_by(.data$Topic, .data$Word) %>%
                dplyr::summarize(Count = dplyr::n())

  res_tibble %>%
    tidyr::spread(key = "Word", value = "Count") -> phi
  phi <- apply(phi, 2, function(x) {ifelse(is.na(x), 0, x)})

  phi <- phi[, 2:ncol(phi)]
  topic_counts <- Matrix::rowSums(phi)
  word_counts <- Matrix::colSums(phi)

  phi <- phi + priors

  if (ncol(phi) == length(vocab)) {
    phi <- phi[, vocab]
  } else {
    # This can happen in `by_strata_TopicWord`, does nothing
  }

  phi <- phi / Matrix::rowSums(phi)
  rownames(phi) <- tnames
  return(list(phi = phi, topic_counts = topic_counts, word_counts = word_counts))
}


#' @noRd
#' @import magrittr
keyATM_output_theta_iter <- function(model, info)
{
  if (model$model %in% c("cov", "ldacov")) {
    posterior_theta <- function(x) {
      Z_table <- model$stored_values$Z_tables[[x]]
      lambda <- model$stored_values$Lambda_iter[[x]]
      Alpha <- exp(model$model_settings$covariates_data_use %*% t(lambda))

      tt <- Z_table + Alpha
      row.names(tt) <- NULL

      return(tt / Matrix::rowSums(tt))
    }
  } else if (model$model %in% c("hmm", "ldahmm")) {
    posterior_theta <- function(x) {
      Z_table <- model$stored_values$Z_tables[[x]]
      R <- model$stored_values$R_iter[[x]] + 1L  # adjust index for R
      R <- R[model$model_settings$time_index]  # retrieve doc level state info

      alphas <- matrix(model$stored_values$alpha_iter[[x]][R],
                       nrow = length(model$Z), ncol = info$allK)

      tt <- Z_table + alphas
      theta <- tt / Matrix::rowSums(tt)
      return(theta)
    }
  } else {
    posterior_theta <- function(x) {
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
#' @importFrom rlang .data
keyATM_output_alpha_iter_base <- function(model, info)
{
  topics <- paste0(1:(info$allK))
  names(model$stored_values$alpha_iter) <- 1:length(model$stored_values$alpha_iter)
  model$stored_values$alpha_iter %>%
    dplyr::bind_rows() %>%
    t() %>%
    tibble::as_tibble(.name_repair = ~topics) %>%
    dplyr::mutate(Iteration = info$used_iter) %>%
    tidyr::gather(key = "Topic", value = "alpha", -"Iteration") %>%
    dplyr::mutate(Topic = as.integer(.data$Topic)) -> alpha_iter
  return(alpha_iter)
}


#' @noRd
#' @import magrittr
#' @importFrom rlang .data
keyATM_output_alpha_iter_hmm <- function(model, info)
{
  topics <- paste0(1:(info$allK))
  model$stored_values$alpha_iter %>%
    purrr::imap_dfr(function(x, i) {
                       x %>%
                         tibble::as_tibble(.name_repair = ~topics) %>%
                         dplyr::mutate(State = 1:(dplyr::n()),
                                       Iteration = info$used_iter[i]) %>%
                         tidyr::gather(key = "Topic", value = "alpha", -"State", -"Iteration")
                     }) %>%
     dplyr::mutate(Topic = as.integer(.data$Topic)) -> alpha_iter
  return(alpha_iter)
}


#' @noRd
#' @export
print.keyATM_output <- function(x, ...)
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
summary.keyATM_output <- function(object, ...)
{
  print(object)
}


#' Save a keyATM_output object
#'
#' @param x a keyATM_output object (see [keyATM()]).
#' @param file file name to create on disk.
#' @seealso [keyATM()], [weightedLDA()], [keyATMvb()]
#' @export
save.keyATM_output <- function(x, file = stop("'file' must be specified"))
{
  saveRDS(x, file = file)
}


#' Show the top words for each topic
#'
#' If \code{show_keyword} is \code{TRUE} then words in their keyword topics
#' are suffixed with a check mark. Words from another keyword topic
#' are labeled with the name of that category.
#'
#' @param x the output (see [keyATM()] and [by_strata_TopicWord()]).
#' @param n integer. The number terms to visualize. Default is \code{10}.
#' @param measure character. The way to sort the terms: \code{probability} (default) or \code{lift}.
#' @param show_keyword logical. If \code{TRUE}, mark keywords. Default is \code{TRUE}.
#'
#' @return An n x k table of the top n words in each topic
#' @export
top_words <- function(x, n = 10, measure = c("probability", "lift"),
                      show_keyword = TRUE)
{
  UseMethod("top_words")
}


#' @noRd
#' @export
top_words.strata_topicword <- function(x, n = 10, measure = c("probability", "lift"),
                                  show_keyword = TRUE)
{
  measure <- match.arg(measure)
  top_words <- lapply(x$phi,  # list of phis
                      function(obj) {
                       top_words_calc(
                         n = n, measure = measure, show_keyword = show_keyword,
                         theta = x$theta, phi = obj$phi,
                         word_counts = obj$word_counts, keywords_raw = x$keywords_raw
                       )
                      })
  return(top_words)
}


#' @noRd
#' @export
top_words.keyATM_output <- function(x, n = 10, measure = c("probability", "lift"),
                                    show_keyword = TRUE)
{
  check_arg_type(x, "keyATM_output")
  modelname <- extract_full_model_name(x)
  measure <- match.arg(measure)

  if (modelname %in% c("lda", "ldacov", "ldahmm"))
     show_keyword <- FALSE

  res <- top_words_calc(n, measure, show_keyword,
                        theta = x$theta, phi = x$phi,
                        word_counts = x$word_counts, keywords_raw = x$keywords_raw)
  return(res)
}


top_words_calc <- function(n, measure, show_keyword,
                           theta, phi, word_counts, keywords_raw)
{
  if (is.null(n))
    n <- nrow(theta)
  if (measure == "probability") {
     measuref <- function(xrow) {
       colnames(phi)[order(xrow, decreasing = TRUE)[1:n]]
     }
  } else if (measure == "lift") {
     wfreq <- word_counts / sum(word_counts)
     measuref <- function(xrow) {
       colnames(phi)[order(xrow / wfreq, decreasing = TRUE)[1:n]]
     }
  }
  res <- apply(phi, 1, measuref)

  if (show_keyword) {
    for (i in 1:ncol(res)) {
      for (j in 1:length(keywords_raw)) {
        # Potential words to check in the top words
        if (i <= length(keywords_raw)) {
          inds_all <- which(res[, i] %in% c(keywords_raw[[i]], keywords_raw[[j]]))
        } else {
          inds_all <- which(res[, i] %in% keywords_raw[[j]])
        }
        display <- res[inds_all, i]

        # If there is a keyword from other topics
        inds <- which(res[inds_all, i] %in% keywords_raw[[j]])
        label <- paste0("[", as.character(j), "]")
        display[inds] <- paste(res[inds_all, i][inds], label)

        # Keywords of topic i, overwriting above
        if (i <= length(keywords_raw)) {
          inds <- which(res[inds_all, i] %in% keywords_raw[[i]])
          display[inds] <- paste(res[inds_all, i][inds], paste0("[", "\U2713" ,"]"))
        }

        # Put it back to the original table
        res[inds_all, i] <- display
      }
    }
  }
  res <- as.data.frame(res, stringsAsFactors = FALSE)

  return(res)
}


#' Show the top topics for each document
#'
#' @param x the output from a keyATM model (see [keyATM()]).
#' @param n integer. The number of topics to show. Default is \code{2}.
#'
#' @return An n x k table of the top n topics in each document.
#' @import magrittr
#' @export
#'
top_topics <- function(x, n = 2)
{
  check_arg_type(x, "keyATM_output")
  check_arg_type(n, "numeric")

  if (n > ncol(x$theta))
    n <- ncol(x$theta)

  measuref <- function(xrow) {
    colnames(x$theta)[order(xrow, decreasing = TRUE)[1:n]]
  }

  res <- t(apply(x$theta, 1, measuref)) %>%
          tibble::as_tibble(.name_repair = ~paste0("Rank", 1:n))
  return(res)
}


#' Show the top documents for each topic
#'
#' @param x the output from a keyATM model (see [keyATM()]).
#' @param n How many documents to show. Default is \code{10}.
#'
#' @return An n x k table of the top n documents for each topic, each number is a document index.
#' @import magrittr
#' @export
top_docs <- function(x, n = 10)
{
  check_arg_type(x, "keyATM_output")
  if (is.null(n))
    n <- nrow(x$theta)

  measuref <- function(xcol) {
    order(xcol, decreasing = TRUE)[1:n]
  }

  res <- apply(x$theta, 2, measuref) %>% as.data.frame()
  return(res)
}


#' Estimate subsetted topic-word distribution
#'
#' @param x the output from a keyATM model (see [keyATM()]).
#' @param keyATM_docs an object generated by [keyATM_read()].
#' @param by a vector whose length is the number of documents.
#'
#' @return strata_topicword object (a list).
#' @import magrittr
#' @export
by_strata_TopicWord <- function(x, keyATM_docs, by)
{
  # Check inputs
  if (!is.vector(by)) {
    stop("`by` should be a vector.")
  }
  if (!"Z" %in% names(x$kept_values)) {
    stop("`Z` and `S` should be in the output. Please check `keep` option in `keyATM()`.")
  }
  if (!"S" %in% names(x$kept_values)) {
    stop("`Z` and `S` should be in the output. Please check `keep` option in `keyATM()`.")
  }
  if (length(keyATM_docs) != length(by)) {
    stop("The length of `by` should be the same as the length of documents.")
  }


  # Get unique values of `by`
  unique_val <- unique(by)
  tnames <- rownames(x$phi)

  # Get phi for each
  obj <- lapply(unique_val,
                function(val) {
                  doc_index <- which(by == val)
                  all_words <- unlist(keyATM_docs[doc_index], use.names = FALSE)
                  all_topics <- as.integer(unlist(x$kept_values$Z[doc_index]), use.names = FALSE)
                  all_s <- as.integer(unlist(x$kept_values$S[doc_index]), use.names = FALSE)
                  pi_estimated <- keyATM_output_pi(x$kept_values$Z[doc_index],
                                                   x$kept_values$S[doc_index],
                                                   x$priors$gamma)
                  vocab <- sort(unique(all_words))
                  phi_obj <- keyATM_output_phi_calc_key(all_words, all_topics, all_s, pi_estimated,
                                                        x$keywords_raw,
                                                        vocab, x$priors, tnames, model = x)
                }
               )
  names(obj) <- unique_val

  res <- list(phi = obj, theta = x$theta, keywords_raw = x$keywords_raw)

  class(res) <- c("strata_topicword", class(res))
  return(res)
}


#' Show covariates information
#'
#' @param x the output from the covariate keyATM model (see [keyATM()]).
#' @export
covariates_info <- function(x) {
  if (x$model != "covariates" | !("keyATM_output" %in% class(x))) {
    stop("This is not an output of the covariate model")
  }
  cat(paste0("Colnames: ", paste(colnames(x$kept_values$model_settings$covariates_data_use), collapse = ", "),
             "\nStandardization: ", as.character(x$kept_values$model_settings$standardize),
             "\nFormula: ", paste(as.character(x$kept_values$model_settings$covariates_formula), collapse = " "), "\n\nPreview:\n"))
  print(utils::head(x$kept_values$model_settings$covariates_data_use))
}


#' Return covariates used in the iteration
#'
#' @param x the output from the covariate keyATM model (see [keyATM()])
#' @export
covariates_get <- function(x) {
  if (x$model != "covariates" | !("keyATM_output" %in% class(x))) {
    stop("This is not an output of the covariate model")
  }
  return(x$kept_values$model_settings$covariates_data_use)
}


#' Estimate document-topic distribution by strata (for covariate models)
#'
#' @param x the output from the covariate keyATM model (see [keyATM()]).
#' @param by_var character. The name of the variable to use.
#' @param labels character. The labels for the values specified in `by_var` (ascending order).
#' @param by_values numeric. Specific values for `by_var`, ordered from small to large. If it is not specified, all values in `by_var` will be used.
#' @param ... other arguments passed on to the [predict()] function.
#' @return strata_topicword object (a list).
#' @import magrittr
#' @importFrom stats predict
#' @export
by_strata_DocTopic <- function(x, by_var, labels, by_values = NULL, ...)
{
  # Check inputs
  variables <- colnames(x$kept_values$model_settings$covariates_data_use)
  if (length(by_var) != 1)
    stop("`by_var` should be a single variable.")
  if (!by_var %in% variables)
    stop(paste0(by_var, " is not in the set of covariates in keyATM model. Check with `covariates_info()`.",
                "Covariates provided are: ",
                paste(colnames(x$kept_values$model_settings$covariates_data_use), collapse=" , ")))

  # Info
  if (is.null(by_values)) {
    by_values <- sort(unique(x$kept_values$model_settings$covariates_data_use[, by_var]))
  }
  if (length(by_values) != length(labels)) {
    stop("Length mismatches. Please check `labels`.")
  }

  set.seed(x$options$seed)
  res <- lapply(1:length(by_values),
                function(i) {
                  value <- by_values[i]
                  new_data <- x$kept_values$model_settings$covariates_data_use
                  new_data[, by_var] <- value
                  obj <- predict(x, new_data, raw_values = TRUE, ...)
                })

  names(res) <- by_values
  obj <- list(theta = tibble::as_tibble(res), by_values = by_values, by_var = by_var, labels = labels)
  class(obj) <- c("strata_doctopic", class(obj))
  return(obj)
}


#' @noRd
#' @export
print.strata_doctopic <- function(x, ...)
{
  cat(paste0("strata_doctopic object for the ", x$by_var, "\n"))
}


#' @noRd
#' @export
summary.strata_doctopic <- function(object, ci = 0.9, method = "hdi", point = "mean", ...)
{
  tables <- lapply(1:length(object$by_values),
                  function(index) {
                     theta <- object$theta[[index]]
                     return(strata_doctopic_CI(theta[, 1:(ncol(theta)-1)],
                                               ci, method, point,
                                               label = object$labels[index]))
                  })
  names(tables) <- object$labels
  return(tables)
}


#' @noRd
#' @importFrom rlang .data
#' @import magrittr
#' @keywords internal
strata_doctopic_CI <- function(theta, ci, method, point, label)
{
  q <- as.data.frame(apply(theta, 2, calc_ci, ci, method, point))
  q$CI <- c("Lower", "Point", "Upper")
  q %>%
    tidyr::gather(key = "Topic", value = "Value", -"CI") %>%
    tidyr::spread(key = "CI", value = "Value") %>%
    dplyr::mutate(TopicId = 1:(dplyr::n())) %>% tibble::as_tibble() -> res

  if (!is.null(label)) {
    res %>% dplyr::mutate(label = label) -> res
  }
  return(res)
}


#' @noRd
#' @keywords internal
calc_ci <- function(vec, ci, method, point)
{
  # Check bayestestR package and Kruschke's book
  if (point == "mean") {
    point <- mean(vec)
  }
  if (point == "median") {
    point <- stats::median(vec)
  }

  if (method == "hdi") {
    sorted_points <- sort.int(vec, method = "quick")
    window_size <- ceiling(ci * length(sorted_points))
    nCIs <-  length(sorted_points) - window_size

    if (window_size < 2 | nCIs < 1) {
      warning("`ci` is too small or interations are not enough, using `eti` option instead.")
      return(calc_ci(vec, ci, method = "eti", point))
    }

    ci_width <- sapply(1:nCIs, function(x) {sorted_points[x + window_size] - sorted_points[x]})
    slice_min <- which.min(ci_width)

    res <- c(sorted_points[slice_min], point, sorted_points[slice_min + window_size])
  }

  if (method == "eti") {
    qua <- stats::quantile(vec, probs = c((1 - ci) / 2, (1 + ci) / 2))
    res <- c(qua[1], point, qua[2])
  }

  names(res) <- c("Lower", "Point", "Upper")
  return(res)
}

