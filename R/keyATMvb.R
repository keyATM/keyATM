#' keyATM with Collapsed Variational Bayes
#'
#' **Experimental feature:** Fit keyATM base with Collapsed Variational Bayes
#'
#' @param docs texts read via [keyATM_read()]
#' @param model keyATM model: \code{base}, \code{covariates}, and \code{dynamic}
#' @param no_keyword_topics the number of regular topics
#' @param keywords a list of keywords
#' @param model_settings a list of model specific settings (details are in the online documentation)
#' @param vb_options a list of settings for Variational Bayes \itemize{
#'           \item \strong{convtol}: the default is \code{1e-4}
#'           \item \strong{init}: \code{mcmc} (default) or \code{random}
#'           }
#' @param priors a list of priors of parameters
#' @param options a list of options same as [keyATM()]. Options are used when initialization method is \code{mcmc}.
#' @param keep a vector of the names of elements you want to keep in output
#' @return A \code{keyATM_output} object
#' @seealso \url{https://keyatm.github.io/keyATM/articles/pkgdown_files/keyATMvb.html}
#' @export
keyATMvb <- function(
  docs,
  model,
  no_keyword_topics,
  keywords = list(),
  model_settings = list(),
  vb_options = list(),
  priors = list(),
  options = list(),
  keep = list()
) {
  cli::cli_alert_info(
    "keyATMvb is an experimental function. DO NOT USE THIS FOR THE FINAL RESULT."
  )
  # Check type
  if (length(keep) != 0) {
    check_arg_type(keep, "character")
  }

  model <- full_model_name(model, type = "keyATM")
  if (is.null(options$iterations)) {
    options$iterations <- 100
  }
  if (is.null(options$seed)) {
    options$seed <- floor(stats::runif(1) * 1e5)
  }
  set.seed(options$seed)

  # Initialize and fit keyATM (if MCMC initialization)
  fitted <- keyATMvb_fit(
    docs,
    model,
    no_keyword_topics,
    keywords,
    model_settings,
    vb_options,
    priors,
    options
  )

  # Get output
  used_iter <- get_used_iter(fitted, resume = FALSE)
  out <- keyATM_output(fitted, keep, used_iter)
  # Add VB options
  out$vb_options <- fitted$vb_options

  return(out)
}


#' Fit a keyATM model with Collapsed Variational Bayes
#'
#' @keywords internal
keyATMvb_fit <- function(
  docs,
  model,
  no_keyword_topics,
  keywords = list(),
  model_settings = list(),
  vb_options = list(),
  priors = list(),
  options = list()
) {
  #
  # Check
  #
  check_arg_type(
    docs,
    "keyATM_docs",
    "Please use `keyATM_read()` to read texts."
  )
  if (!is.integer(no_keyword_topics) & !is.numeric(no_keyword_topics)) {
    cli::cli_abort("`no_keyword_topics` is neigher numeric nor integer.")
  }

  no_keyword_topics <- as.integer(no_keyword_topics)

  if (!model %in% c("base", "cov", "hmm", "lda", "ldacov", "ldahmm")) {
    cli::cli_abort("Please select a correct model.")
  }

  info <- list(
    models_keyATM = c("base", "cov", "hmm"),
    models_lda = c("lda", "ldacov", "ldahmm")
  )
  keywords <- check_arg(keywords, "keywords", model, info)

  # Get Info
  if (is.null(docs$doc_index)) {
    info$use_doc_index <- get_doc_index(docs$W_raw)
  } else {
    info$use_doc_index <- docs$doc_index
    if (length(docs$doc_index) != length(docs$W_raw)) {
      cli::cli_warn(
        "Some documents have 0 length. Please review the preprocessing steps."
      )
    }
  }
  docs$W_raw <- docs$W_raw[info$use_doc_index]
  info$num_doc <- length(docs$W_raw)
  info$keyword_k <- length(keywords)
  info$total_k <- length(keywords) + no_keyword_topics
  info$num_core <- max(
    1,
    parallel::detectCores(all.tests = FALSE, logical = TRUE) - 2L
  )

  # Set default values
  model_settings <- check_arg(model_settings, "model_settings", model, info)
  priors <- check_arg(priors, "priors", model, info)
  options <- check_arg(options, "options", model, info)
  vb_options <- check_arg(vb_options, "vb_options", model, info)
  info$parallel_init <- options$parallel_init

  #
  # Initialization
  #
  set.seed(options$seed)

  ## Initialize W, Z, S
  if (is.null(docs$wd_names)) {
    info$wd_names <- unique(unlist(
      docs$W_raw,
      use.names = FALSE,
      recursive = FALSE
    ))
    check_vocabulary(info$wd_names)
  } else {
    info$wd_names <- docs$wd_names
  }

  # Check keywords
  keywords <- check_keywords(info$wd_names, keywords, options$prune)
  keywords_raw <- keywords # keep raw keywords (not word_id)
  info$keywords_raw <- keywords_raw

  # Create empty W (placeholder)
  initialized <- list()
  initialized$W <- lapply(docs$W_raw, function(x) {
    rep(-1L, length(x))
  })
  initialized$Z <- rlang::duplicate(initialized$W)
  initialized$keywords_id <- lapply(keywords, function(x) {
    rep(-1L, length(x))
  })
  initialized$S <- rlang::duplicate(initialized$W)

  if (model %in% info$models_keyATM) {
    initialized$model_key <- 1L
  } else {
    initialized$model_key <- 0L
  }

  initialized <- make_wsz_cpp(docs$W_raw, info, initialized)
  W <- initialized$W
  Z <- initialized$Z
  keywords_id <- initialized$keywords_id
  S <- initialized$S

  # Organize
  stored_values <- list(
    vocab_weights = rep(-1, length(info$wd_names)),
    doc_index = info$use_doc_index,
    keyATMdoc_meta = docs[-c(1, 3)]
  )

  if (model %in% c("base")) {
    if (options$estimate_alpha) {
      stored_values$alpha_iter <- list()
    }
  }

  if (model %in% c("hmm")) {
    options$estimate_alpha <- 1
    stored_values$alpha_iter <- list()
  }

  if (model %in% c("cov")) {
    stored_values$Lambda_iter <- list()
  }

  if (model %in% c("hmm")) {
    stored_values$R_iter <- list()

    if (options$store_transition_matrix) {
      stored_values$P_iter <- list()
    }
  }

  if (model %in% info$models_keyATM) {
    if (options$store_pi) {
      stored_values$pi_vectors <- list()
    }
  }

  if (options$store_theta) {
    stored_values$Z_tables <- list()
  }

  key_model <- list(
    W = W,
    Z = Z,
    S = S,
    model = abb_model_name(model),
    keywords = keywords_id,
    keywords_raw = keywords_raw,
    no_keyword_topics = no_keyword_topics,
    keyword_k = length(keywords_raw),
    vocab = info$wd_names,
    model_settings = model_settings,
    priors = priors,
    options = options,
    vb_options = vb_options,
    stored_values = stored_values,
    model_fit = list(),
    call = match.call()
  )
  rm(info)
  class(key_model) <- c("keyATM_model", model, class(key_model))

  keyATM_initialized <- list(
    model = key_model,
    model_name = model,
    options = options
  )
  class(keyATM_initialized) <- c(
    "keyATM_initialized",
    class(keyATM_initialized)
  )

  # MCMC initialization
  if (vb_options$init == "mcmc") {
    cli::cli_inform("Initializing with MCMC..")
    key_model <- keyATM_fit(keyATM_initialized)
  } else {
    # Do nothing
    keyATM_initialized$model$options$iterations <- 1L
    keyATM_initialized$model$options$iter_new <- 1L
    key_model <- keyATM_fit(keyATM_initialized)
  }

  # Fit VB
  set.seed(key_model$options$seed)
  cli::cli_inform("Fitting Variational Bayes...")
  key_model$vb_options$Perplexity_VB <- list(value = list(), iter = list())
  key_model <- keyATMvb_call(key_model)
  class(key_model) <- c("keyATM_fitted_VB", class(key_model))

  return(key_model)
}


check_arg_vboptions <- function(obj, model, info) {
  check_arg_type(obj, "list")
  allowed_arguments <- c("convtol", "init")

  if (is.null(obj$convtol)) {
    obj$convtol <- 1e-4
  } else {
    if (!is.numeric(obj$convtol)) {
      cli::cli_abort("`vb_options$convtol` should be a numeric value.")
    }
    if (obj$convtol <= 0) {
      cli::cli_abort("`vb_options$convtol` should be a positive value.")
    }
  }

  if (is.null(obj$init)) {
    obj$init <- "mcmc"
    if (!obj$init %in% c("mcmc", "random")) {
      cli::cli_abort("`vb_options$init` should be `mcmc` or `random`.")
    }
  }

  show_unused_arguments(obj, "`vb_options`", allowed_arguments)
  return(obj)
}
