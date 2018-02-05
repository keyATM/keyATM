#' Get posterior quantities from model output
#'
#' Constructs a (N x K) matrix \code{theta} and (K x V) matrix \code{beta}
#' plus their margins from the sample of Z and W in \code{model}.
#' These statistics implicitly marginalize over X.
#'
#' @param model a topicdict model, fitted or simply initialized
#'
#' @return a list with elements
#'   \itemize{
#'     \item{seed_K}{ Number of seeded topics}
#'     \item{extra_K}{ Number of regular unseeded topics}
#'     \item{V}{ Number of word types}
#'     \item{N}{ Number of documents}
#'     \item{theta}{ Normalized tpoic proportions for each document}
#'     \item{beta}{ Normalized topic specific word generation probabilities}
#'     \item{topic_counts}{ Number of tokens assigned to each topic}
#'     \item{word_counts}{ Number of times each word type appears}
#'     \item{doc_lens}{ Length of each document in tokens}
#'     \item{vocab}{ Words in the vocabulary}
#'   }
#' @export
posterior <- function(model){
  check_arg_type(model, "topicdict")
  allK <- model$extra_k + length(model$dict)
  V <- length(model$vocab)
  tnames <- c(names(model$seeds), paste0("T_", 1:model$extra_k))

  tNZ <- do.call(rbind,
                 lapply(model$Z, function(x){ table(factor(x, levels = 1:allK - 1)) }))
  rownames(tNZ) <- basename(model$files)
  doc_lens <- rowSums(tNZ)
  tNZ <- tNZ / doc_lens
  colnames(tNZ) <- tnames # label seeded topics

  tZW <- Reduce(`+`,
                 mapply(function(a, b){ table(factor(a, levels = 1:allK - 1),
                                              factor(b, levels = 1:V - 1)) },
                        model$Z, model$W, SIMPLIFY = FALSE))
  colnames(tZW) <- model$vocab
  word_counts <- colSums(tZW)
  topic_counts <- rowSums(tZW)
  tZW <- tZW / topic_counts
  rownames(tZW) <- tnames

  ## TODO fix this naming nonsense
  dict <- model$dict
  names(dict) <- names(model$seeds)

  ll <- list(seed_K = length(model$dict), extra_K = model$extra_k,
             V = ncol(tZW), N = nrow(tNZ),
             theta = tNZ, beta = as.matrix(as.data.frame.matrix(tZW)),
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = doc_lens, vocab = model$vocab,
             dict = dict)
  class(ll) <- c("topicdict_posterior", class(ll))
  ll
}

# a more than usually informative error message for handing in the
# wrong type to a function
check_arg_type <- function(arg, typename){
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename))
    stop(paste("'", argname, '" is not a ', typename))
}

#' Suggest composite names for each topic
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param measure Method to find topics for new names. See \code{top_terms}
#' @param n How many topic terms to use in the name: default 2
#'
#' @return A vector of new topic names constructed from top terms
#' @export
suggest_topic_names <- function(x, measure = c("probability", "lift"), n = 3){
  check_arg_type(x, "topicdict_posterior")
  meas <- match.arg(measure)
  tt <- top_terms(x, measure = meas, n = n)
  apply(tt, 2, function(x){ paste(x, collapse = "-") })
}

#' Set topic names
#'
#' @param x Posterior from a seededlda model (see \code{posterior})
#' @param topic_names new names for topics
#'
#' @return a posterior object with new topic names in its components
#' @export
#'
set_topic_names <- function(x, topic_names){
  check_arg_type(x, "topicdict_posterior")
  colnames(x$theta) <- topic_names
  names(x$topic_counts) <- topic_names
  rownames(x$beta) <- topic_names
  x
}

#' Set document names
#'
#' @param x Posterior from a seededlda model (see \code{posterior})
#' @param doc_names new names for documents
#'
#' @return a posterior object with new document names in its components
#' @export
#'
set_doc_names <- function(x, doc_names){
  check_arg_type(x, "topicdict_posterior")
  rownames(x$theta) <- doc_names
  names(x$doc_lens) <- doc_names
  x
}

#' Show the top terms for each topic
#'
#' If \code{show_seed} is true then words in their seeded categories
#' are suffixed with a check mark. Words from another seeded category
#' are labeled with the name of that category.
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#' @param n How many terms to show. Default: NULL, which shows all
#' @param show_seed Mark seeded vocabulary. See below for details.
#'
#' @return An n x k table of the top n words in each topic
#' @export
#'
top_terms <- function(x, measure = c("probability", "lift"), n = 10,
                      show_seed = FALSE){
  check_arg_type(x, "topicdict_posterior")
  if (is.null(n))
    n <- nrow(x$theta)
  measure <- match.arg(measure)
  if (measure == "probability") {
     measuref <- function(xrow){
       colnames(x$beta)[order(xrow, decreasing = TRUE)[1:n]]
     }
  } else if (measure == "lift") {
     wfreq <- x$word_counts / sum(x$word_counts)
     measuref <- function(xrow){
       colnames(x$beta)[order(xrow / wfreq, decreasing = TRUE)[1:n]]
     }
  }
  res <- apply(x$beta, 1, measuref)
  if (show_seed) {
    for (i in 1:ncol(res)) {
      for (j in 1:length(x$dict)) {
         inds <- which(res[,i] %in% x$dict[[j]])
         label <- ifelse(i == j,
                         paste0("[", "\U2713" ,"]"),
                         paste0("[", names(x$dict)[j], "]"))
         res[inds, i] <- paste(res[inds, i], label)
      }
    }
  }
  res
}

#' Show the top topics for each document
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param measure How to sort the topics: 'probability' (default) or 'lift'
#' @param n How many topics to show. Default: 2
#'
#' @return An n x k table of the top n topics in each document
#' @export
#'
top_topics <- function(x, measure = c("probability", "lift"), n = 2){
  check_arg_type(x, "topicdict_posterior")
  if (is.null(n))
    n <- nrow(x$theta)

  measure <- match.arg(measure)
  if (measure == "probability") {
    measuref <- function(xrow){
      colnames(x$theta)[order(xrow, decreasing = TRUE)[1:n]]
    }
  } else if (measure == "lift"){
    wfreq <- x$topic_counts / sum(x$topic_counts)
    measuref <- function(xrow){
      colnames(x$theta)[order(xrow / wfreq, decreasing = TRUE)[1:n]]
    }
  }
  t(apply(x$theta, 1, measuref))
}

#' Show the top documents for each topic
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#' @param n How many documents to show. Default: 10
#'
#' @return An n x k table of the top n documents for each topic
#' @export
top_docs <- function(x, measure = c("probability", "lift"), n = 10){
  check_arg_type(x, "topicdict_posterior")
  if (is.null(n))
    n <- nrow(x$theta)

  measure <- match.arg(measure)
  if (measure == "probability"){
    measuref <- function(xcol){
      rownames(x$theta)[order(xcol, decreasing = TRUE)[1:n]]
    }
    apply(x$theta, 2, measuref)
  } else if (measure == "lift"){
    tfreq <- x$topic_counts / sum(x$topic_counts)
    measuref <- function(xcol){
      rownames(x$theta)[order(xcol, decreasing = TRUE)[1:n]]
    }
    apply(x$theta / outer(1:x$N, tfreq), 2, measuref)
  }
}

