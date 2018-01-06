#' Get posterior quantities from model output
#'
#' Constructs a (N x K) matrix \code{theta} and (K x V) matrix \code{beta}
#' plus their margins from the sample of Z and W in \code{model}.
#' These statistics implicitly marginalize over X.
#'
#' @param model a model, fitted or simply initialized
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
  allK <- model$extra_k + length(model$dict)
  V <- length(model$vocab)

  tNZ <- do.call(rbind,
                 lapply(model$Z, function(x){ table(factor(x, levels = 1:allK - 1)) }))
  rownames(tNZ) <- basename(model$files)
  doc_lens <- rowSums(tNZ)
  tNZ <- tNZ / doc_lens

  tZW <- Reduce(`+`,
                 mapply(function(a, b){ table(factor(a, levels = 1:allK - 1),
                                              factor(b, levels = 1:V - 1)) },
                        model$Z, model$W, SIMPLIFY = FALSE))
  colnames(tZW) <- model$vocab
  word_counts <- colSums(tZW)
  topic_counts <- rowSums(tZW)
  tZW <- tZW / topic_counts

  ll <- list(seed_K = length(model$dict), extra_K = model$extra_k,
             V = ncol(tZW), N = nrow(tNZ),
             theta = tNZ, beta = as.matrix(as.data.frame.matrix(tZW)),
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = doc_lens, vocab = model$vocab)
  class(ll) <- c("seededlda_posterior", class(ll))
  ll
}

#' Suggest composite names for each topic
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param n How many topic terms to use in the name: default 2
#' @param measure Method to find topics for new names. See \code{top_terms}
#'
#' @return A vector of new topic names constructed from top terms
#' @export
suggest_topic_names <- function(x, n = 3,
                                measure = c("probability", "lift")){
  tt <- top_terms(x, n, measure)
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
  rownames(x$theta) <- doc_names
  names(x$doc_lens) <- doc_names
  x
}

#' Show the top terms for each topic
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param n How many terms to show: default NULL shows all
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n words in each topic
#' @export
#'
top_terms <- function(x, n = 10,
                      measure = c("probability", "lift")){
  if (is.null(n))
    n <- nrow(x$theta)
  measure <- match.arg(measure)
  if (measure == "probability"){
     measuref <- function(xrow){
       colnames(x$beta)[order(xrow, decreasing = TRUE)[1:n]]
     }
  } else if (measure == "lift"){
     wfreq <- x$word_counts / sum(x$word_counts)
     measuref <- function(xrow){
       colnames(x$beta)[order(xrow / wfreq, decreasing = TRUE)[1:n]]
     }
  }
  apply(x$beta, 1, measuref)
}

#' Show the top topics for each document
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param n How many topics to show: default 2
#' @param measure How to sort the topics: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n topics in each document
#' @export
#'
top_topics <- function(x, n = 2,
                       measure = c("probability", "lift")){
  if (is.null(n))
    n <- nrow(x$theta)

  measure <- match.arg(measure)
  if (measure == "probability"){
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
#' @param n How many documents to show: default 10
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n documents for each topic
#' @export
top_docs <- function(x, n = 10, measure = c("probability", "lift")){
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

