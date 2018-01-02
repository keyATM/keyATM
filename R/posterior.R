#' Get posterior quantities from model output
#'
#' @param model_results The result of running \code{seededlda}
#' @param docnames Names of documents (not stored in the C++
#'                 seeededlda yet)
#'
#' @return a list with elements
#'   \itemize{
#'     \item{topics}{An N x K matrix of topic proportions}
#'     \item{terms}{A K x V matrix of topic-specific term generation
#'                  probabilities}
#'     \item{doclens}{A vector of document lengths in words, with which
#'                    Z can be reinflated from theta}
#'   }
#' @export
posterior <- function(x, docnames = NULL){
  K <- Reduce(max, Map(max, x$Z)) # 0:K
  V <- Reduce(max, Map(max, x$WordIDs)) # 0:V
  N <- length(x$WordIDs)
  vocab <- unique(data.frame(word = unlist(x$RawWords),
                             id = unlist(x$WordIDs)))
  vocab <- vocab[order(vocab$id),] # sorted by id
  vocab_freq <- as.numeric(table(unlist(x$WordIDs)))

  Z <- do.call(rbind,
               lapply(x$Z, function(x){ table(factor(x, levels = 0:K)) }))
  theta <- Z / rowSums(Z)
  Ns <- unlist(lapply(x$WordIDs, length))
  if (!is.null(docnames))
    rownames(theta) <- names(Ns) <- docnames
  else {
    rownames(theta) <- 1:N
    names(Ns) <- 1:N
  }
  z_freq <- colSums(Z)
  beta <- Reduce(`+`,
                 mapply(function(a, b){ table(factor(a, levels = 0:K),
                                              factor(b, levels = 0:V)) },
                        x$Z, x$WordIDs, SIMPLIFY = FALSE))
  beta <- beta / rowSums(beta)
  colnames(beta) <- vocab$word

  ll <- list(K = K,
             V = V,
             N = N,
             theta = theta,
             z_freq = z_freq,
             beta = beta,
             doclens = Ns,
             vocab = vocab,
             vocab_freq = vocab_freq)
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
  apply(tt, 2, function(x){ paste(x, collapse="-") })
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
  names(x$z_freq) <- topic_names
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
  names(x$doclens) <- doc_names
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
  measure <- match.arg(measure)
  if (is.null(n))
    n <- nrow(x$theta)

  if (measure == "probability"){
     measuref <- function(xrow){
       colnames(x$beta)[order(xrow, decreasing = TRUE)[1:n]]
     }
  } else if (measure == "lift"){
     wfreq <- x$vocab_freq / sum(x$vocab_freq)
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
  measure <- match.arg(measure)
  if (is.null(n))
    n <- nrow(x$theta)

  if (measure == "probability"){
    measuref <- function(xrow){
      colnames(x$theta)[order(xrow, decreasing = TRUE)[1:n]]
    }
  } else if (measure == "lift"){
    wfreq <- x$topic_freq / sum(x$topic_freq)
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
  measure <- match.arg(measure)
  if (is.null(n))
    n <- nrow(x$theta)

  if (measure == "probability"){
    measuref <- function(xcol){
      rownames(x$theta)[order(xcol, decreasing = TRUE)[1:n]]
    }
    apply(x$theta, 2, measuref)
  } else if (measure == "lift"){
    tfreq <- x$z_freq / sum(x$z_freq)
    measuref <- function(xcol){
      rownames(x$theta)[order(xcol, decreasing = TRUE)[1:n]]
    }
    apply(x$theta / outer(1:x$N, tfreq), 2, measuref)
  }
}

