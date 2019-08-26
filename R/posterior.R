#' Get posterior (deprecated)
#'
#' @export
posterior <- function(...){
  message("Warning: `posterior` is deprecated, please use `keyATM_posterior` instead.")
  return(keyATM_posterior(...))
}


#' Get posterior quantities from model output
#'
#' Constructs a (N x K) matrix \code{theta} and (K x V) matrix \code{beta}
#' plus their margins from the sample of Z and W in \code{model}.
#' These statistics implicitly marginalize over X.
#'
#' @param model a fitted keyATM model
#'
#' @return A list containing:
#'   \describe{
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
#'     \item{alpha}{ \code{alpha} during the iteration}
#'     \item{modelfit}{ Perplexity and log-likelihood}
#'     \item{p}{ Estimated p}
#'   }
#' @export
keyATM_posterior <- function(model){
  message("Creating a posterior object. It may take time...")

  check_arg_type(model, "keyATM")
  allK <- model$extra_k + length(model$dict)
  V <- length(model$vocab)
  N = length(model$W)
  doc_lens <- sapply(model$W, length)

  if(model$extra_k > 0){
    tnames <- c(names(model$keywords), paste0("T_", 1:model$extra_k))
  }else{
    tnames <- c(names(model$keywords))
  }

  posterior_z <- function(zvec){
    tt <- table(factor(zvec, levels = 1:allK - 1))
    # (tt + model$alpha) / (sum(tt) + sum(model$alpha)) # posterior mean
    (tt) / (sum(tt)) # posterior mean
  }
  theta <- do.call(rbind, lapply(model$Z, posterior_z))
  rownames(theta) <- basename(model$files)
  colnames(theta) <- tnames # label seeded topics


  all_words <- model$vocab[as.integer(unlist(model$W)) + 1]
  all_topics <- as.integer(unlist(model$Z))
  
  res_tibble <- tibble(
                        Word = all_words,
                        Topic = all_topics
                       ) %>%
                group_by(Topic, Word) %>%
                summarize(Count = n())
  
  res_tibble %>%
    tidyr::spread(key=Word, value=Count)  -> beta
  beta <- apply(beta, 2, function(x){ifelse(is.na(x), 0, x)})
  beta <- beta[, 2:ncol(beta)]
  beta <- beta[, model$vocab]

  topic_counts <- Matrix::rowSums(beta)
  word_counts <- Matrix::colSums(beta)

  tZW <- beta / topic_counts
  rownames(tZW) <- tnames

  # alpha
  res_alpha <- data.frame(model$alpha_iter)
  colnames(res_alpha) <- NULL
  res_alpha <- data.frame(t(res_alpha))
  if(nrow(res_alpha) > 0){
    colnames(res_alpha) <- paste0("EstTopic", 1:ncol(res_alpha))
    res_alpha$iter <- 1:nrow(res_alpha)
  }

  # model fit
  modelfit <- data.frame(model$model_fit)
  colnames(modelfit) <- NULL
  if(nrow(modelfit) > 0){
    modelfit <- data.frame(t(modelfit))
    colnames(modelfit) <-  c("Iteration", "Log Likelihood", "Perplexity")
  }

  # p
  collapse <- function(obj){
  temp <- unlist(obj)
  names(temp) <- NULL
  return(temp)
  }

  data <- data.frame(Z=collapse(model$Z), X=collapse(model$X))
  data %>%
    mutate_(Topic='Z+1') %>%
    select(-starts_with("Z")) %>%
    group_by_('Topic') %>%
    summarize_(count = 'n()', sumx='sum(X)') %>%
    ungroup() %>%
    mutate_(Proportion='round(sumx/count*100, 3)') -> p_estimated

  ## TODO fix this naming nonsense
  dict <- model$dict
  names(dict) <- names(model$keywords)

  ll <- list(seed_K = length(model$dict), extra_K = model$extra_k,
             V = V, N = N,
             model=model$mode,
             theta = theta, beta = tZW, # as.matrix(as.data.frame.matrix(tZW)),
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = doc_lens, vocab = model$vocab,
             dict = dict,
             alpha=res_alpha, modelfit=modelfit, p=p_estimated)
  class(ll) <- c("keyATM_posterior", class(ll))
  ll
}

# a more than usually informative error message for handing in the
# wrong type to a function
check_arg_type <- function(arg, typename){
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename))
    stop(paste("'", argname, '" is not a ', typename))
}



#' Set topic names
#'
#' @param x Posterior from a keyATM model (see \code{keyATM_posterior})
#' @param topic_names new names for topics
#'
#' @return a posterior object with new topic names in its components
#' @export
#'
set_topic_names <- function(x, topic_names){
  check_arg_type(x, "keyATM_posterior")
  colnames(x$theta) <- topic_names
  names(x$topic_counts) <- topic_names
  rownames(x$beta) <- topic_names
  x
}

#' Set document names
#'
#' @param x Posterior from a keyATM model (see \code{keyATM_posterior})
#' @param doc_names new names for documents
#'
#' @return a posterior object with new document names in its components
#' @export
#'
set_doc_names <- function(x, doc_names){
  check_arg_type(x, "keyATM_posterior")
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
#' @param n How many terms to show. Default: NULL, which shows all
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#' @param show_seed Mark seeded vocabulary. See below for details (Default: TRUE)
#'
#' @return An n x k table of the top n words in each topic
#' @export
#'
top_terms <- function(x, n = 10, measure = c("probability", "lift"),
                      show_seed = TRUE){
  check_arg_type(x, "keyATM_posterior")
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
#' @param n How many topics to show. Default: 2
#' @param measure How to sort the topics: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n topics in each document
#' @export
#'
top_topics <- function(x, n = 2, measure = c("probability", "lift")){
  check_arg_type(x, "keyATM_posterior")
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
#' @param n How many documents to show. Default: 10
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n documents for each topic
#' @export
top_docs <- function(x, n = 10, measure = c("probability", "lift")){
  check_arg_type(x, "keyATM_posterior")
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


#' Show a diagnosis plot of alpha
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param start Slice iteration
#' @param show_topic a vector to specify topic indexes to show
#' @param true_vec a vector to visualize true values of alpha
#' @param scale a parameter to control the scale of y-axis: 'free' adjusts y-axis for parameters
#'
#' @return ggplot2 object
#' @importFrom stats as.formula
#' @import ggplot2
#' @export
diagnosis_alpha <- function(x, start = NULL, show_topic = NULL, true_vec = NULL,
                            scale = ""){


  if("keyATM" %in% class(x)){
    num_topic <-  length(x$dict) + x$extra_k

    res_alpha <- data.frame(x$alpha_iter)
    colnames(res_alpha) <- NULL
    res_alpha <- data.frame(t(res_alpha))
    if(nrow(res_alpha) > 0){
      colnames(res_alpha) <- paste0("EstTopic", 1:ncol(res_alpha))
      res_alpha$iter <- 1:nrow(res_alpha)
    }
    
  }else if("keyATM_posterior" %in% class(x)){
    num_topic <-  x$seed_K + x$extra_k
    res_alpha <- x$alpha  
  }

  if(!is.null(show_topic)){
    # show topic is a vector of column index e.g., c(1,3,5)
    res_alpha <- res_alpha[, show_topic]
  }
  res_alpha$iter <- 1:nrow(res_alpha)

  if(!is.null(start)){
    res_alpha <- res_alpha[start:nrow(res_alpha), ]
  }


  parameters <- tidyr::gather(res_alpha, key = "parameter", value = "value",
                              -"iter")

  p <- ggplot(data=parameters, aes_string(x = 'iter', y = 'value',
                                          group = 'parameter', color = 'parameter')) +
     geom_line() +
     geom_point(size = 0.3)

  if(scale == ""){
    p <- p + facet_wrap(as.formula(paste("~", "parameter")), ncol = 2)
  } else if(scale == "free"){
    p <- p + facet_wrap(as.formula(paste("~", "parameter")), ncol = 2,
                        scales = "free")
  }

  if(!is.null(true_vec)){
    true <- data.frame(
       parameter = paste0("EstTopic", 1:length(true_vec)),
       value = true_vec
       )
    if (!is.null(show_topic)){
      true <- true[show_topic,]
    }

    p <- p + geom_hline(data = true, aes_string(yintercept = 'value'), color="black")
  }

  p <- p + ylab("Value") +
    ggtitle("Estimated Alpha") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}

#' Show a diagnosis plot of log-likelihood and perplexity
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param start Slice iteration
#'
#' @return ggplot2 object
#' @import ggplot2
#' @importFrom stats as.formula
#' @export
diagnosis_model_fit <- function(x, start=NULL){

  if("keyATM_posterior" %in% class(x)){
    modelfit <- x$modelfit
  }else if("keyATM" %in% class(x)){
    modelfit <- data.frame(x$model_fit)
    colnames(modelfit) <- NULL
    if(nrow(modelfit) > 0){
      modelfit <- data.frame(t(modelfit))
      colnames(modelfit) <-  c("Iteration", "Log Likelihood", "Perplexity")
    }  
  }

  if(!is.null(start)){
    modelfit <- modelfit[ modelfit$Iteration >= start, ]
  }

  modelfit <- tidyr::gather(modelfit, key="Measures", value="value", -"Iteration")

  p <- ggplot(data=modelfit, aes_string(x='Iteration', y='value',
                                        group='Measures', color='Measures')) +
     geom_line(show.legend = F) +
     geom_point(size=0.3, show.legend = F) +
     facet_wrap(as.formula(paste("~", "Measures")), ncol=2, scales = "free") +
     ylab("Value")

  p <- p + ggtitle("Model Fit") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


#' Show a diagnosis plot of p
#'
#' @param x The posterior from a fitted model (see \code{posterior})
#' @param topicvec A topic vector to reorder
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import dplyr
#' @export
diagnosis_p <- function(x, topicvec=c()){

  num <- length(unique(x$p$Topic))
  if(is.null(topicvec)){
    topicvec <- 1:num
  }else if(length(topicvec) != num){
    message("Topicvec length does not match with the topic number")
    topicvec <- 1:num
  }

  temp <- x$p
  temp$Topic <- paste0("EstTopic", temp$Topic)
  g  <- ggplot(temp, aes_string(x='Topic', y='Proportion')) +
      geom_bar(stat="identity") +
      theme_bw() +
      scale_x_discrete(limits = paste0("EstTopic", get("topicvec"))) +
      ylab("Proportion (%)") +
      xlab("Topic") +
      ggtitle("Proportion of words drawn from seed topic-word distribution") +
      theme(plot.title = element_text(hjust = 0.5))

  return(g)
}








