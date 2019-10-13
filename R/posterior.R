#' Get posterior quantities from model output
#'
#' Constructs a (N x K) matrix \code{theta} and (K x V) matrix \code{beta}
#' plus their margins from the sample of Z and W in \code{model}.
#' \code{model} is an output of \code{keyATM_run()}.
#' These statistics implicitly marginalize over X.
#'
#' @param model a fitted keyATM model
#'
#' @return A list containing:
#'   \describe{
#'     \item{seed_K}{ Number of keyword topics}
#'     \item{regular_k}{ Number of regular unseeded topics}
#'     \item{V}{ Number of word types}
#'     \item{N}{ Number of documents}
#'     \item{theta}{ Normalized tpoic proportions for each document}
#'     \item{phi}{ Normalized topic specific word generation probabilities}
#'     \item{topic_counts}{ Number of tokens assigned to each topic}
#'     \item{word_counts}{ Number of times each word type appears}
#'     \item{doc_lens}{ Length of each document in tokens}
#'     \item{vocab}{ Words in the vocabulary}
#'     \item{alpha}{ \code{alpha} during the iteration}
#'     \item{modelfit}{ Perplexity and log-likelihood}
#'     \item{p}{ Estimated p}
#'     \item{options}{ Options used in the \code{model}}
#'   }
#' @export
keyATM_output <- function(model){
  message("Creating an output object. It may take time...")

  check_arg_type(model, "keyATM_fitted")
  values_iter <- list()  # store values by iteration

  # Make info
  info <- list()
  info$allK <- model$regular_k + length(model$keywords)
  info$V <- length(model$vocab)
  info$N <- length(model$W)
  info$doc_lens <- sapply(model$W, length)

  if(model$regular_k > 0 & length(model$keywords) != 0){
    info$tnames <- c(paste0("", 1:length(model$keywords)), paste0("T_", 1:model$regular_k))
  }else if(model$regular_k > 0 & length(model$keywords) == 0) {
    # No keywords (= lda models)
    info$tnames <- paste0("T_", 1:model$regular_k)
  }else{
    # Keywords only
    info$tnames <- c(paste0("", 1:length(model$keywords)))
  }


  # theta (document-topic distribution)
  theta <- keyATM_output_theta(model, info)

  # theta iter
  if(model$options$store_theta){
    values_iter$theta_iter <- keyATM_output_theta_iter(model, info)  
  }

  # Phi (topic-word distribution)
  phi <- keyATM_output_phi(model, info)
  

  # alpha
  if(model$model %in% c("hmm", "ldahmm")){
    res_alpha <- list()

    for (i in 1:model$model_settings$num_states) {
      res_alpha[[i]] <- do.call(rbind, lapply(model$stored_values$alpha_iter, function(x){ x[i, ] }))
    }
    res_alpha <- lapply(res_alpha, 
      function(x){colnames(x) <- paste0("EstTopic", 1:ncol(x)); y <- data.frame(x); y$iter <- 1:nrow(x); return(y)})

  }else{
        res_alpha <- data.frame(model$stored_values$alpha_iter)
        colnames(res_alpha) <- NULL
        res_alpha <- data.frame(t(res_alpha))
        if(nrow(res_alpha) > 0){
          colnames(res_alpha) <- paste0("EstTopic", 1:ncol(res_alpha))
          res_alpha$iter <- 1:nrow(res_alpha)
        }
  }

  # model fit
	modelfit <- NULL
	if(length(modelfit) > 0){
		model$model_fit %>%
			purrr::set_names(1:length(.)) %>%
			dplyr::bind_rows() %>%
			t() %>%
			as_tibble(., .name_repair = ~c("Iteration", "Log Likelihood", "Perplexity")) -> modelfit
	}

  # p
  data <- tibble::tibble(Z=unlist(model$Z, use.names=F),
                         X=unlist(model$X, use.names=F))
  data %>%
    dplyr::mutate(Topic=Z+1L) %>%
    dplyr::select(-starts_with("Z")) %>%
    dplyr::group_by(Topic) %>%
    dplyr::summarize(count = n(), sumx=sum(X)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Proportion=round(sumx/count*100, 3)) -> p_estimated

  # Make an object to return
  ll <- list(keyword_K = length(model$keywords), regular_k = model$regular_k,
             V = V, N = N,
             model=model$model,
             theta = theta, phi = phi,
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = doc_lens, vocab = model$vocab,
             keywords_raw = model$keywords_raw,
             alpha=res_alpha, modelfit=modelfit, p=p_estimated,
             values_iter=values_iter)
  class(ll) <- c("keyATM_output", class(ll))
  ll
}


keyATM_output_theta <- function(model, info)
{
  # Theta
  if(model$model %in% c("cov")){
    Alpha <- exp(model$model_settings$covariates_data %*% t(model$stored_values$Lambda_iter[[length(model$stored_values$Lambda_iter)]]))

    posterior_z <- function(docid){
      zvec <- model$Z[[docid]]
      alpha <- Alpha[docid, ]
      tt <- table(factor(zvec, levels = 1:(info$allK) - 1L))
      (tt + alpha) / (sum(tt) + sum(alpha)) # posterior mean
    }

    theta <- do.call(dplyr::bind_rows, lapply(1:length(model$Z), posterior_z))

  }else if(model$model %in% c("basic", "lda")){
    alpha <- model$stored_values$alpha_iter[[length(model$stored_values$alpha_iter)]]  

    posterior_z <- function(zvec){
      tt <- table(factor(zvec, levels = 1:(info$allK) - 1L))
      (tt + alpha) / (sum(tt) + sum(alpha)) # posterior mean
    }  

    theta <- do.call(dplyr::bind_rows, lapply(model$Z, posterior_z))

  }else if(model$model %in% c("hmm", "ldahmm")){
    S <- model$stored_values$S_iter[[length(model$stored_values$S_iter)]] + 1  # adjust index for R
    S <- S[model$model_settings$time_index]  # retrieve doc level state info
    alphas <- matrix(model$stored_values$alpha_iter[[length(model$stored_values$alpha_iter)]][S],
                     nrow=length(model$W), ncol=info$allK)

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


keyATM_output_phi <- function(model, info)
{
  all_words <- model$vocab[as.integer(unlist(model$W)) + 1L]
  all_topics <- as.integer(unlist(model$Z))
  
  res_tibble <- data.frame(
                        Word = all_words,
                        Topic = all_topics
                       ) %>%
                dplyr::group_by(Topic, Word) %>%
                dplyr::summarize(Count = dplyr::n())
  
  res_tibble %>%
    tidyr::spread(key=Word, value=Count)  -> beta
  beta <- apply(beta, 2, function(x){ifelse(is.na(x), 0, x)})
  beta <- beta[, 2:ncol(beta)] + model$priors$beta
  beta <- beta[, model$vocab]

  topic_counts <- Matrix::rowSums(beta)
  word_counts <- Matrix::colSums(beta)

  phi <- beta / topic_counts
  rownames(phi) <- tnames

  return(phi)
}


keyATM_output_theta_iter <- function(model, info)
{
  if(model$model %in% c("cov")){
    posterior_theta <- function(x){
      Z_table <- model$stored_values$Z_tables[[x]]
      lambda <- model$stored_values$Lambda_iter[[x]]
      Alpha <- exp(model$model_settings$covariates_data %*% t(lambda))

      tt <- Z_table + Alpha
      row.names(tt) <- NULL

      return(tt / Matrix::rowSums(tt))
    }
  }else if(model$model %in% c("hmm", "ldahmm")){
    posterior_theta <- function(x){
      Z_table <- model$stored_values$Z_tables[[x]]
      S <- model$stored_values$S_iter[[x]] + 1L  # adjust index for R
      S <- S[model$model_settings$time_index]  # retrieve doc level state info

      alphas <- matrix(model$stored_values$alpha_iter[[x]][S],
                       nrow=length(model$W), ncol=info$allK)
    
      tt <- Z_table + alphas
      theta <- tt / Matrix::rowSums(tt)
      return(theta)
    }
  }else{
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
#' @export
print.keyATM_output <- function(x){
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
summary.keyATM_output <- function(x){
  cat(
      paste0(
             "keyATM_output object for the ",
             x$model,
             " model. ",
             "\n"
      )
     )
}



# a more than usually informative error message for handing in the
# wrong type to a function
check_arg_type <- function(arg, typename, message=NULL){
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename)){
    if(is.null(message))
      stop(paste('"', argname, '" is not a ', typename))
    else
      stop(message)
  }
}



#' Set topic names
#'
#' @param x the output from a keyATM model (see \code{keyATM_output})
#' @param topic_names new names for topics
#'
#' @return an output object with new topic names in its components
#' @export
#'
set_topic_names <- function(x, topic_names){
  check_arg_type(x, "keyATM_output")
  colnames(x$theta) <- topic_names
  names(x$topic_counts) <- topic_names
  rownames(x$beta) <- topic_names
  x
}

#' Set document names
#'
#' @param x the output from a keyATM model (see \code{keyATM_output})
#' @param doc_names new names for documents
#'
#' @return an output object with new document names in its components
#' @export
#'
set_doc_names <- function(x, doc_names){
  check_arg_type(x, "keyATM_output")
  rownames(x$theta) <- doc_names
  names(x$doc_lens) <- doc_names
  x
}



#' Show the top words for each topic
#'
#' If \code{show_keyword} is true then words in their seeded categories
#' are suffixed with a check mark. Words from another seeded category
#' are labeled with the name of that category.
#'
#' @param x the output from a keyATM model (see \code{keyATM_output})
#' @param n How many terms to show. Default: NULL, which shows all
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#' @param show_keyword Mark keywords. (default: TRUE)
#'
#' @return An n x k table of the top n words in each topic
#' @export
#'
top_words <- function(x, n = 10, measure = c("probability", "lift"),
                      show_keyword = TRUE){
  check_arg_type(x, "keyATM_output")

  if(x$model %in% c("lda", "ldacov", "ldahmm"))
     show_keyword <- FALSE

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
  if (show_keyword) {
    for (i in 1:ncol(res)) {
      for (j in 1:length(x$keywords_raw)) {
         inds <- which(res[,i] %in% x$keywords_raw[[j]])
         label <- ifelse(i == j,
                         paste0("[", "\U2713" ,"]"),
                         paste0("[", names(x$keywords_raw)[j], "]"))
         res[inds, i] <- paste(res[inds, i], label)
      }
    }
  }
  res
}

#' Show the top topics for each document
#'
#' @param x the output from a keyATM model (see \code{keyATM_output})
#' @param n How many topics to show. Default: 2
#' @param measure How to sort the topics: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n topics in each document
#' @export
#'
top_topics <- function(x, n = 2, measure = c("probability", "lift")){
  check_arg_type(x, "keyATM_output")
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
#' @param x the output from a keyATM model (see \code{keyATM_output})
#' @param n How many documents to show. Default: 10
#' @param measure How to sort the terms: 'probability' (default) or 'lift'
#'
#' @return An n x k table of the top n documents for each topic
#' @export
top_docs <- function(x, n = 10, measure = c("probability", "lift")){
  check_arg_type(x, "keyATM_output")
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
#' @param x the output from a keyATM model (see \code{keyATM_output})
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

  check_arg_type(x, "keyATM_output")
  if("keyATM" %in% class(x)){
    num_topic <-  length(x$keywords_raw) + x$regular_k

    res_alpha <- data.frame(x$alpha_iter)
    colnames(res_alpha) <- NULL
    res_alpha <- data.frame(t(res_alpha))
    if(nrow(res_alpha) > 0){
      colnames(res_alpha) <- paste0("EstTopic", 1:ncol(res_alpha))
      res_alpha$iter <- 1:nrow(res_alpha)
    }
    
  }else if("keyATM_output" %in% class(x)){
    num_topic <-  x$seed_K + x$regular_k
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
#' @param x the output from a keyATM model (see \code{keyATM_output})
#' @param start Slice iteration
#'
#' @return ggplot2 object
#' @import ggplot2
#' @importFrom stats as.formula
#' @export
diagnosis_model_fit <- function(x, start=NULL){

  if("keyATM_output" %in% class(x)){
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
#' @param x the output from a keyATM model (see \code{keyATM_output})
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








