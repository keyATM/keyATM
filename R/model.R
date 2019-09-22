

#' Initialize a model (deprecated)
#'
#' @export
topicdict_model <- function(...){
  message("Warning: `topicdict_model` is deprecated, please use `keyATM_read` instead. keyATM does not take `quanteda` options.")
  return(keyATM_model(...))
}


#' Initialize a keyATM model
#'
#' This function reads a text input and initialize a keyATM model.
#' Internally, \code{keyATM_read()} calls \code{keyATM_model()} to initialize a model.
#'
#' @param texts Inputs. It can take quanteda dfm, data.frame, tibble, and a vector of characters.
#' @param keywords a quanteda dictionary or a list of character vectors
#' @param mode "basic", "cov", "tot", "totcov", and "ldaweight"
#' @param iteration number of iteration
#' @param extra_k number of regular topics in addition to the keyword topics by
#'                \code{keywords}
#' @param covariates_data covariate
#' @param covariates_formula formula applied to covariate data
#' @param timestamps timestamps
#' @param options options are seed, use_weights,
#'           visualize_keywords, and output_per.
#'
#' @return keyATM object, which is a list containing \describe{
#'         \item{W}{a list of vectors of word indexes}
#'         \item{Z}{a list of vectors of topic indicators isomorphic to W}
#'         \item{C}{a covariate matrix is there is an input}
#'         \item{X}{a list of vectors of seed indicators (0/1) isomorphic to W}
#'         \item{vocab}{a vector of vocabulary items}
#'         \item{files}{a vector of document filenames}
#'         \item{dict}{a tokenized version of the keyword dictionary}
#'         \item{keywords}{a list of keywords in dict, named by dictionary category}
#'         \item{extra_k}{how many extra non-seeded topics are required}
#'         \item{alpha}{a vector of topic proportion hyperparameters. If you use the model with covariates, it is not used.}
#'         \item{alpha_iter}{a list to store topic proportion hyperparameters}
#'         \item{Lambda_iter}{a list to store coefficients of the covariates}
#'         \item{S_iter}{a list to store states sampled in HMM}
#'         \item{tot_beta}{a list to store sampled beta parameters in topic-over-time}
#'         \item{sampling_info}{information related to sampling}
#'         \item{model_fit}{a list to store perplexity and log-likelihood}
#'         \item{gamma1}{First prior probability parameter for X (currently the same for all topics)}
#'         \item{gamma2}{Second prior probability parameter for X (currently the same for all topics)}
#'         \item{beta}{prior parameter for the non-seeded word generation probabilities}
#'         \item{beta_s}{prior parameter for the seeded word generation probabilities}
#'         \item{use_cov}{boolean, whether or not use the covariate}
#'         \item{num_states}{number of states in HMM}
#'         \item{timestamps}{time stamp for topics-over-time model}
#'         \item{visualize_keywords}{ggplot2 object}
#'         \item{call}{details of the function call}
#'         }.
#'
#' @export
keyATM_read <- function(texts, keywords, mode, extra_k,
                        iteration=1000,
                        covariates_data=NULL, covariates_formula= ~.+0,
                        timestamps=NULL,
                        options=list(
                                     seed=225,
                                     output_per=10,
                                     use_weights=TRUE,
                                     visualize_keywords=TRUE,
                                     x_prior=NULL
                                    )
                       )
{
  # Set option
  if(is.null(options$output_per))
    options$output_per <- 10

  if(is.null(options$visualize_keywords))
    options$visualize_keywords <- TRUE

  # Set random seed
  if(is.null(options$seed))
    options$seed <- 225
  set.seed(options$seed)

  # Detect input
  if("data.frame" %in% class(texts)){
    text_dfm <- NULL
    files <- NULL
    text_df <- texts
  }else if(class(texts) == "dfm"){
    text_dfm <- texts
    files <- NULL
    text_df <- NULL
  }else if(class(texts) == "character"){
    text_dfm <- NULL
    files <- texts
    text_df <- NULL
  }else{
    stop("Check `texts` argument.\n
         It can take quanteda dfm, data.frame, tibble, and a vector of characters.")  
  }

  # Reformat keywords
  if(class(keywords) != "list"){
      stop("`keywords` should be a list of character vectors")
  }

  # Initialize model
  message("Initializing keyATM...")
  model <- keyATM_model(
                          files=files, text_df=text_df, text_dfm=text_dfm,
                          extra_k=extra_k,
                          dict=keywords,
                          mode=mode,
                          covariates_data=covariates_data, covariates_formula=covariates_formula,
                          timestamps=timestamps,
                          options=options
                        )

  return(model)
}


#' Fit a keyATM model
#'
#' Fit keyATM model.
#'
#' @param model keyATM object (the output of \code{keyATM_read()})
#' @param iteration number of iteration
#'
#' @return keyATM object, which is a list containing \describe{
#'         \item{W}{a list of vectors of word indexes}
#'         \item{Z}{a list of vectors of topic indicators isomorphic to W}
#'         \item{C}{a covariate matrix is there is an input}
#'         \item{X}{a list of vectors of seed indicators (0/1) isomorphic to W}
#'         \item{vocab}{a vector of vocabulary items}
#'         \item{files}{a vector of document filenames}
#'         \item{dict}{a tokenized version of the keyword dictionary}
#'         \item{keywords}{a list of keywords in dict, named by dictionary category}
#'         \item{extra_k}{how many extra non-seeded topics are required}
#'         \item{alpha}{a vector of topic proportion hyperparameters. If you use the model with covariates, it is not used.}
#'         \item{alpha_iter}{a list to store topic proportion hyperparameters}
#'         \item{Lambda_iter}{a list to store coefficients of the covariates}
#'         \item{S_iter}{a list to store states sampled in HMM}
#'         \item{tot_beta}{a list to store sampled beta parameters in topic-over-time}
#'         \item{sampling_info}{information related to sampling}
#'         \item{model_fit}{a list to store perplexity and log-likelihood}
#'         \item{gamma1}{First prior probability parameter for X (currently the same for all topics)}
#'         \item{gamma2}{Second prior probability parameter for X (currently the same for all topics)}
#'         \item{beta}{prior parameter for the non-seeded word generation probabilities}
#'         \item{beta_s}{prior parameter for the seeded word generation probabilities}
#'         \item{use_cov}{boolean, whether or not use the covariate}
#'         \item{num_states}{number of states in HMM}
#'         \item{timestamps}{time stamp for topics-over-time model}
#'         \item{visualize_keywords}{ggplot2 object}
#'         \item{call}{details of the function call}
#'         }.
#'
#' @export
keyATM_fit <- function(model, iteration=1000){

  argname <- deparse(match.call()[['model']])
  if (!inherits(model, "keyATM"))
    stop(paste("'", argname, '" is not a keyATM object. Please use `keyATM_read()`'))


  # Prepare for store_theta
  if(model$options$store_theta){
    # We need matrices to store theta  
    model$options$Z_tables <- list()
  }  


  mode <- model$mode
  set.seed(model$options$seed)

  message("Start fitting keyATM...")  
  if(mode == "basic"){
    res <- keyATM_train(model, iter=iteration, output_per=model$options$output_per)
  }else if(mode == "cov"){
    res <- keyATM_train_cov(model, iter=iteration, output_per=model$options$output_per)
  }else if(mode == "tot"){
    res <- keyATM_train_tot(model, iter=iteration, output_per=model$options$output_per)
  }else if(mode == "totcov"){
    res <- keyATM_train_totcov(model, iter=iteration, output_per=model$options$output_per)
  }else if(mode == "ldaweight"){
    res <- LDA_weight(model, iter=iteration, output_per=model$options$output_per)  
  }else{
    stop("Please check `mode`.")  
  }

  return(res)

}


#' Fit topic model (deprecated)
#'
#' `topicdict_train()` is deprecated. Use `keyATM_fit()`
#'
#' @export
topicdict_train <- function(...){
  message("Warning: `topicdict_train` is deprecated, please use `keyATM_fit` instead.")
  return(keyATM_train(...))
}

#' Fit topic model (deprecated)
#'
#' `topicdict_train_cov()` is deprecated. Use `keyATM_fit()`
#'
#' @export
topicdict_train_cov <- function(...){
  message("Warning: `topicdict_train_cov` is deprecated, please use `keyATM_fit` instead.")
  return(keyATM_train_cov(...))
}

#' Fit topic model (deprecated)
#'
#' `topicdict_train_tot()` is deprecated. Use `keyATM_fit()`
#'
#' @export
topicdict_train_tot <- function(...){
  message("Warning: `topicdict_train_tot` is deprecated, please use `keyATM_fit` instead.")
  return(keyATM_train_tot(...))
}

#' Fit topic model (deprecated)
#'
#' `topicdict_train_totcov()` is deprecated. Use `keyATM_fit()`
#'
#' @export
topicdict_train_totcov <- function(...){
  message("Warning: `topicdict_train_totcov` is deprecated, please use `keyATM_fit` instead.")
  return(keyATM_train_totcov(...))
}



keyATM_model <- function(files=NULL, dict=NULL, text_df=NULL, text_dfm=NULL,
                         mode="",
                         extra_k = 1,
                         covariates_data=NULL, covariates_formula=NULL,
                         num_states=NULL, timestamps=NULL,
                         alpha = 50/(length(dict) + extra_k),
                         beta = 0.01, beta_s = 0.1,
                         options = list()
                        )
{
  cl <- match.call()

  ##
  ## Check format
  ##
  if(mode %in% c("basic", "cov", "hmm", "tot", "totcov", "ldaweight")){
  }else{
    stop(paste0("Unknown model:", mode))  
  }


  if(!is.null(covariates_data) & !is.null(covariates_formula) & (mode != "cov" & mode != "totcov")){
    stop("Covariates information provided, specify the model.")  
  }

  if(is.null(num_states) & mode == "hmm"){
    stop("Provide the number of states.")  
  }

  if(!is.null(text_df)){
    if("text" %in% names(text_df)){
      text_df <- text_df["text"]
    }else{
      stop("text_df should have a 'text' colum that has documents.")
    }  
  }

  if(mode == "tot" | mode == "totcov"){
    if(is.null(timestamps)){
      stop("Please provide time stamps.")  
    }else if(min(timestamps) < 0 | max(timestamps) >= 1){
      stop("Time stamps sholud be between 0 and 1. Please use `make_timestamps()` function to format time stamps.")  
    }
  }


  ##
  ## Check length
  ##

  proper_len <- length(dict) + extra_k

  if(mode == "cov" | mode == "totcov"){
    # make sure covariates are provided for all documents  
    doc_num <- ifelse(!is.null(files), length(files),
                      ifelse(!is.null(text_df), nrow(text_df),
                             ifelse(!is.null(text_dfm), nrow(text_dfm),
                                    0)
                             )
                      )
    if( nrow(covariates_data) != doc_num ){
      stop("Covariates dimension does not match with the number of documents.")  
    }
  }


  ##
  ## Set options
  ##
  if(is.null(options$use_weights)){
    options$use_weights <- 1  
  }else{
    options$use_weights <- as.numeric(options$use_weights)
  }
  if(!is.null(options$logsumexp_approx)){
    message("Warning: `logsumexp_approx` option will be deprecated.")
  }
  if(is.null(options$logsumexp_approx)){
    # 0: do not approximate logsumexp
    options$logsumexp_approx <- 0
  }
  if(is.null(options$slice_shape)){
    # parameter for slice sampling
    options$slice_shape <- 1.2
  }
  if(is.null(options$use_mom)){
    # Method of Moments in TOT
    options$use_mom <- 0
  }
  if(!is.null(options$alpha)){
    # If alpha value is overwritten
    alpha = options$alpha
  }
  if(is.null(options$store_theta)){
    options$store_theta <- 0
  }else{
    options$store_theta <- as.numeric(options$store_theta)  
  }
  if(is.null(options$x_prior)){
    options$x_prior <- matrix(1.0, nrow=proper_len, ncol=2)  
  }
  if(!is.null(options$x_prior)){
    if(dim(options$x_prior)[1] != proper_len)  
      stop("Check the dimension of `options$x_prior`")
    if(dim(options$x_prior)[2] != 2)  
      stop("Check the dimension of `options$x_prior`")
  }


  ## Create alpha 
  if(is.null(covariates_data) || is.null(covariates_formula)){
    # If it doesn't use covariates, make alpha inside
    if (length(alpha) == 1){
      message("All ", proper_len, " values for alpha starting as ", alpha)
      alpha = rep(alpha, proper_len)
    } else if (length(alpha) != proper_len)
      stop("Starting alpha must be a scalar or a vector of length ", proper_len)
  }


  ##
  ## Text preprocessing
  ##

  # If you have quanteda object
  if(!is.null(text_dfm)){
    message("Use quanteda dfm.")  

    vocabulary <- colnames(text_dfm)
    texts <- apply(text_dfm, 1,
                   function(x){
                       single_text <- paste(rep(vocabulary, x), collapse=" ")
                       return(single_text)
                   })

    text_df <- data.frame(text = texts)
    text_df$text <- as.character(text_df$text)
  }else{
    ## preprocess each text
    # Use files <- list.files(doc_folder, pattern="txt", full.names=T) when you pass
    if(is.null(text_df)){
      text_df <- data.frame(text = unlist(lapply(files, 
                                                 function(x)
                                                 { 
                                                     paste0(readLines(x, encoding = encoding),
                                                           collapse = "\n") 
                                                 })),
                            stringsAsFactors = FALSE)
    }
    
    text_df$doc_id <- paste0("text", 1:nrow(text_df))
  }

  # args$x <- corpus(readtext(file_pattern, encoding = encoding))
  # for new version quanteda, you need this
  args$x$documents$doc_id <- paste0("text", 1:quanteda::ndoc(args$x))
  doc_names <- quanteda::docvars(args$x, "doc_id") # docnames
  toks <- do.call(quanteda::tokens, args = args)
  # if (lowercase)
    # toks <- tokens_tolower(toks)
  if (!is.null(stopwords))
    toks <- tokens_remove(toks, stopwords)
  if (!is.null(stem_language))
    toks <- tokens_wordstem(toks, language = stem_language)

  ## apply the same preprocessing to the seed words
  args$x <- do.call(rbind, lapply(as.list(dict), paste0, collapse = " "))
  dtoks <- do.call(quanteda::tokens, args = args)
  # if (lowercase)
    # dtoks <- tokens_tolower(dtoks)
  if (!is.null(stopwords))
    dtoks <- tokens_remove(dtoks, stopwords)
  if (!is.null(stem_language))
    dtoks <- tokens_wordstem(dtoks, language = stem_language)
  K <- length(dtoks) # number of seeded categories a.k.a. size of dictionary


  ##
  ## Visualize keywords
  ##
  if(options$visualize_keywords){
    dfm_ <- quanteda::dfm(toks, tolower=F)  
    data <- tidytext::tidy(dfm_)
    totalwords <- sum(data$count)

    data %>%
      rename(Word=term) %>%
      group_by(Word) %>%
      summarize(WordCount = sum(count)) %>%
      ungroup() %>%
      mutate(`Proportion(%)` = round(WordCount/totalwords*100, 3)) %>%
      arrange(desc(WordCount)) %>%
      mutate(Ranking = 1:n()) -> data


    seed_list <- dict

    names(seed_list) <- paste0("Topic", 1:length(seed_list))
    seeds <- lapply(seed_list, function(x){unlist(strsplit(x," "))})
    ext_k <- length(seeds)
    max_num_words <- max(unlist(lapply(seeds, function(x){length(x)})))

    seeds_df <- data.frame(Topic=1, Word=1)
    for(k in 1:ext_k){
      words <- seeds[[k]]
      numwords <- length(words)
      topicname <- paste0("Topic", k)
      for(w in 1:numwords){
        seeds_df <- rbind(seeds_df, data.frame(Topic=topicname, Word=words[w]))
      }
    }
    seeds_df <- seeds_df[2:nrow(seeds_df), ]

    dplyr::inner_join(data, seeds_df, by="Word") %>%
      group_by(Topic) %>%
      mutate(Ranking = 1:n()) -> temp

    visualize_keywords <- 
      ggplot(temp, aes(x=Ranking, y=`Proportion(%)`, colour=Topic)) +
        geom_line() +
        geom_point() +
        geom_label_repel(aes(label = Word), size=2.8,
                         box.padding = 0.20, label.padding = 0.12,
                         arrow=arrow(angle=10, length = unit(0.10, "inches"), ends = "last", type = "closed"),
                         show.legend = F) +
        scale_x_continuous(breaks=1:max_num_words) +
        ylab("Proportion (%)") +
        theme_bw()
  
  }else{
    visualize_keywords <- NULL  
  }


  ##
  ## Initialization
  ##

  ## construct W and a vocab list (W elements are 0 based ids)
  wd_names <- attr(toks, "types") # vocab
  wd_map <- hashmap::hashmap(wd_names, as.integer(1:length(wd_names) - 1))
  W <- lapply(toks, function(x){ wd_map[[x]] })

  # zx_assigner maps seed words to category ids
  seed_wdids <- unlist(lapply(dtoks, function(x){ wd_map$find(x) }))
  cat_ids <- rep(1:K - 1, unlist(lapply(dtoks, length)))
  zx_assigner <- hashmap(as.integer(seed_wdids), as.integer(cat_ids))

  ## xx indicates whether the word comes from a seed topic-word distribution or not
  make_x <- function(x){
    seeded <- as.numeric(zx_assigner$has_keys(x)) # 1 if they're a seed
    # Use x structure
    x[seeded == 0] <- 0 # non-seeded words have x=0
    x[seeded == 1] <- sample(0:1, length(x[seeded == 1]), prob = c(0.3, 0.7), replace = TRUE)
      # seeded words have x=1 probabilistically
    x
  }

  X <- lapply(W, make_x)

  # if the word is a seed, assign the appropriate (0 start) Z, else a random Z
  make_z <- function(x){
    zz <- zx_assigner[[x]] # if it is a seed word, we already know the topic
    zz[is.na(zz)] <- sample(1:(K + extra_k) - 1,
                            sum(as.numeric(is.na(zz))),
                            replace = TRUE)
    zz
  }
  Z <- lapply(W, make_z)

  # dictionary category names -> vector of word_id.
  # (Later processes ignore names)
  keywords <- lapply(dtoks, function(x){ wd_map$find(x) })
  names(keywords) <- names(dict)

  # Covariate
  if(is.null(covariates_data) || is.null(covariates_formula)){
    C <- matrix()
    use_cov <- FALSE
  }else{
    C <- model.matrix(covariates_formula, covariates_data)
    use_cov <- TRUE
    if(sum(is.na(C)) != 0){
      stop("Covariate data should not contain missing values.")
    }
  }


  ll <- list(W = W, Z = Z, X = X, vocab = wd_names, mode=mode,
             files = doc_names, dict = dtoks, keywords = keywords, extra_k = extra_k,
             alpha = alpha,
             beta = beta, beta_s = beta_s,
             alpha_iter = list(), Lambda_iter = list(), S_iter = list(),
             model_fit = list(), sampling_info = list(), tot_beta = list(),
             C=C, use_cov=use_cov,
             num_states=num_states,
             timestamps=timestamps,
             options=options,
             visualize_keywords=visualize_keywords,
             call = cl)

  class(ll) <- c("keyATM", class(ll))

  return(ll)
}

