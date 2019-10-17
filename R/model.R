#' Read Texts
#'
#' This function read texts and creates a tibble. 
#'
#' @param texts Input. keyATM takes dfm, data.frame, tibble, and a vector of file paths.
#'
#' @return A list whose elements are splitted texts.
#'
#' @examples
#' \dontrun{
#'  # Use quanteda dfm 
#'  keyATM_docs <- keyATM_read(quanteda_dfm) 
#'   
#'  # Use data.frame or tibble (texts should be stored in a column named `text`)
#'  keyATM_docs <- keyATM_read(data_frame_object) 
#'  keyATM_docs <- keyATM_read(tibble_object) 
#' 
#'  # Use a vector that stores full paths to the files  
#'  files <- list.files(doc_folder, pattern = "*.txt", full.names = T) 
#'  keyATM_docs <- keyATM_read(files) 
#' }
#'
#' @import magrittr
#' @export
keyATM_read <- function(texts, encoding = "UTF-8")
{

  # Detect input
  if("tbl" %in% class(texts)){
    message("Using tibble.")
    text_dfm <- NULL
    files <- NULL
    text_df <- texts
  }else if("data.frame" %in% class(texts)){
    message("Using data.frame.")
    text_dfm <- NULL
    files <- NULL
    text_df <- tibble::as_tibble(texts)
  }else if(class(texts) == "dfm"){
    message("Using quanteda dfm.")
    text_dfm <- texts
    files <- NULL
    text_df <- NULL
  }else if(class(texts) == "character"){
    message("Reading from files.")
    text_dfm <- NULL
    files <- texts
    text_df <- NULL
    message(paste0("Encoding: ", encoding))
  }else{
    stop("Check `texts` argument.\n
         It can take quanteda dfm, data.frame, tibble, and a vector of characters.")  
  }


  # Read texts

  # If you have quanteda object
  if(!is.null(text_dfm)){
    vocabulary <- colnames(text_dfm)
    text_df <- tibble::tibble(
                              text_split = 
                                apply(text_dfm, 1,
                                       function(x){
                                        return(rep(vocabulary, x))
                                       }
                                    )
                             )
    names(text_df$text_split) <-NULL
  }else{
    ## preprocess each text
    # Use files <- list.files(doc_folder, pattern = "txt", full.names = T) when you pass
    if(is.null(text_df)){
      text_df <- tibble::tibble(text = unlist(lapply(files,
                                                 function(x)
                                                 { 
                                                     paste0(readLines(x, encoding = encoding),
                                                           collapse = "\n") 
                                                 })))
    }
    text_df <- text_df %>% dplyr::mutate(text_split = stringr::str_split(text, pattern = " "))
  }

  W_raw <- text_df %>% dplyr::pull(text_split)
  class(W_raw) <- c("keyATM_docs", class(W_raw))

  return(W_raw)
}


#' @noRd
#' @export
print.keyATM_docs <- function(x)
{
  cat(paste0("keyATM_docs object of ",
                 length(x), " documents",
                 ".\n"
                )
      )
}


#' @noRd
#' @export
summary.keyATM_docs <- function(x)
{
  doc_len <- sapply(x, length)
  cat(paste0("keyATM_docs object of: ",
              length(x), " documents",
              ".\n",
              "Length of documents:",
              "\n  Avg: ", round(mean(doc_len),3),
              "\n  Min: ", round(min(doc_len),3),
              "\n  Max: ", round(max(doc_len),3),
              "\n   SD: ", round(sd(doc_len),3),
              "\n"
             )  
         )
}



#' Visualize texts
#'
#' This function visualizes the proportion of keywords in the documents.
#'
#' @param keyATM_docs A list of texts read via \code{keyATM_read()} function
#' @param keywords A list of keywords
#' @param label_size The size of the keyword labels
#'
#' @return A list containing \describe{
#'    \item{figure}{a ggplot2 object}
#'    \item{values}{a tibble object that stores values}
#' }
#'
#' @examples
#' \dontrun{
#'  # Prepare a keyATM_docs object
#'  keyATM_docs <- keyATM_read(input) 
#'   
#'  # Keywords are in a list  
#'  keywords <- list(
#'                    c("education", "child", "student"),  # Education
#'                    c("public", "health", "program"),  # Health
#'                  )
#'
#'  # Visualize keywords
#'  keyATM_viz <- visualize_keywords(keyATM_docs, keywords)
#'
#'  # View a figure
#'  keyATM_viz
#'    # Or: `keyATM_viz$figure`
#' 
#'  # Save a figure 
#'  ggsave(filename, keyATM_viz$figure)
#'
#' }
#'
#' @import magrittr
#' @import ggplot2
#' @export
visualize_keywords <- function(keyATM_docs, keywords, label_size = 3.2)
{
  # Check type
  check_arg_type(keyATM_docs, "keyATM_docs", "Please use `keyATM_read()` to read texts.")
  check_arg_type(keywords, "list")
  c <- lapply(keywords, function(x){check_arg_type(x, "character")})

  unnested_data <- tibble::tibble(text_split = unlist(keyATM_docs,
                                                      recursive = FALSE, use.names = FALSE))

  # Organize data
  totalwords <- nrow(unnested_data)

  unnested_data %>%
    dplyr::rename(Word = text_split) %>%
    dplyr::group_by(Word) %>%
    dplyr::summarize(WordCount = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`Proportion(%)` = round(WordCount / totalwords * 100, 3)) %>%
    dplyr::arrange(desc(WordCount)) %>%
    dplyr::mutate(Ranking = 1:(dplyr::n())) -> data

  names(keywords) <- paste0("Topic", 1:length(keywords))
  keywords <- lapply(keywords, function(x){unlist(strsplit(x," "))})
  ext_k <- length(keywords)
  max_num_words <- max(unlist(lapply(keywords, function(x){length(x)}), use.names = F))

  # Make keywords_df
  keywords_df <- data.frame(Topic = 1, Word = 1)
  for(k in 1:ext_k){
    words <- keywords[[k]]
    numwords <- length(words)
    topicname <- paste0("Topic", k)
    for(w in 1:numwords){
      keywords_df <- rbind(keywords_df, data.frame(Topic = topicname, Word = words[w]))
    }
  }
  keywords_df <- keywords_df[2:nrow(keywords_df), ]

  dplyr::right_join(data, keywords_df, by = "Word") %>%
    dplyr::group_by(Topic) %>%
    dplyr::arrange(desc(`Proportion(%)`)) %>%
    dplyr::mutate(Ranking  =  1:(dplyr::n())) %>%
    dplyr::arrange(Topic, Ranking) -> temp


  # Check keywords existence
  temp %>%
    dplyr::filter(is.na(WordCount)) %>%
    dplyr::pull(Word) -> non_appearence

  if(length(non_appearence) != 0){
    if(length(non_appearence) == 1){
      stop("A keyword not found in texts: ", paste(non_appearence, collapse = ", "))
    }else{
      stop("Keywords not found in texts: ", paste(non_appearence, collapse = ", "))
    }
  }


  visualize_keywords <- 
    ggplot(temp, aes(x = Ranking, y=`Proportion(%)`, colour = Topic)) +
      geom_line() +
      geom_point() +
      ggrepel::geom_label_repel(aes(label = Word), size = label_size,
                       box.padding = 0.20, label.padding = 0.12,
                       arrow = arrow(angle = 10, length = unit(0.10, "inches"),
                                   ends = "last", type = "closed"),
                       show.legend = F) +
      scale_x_continuous(breaks = 1:max_num_words) +
      ylab("Proportion (%)") +
      theme_bw()

  keyATM_viz <- list(figure = visualize_keywords, values = temp)
  class(keyATM_viz) <- c("keyATM_viz", class(keyATM_viz))
  
  return(keyATM_viz)

}


#' @noRd
#' @export
print.keyATM_viz <- function(x)
{
  print(x$figure)  
}


#' @noRd
#' @export
summary.keyATM_viz <- function(x)
{
  return(x$values)  
}


#' @noRd
#' @export
save.keyATM_viz <- function(x, file = stop("'file' must be specified"))
{
  save(x, file = file, compress="xz", compression_level = 3)
}



#' Fit keyATM model
#'
#' Select and specify one of the keyATM models and fit the model.
#'
#'
#' @param keyATM_docs texts read via \code{keyATM_read()}
#' @param model keyATM model: "basic", "cov", "hmm", "lda", "ldacov" and "ldahmm"
#' @param regular_k the number of regular topics
#' @param keywords a list of keywords
#' @param model_settings a list of model specific settings
#' @param priors a list of priors of parameters
#' @param options a list of options
#'
#' @return keyATM_model object, which is a list containing \describe{
#'   \item{W}{a list of vectors of word indexes}
#'   \item{Z}{a list of vectors of topic indicators isomorphic to W}
#'   \item{X}{a list}
#'   \item{model}{the name of the model}
#'   \item{keywords}{}
#'   \item{keywords_raw}{}
#'   \item{regular_k}{the number of regular topics}
#'   \item{model_settings}{a list of settings}
#'   \item{priors}{a list of priors}
#'   \item{options}{a list of options}
#'   \item{stored_values}{a list of stored_values}
#'   \item{call}{details of the function call}
#' } 
#'
#' @examples
#' \dontrun{
#'   # keyATM Basic
#'   fitted <- keyATM_fit(
#'                        keyATM_docs, model = "basic", regular_k = 5, keywords = keywords_list
#'                       )
#'
#'   # keyATM Cov
#'   fitted <- keyATM_fit(
#'                        keyATM_docs, model = "cov", regular_k = 5, keywords = keywords_list,
#'                        model_settings(covariates_data = cov)
#'                       )
#'
#'   # keyATM HMM
#'   fitted <- keyATM_fit(
#'                        keyATM_docs, model = "hmm", regular_k = 5, keywords = keywords_list,
#'                        model_settings(time_index = time_index_vec, num_states = 5)
#'                       )
#'
#'   # Weighted LDA
#'   fitted <- keyATM_fit(
#'                        keyATM_docs, model = "lda", regular_k = 5
#'                       )
#'
#'   # Weighted LDA Cov
#'   fitted <- keyATM_fit(
#'                        keyATM_docs, model = "ldacov", regular_k = 5,
#'                        model_settings(covariates_data = cov)
#'                       )                   
#'
#'   # Weighted LDA HMM
#'   fitted <- keyATM_fit(
#'                        keyATM_docs, model = "ldahmm", regular_k = 5,
#'                        model_settings(time_index = time_index_vec, num_states = 5)
#'                       )
#'
#' }
#'
#' @export
keyATM_fit <- function(keyATM_docs, model, regular_k,
                       keywords = list(), model_settings = list(),
                       priors = list(), options = list()) 
{
  ##
  ## Check
  ##

  # Check type
  check_arg_type(keyATM_docs, "keyATM_docs", "Please use `keyATM_read()` to read texts.")
  if(!is.integer(regular_k) & !is.numeric(regular_k))
    stop("`regular_k` is neigher numeric nor integer.")

  regular_k <- as.integer(regular_k)

  if(!model %in% c("basic", "cov", "hmm", "lda", "ldacov", "ldahmm")){
    stop("Please select a correct model.")  
  }

  info <- list(
                models_keyATM = c("basic", "cov", "hmm"),
                models_lda = c("lda", "ldacov", "ldahmm")
              )
  keywords <- check_arg(keywords, "keywords", model, info)

  # Get Info
  info$num_doc <- length(keyATM_docs)
  info$keyword_k <- length(keywords)
  info$total_k <- length(keywords) + regular_k

  # Set default values
  model_settings <- check_arg(model_settings, "model_settings", model, info)
  priors <- check_arg(priors, "priors", model, info)
  options <- check_arg(options, "options", model, info)

  ##
  ## Initialization
  ##
  message("Initializing the model...")
  set.seed(options$seed)

  # W
  info$wd_names <- unique(unlist(keyATM_docs, use.names = F))
  if(" " %in% info$wd_names){
    stop("A space is recognized as a vocabulary.
          Please remove an empty document or consider using quanteda::dfm.")  
  }

  info$wd_map <- hashmap::hashmap(info$wd_names, as.integer(1:length(info$wd_names) - 1L))
  W <- lapply(keyATM_docs, function(x){ info$wd_map[[x]] })


  # Check keywords
  c <- sapply(unlist(keywords), 
         function(x){if(! x %in% info$wd_names)
           stop(paste0('"', x, '"', " does not appear in texts. Please check keywords."))})

  keywords_raw <- keywords  # keep raw keywords (not word_id)
  keywords_id <- lapply(keywords, function(x){ as.integer(info$wd_map$find(x)) })

  # Assign X and Z
  if(model %in% info$models_keyATM){
    res <- make_xz_key(W, keywords, info)
    X <- res$X
    Z <- res$Z
  }else{
    # LDA based models
    res <- make_xz_lda(W, info)
    X <- res$X
    Z <- res$Z
  }
  rm(res)

  # Organize
  stored_values <- list()

  if(model %in% c("basic", "lda")){
    if(options$estimate_alpha)
      stored_values$alpha_iter <- list()  
  }

  if(model %in% c("hmm", "ldahmm")){
    options$estimate_alpha <- 1
    stored_values$alpha_iter <- list()  
  }


  if(model %in% c("cov", "ldacov")){
    stored_values$Lambda_iter <- list()
  }

  if(model %in% c("hmm", "ldahmm")){
    stored_values$S_iter <- list()

    if(options$store_transition_matrix){
      stored_values$P_iter <- list()  
    }
  }

  if(options$store_theta)
    stored_values$Z_tables <- list()

  key_model <- list(
                    W = W, Z = Z, X = X,
                    model = model,
                    keywords = keywords_id, keywords_raw = keywords_raw,
                    regular_k = regular_k,
                    vocab = info$wd_names,
                    model_settings = model_settings,
                    priors = priors,
                    options = options,
                    stored_values = stored_values,
                    model_fit = list(),
                    call = match.call()
                   )

  rm(info)
  class(key_model) <- c("keyATM_model", class(key_model))

  if(options$iterations == 0){
    message("`options$iterations` is 0. keyATM returns an initialized object.")  
    return(key_model)
  }


  ##
  ## Fitting
  ##
  message(paste0("Fitting the model. ", options$iterations, " iterations..."))
  set.seed(options$seed)

  if(model == "basic"){
    key_model <- keyATM_train(key_model, iter = options$iterations, output_per = options$output_per)
  }else if(model == "cov"){
    key_model <- keyATM_train_cov(key_model, iter = options$iteration, output_per = options$output_per)
  }else if(model == "lda"){
    key_model <- LDA_weight(key_model, iter = options$iteration, output_per = options$output_per)  
  }else if(model == "hmm"){
    key_model <- keyATM_train_HMM(key_model, iter = options$iteration, output_per = options$output_per)  
  }else if(model == "ldahmm"){
    key_model <- keyATM_train_LDAHMM(key_model, iter = options$iteration, output_per = options$output_per)  
  }else{
    stop("Please check `mode`.")  
  }

  class(key_model) <- c("keyATM_fitted", class(key_model))
  return(key_model)
}


#' @noRd
#' @export
print.keyATM_model <- function(x)
{
  cat(
      paste0(
             "keyATM_model object for the ",
             x$model,
             " model.",
             "\n"
            )
     )
}


#' @noRd
#' @export
summary.keyATM_model <- function(x)
{
  cat(
      paste0(
             "keyATM_model object for the ",
             x$model,
             " model.",
             "\n"
            )
     )
}


#' @noRd
#' @export
save.keyATM_model <- function(x, file = stop("'file' must be specified"))
{
  save(x, file = file, compress="xz", compression_level = 3)
}


#' @noRd
#' @export
print.keyATM_fitted <- function(x)
{
  cat(
      paste0(
             "keyATM_model object for the ",
             x$model,
             " model. ",
             x$options$iterations, " iterations.\n",
             length(x$W), " documents | ",
             length(x$keywords), " keyword topics",
             "\n"
      )
     )
}


#' @noRd
#' @export
summary.keyATM_fitted <- function(x)
{
  cat(
      paste0(
             "keyATM_model object for the ",
             x$model,
             " model. ",
             x$options$iterations, " iterations.\n",
             length(x$W), " documents | ",
             length(x$keywords), " keyword topics",
             "\n"
      )
     )
}


#' @noRd
#' @export
save.keyATM_fitted <- function(x, file = stop("'file' must be specified"))
{
  save(x, file = file, compress="xz", compression_level = 3)
}


check_arg <- function(obj, name, model, info = list())
{
  if(name == "keywords"){
    return(check_arg_keywords(obj, model, info))
  }

  if(name == "model_settings"){
    return(check_arg_model_settings(obj, model, info))  
  }

  if(name == "priors"){
    return(check_arg_priors(obj, model, info))  
  }

  if(name == "options"){
    return(check_arg_options(obj, model, info))  
  }
}


check_arg_keywords <- function(keywords, model, info)
{
  check_arg_type(keywords, "list")

  if(length(keywords) == 0 & model %in% info$models_keyATM){
    stop("Please provide keywords.")  
  }

  if(length(keywords) != 0 & model %in% info$models_lda){
    stop("This model does not take keywords.")  
  }


  # Name of keywords topic
  if(model %in% info$models_keyATM){
    c <- lapply(keywords, function(x){check_arg_type(x, "character")})
  
    if(is.null(names(keywords))){
      names(keywords)  <- paste0(1:length(keywords))
    }else{
      names(keywords)  <- paste0(1:length(keywords), "_", names(keywords))
    }

  }

  return(keywords)
}

show_unused_arguments <- function(obj, name, allowed_arguments)
{
  unused_input <- names(obj)[! names(obj) %in% allowed_arguments]
  if(length(unused_input) != 0)
    stop(paste0(
                "keyATM doesn't recognize some of the arguments ",
                "in ", name, ": ",
                paste(unused_input, collapse=", ")
               )
        )
}


check_arg_model_settings <- function(obj, model, info)
{
  check_arg_type(obj, "list")
  allowed_arguments <- c()

  if(model %in% c("cov", "ldacov")){
     if(is.null(obj$covariates_data)){
      stop("Please provide `obj$covariates_data`.")  
    }

    check_arg_type(obj$covariates_data, "matrix")

    if(nrow(obj$covariates_data) != info$num_doc){
      stop("The row of `model_settings$covariates_data` should be the same as the number of documents.")  
    }

    if(sum(is.na(obj$covariates_data)) != 0){
      stop("Covariate data should not contain missing values.")
    }

    if(is.null(obj$covariates_formula)){
      obj$covariates_formula <- NULL
    }else if(is.formula(obj$covariates_formula)){
      message("Convert covariates data using `obj$covariates_formula`.")
      obj$covariates_data <- stats::model.matrix(obj$covariates_formula, obj$covariates_data)
    }else{
      stop("Check `model_settings$covariates_formula`.")  
    }

    allowed_arguments <- c(allowed_arguments, "covariates_data", "covariates_formula")
  }


  if(model %in% c("hmm", "ldahmm")){
    if(is.null(obj$num_states)){
      stop("`model_settings$num_states` is not provided.")  
    }

    if(is.null(obj$time_index)){
      stop("`model_settings$time_index` is not provided.")
    }

    if(length(obj$time_index) != info$num_doc){
      stop("The length of the `model_settings$time_index` does not match with the number of documents.")  
    }
    
    if(min(obj$time_index) != 1 | max(obj$time_index) > info$num_doc){
      stop("`model_settings$time_index` should start from 1 and not exceed the number of documents.")
    }

    if(max(obj$time_index) < obj$num_states)
      stop("`model_settings$num_states` should not exceed the maximum of `model_settings$time_index`.")

    check <- unique(obj$time_index[2:length(obj$time_index)] - lag(obj$time_index)[2:length(obj$time_index)])
    if(sum(!unique(check) %in% c(0,1)) != 0)
      stop("`model_settings$num_states` does not increment by 1.")

    obj$time_index <- as.integer(obj$time_index)

    allowed_arguments <- c(allowed_arguments, "num_states", "time_index")
    
  }

  show_unused_arguments(obj, "`model_settings`", allowed_arguments)

  return(obj)
}


check_arg_priors <- function(obj, model, info)
{
  check_arg_type(obj, "list")
  # Base arguments
  allowed_arguments <- c("beta")

  # prior of pi
  if(model %in% info$models_keyATM){
    if(is.null(obj$gamma)){
      obj$gamma <- matrix(1.0, nrow = info$total_k, ncol = 2)  
    }

    if(!is.null(obj$gamma)){
      if(dim(obj$gamma)[1] != info$total_k)  
        stop("Check the dimension of `priors$gamma`")
      if(dim(obj$gamma)[2] != 2)  
        stop("Check the dimension of `priors$gamma`")
    }


    if(info$keyword_k < info$total_k){
      # Regular topics are used in keyATM models
      # Priors of regular topics should be 0
      if(sum(obj$gamma[(info$keyword_k+1):info$total_k, ]) != 0){
        obj$gamma[(info$keyword_k+1):info$total_k, ] <- 0
      }
    }

    allowed_arguments <- c(allowed_arguments, "gamma")
  }


  # beta
  if(is.null(obj$beta)){
    obj$beta <- 0.01  
  }

  if(model %in% info$models_keyATM){
    if(is.null(obj$beta_s)){
      obj$beta_s <- 0.1  
    }  
    allowed_arguments <- c(allowed_arguments, "beta_s")
  }


  if(model %in% c("basic", "lda")){
    # alpha
    if(is.null(obj$alpha)){
      obj$alpha <- rep(1/info$total_k, info$total_k)
    }
    if(length(obj$alpha) != info$total_k){
      stop("Starting alpha must be a vector of length ", info$total_k)
    }
    allowed_arguments <- c(allowed_arguments, "alpha")
  
  }

  show_unused_arguments(obj, "`priors`", allowed_arguments)

  return(obj)
}


check_arg_options <- function(obj, model, info)
{
  check_arg_type(obj, "list")
  allowed_arguments <- c("seed", "output_per", "thinning",
                         "iterations",
                         "use_weights",
                         "store_theta", "slice_shape")

  # Output per
  if(is.null(obj$output_per))
    obj$output_per <- 10L

  if(!is.numeric(obj$output_per) | obj$output_per < 0 | obj$output_per%%1!=0){
      stop("An invalid value in `options$output_per`")  
  }

  # thinning
  if(is.null(obj$thinning))
    obj$thinning <- 1L

  if(!is.numeric(obj$thinning) | obj$thinning < 0| obj$thinning%%1!=0){
      stop("An invalid value in `options$thinning`")  
  }

  # seed
  if(is.null(obj$seed))
    obj$seed <- floor(runif(1)*1e5)

  # iterations
  if(is.null(obj$iterations))
    obj$iterations <- 1500L
  if(!is.numeric(obj$iterations) | obj$iterations < 0| obj$iterations%%1!=0){
      stop("An invalid value in `options$iterations`")  
  }

  # Store theta
  if(is.null(obj$store_theta)){
    obj$store_theta <- 0L
  }else{
    obj$store_theta <- as.integer(obj$store_theta)  
    if(!obj$store_theta %in% c(0, 1)){
      stop("An invalid value in `options$store_theta`")  
    }
  }

  # Estimate alpha
  if(model %in% c("basic", "lda")){
    if(is.null(obj$estimate_alpha)){
      obj$estimate_alpha <- 1L
    }else{
      obj$estimate_alpha <- as.integer(obj$estimate_alpha)  
      if(!obj$estimate_alpha %in% c(0, 1)){
        stop("An invalid value in `options$estimate_alpha`")  
      }

    }
    allowed_arguments <- c(allowed_arguments, "estimate_alpha")
  }
  
  # Slice shape
  if(is.null(obj$slice_shape)){
    # parameter for slice sampling
    obj$slice_shape <- 1.2
  }
  if(!is.numeric(obj$slice_shape) | obj$slice_shape < 0){
      stop("An invalid value in `options$slice_shape`")  
  }

  # Use weights
  if(is.null(obj$use_weights)){
    obj$use_weights <- 1L 
  }else{
    obj$use_weights <- as.integer(obj$use_weights)
    if(!obj$use_weights %in% c(0, 1)){
      stop("An invalid value in `options$use_weights`")  
    }
  }

  if(model %in% c("hmm", "ldahmm")){
    if(is.null(obj$store_transition_matrix)){
      obj$store_transition_matrix <- 0L  
    }
    if(!obj$store_transition_matrix %in% c(0, 1)){
      stop("An invalid value in `options$store_transition_matrix`")  
    }
    allowed_arguments <- c(allowed_arguments, "store_transition_matrix")
  }

  # Check unused arguments
  show_unused_arguments(obj, "`options`", allowed_arguments)
  return(obj)
}


make_xz_key <- function(W, keywords, info)
{
  # zx_assigner maps keywords to category ids
  key_wdids <- unlist(lapply(keywords, function(x){ info$wd_map$find(x) }))
  cat_ids <- rep(1:(info$keyword_k) - 1L, unlist(lapply(keywords, length)))

  if(length(key_wdids) == length(unique(key_wdids))){
    #
    # No keyword appears more than once
    #
    zx_assigner <- hashmap::hashmap(as.integer(key_wdids), as.integer(cat_ids))

    # if the word is a keyword, assign the appropriate (0 start) Z, else a random Z
    topicvec <- 1:(info$total_k) - 1L
    make_z <- function(x, topicvec){
      zz <- zx_assigner[[x]] # if it is a keyword word, we already know the topic
      zz[is.na(zz)] <- sample(topicvec,
                              sum(is.na(zz)),
                              replace = TRUE)
      return(zz)
    }

  }else{
    #
    # Some keywords appear multiple times
    #
    keys_df <- data.frame(wid = key_wdids, cat = cat_ids)
    keys_char <- sapply(unique(key_wdids),
                        function(x){
                          paste(as.character(keys_df[keys_df$wid==x, "cat"]), collapse=",")
                        })
    zx_hashtable <- hashmap::hashmap(as.integer(unique(key_wdids)), keys_char)

    zx_assigner <- function(x){
      topic <- zx_hashtable[[x]]
      topic <- strsplit(topic, split=",")
      topic <- lapply(topic, sample, 1)
      topic <- as.integer(unlist(topic))
      return(topic)
    }

    # if the word is a seed, assign the appropriate (0 start) Z, else a random Z
    topicvec <- 1:(info$total_k) - 1L
    make_z <- function(x, topicvec){
      zz <- zx_assigner(x) # if it is a seed word, we already know the topic
      zz[is.na(zz)] <- sample(topicvec,
                              sum(is.na(zz)),
                              replace = TRUE)
      return(zz)
    }
  }


  ## xx indicates whether the word comes from a seed topic-word distribution or not
  make_x <- function(x){
    key <- as.numeric(x %in% key_wdids) # 1 if they're a seed
    # Use x structure
    x[key == 0] <- 0L # non-keyword words have x = 0
    x[key == 1] <- sample(0:1, length(x[key == 1]), prob = c(0.3, 0.7), replace = TRUE)
      # keywords have x = 1 probabilistically
    return(x)
  }

 X <- lapply(W, make_x)
 Z <- lapply(W, make_z, topicvec)

 return(list(X = X, Z = Z))
}


make_xz_lda <- function(W, info)
{
  topicvec <- 1:(info$total_k) - 1L
  make_z <- function(x, topicvec){
    zz <- sample(topicvec,
                 length(x),
                 replace = TRUE)
    return(as.integer(zz))
  }  

  make_x <- function(x){
    return(rep(0L, length(x)))  
  }

 X <- lapply(W, make_x)
 Z <- lapply(W, make_z, topicvec)

 return(list(X = X, Z = Z))
}


#' keyATM Main function
#'
#' This is a wrapper function of \code{keyATM_fit()} and \code{keyATM_output()}.
#'
#'
#' @param keyATM_docs texts read via \code{keyATM_read()}
#' @param model keyATM model: "basic", "cov", "hmm", "lda", "ldacov" and "ldahmm"
#' @param regular_k the number of regular topics
#' @param keywords a list of keywords
#' @param model_settings a list of model specific settings
#' @param priors a list of priors of parameters
#' @param options a list of options
#' @param keep a vector of the names of elements you want to keep from \code{keyATM_fit()} output
#'
#' @return A keyATM_output object containing:
#'   \describe{
#'     \item{keyword_k}{Number of keyword topics}
#'     \item{regular_k}{Number of regular unseeded topics}
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
#'     \item{kept_values}{Output from \code{keyATM_fit()} you specified to store.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # keyATM Basic
#'   out <- keyATM(
#'                 keyATM_docs, model = "basic", regular_k = 5, keywords = keywords_list
#'                )
#'
#'   # keyATM Cov
#'   out <- keyATM(
#'                 keyATM_docs, model = "cov", regular_k = 5, keywords = keywords_list,
#'                 model_settings(covariates_data = cov)
#'                )
#'
#'   # keyATM HMM
#'   out <- keyATM(
#'                    keyATM_docs, model = "hmm", regular_k = 5, keywords = keywords_list,
#'                    model_settings(time_index = time_index_vec, num_states = 5)
#'                   )
#'
#'   # Weighted LDA
#'   out <- keyATM(
#'                 keyATM_docs, model = "lda", regular_k = 5
#'                )
#'
#'   # Weighted LDA Cov
#'   out <- keyATM(
#'                 keyATM_docs, model = "ldacov", regular_k = 5,
#'                 model_settings(covariates_data = cov)
#'                )                   
#'
#'   # Weighted LDA HMM
#'   out <- keyATM(
#'                 keyATM_docs, model = "ldahmm", regular_k = 5,
#'                 model_settings(time_index = time_index_vec, num_states = 5)
#'                )
#'
#' }
#'
#' @export
keyATM <- function(keyATM_docs, model, regular_k,
                   keywords = list(), model_settings = list(),
                   priors = list(), options = list(), keep = c())
{
  # Check type
  if(length(keep) != 0)
    check_arg_type(keep, "character")

  # Fit keyATM
  fitted <- keyATM_fit(
                       keyATM_docs, model, regular_k,
                       keywords, model_settings, priors, options
                      )

  # Get output
  out <- keyATM_output(fitted)

  # Keep some objects if specified
  if(length(keep) != 0){
    kept_values <- list()
    use_elements <- keep[keep %in% names(fitted)]
    for(i in 1:length(use_elements)){
      kept_values[use_elements[i]]  <- fitted[use_elements[i]]
    }
    out$kept_values <- kept_values
  }

  return(out)
}


