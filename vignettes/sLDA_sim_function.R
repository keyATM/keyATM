#####################################################################################
#####################################################################################
## functions to run simulation
## if you run code in vignette, you can use the following function by 
## >> source("./sLDA_sim_function")
#####################################################################################
#####################################################################################

#######################################
#######################################
## functions to create topicdict model
#######################################
#######################################

create_model <- function(folder, extra_k, 
                                     seed_file_name="", seed_list = NULL,
                                     cull_seed=NULL){

  # cull_seed: Which row (=topic) and column (=number of seed) to use

  # folder <- paste0(data_folder, data_folder_name, "/W")
  docs <- list.files(folder, pattern = "*.txt", full.names = TRUE)

  if(seed_file_name != ""){
    seed_file <- paste0(data_folder, data_folder_name, "/", seed_file_name)
    seeddata <- invisible(readr::read_csv(seed_file))
  }

  if(!is.null(cull_seed)){
    seeddata <- seeddata[1:cull_seed[1], 1:cull_seed[2]]
  }

  if(is.null(seed_list)){
    seed_list <- as.list( apply(seeddata, 1, function(x){return(paste0(x, collapse=" "))}) )
  }
  names(seed_list) <- 1:length(seed_list)
  dict <- quanteda::dictionary(seed_list)

  model <- topicdict_model(docs,
               dict = dict, extra_k = extra_k,
               remove_numbers = FALSE, # For simulation, make it false
               remove_punct = TRUE,
               remove_symbols = TRUE,
               remove_separators = TRUE)

  return(model)
}

#######################################
#######################################
## functions to create simultion data
#######################################
#######################################
rbern <- function(n, prob){
  return(rbinom(n, 1, prob))
}
rcat <- function(n, p){
  if (is.vector(p)) {
    x <- as.vector(which(rmultinom(n, size=1, prob=p) == 1, arr.ind=TRUE)[, "row"])
  } else {
    d <- dim(p)
    n <- d[1]
    k <- d[2]
    lev <- dimnames(p)[[2]]
    if (!length(lev)) lev <- 1:k
    z <- colSums(p)
    U <- apply(p, 1, cumsum)
    U[,k] <- 1
    un <- rep(runif(n), rep(k,n))
    x <- lev[1 + colSums(un > U)]}
    return(x)
}
## From MCMCpack
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l*n,alpha), ncol=l, byrow=TRUE)
  sm <- x %*% rep(1,l)
  return(x / as.vector(sm))
}

CreateDir_wzx <- function(saveDir){
  dir1 <- paste0(saveDir, "W")
  dir.create(dir1)
  dir2 <- paste0(saveDir, "Z")
  dir.create(dir2)
  dir3 <- paste0(saveDir, "X")
  dir.create(dir3)
  dir4 <- paste0(saveDir, "Wraw")
  dir.create(dir4)
  return(list(dir1, dir2, dir3, dir4))
}

Gen_ND <- function(document, lambda){
  nd <- numeric(document)
  ## generate document length
  nd <- rpois(document, lambda)
  return(nd)
}

Gen_theta <- function(doc_len, alpha_vec){
  theta <- rdirichlet(doc_len, alpha_vec)
  return(theta)
}

Gen_z <- function(theta, d, doc_len){
  z <- rcat(doc_len, theta[d,])
  return(z)
}

Gen_x <- function(topic_num, probX){
  x <- rbern(1, probX[topic_num])
  return(x)
}

Gen_w <- function(phi, topic){
  tmp <- rcat(1, phi[topic,])
  words <- colnames(phi)
  return(words[tmp])
}

Gen_w_seeds <- function(index, phiR, phiS, z, x){
  topic <- z[index]
  indicator <- x[index]

  if (indicator == 0){
    # Regular words
    word <- Gen_w(phiR, topic)
    # word <- original_vocabulary[as.numeric(word)] # original_vocabulary is extracted from the trained model outside
    word <- paste0("W", as.character(word), "t", as.character(topic))
  } else {
    # Seed words
    word <- Gen_w(phiS, topic)
    # word <- original_vocabulary[as.numeric(word)]
    word <- paste0("W", as.character(word), "t", as.character(topic))
  }
  return(word)
}

Gen_wzx <- function(doc_id, doc_len, phiR, phiS, p_vec, theta){
  z <- Gen_z(theta, doc_id, doc_len)
  x <- sapply(z, Gen_x, probX=p_vec)
  w <- sapply(1:doc_len, Gen_w_seeds, phiR=phiR, phiS=phiS, z=z, x=x)
  # add functions to take away the topic id on words
  # print(w)
  w_raw <- sub("t\\d", "", w)
  # print(w_raw)
  return(list(doc_id=rep(doc_id, doc_len), w=w,z=z,x=x, w_raw = w_raw))
}

Write_wzx <- function(content, name, saveDir){
  txt <- paste(as.character(content), collapse = " ")
  name <- paste0(saveDir, name, ".txt")
  write(txt, name)
}

Save_wzx <- function(data, i, saveDir){
  ## save W for document i
  Write_wzx(data[["w"]], paste0("text", i), paste0(saveDir,"W/"))

  ## save Z for document i
  Write_wzx(data[["z"]], paste0("z_", i), paste0(saveDir,"Z/"))

  ## save X for document i
  Write_wzx(data[["x"]], paste0("x_", i), paste0(saveDir,"X/"))

  ## save raw W (without topic indicator) for docuent i
  Write_wzx(data[["w_raw"]], paste0("w_raw_", i), paste0(saveDir,"Wraw/"))
}

create_sim_data_from_fitted <- function(
  saveDir, alpha, phiR, phiS, p,
  D=200, lambda=300){
  # Set Seed
  # set.seed(rand_seed)

  # Prepare Directory
  dir.create(saveDir)
  Dirs <- CreateDir_wzx(saveDir)

  # Set the number of topics
  K <- nrow(phiS)
  phiR <- phiR[1:K, ]
  alpha_vec <- alpha[1:K]
  V <- ncol(phiR)
  p_vec <- p[1:K]

  # Length of the documents
  nd <- Gen_ND(D, lambda)

  # Get theta
  theta <- Gen_theta(D, alpha_vec)

  ## Generate topic and word for each word in each document
  doc_list <- vector("list", D)
  
  for (i in 1:D){
    # the length of the document i
    tmp_nd <- nd[i]
  
    # ## a list of vector to contain the results
    # tmp <- WZX_Vec(tmp_nd)
  
    ## Generate W, Z, and X
    wzx_list <- Gen_wzx(i, tmp_nd, phiR, phiS, p_vec, theta)
  
    ## save W, Z, X for document i
    Save_wzx(wzx_list, i, saveDir)
  
    ## Store Created Docs
    doc_list[[i]] <- wzx_list
  
  }
}

show_true_seeds <- function(phiS){
  n <- nrow(phiS)
  seeds_ <- c()
  for(i in 1:n){
    temp <- phiS[i, ]
    seeds_ <- c(seeds_, names(temp[temp > 1e-8]))
  }
  return(seeds_)
}


#######################################
#######################################
## for diagnosis
#######################################
#######################################


compare_seededLDA <- function(seed_list, total_k, folder_name="Sim1"){
  extra_k <- total_k - length(seed_list)
  model <- create_model(folder_name, seed_list=seed_list, extra_k=extra_k)
  res <- topicdict_train(model, iter = iter_num)
  post <- topicdict::posterior(res)
  diagnosis_topic_recovery_heatmap(post, 25)
}

## return a plot which shows the proportion of match between true topic and estimated topic among the top words
## you can use this for text data with topic indicator "w1t1", not just word
## post: posteroior object
diagnosis_topic_recovery_heatmap <- function(post, n=25, 
                        topicvec=c(), merge=list()){
  topwords <- top_terms(post, n=n)
  topwords <- data.frame(topwords)
  colnames(topwords) <- paste0("EstTopic", 1:ncol(topwords))

  topwords <- tidyr::gather(topwords, key=EstTopic, value=Word)

  topwords %>%
    mutate(RawWord = Word) %>%
    tidyr::separate(Word,
        into=c("word_id", "TrueTopic"),
        sep="t") %>%
    mutate(TrueTopic = paste0("True", as.character(TrueTopic))) -> res_

  res_$TrueTopic <- sub(" \\[✓\\]", "", res_$TrueTopic)

  merge_length <- length(merge)
  if(merge_length != 0){
    # Merge Topics
    for(i in 1:merge_length){
      m <- merge[[i]]
      mt <- paste0("True", m)

      res_ %>%
        mutate(TrueTopic=replace(TrueTopic, TrueTopic==mt[1], mt[3])) %>%
        mutate(TrueTopic=replace(TrueTopic, TrueTopic==mt[2], mt[3])) -> res_
    }
  }

  res_ %>%
    group_by(EstTopic, TrueTopic) %>%
    summarise(counts = n()) %>%
    ungroup() %>%
    group_by(EstTopic) %>% 
    mutate(topicsum = sum(counts)) %>%
    ungroup() %>%
    mutate(Proportion = counts / topicsum * 100) -> res_

  num <- length(unique(res_$EstTopic))
  # if(is.null(topicvec)){
    # res_ %>%
    #   group_by(EstTopic) %>%
    #   top_n(1, Proportion) %>%
    #   arrange(TrueTopic) %>%
    #   select(EstTopic) -> topicvec 
    # topicvec <- unique(as.integer(gsub("EstTopic", "", topicvec$EstTopic)))
  # }else if(length(topicvec) != num){
    # message("topicvec length does not match")
    # topicvec <- 1:num
  # }
  topicvec <- 1:num
  truenum <- length(unique(res_$TrueTopic))

  title <- paste0("Seeded LDA: Top ", as.character(n), " words")

  g <- ggplot(res_, aes(EstTopic, TrueTopic)) +
        geom_tile(aes(fill=Proportion)) + 
        scale_fill_gradient(limits=c(0, 100), low="#e8e8e8", high="#0072B2", name = "Proportion") +
        scale_x_discrete(limits = rev(paste0("EstTopic", topicvec))) +
        coord_flip() +
        # scale_y_discrete(limits = paste0("True", 1:truenum)) +
        xlab("Estimated Topics") + ylab("True Topic") + theme_bw(base_size=13) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))

  return(g)
}

## function to create semi-auto dictionary based on term-frequency
## eobj should be an outut of explore function
semiauto_dictionary <- function(eobj, num_seeds, 
                            prop_center=rep(0.3, length(num_seeds))) {

  if(length(num_seeds) != length(prop_center) | length(num_seeds) != length(window) ){
    stop("Length of prop_center and num_seeds should be equal to the number of seeded topics.")
  }

  seed_list <- list()
  for(i in 1:length(num_seeds)){
    num <- num_seeds[[i]]
    min <- prop_center[i] - window[i]
    max <- prop_center[i] + window[i]
    eobj$data %>%
      filter(`Proportion(%)` > get("min") & `Proportion(%)` < get("max")) -> temp

    if(nrow(temp) == 0){
      stop(paste0("There is no word that meets the specified condition. Recheck topic "), i, ".")
    } else if(nrow(temp) < num){
      warning(paste0("The number of words that match the condition is fewer than the number you specidied. Recheck topic "), i, ".")
      seed_list <- c(seed_list, list(paste(temp$Word, collapse=" ")))
    }else{
      select <- sample(1:nrow(temp), num)
      seed_list <- c(seed_list, list(paste(temp$Word[select], collapse=" ")))
    }
  }

  return(seed_list)
}

## eobj: output of explore object
## topic_num: the number of topics to be assigned keywords
## keyword_num: the number of keywords to be assigned to each topic
## top: choose top frequent words as keywords or not
auto_dictionary <- function(eobj, topic_num = 3, 
                            keyword_num = 3, 
                            prop = 3, top = TRUE){
  seed_list <- list()
  ## if the length of key_word_num does not match topic_num
  ## forced to choose the first element of keyword_num and prop
  if(length(keyword_num) != topic_num){keyword_num <- rep(keyword_num[1], topic_num)}
  ## if you want to choose keywords from top words
  if (top == TRUE){
    for (i in 1:topic_num){
      ## set topic indicator as topic_id (i.g. "t1$")
      topic_id <- paste0("t", i, "$")
      ## set topic indicator as topic_id (i.g. "t1$")
      rows <- head(grep(topic_id, eobj$data$Word), keyword_num[i])
      keywords <- gsub("t\\d*$", "", eobj$data$Word[rows])
      seed_list[[i]] <- paste(keywords, collapse = " ")
    }
  } else {
    for (i in 1:topic_num){
      if(length(prop) != topic_num){prop <- rep(prop[1], topic_num)}
      cut_rows <- which.min(abs(eobj$data$`Proportion(%)` - prop[i])) - 1 # rows ignored
      if (cut_rows != 0) {
          eobj_dw_cutted <- eobj$data$Word[-c(1:cut_rows)]
        } else {
          eobj_dw_cutted <- eobj$data$Word
        }      
      topic_id <- paste0("t", i, "$")
      rows <- head(grep(topic_id, eobj_dw_cutted), keyword_num[i])
      keywords <- gsub("t\\d*$", "", eobj_dw_cutted[rows])
      seed_list[[i]] <- paste(keywords, collapse = " ")
    }
  }
  return(seed_list)
}


## function to create a plot which shows the accuracy of LDA result
## data_path should be a file that contains the txt files ("/W")
get_lda_result <- function(data_path, iter, k, topicvec=1:k, show_n=25){
  # folder <- paste0(data_folder, data_folder_name, "/W")
  folder <- data_path

  # Prepare Data
  corpus <- Corpus(DirSource(folder))
  strsplit_space_tokenizer <- function(x)
      unlist(strsplit(as.character(x), "[[:space:]]+"))

  dtm <- DocumentTermMatrix(corpus,
                           control = list(tokenize=strsplit_space_tokenizer, 
                           stopwords = F, tolower = F, 
                           stemming = F, wordLengths = c(1, Inf)))

  lda <- LDA(dtm, k = k, control = list(seed = 225, iter=iter), method="Gibbs")

  assign <-tidytext::augment(lda, dtm)
  assign %>% 
        mutate(Word=term) %>%
        tidyr::separate(term,
            into=c("word_id", "TrueTopic"),
            sep="t") %>%
        mutate(TrueTopic = paste0("True", as.character(TrueTopic))) %>%
        group_by(.topic, TrueTopic, Word) %>%
        summarise(counts = sum(count)) %>%
        ungroup() %>%
        group_by(.topic) %>%
        top_n(show_n, counts) %>%
        ungroup() %>%
        group_by(.topic, TrueTopic) %>%
        summarise(counts = sum(counts)) %>%
        ungroup() %>%
        group_by(.topic) %>%
        mutate(topicsum = sum(counts)) %>%
        ungroup() %>%
        mutate(Proportion = counts / topicsum * 100,
               EstTopic = paste0("EstTopic", .topic)) -> res

  res %>%
    group_by(EstTopic) %>%
    top_n(1, Proportion) %>%
    arrange(TrueTopic) %>%
    select(EstTopic) -> topicvec 
  topicvec <- unique(as.integer(gsub("EstTopic", "", topicvec$EstTopic)))

  num <- length(unique(res$EstTopic))
  truenum <- length(unique(res$TrueTopic))
  title <- paste0("LDA: Top ", as.character(show_n), " words")

  g <- ggplot(res, aes(EstTopic, TrueTopic)) +
        geom_tile(aes(fill=Proportion)) + 
        scale_fill_gradient(limits=c(0, 100), low="#e8e8e8", high="#0072B2") +
        scale_x_discrete(limits = rev(paste0("EstTopic", topicvec))) +
        coord_flip() +
        scale_y_discrete(limits = paste0("True", 1:truenum)) +
        xlab("Estimated Topics") + ylab("True Topic") + theme_bw(base_size=13) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))

  return(g)
}

## function to return the result of LDA
## data_path should be a file that contains the txt files ("/W")
get_lda_result2 <- function(data_path, iter, k){
  # folder <- paste0(data_folder, data_folder_name, "/W")
  folder <- data_path

  # Prepare Data
  corpus <- Corpus(DirSource(folder))
  strsplit_space_tokenizer <- function(x)
      unlist(strsplit(as.character(x), "[[:space:]]+"))

  dtm <- DocumentTermMatrix(corpus,
                           control = list(tokenize=strsplit_space_tokenizer, 
                           stopwords = F, tolower = F, 
                           stemming = F, wordLengths = c(1, Inf)))

  lda <- LDA(dtm, k = k, control = list(seed = 225, iter=iter), method="Gibbs")
  return (lda)
}


## function to return the result of LDA, given the result of LDA
## data_path should be a file that contains the txt files ("/W")
plot_lda_result <- function(data_path, lda_result, topicvec=1:k, show_n=25){
  folder <- data_path

  # Prepare Data
  corpus <- Corpus(DirSource(folder))
  strsplit_space_tokenizer <- function(x)
      unlist(strsplit(as.character(x), "[[:space:]]+"))

  dtm <- DocumentTermMatrix(corpus,
                           control = list(tokenize=strsplit_space_tokenizer, 
                           stopwords = F, tolower = F, 
                           stemming = F, wordLengths = c(1, Inf)))

  assign <-tidytext::augment(lda_result, dtm)
  assign %>% 
        mutate(Word=term) %>%
        tidyr::separate(term,
            into=c("word_id", "TrueTopic"),
            sep="t") %>%
        mutate(TrueTopic = paste0("True", as.character(TrueTopic))) %>%
        group_by(.topic, TrueTopic, Word) %>%
        summarise(counts = sum(count)) %>%
        ungroup() %>%
        group_by(.topic) %>%
        top_n(show_n, counts) %>%
        ungroup() %>%
        group_by(.topic, TrueTopic) %>%
        summarise(counts = sum(counts)) %>%
        ungroup() %>%
        group_by(.topic) %>%
        mutate(topicsum = sum(counts)) %>%
        ungroup() %>%
        mutate(Proportion = counts / topicsum * 100,
               EstTopic = paste0("EstTopic", .topic)) -> res

  res %>%
    group_by(EstTopic) %>%
    top_n(1, Proportion) %>%
    arrange(TrueTopic) %>%
    select(EstTopic) -> topicvec 
  topicvec <- unique(as.integer(gsub("EstTopic", "", topicvec$EstTopic)))

  num <- length(unique(res$EstTopic))
  truenum <- length(unique(res$TrueTopic))
  title <- paste0("LDA: Top ", as.character(show_n), " words")

  g <- ggplot(res, aes(EstTopic, TrueTopic)) +
        geom_tile(aes(fill=Proportion)) + 
        scale_fill_gradient(limits=c(0, 100), low="#e8e8e8", high="#0072B2") +
        scale_x_discrete(limits = rev(paste0("EstTopic", topicvec))) +
        coord_flip() +
        scale_y_discrete(limits = paste0("True", 1:truenum)) +
        xlab("Estimated Topics") + ylab("True Topic") + theme_bw(base_size=13) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))

  return(g)

}

## function to compare the estimated results VS true W
## in this function, put the words with topic indicator, "w1t1"
read_true_word <- function(folder, encoding = "UTF-8"){
  files <- list.files(folder, pattern = "*.txt", full.names = TRUE)
  # df <- lapply(files, function(x){readLines(x, encoding = encoding)})
  
  df <- lapply(files, function(x){ paste0(readLines(x, encoding = encoding), collapse = "\n") })
  df <- sapply(df, function(x){ strsplit(x, split = " ")})
  return (df)
}


## return plot
## post: posterior object
## res: a list which contains true words for all the docuemnts (term with topic id), topicdict object
## true_word: a list which contains teh true word vectors, output of the read_true_word
## n: top n terms you use to calculate
## num_true_topic: the number of true topics (maximum number followed by "t")
match_truth_Z <- function(post, res, true_word, n = 10, num_true_topic = 6){
  if ("plyr" %in% (.packages())){detach(package:plyr)}
  # check_arg_type(post, "topicdict_posterior")
  topterms <- top_terms(post, n)
  topterms <- data.frame(topterms)
  colnames(topterms) <- paste0("EstTopic", 1:ncol(topterms))

  ## output is a dataframe
  ## first column: "EstTopic*", which is an estimated topic number each
  ## second column: terms 
  topterms <- tidyr::gather(topterms, key = EstTopic, value = term)

  ## take away the checkmark which indicates keywords
  topterms$term <- gsub(" \\[✓\\]", "", topterms$term)

  ## "W" and "w" is mixed so use "W"
  topterms$term <- sub("w", "W", topterms$term)

  ## new column which contains the estimated topic number
  topterms$EstTopicNum <- gsub("EstTopic", "", topterms$EstTopic)

  ## mid is a matrix. Each row is a term and each column shows counts of the true topic for the term (given the term XXX, what is the true topic of XXXs?)
  ## the number of estimated topic \times (n \times estimated topic)
  mid <- apply(topterms, 1, function(x, res, true_word, num_true_topic){count_Z(res, true_word, x[2], x[3], num_true_topic)}, res, true_word, num_true_topic)

  ## bind with dataframe which contains the information of EstTopic* and terms
  top <- cbind(topterms, t(mid))
  colnames(top)[4:(4 + num_true_topic - 1)] <- paste0("True", 1:num_true_topic)

  tt <- top[-c(2,3)]

  ## Calculate the proportion of each true topic for each term
  if ("plyr" %in% (.packages())){detach(package:plyr)}
  tt %>% 
    tidyr::gather(key=True, value = term_freq, -EstTopic) %>% 
    group_by(EstTopic, True) %>%
    mutate(freq = sum(term_freq)) %>%
    ungroup(EstTopic, True) %>%
    select(-c(term_freq)) %>%
    distinct(EstTopic, True, .keep_all = TRUE) %>%
    group_by(EstTopic) %>% 
    mutate(Proportion = freq/sum(freq)) -> tt_new ## tt_new: a tibble

  ## vectors used for the label of axis
  est_topicvec <- 1:length(unique(tt_new$EstTopic))
  true_topicvec <- 1:num_true_topic
  g <- ggplot(tt_new, aes(EstTopic, True)) +
      geom_tile(aes(fill=Proportion)) + 
      scale_fill_gradient(limits=c(0, 1), low="white", high="#0072B2") +
      scale_x_discrete(limits = rev(paste0("EstTopic", est_topicvec))) +
      coord_flip() +
      scale_y_discrete(limits = paste0("True", true_topicvec)) +
      xlab("Estimated Topics") + ylab("True Topic") + theme_bw(base_size=13) +
      theme(plot.title = element_text(hjust = 0.5)) 
  return(g)
  # return(tt_new)
}

## find a match between estimated topic for a term and the actual true topic assigned to that term, and then count the matches
## return a scaler, numeric
## x and y should be a vector
## x: res$Z[[i]], which is estimated topic
## y: a vector which contains true words for one document (each document, term with topic id) 
## term: a term (string) you are lookng at with topic id ("w1t1")
## lookup_id: id which you are looking at
find_match <- function(x, y, term, lookup_id){
  id <- which(y %in% term) # if there is no match return NA
  if (identical(id, integer(0)) != TRUE){
    est_topic <- x[id]
    est_topic <- est_topic + 1 # because topic id starts with 0
    return (sum( as.numeric(est_topic == lookup_id)) ) ## count the number of match
  }
}

## given term and estimated topic, this is to count the frequency
## return vector
## post: posterior object
## true_word: a list which contains true words for all the docuemnts (term with topic id) 
## term: a term (string) you are lookng at ("w1", without topic id)
## lookup_id: id which you are looking at
## num_true_topic: the number of true topics (maximum number followed by "t")
count_Z <- function(res, true_word, term, lookup_id, num_true_topic = 6){
  out <- c()
  for (i in 1:num_true_topic){
    term_new <- paste0(term, "t", i)
    out[i] <- sum( unlist( mapply(find_match, res$Z, true_word, term_new, lookup_id) ) )
  }
  return(out)
}

## return plot
## post: posterior object
## res: a list which contains true words for all the docuemnts (term with topic id), topicdict object
## true_word: a list which contains teh true word vectors, output of the read_true_word
## n: top n terms you use to calculate
## num_true_topic: the number of true topics (maximum number followed by "t")
match_truth_Z_part <- function(post, res, true_word, n = 10, num_true_topic = 6){
  
  if ("plyr" %in% (.packages())){detach(package:plyr)}
  # check_arg_type(post, "topicdict_posterior")
  topterms <- top_terms(post, n)
  topterms <- data.frame(topterms)
  colnames(topterms) <- paste0("EstTopic", 1:ncol(topterms))

  ## output is a dataframe
  ## first column: "EstTopic*", which is an estimated topic number each
  ## second column: terms 
  topterms <- tidyr::gather(topterms, key = EstTopic, value = term)

  ## take away the checkmark which indicates keywords
  topterms$term <- gsub(" \\[✓\\]", "", topterms$term)

  ## "W" and "w" is mixed so use "W"
  topterms$term <- sub("w", "W", topterms$term)

  ## new column which contains the estimated topic number
  topterms$EstTopicNum <- gsub("EstTopic", "", topterms$EstTopic)

  ## mid is a matrix. Each row is a term and each column shows counts of the true topic for the term (given the term XXX, what is the true topic of XXXs?)
  ## the number of estimated topic \times (n \times estimated topic)
  mid <- apply(topterms, 1, function(x, res, true_word, num_true_topic){count_Z(res, true_word, x[2], x[3], num_true_topic)}, res, true_word, num_true_topic)

  ## bind with dataframe which contains the information of EstTopic* and terms
  top <- cbind(topterms, t(mid))
  colnames(top)[4:(4 + num_true_topic - 1)] <- paste0("True", 1:num_true_topic)

  tt <- top[-c(2,3)]

  ## Calculate the proportion of each true topic for each term
  if ("plyr" %in% (.packages())){detach(package:plyr)}
  tt %>% 
    tidyr::gather(key=True, value = term_freq, -EstTopic) %>% 
    group_by(EstTopic, True) %>%
    mutate(freq = sum(term_freq)) %>%
    ungroup(EstTopic, True) %>%
    select(-c(term_freq)) %>%
    distinct(EstTopic, True, .keep_all = TRUE) %>%
    group_by(EstTopic) %>% 
    mutate(Proportion = freq/sum(freq)) -> tt_new ## tt_new: a tibble

  return(tt_new)
}




## Calculate KL-divergence
## compare KL divergence between LDA results and seeded lda results
## extract "true" beta
## true_word: raw word, with topic indicator, output of the read_true_word
## return list which contains 
## [[1]]: term-topic matrix for slda results
## [[2]]: term-topic matrix for lda results
## [[3]]: term-topic matrix for truth
term_topic_mat <- function(true_word, post, lda_res, true_K){
  ## slda
  slda.mm <- post$beta
  slda.mm <- slda.mm[, order(colnames(slda.mm))]
  colnames(slda.mm) <- gsub("w", "W", colnames(slda.mm))

  ## extract term-topic distirbution from lda results
  lda.mm <- as.matrix(topicmodels::posterior(lda_res)$terms)
  ## sort colmanes: note that since lda and slda use the same document, 
  ## order function should result in the same result
  lda.mm <- lda.mm[, order(colnames(lda.mm))]

  ## extract unique words
  unique_term <- colnames(slda.mm)

  ## create matrix to contain the "true" topic assignment
  mat <- matrix(0, nrow =  true_K, ncol = length(unique_term))
  colnames(mat) <- unique_term
  ## read topic assignment and put into matrix mat
  true_count <- table( unlist(true_word) )
  terms <- unlist(attr(true_count, "dimnames"))
  for (i in 1:length(true_count)){
    topic_id <- as.numeric(  gsub("W\\d+t","", terms[i]) )
    term <- gsub("t\\d+","", terms[i])
    mat[topic_id, term] <- true_count[i]
  }
  mat <- t(apply(mat, 1, function(x){x/sum(x)}))
  ret <- list(slda.mm, lda.mm, mat)
  return(ret)
}

## px and qx should be vectors
cal_KL <- function(px, qx){
  kl <- (px + 1e-10) %*% log(px + 1e-10) - (px + 1e-10) %*% log(qx + 1e-10)
  return(kl)
}

posterior_simulation <- function(model){
  check_arg_type(model, "topicdict")
  allK <- model$extra_k + length(model$dict)
  V <- length(model$vocab)
  N = length(model$W)
  doc_lens <- sapply(model$W, length)

  if(model$extra_k > 0){
    tnames <- c(names(model$seeds), paste0("T_", 1:model$extra_k))
  }else{
    tnames <- c(names(model$seeds))
  }

  posterior_z <- function(zvec){
    tt <- table(factor(zvec, levels = 1:allK - 1))
    (tt + model$alpha) / (sum(tt) + sum(model$alpha)) # posterior mean
  }
  theta <- do.call(rbind, lapply(model$Z, posterior_z))
  rownames(theta) <- basename(model$files)
  colnames(theta) <- tnames # label seeded topics

  tZW <- Reduce(`+`,
                 mapply(function(z, w){ table(factor(z, levels = 1:allK - 1),
                                              factor(w, levels = 1:V - 1)) },
                        model$Z, model$W, SIMPLIFY = FALSE))
  word_counts <- colSums(tZW)

  colnames(tZW) <- model$vocab
  topic_counts <- rowSums(tZW)
  tZW <- tZW / topic_counts
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
    colnames(modelfit) <- c("Iteration", "Log Likelihood", "Perplexity")
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
  names(dict) <- names(model$seeds)

  ll <- list(seed_K = length(model$dict), extra_K = model$extra_k,
             V = V, N = N, Z=model$Z
             theta = theta, beta = as.matrix(as.data.frame.matrix(tZW)),
             topic_counts = topic_counts, word_counts = word_counts,
             doc_lens = doc_lens, vocab = model$vocab,
             dict = dict,
             alpha=res_alpha, modelfit=modelfit, p=p_estimated)
  class(ll) <- c("topicdict_posterior", class(ll))
  ll
}