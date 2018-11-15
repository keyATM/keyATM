library(grid)
library(gridExtra)
library(tm)

multiple_simulations <- function(trueK,
                                 estimatedK,
                                 seed_only=F,
                                 seed_contamination=0, seed_p_base=0.15){

  # Run Simulation
  run_simulations(trueK, estimatedK, seed_only=seed_only, seed_contamination=seed_contamination,
									seed_p_base=seed_p_base)

  # Create Figure
  create_simulation_figure(trueK, estimatedK)
}


run_simulations <- function(trueK, estimatedK, seeds_len=6,
                            seed_only=F, seed_contamination=0, seed_p_base=0.15,
														show_n=15){
	# show_n: number of topwords to use
	
  # Create Combinations
  combinations <- expand.grid(trueK, estimatedK)
  num_combinations <- nrow(combinations)

  # Run Simulations
  for(s in 1:num_combinations){
    trueK_ <- combinations[s, 1]
    estimatedK_ <- combinations[s, 2]

    # Create Data
    set.seed(225)
    data_folder <- tempfile()
    seed_list <- create_sim_data(saveDir=paste0(data_folder, "Sim1"), 
                                 D=1000, K=trueK_, TotalV=3000, alpha=0.1, 
                                 beta_r=0.1, beta_s=0.1, 
                                 p=rep(seed_p_base, trueK_) + rnorm(trueK_, mean=0, sd=0.04), 
                                 lambda=200, seeds_len=seeds_len)
    seed_list <- lapply(seed_list, function(x){tolower(x)})

    # Seed contamination
    if(seed_contamination != 0){
      for(i in 1:seed_contamination){
        seed_list <- lapply(seed_list, function(x){
                      x[sample(1:seeds_len, 1)] <- seed_list[[sample(1:trueK_, 1)]][sample(1:seeds_len, 1)] 
                      return(x)
                      })
      }
    }

    ## Fit Seeded LDA
    extra_k_ <- estimatedK_ - length(seed_list)
    if(extra_k_ < 0){
      message("extra_k is negative, setting it to 0")
      extra_k_ <- 0
      seed_list <- seed_list[1:estimatedK_]
    }
    doc_folder <- paste0(data_folder, "Sim1", "/W")
    docs <- list.files(doc_folder, pattern = "*.txt", full.names = TRUE)
    model <- create_model(docs, seed_list, extra_k=extra_k_)
    res <- topicdict_train(model, iter = iter_num)
    post <- topicdict::posterior(res)

		# Keep seed words
    if(seed_only){
      g <- diagnosis_topic_recovery_heatmap(post, show_n, title_=F, seed_list=seed_list)
    }else{
      g <- diagnosis_topic_recovery_heatmap(post, show_n, title_=F)
    }

    # Save
    saveRDS(g, file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "SeededFig_T", trueK_, "_E", estimatedK_, ".obj"))
    saveRDS(post, file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "SeededPost_T", trueK_, "_E", estimatedK_, ".obj"))

		# Remove seed words
    if(seed_only){
      g <- diagnosis_topic_recovery_heatmap(post, show_n, title_=F, seed_list=seed_list, remove_seed_words=T)
    }else{
      g <- diagnosis_topic_recovery_heatmap(post, show_n, title_=F, remove_seed_words=T, seed_list=seed_list)
    }

    # Save
    saveRDS(g, file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "SeededFig_T_NS", trueK_, "_E", estimatedK_, ".obj"))


		## Fit Standard LDA
		# Prepare Data
		corpus <- Corpus(DirSource(doc_folder))
		strsplit_space_tokenizer <- function(x)
				unlist(strsplit(as.character(x), "[[:space:]]+"))

		dtm <- DocumentTermMatrix(corpus,
														 control = list(tokenize=strsplit_space_tokenizer, 
														 stopwords = F, tolower = F, 
														 stemming = F, wordLengths = c(1, Inf)))

		lda <- LDA(dtm, k = estimatedK_, control = list(seed = 225, iter=iter_num), method="Gibbs")

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

		num <- length(unique(res$EstTopic))
		res %>%
			group_by(EstTopic) %>%
			top_n(1, Proportion) %>%
			# mutate(forranking = as.integer(gsub("EstTopic", "", EstTopic))) %>%
			# mutate(forranking = as.integer(gsub("True", "", TrueTopic))) %>%
			arrange(TrueTopic) %>%
			ungroup() %>%
			select(EstTopic) %>% as.matrix() %>% as.vector() -> topicvec 

		num <- length(unique(res$EstTopic))
		truenum <- length(unique(res$TrueTopic))

		g <- ggplot(res, aes(EstTopic, TrueTopic)) +
					geom_tile(aes(fill=Proportion)) + 
					scale_fill_gradient(limits=c(0, 100), low="#e8e8e8", high="#0072B2") +
					scale_x_discrete(limits = (topicvec)) +
					coord_flip() +
					# scale_y_discrete(limits = paste0("True", 1:truenum)) +
					xlab("Estimated Topics") + ylab("True Topic") + theme_bw(base_size=13) 

    # Save
    saveRDS(g, file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "LDAFig_T", trueK_, "_E", estimatedK_, ".obj"))
    saveRDS(post, file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "LDAPost_T", trueK_, "_E", estimatedK_, ".obj"))

    message(paste0("Done: ", s, "/", num_combinations))
  }
}

create_simulation_figure <- function(trueK, estimatedK, title_="Simulation Results"){
  # Create Combinations
  combinations <- expand.grid(trueK, estimatedK) %>%
                    arrange(rev(Var2)) 
  num_combinations <- nrow(combinations)

	#### Seeded LDA with keywords
  # Load Data
  figures <- list()
  for(s in 1:num_combinations){
    trueK_ <- combinations[s, 1]
    estimatedK_ <- combinations[s, 2]

    figures[[s]] <- readRDS(file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "SeededFig_T", trueK_, "_E", estimatedK_, ".obj"))
  }

  ## Create Figure
  # Get Information
  g1 <- ggplotGrob(figures[[1]])
  id.legend <- grep("guide", g1$layout$name)
  legend <- g1[["grobs"]][[id.legend]]


  # Edit Figure
  edit_figure <- theme(legend.position="none",
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank())
  figures <- lapply(figures, function(x){x + edit_figure})

  # New Pictures cf. https://stackoverflow.com/a/11093069/4357279
  ga <- arrangeGrob(grobs=figures, 
               nrow = length(estimatedK),
               right = legend,
               top = textGrob(paste0(title_, ": Seeded LDA")),
               left = textGrob("Estimated Topic", rot = 90, vjust = 1),
               bottom = textGrob("True Topic", vjust = -0.1))

	grid.newpage()
  grid.draw(ga, recording=F) # Show plot

	#### Seeded LDA without keywords
  # Load Data
  figures <- list()
  for(s in 1:num_combinations){
    trueK_ <- combinations[s, 1]
    estimatedK_ <- combinations[s, 2]

    figures[[s]] <- readRDS(file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "SeededFig_T_NS", trueK_, "_E", estimatedK_, ".obj"))
  }

  ## Create Figure
  # Get Information
  g1 <- ggplotGrob(figures[[1]])
  id.legend <- grep("guide", g1$layout$name)
  legend <- g1[["grobs"]][[id.legend]]


  # Edit Figure
  edit_figure <- theme(legend.position="none",
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank())
  figures <- lapply(figures, function(x){x + edit_figure})

  # New Pictures cf. https://stackoverflow.com/a/11093069/4357279
  gb <- arrangeGrob(grobs=figures, 
               nrow = length(estimatedK),
               right = legend,
               top = textGrob(paste0(title_, ": Seeded LDA without keywords")),
               left = textGrob("Estimated Topic", rot = 90, vjust = 1),
               bottom = textGrob("True Topic", vjust = -0.1))

	grid.newpage()
  grid.draw(gb, recording=F) # Show plot

	#### LDA
  # Load Data
  figures <- list()
  for(s in 1:num_combinations){
    trueK_ <- combinations[s, 1]
    estimatedK_ <- combinations[s, 2]

    figures[[s]] <- readRDS(file = paste0("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict/vignettes/obj/", 
                            "LDAFig_T", trueK_, "_E", estimatedK_, ".obj"))
  }

  ## Create Figure
  # Get Information
  g1 <- ggplotGrob(figures[[1]])
  id.legend <- grep("guide", g1$layout$name)
  legend <- g1[["grobs"]][[id.legend]]


  # Edit Figure
  edit_figure <- theme(legend.position="none",
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank())
  figures <- lapply(figures, function(x){x + edit_figure})

  # New Pictures cf. https://stackoverflow.com/a/11093069/4357279
  gc <- arrangeGrob(grobs=figures, 
               nrow = length(estimatedK),
               right = legend,
               top = textGrob(paste0(title_, ": Standard LDA")),
               left = textGrob("Estimated Topic", rot = 90, vjust = 1),
               bottom = textGrob("True Topic", vjust = -0.1))

	grid.newpage()
  grid.draw(gc) # Show plot
}



diagnosis_topic_recovery_heatmap <- function(post, n=15, title_=F,
                        seed_list=NULL, seed_only=F, # seed_only: show only seed words
                        topicvec=c(), merge=list(), remove_seed_words=F){
  topwords <- top_terms(post, n=n)
  topwords <- data.frame(topwords)
  colnames(topwords) <- paste0("EstTopic", 1:ncol(topwords))

  topwords <- tidyr::gather(topwords, key=EstTopic, value=Word) %>%
                mutate(Word = gsub("\\s.*$", "", Word)) # remove unnecessary brackets

	if(remove_seed_words){
		# Remove keywords from the topwords
		keywords <- unlist(seed_list)
		topwords <- topwords[ apply(topwords, 2, function(x){ return(!(x %in% keywords)) })[, "Word"], ] # remove words not in seed_list
	}

  topwords %>%
    mutate(RawWord = Word) %>%
    tidyr::separate(Word,
        into=c("word_id", "TrueTopic"),
        sep="t") %>%
    mutate(TrueTopic = paste0("True", as.character(TrueTopic))) -> res_

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

  if(seed_only){
    # Use only topics with keywords
    num <- length(seed_list)
    seed_list_name <- paste0("EstTopic", 1:num)
    res_ %>%
      filter(EstTopic %in% get("seed_list_name")) -> res_
  }

  num <- length(unique(res_$EstTopic))
  if(is.null(topicvec)){
    res_ %>%
      group_by(EstTopic) %>%
      top_n(1, Proportion) %>%
      mutate(forranking = as.integer(gsub("EstTopic", "", EstTopic))) %>%
      arrange(forranking) %>%
      select(EstTopic) -> topicvec 
    topicvec <-  unique(as.integer(gsub("EstTopic", "", topicvec$EstTopic)))
  }else if(length(topicvec) != num){
    message("topicvec length does not match")
    topicvec <- 1:num
  }

  truenum <- length(unique(res_$TrueTopic))

  title <- paste0("Seeded LDA: Top ", as.character(n), " words")

  g <- ggplot(res_, aes(EstTopic, TrueTopic)) +
        geom_tile(aes(fill=Proportion)) + 
        scale_fill_gradient(limits=c(0, 100), low="#e8e8e8", high="#0072B2", name = "Proportion") +
        scale_x_discrete(limits = rev(paste0("EstTopic", topicvec))) +
        coord_flip() +
        scale_y_discrete(limits = paste0("True", 1:truenum)) +
        xlab("Estimated Topics") + ylab("True Topic") + theme_bw(base_size=13)
  
  if(title_){
    g <- g + ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5))
  }


  return(g)
}


