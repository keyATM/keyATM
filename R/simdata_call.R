#' @import repurrrsive dplyr tibble purrr
NULL

#' Create Simulation Data
#' 
#' @param saveDir data storage (set full path)
#' @param D the number of documents
#' @param K the number of topics
#' @param TotalV maximum number of unique words
#' @param alpha a parameter of dirichlet of theta, a vector or a scalar
#' @param beta_r a parameter of dirichlet of phi_regular
#' @param beta_s a parameter of dirichlet of phi_seed
#' @param gamma1 a parameter of beta (set either gamma or p)
#' @param gamma2 a parameter of beta (set either gamma or p)
#' @param p a probability vector of using seed (set either gamma or p)
#' @param lambda a parameter of poisson (average length of documents)
#' @param seeds_len the number of seeds words per each topic
#' @param prop_seeds the proportion of seed words 
#' @examples
#' create_sim_data(saveDir="/Users/Shusei/Desktop/temp/SimulationData/", D=200, K=10, TotalV=8000, alpha=0.1, beta_r=0.1, beta_s=0.1, gamma1=1, gamma2=20, lambda=300, seeds_len=5)
#' create_sim_data(saveDir="/Users/Shusei/Desktop/temp/SimulationData/", D=200, K=10, TotalV=8000, alpha=0.1, beta_r=0.1, beta_s=0.1, p=rep(0.2, 10), lambda=300, seeds_len=5)
#'
#' @export
create_sim_data <- function(saveDir, D=200, K=10, TotalV=8000, alpha=0.1, beta_r=0.1, beta_s=0.1, p=NULL, gamma1=NULL, gamma2=NULL, lambda=300, seeds_len=5, prop_seeds=NULL, rand_seed=123){

	# Check values
	if(is.null(p) && is.null(gamma1) && is.null(gamma2))
		stop("Set p, or gamma1 and gamma2.")

	if( is.null(p) && ( is.null(gamma1) | is.null(gamma2)))
		stop("Set both gamma1 and gamma2.")

	if( !is.null(p) && ( !is.null(gamma1) | !is.null(gamma2)))
		stop("You cannot set p and gamma at the same time.")

	if( is.vector(p)){ # specify the each topic's proportion, length should be equal to the number of topics
		if(length(p) != K){
			stop("p: p should be given as a vector. Length of p does not match with the number of topic. You need to specify each topic.")
		}

	if(is.vector(alpha) && length(alpha) > 1){ # specify the each topic's proportion, length should be equal to the number of topics
		if(length(alpha) != K){
			stop("alpha: alpha should be given as a vector. Length of alpha does not match with the number of topic. You need to specify each topic.")
		}
	}

		flag <- 0
		for(i in length(p):1){
			if(p[i] != 0){
				flag <- 1
			}
			if(flag == 1 && p[i]==0){
				message("p: The output looks prettier if 0 (0s) comes last")
			}
		}
	}

	if( !is.null(prop_seeds)){
		if( length(prop_seeds) != seeds_len){
			stop("Length of prop_seeds should be the same as the number of seed words (seeds_len)")
		}
	}

	# Modify path
	lastcharacter <- stringi::stri_sub(saveDir,-1,-1)
	if (lastcharacter != "/"){
		saveDir <- paste0(saveDir, "/")
	}

	dir.create(saveDir)
	
	## if you want to save Z and X as well, create folder for each
	Dirs <- CreateDir_wzx(saveDir)
	
	## Data generating process
	## parameters and values
	# D <- 200 # the number of document
	# K <- 10 # the number of topic
	# TotalV <- 8000
	V <- as.integer(TotalV/K) # the number of unique vocabulary for each topic
	# # V_each <- V/10 # the number of unique vocabulary, divided by 10
	# alpha <- 0.1 # parameter of dirichlet of theta
	# beta_r <- 0.1 # parameter of dirichlet of phi_regular
	# beta_s <- 0.1 # parameter of dirichlet of phi_seed
	# gamma1 <- 1 # shape parameter of beta
	# gamma2 <- 20 # shape parameter v beta
	# lambda <- 300 # paraemter of poisson
	# seeds_len <- 5 # number of seeds words per each topic
	
	## random
	set.seed(rand_seed)
	start_time <- proc.time()
	
	## A. generate topic word distribution
	## 1. generate regular topic word distribution
	## phiR: K \times V
	## rdirichlet: each row represent a vector
	phiR <- Gen_phiR(K, rep(beta_r, V))
	
	## 2. generate seed topic word distribution
	if (is.null(prop_seeds))
		prop_seeds <- rep(1, seeds_len)
	
	prob <- prop_seeds / sum(prop_seeds)
	phiS <- manual_phiS(K, prob)
	# phiS <- Gen_phiS(K, rep(beta_s, seeds_len))
	seed_words <- Gen_seedwords(K, seeds_len, V)
	save_seed_words(seed_words, p, saveDir)
	
	
	## 3. generate \pi which control the draw from seeds words
	if(is.null(p)){
		p_use_seed <- Gen_p(K, gamma1, gamma2)
	}else{
		if(is.vector(p)){
			p_use_seed <- p
		} else {
			# Proportion is the same for all topics
			p_use_seed <- rep(p, K)
		}
	}
	
	## B. generate document
	# 0. set the length of each document
	## nd is a matrix, D \times length
	nd <- Gen_ND(D, lambda)
	
	## 1. generate \theta for each document
	if( is.vector(alpha) && length(alpha) > 1 ){
		alpha_vec <- alpha
	} else {
		alpha_vec <- Gen_alpha(alpha, K)
	}
	theta <- Gen_theta(D, alpha_vec)
	
	## 2. generate topic and word for each word in each document
	doc_list <- vector("list", D)
	
	for (i in 1:D){
		# the length of the document i
		tmp_nd <- nd[i]
	
		# ## a list of vector to contain the results
		# tmp <- WZX_Vec(tmp_nd)
	
		## Generate W, Z, and X
		wzx_list <- Gen_wzx(i, tmp_nd, phiR, phiS, p_use_seed, theta, V, seed_words)
	
		## save W, Z, X for document i
		Save_wzx(wzx_list, i, saveDir)
	
		## Store Created Docs
		doc_list[[i]] <- wzx_list
	
	}
	
	# Get Number of Unique Vocabulary
	unique_V <- length(unique(unlist(map(doc_list, "w"))))
	
	## topic word assignment across documents
	doc_list %>% map_df(`[`, c("w", "z")) %>%
		group_by(w,z) %>%
		summarize(count=n()) %>%
		tidyr::spread(key=w,  value=count) -> TopicVoc
	
	## document_vocabulary assignment
	doc_list %>% map_df(`[`, c("w", "doc_id")) %>%
		group_by(w, doc_id) %>%
		summarize(count=n()) %>%
		tidyr::spread(key=w,  value=count) -> DocVoc
	
	
	## document_topic assignment
	doc_list %>% map_df(`[`, c("doc_id", "z")) %>%
		group_by(doc_id, z) %>%
		summarize(count=n()) %>%
		tidyr::spread(key=z,  value=count) -> DocTopic
	
	## save topic-word assignment, word assignment, and topic assignment
	
	## save topic-word assignment etc
	Save_info(TopicVoc, DocVoc, DocTopic, phiS, p, saveDir)
	## Save parameters
	Save_param(D, unique_V, K, alpha, beta_r, beta_s, gamma1, gamma2, lambda, seeds_len, saveDir)
	
	## save word frequency per topic
	Save_Index_Freq(doc_list, 20, saveDir)
	
	## Save assignments of seed words
	Save_seeds_len(doc_list, saveDir)
	
	end_time <- proc.time()
	
	time <- end_time - start_time
	print("Finished: ")
	print(time)


}


