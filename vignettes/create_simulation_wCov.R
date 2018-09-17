library(tidyverse)

# Generate Simulation Data
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

## generate regular topic word distribution
## beta_vec: a vector with length V
## topic_len: a integer K, the number of dirichlet draws you want to create
Gen_phiR <- function(topic_len, beta_vec){
	phiR <- rdirichlet(topic_len, beta_vec)
	return(phiR)
}

## generate seeds topic word distribution
Gen_seedwords <- function(K, seeds_len, V){
	seed_words <- matrix(NA, K, seeds_len)
	for(k in 1:K){
		# seed_words[k, ] <- paste0(as.character(k), "T", as.character(sample(1:V, seeds_len)))
		# seed_words[k, ] <- paste0("W", as.character(sample(1:V, seeds_len)), "t", as.character(k))
		seed_words[k, ] <- paste0("W", as.character(sample(1:V, seeds_len)))
	}
	return(seed_words)	
}

manual_phiS <- function(topic_len, prob){
	prob <- prob / sum(prob)
	tmp <- t( matrix(rep(prob, topic_len), length(prob), topic_len) )
	return(tmp)
}

Gen_phiS <- function(topic_len, beta_vec){
	phiS <- rdirichlet(topic_len, beta_vec)
	return(phiS)
}

## generate probabiliity of choosing word from seed words
Gen_p <- function(topic_len, shape1, shape2){
	p <- rbeta(topic_len, shape1, shape2)
	return(p)
}

## generate document length for each document
Gen_ND <- function(document, lambda){
	nd <- numeric(document)
	## generate document length
	nd <- rpois(document, lambda)
	return(nd)
}

## generate document-topic distribution (theta)
## alpha_vec: a vector with length K
## topic_len: a integer D, the number of dirichlet draws you want to create
Gen_alpha <- function(alpha, K){
	return (rep(alpha, K))
}

Gen_theta <- function(doc_len, Alpha){
	theta <- t(apply(Alpha, 1, rdirichlet, n=1))

	# Slightly modify
	theta <- theta + 0.0000001
	theta <- theta / rowSums(theta)
	return(theta)
}

## generate topic (z) given topic-document distribution
Gen_z <- function(theta, d, doc_len){
	z <- rcat(doc_len, theta[d,])
	return(z)
}

## generate indicator (x) whether draw from seeds or usual topic-word distribution
Gen_x <- function(topic_num, probX){
	x <- rbern(1, probX[topic_num])
	return(x)
}

## generate word (w) given topic and topic-word distribution
Gen_w <- function(phi, topic){
	tmp <- rcat(1, phi[topic,])
	return(tmp)
}

## generate word (w) given seeds, uniform 
Gen_w_seeds <- function(index, phiR, phiS, z, x, seed_words){
	topic <- z[index]
	indicator <- x[index]

	if (indicator == 0){
		# Regular words
		word <- Gen_w(phiR, topic)
		# word <- paste0("W", as.character(word), "t", as.character(topic)) # no overlap of words between topics
		word <- paste0("W", as.character(word)) # allows overlap of words between topics
	} else {
		# Seed words
		seed_index <- rcat(1, phiS[topic, ])
		word <- seed_words[topic, seed_index]
	}
	return(word)
}


## generate each word
Gen_wzx <- function(doc_id, doc_len, phiR, phiS, p, theta, V, seed_words, prob=NULL){
	z <- Gen_z(theta, doc_id, doc_len)
	x <- sapply(z, Gen_x, probX=p)
	w <- sapply(1:doc_len, Gen_w_seeds, phiR=phiR, phiS=phiS, z=z, x=x, seed_words=seed_words)
	# print(w)
	return(list(doc_id=rep(doc_id, doc_len), w=w,z=z,x=x))
}

## generate vectors to contain each document
WZX_Vec <- function(doc_len){
	tmp_1 <- numeric(doc_len)
	tmp_2 <- numeric(doc_len)
	tmp_3 <- numeric(doc_len)
	return(list(tmp_1, tmp_2, tmp_3))
}


## create directory to save word assignment, topic assignment, and indicator assignment
CreateDir_wzx <- function(saveDir){
	dir1 <- paste0(saveDir, "W")
	dir.create(dir1)
	dir2 <- paste0(saveDir, "Z")
	dir.create(dir2)
	dir3 <- paste0(saveDir, "X")
	dir.create(dir3)
	return(list(dir1, dir2, dir3))
}

## return index of frequent words within each topic
Order_Index <- function(vec, num){
	ord <- order(vec,decreasing=T)[1:num]
	return(ord) # return vector
}

## create matrix with all elements zero
CreateZeroMat <- function(row, col){
	return( matrix(0, nrow = row, ncol = col) )
}


## return the number of each frequent word within each topic given index
Index_to_Freq <- function(vec, index){
	freq <- vec[index]
	return(freq) # return vector
}

## create a ".txt" which contains the frequent words within each topic
Save_Index_Freq <- function(doc_list, show_num, saveDir){
	doc_list %>% map_df(`[`, c("w", "z")) %>%
		group_by(w,z) %>%
		summarize(count=n()) %>%
		ungroup() %>% group_by(z) %>%
		mutate(ranking=row_number(desc(count))) %>%
		filter(ranking <= get("show_num")) %>%
		arrange(z, ranking) %>% select(-ranking) %>%
		select(z, w, count) -> TopicVocFreq

		write.csv(TopicVocFreq, paste0(saveDir,"topic_word_freq.csv"), na="0", row.names=FALSE)
}


## Save

save_seed_words <- function(seed_words, p, saveDir){
	if(!is.vector(p)){
		write.csv(data.frame(seed_words), paste0(saveDir,"SeedWords.csv"), na="0", row.names=FALSE)
	}else{
		write.csv(data.frame(seed_words)[p != 0, ], paste0(saveDir,"SeedWords.csv"), na="0", row.names=FALSE)
	}
}

## content is a vector
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
}

## save topic and word assignments
Save_info <- function(topic_voc, doc_voc, doc_topic, phi, p, C, Lambda,
											saveDir){
	## save topic-word assignment
	write.csv(topic_voc, paste0(saveDir,"topic_word.csv"), na="0", row.names=FALSE)
	write.csv(doc_voc, paste0(saveDir,"doc_voc.csv"), na="0", row.names=FALSE)
	write.csv(doc_topic, paste0(saveDir,"doc_topic.csv"), na="0", row.names=FALSE)
	write.csv(C, paste0(saveDir,"doc_covariates.csv"), row.names=FALSE)
	write.csv(Lambda, paste0(saveDir,"doc_Lambda.csv"), row.names=FALSE)
	Save_phiS(phi, saveDir)
	write(p, paste0(saveDir,"p.txt"))

	## Recover Topic Proportion
	total_words <- sum(topic_voc[, 2:ncol(topic_voc)], na.rm=TRUE)
	num_topicwords <- apply(topic_voc[, 2:ncol(topic_voc)], 1, sum, na.rm=TRUE)
	prop_topicwords <- num_topicwords / total_words
	num_seedwords <- num_topicwords * p
	prop_seedwords <- num_seedwords / total_words

	word_info <- paste("Proportion of Topic in Total (observed):",
				paste(as.character(round(prop_topicwords,3)), collapse=", "),
				"\nProportion of Seed words (with respect to toal words / theoretical value calculated from p)",
				paste(as.character(round(prop_seedwords,3)), collapse=", "),
				"\nNumber of Seed words (with respect to toal words, theoretical value calculated from p)",
				paste(as.character(round(num_seedwords,3)), collapse=", "),
				"\nTotal number of words:",
				paste(as.character(total_words,3), collapse=", "),
				sep="\n")
	write(word_info, paste0(saveDir,"info_topic_seed.txt"))
}

## save phiS
Save_phiS <- function(phiS, saveDir){
	tmp <- data.frame(topic_id=1:nrow(phiS))
	tmp <- cbind(tmp, phiS)	
	write.csv(tmp, paste0(saveDir,"phiS.csv"), na="0", row.names=FALSE)
}

## Save paraemters
Save_param <- function(D, V, K, alpha, beta_r, beta_s, gamma1, gamma2, lambda, seeds_len,
											 TotalV, num_covariates, Lambda_sigma,
											 saveDir){
	txt <- paste( paste("D", D, sep = ":"),
		paste("K", K, sep = ":"),
		paste("V (used)", V, sep = ":"),
		paste("beta_r", beta_r, sep = ":"),
		paste("beta_s", beta_s, sep = ":"),
		paste("gamam1", gamma1, sep = ":"),
		paste("gamma2", gamma2, sep = ":"),
		paste("lambda", lambda, sep = ":"),
		paste("TotalV", TotalV, sep = ":"),
		paste("num_covariates", num_covariates, sep = ":"),
		paste("Lambda_sigma", lambda, sep = ":"),
		paste("seeds_len", seeds_len, sep = ":"), sep = "\n")
	write(txt, paste0(saveDir,"parameters.txt"))
}

## Save assignments of seed words
Save_seeds_len <- function(doc_list, saveDir){

	doc_list %>% map_df(`[`, c("z", "w", "x")) %>%
		filter(x==1) %>%
		select(-x) %>%
		group_by(z, w) %>%
		summarize(count=n()) %>%
		arrange(z, desc(count)) %>%
		tidyr::spread(key=w,  value=count) -> tmp

	write.csv(tmp, paste0(saveDir,"seed_assignment_truth.csv"), na="0", row.names=FALSE)
}


## write w, z, x across document each into one txt file
## content is a list each of which contains integer
Write_wxz_all <- function(content, name, saveDir){
	len <- length(content)
	txt <- paste(as.character(content[[1]]), collapse = " ")
	for (i in 2:len){
		tmp <- as.character(content[[i]])
		tmp <- paste(res_tmp, collapse=" " )
		txt <- paste(txt, tmp, sep = "\n") 	
	}
	name <- paste0(saveDir, name, ".txt")
	write(txt, name)
}

## check simualation data and results
## to calculate p from the results of seeded LDA
## Create simulation data

## find index that match the value
whichVec <- function(vector, word){
	result <- which(vector == word)
	return(result) # return vector
}

## find index that match the value for each element of lists
whichVec_list <- function(data, id){
	result <- lapply(data, whichVec, id)
	return(result) # return list contains a vector
}

## find match of three vectors
intersect3 <- function(x, y, z){
	return ( intersect( intersect(x, y), z) )
}

## find matches across two lists
Match_wz <- function(w, z){
	match <- mapply(intersect, w, z) # find the match of topic, word and X
	match <- length(unlist(match))
	return(match)
}

## find matches across three lists
Match_wzx <- function(w,z,x){
	match <- mapply(intersect3, w,z,x) # find the match of topic, word and X
	match <- length( unlist(match) )
	return(match)
}

## check observed p given word and z
Obs_p <- function(w,z,x, word, z_id){
	w_tmp <- whichVec_list(w, word)
	z_tmp <- whichVec_list(z, z_id)
	x_tmp <- whichVec_list(x, 1)

	enum <- Match_wzx(w_tmp, z_tmp, x_tmp)
	denom <- Match_wz(w_tmp, z_tmp)
	p <- enum/denom
	# p <- enum/z_tmp
	return(c(p, denom, enum ))
}

Obs_p_each <- function(w,z,x, word, z_id){
	w_tmp <- whichVec_list(w, word)
	z_tmp <- whichVec_list(z, z_id)
	x_tmp <- whichVec_list(x, 1)

	enum <- Match_wzx(w_tmp, z_tmp, x_tmp)
	# denom <- Match_wz(w_tmp, z_tmp)
	# p <- enum/denom
	# p <- enum/z_tmp
	# return(c(p, denom, enum ))
	return(enum)
}

Obs_p_all <- function(w,z,x,seeds_len,topic_num){
	p <- numeric(topic_num)

	for(i in 0:(topic_num - 1)){
		J <- seq(1, seeds_len) + i * seeds_len
		enum_vec <- numeric(seed_num)
		for (j in 1:5){
			enum_vec[j] <- Obs_p_each(res$RawWords, res$Z, res$X, J[j], i)
		# print(enum_vec[j])
		}
		sum_enum_vec <- sum(enum_vec)
		z_tmp <- whichVec_list(res$Z, i)
		z_tmp <- Reduce("+", purrr::map(z_tmp, length))

		p[i + 1] <- sum_enum_vec/z_tmp
		print(p[i + 1])
	}
	return(p)
}


create_sim_data_cov <- function(saveDir, D=200, K=10, TotalV=1500, 
														# dim = 3, 
														beta_r=0.1, beta_s=0.1, p=NULL, 
														gamma1=NULL, gamma2=NULL, lambda=300, 
														num_covariates=3, Lambda_sigma=1,
														seeds_len=5, prop_seeds=NULL, rand_seed=123){

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
	V <- TotalV #as.integer(TotalV/K) # the number of unique vocabulary for each topic
	# # V_each <- V/10 # the number of unique vocabulary, divided by 10
	# alpha <- 0.1 # parameter of dirichlet of theta
	# beta_r <- 0.1 # parameter of dirichlet of phi_regular
	# beta_s <- 0.1 # parameter of dirichlet of phi_seed
	# gamma1 <- 1 # shape parameter of beta
	# gamma2 <- 20 # shape parameter v beta
	# lambda <- 300 # paraemter of poisson
	# seeds_len <- 5 # number of seeds words per each topic

	# Prepare Covariates and alpha
	C <- MASS::mvrnorm(n=D, mu=rep(0, num_covariates), Sigma=diag(0.5, num_covariates))
	# C[, 1] <- 1  # intercept
	Lambda <- MASS::mvrnorm(n=K, mu=rep(0, num_covariates), Sigma=diag(Lambda_sigma, num_covariates))
	# Lambda[, 1] <- 1  # fix intercept // just a try
	Alpha <- exp(C %*% t(Lambda))
	
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
	theta <- Gen_theta(D, Alpha)
	
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
	unique_V <- length(unique(unlist(purrr::map(doc_list, "w"))))
	
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
	CName <- paste0("C_", 1:D)
	C_wName <- cbind(CName, data.frame(C))
	Save_info(TopicVoc, DocVoc, DocTopic, phiS, p,
						C_wName, Lambda,
						saveDir)

	## Save parameters
	Save_param(D, unique_V, K, alpha, beta_r, beta_s, gamma1, gamma2, lambda, seeds_len,
						 TotalV, num_covariates, Lambda_sigma,
						 saveDir)
	
	## save word frequency per topic
	Save_Index_Freq(doc_list, 20, saveDir)
	
	## Save assignments of seed words
	Save_seeds_len(doc_list, saveDir)
	
	end_time <- proc.time()
	
	time <- end_time - start_time
	print("Finished: ")
	print(time)

	seed_list <- list()
	for(r in 1:nrow(seed_words)){
		seed_list[[r]] <- seed_words[r, ]
	}

	# lapply(seed_list, function(x){print(paste(x, collapse=" "))})

	return(seed_list)

}

# create_sim_data_cov(saveDir="/Users/Shusei/Desktop/temp/SimulationData/SeededCov",
# 														D=800, K=4, TotalV=1000, 
# 														# dim=2, 
# 														p=rep(0.3, 4),
# 														beta_r=0.01, beta_s=0.01,  lambda=350, 
# 														num_covariates=3, Lambda_sigma=0.8,
# 														seeds_len=5, rand_seed=123)
