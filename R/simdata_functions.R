#' @import repurrrsive dplyr tibble purrr
NULL


## Create simulation data
########################################################################################

# library(MCMCpack)
# library(LaplacesDemon)
# library(repurrrsive)
# library(dplyr)
# library(tibble)
# library(purrr)

## initial parameters and values
# D: the number of document
# V: the number of unique vocabulary
# K: the number of topic
# alpha: parameter of dirichlet of theta
# beta: parameter of dirichlet of phi
# gamma1: shape parameter of beta
# gamma2: shape parameter v beta
# lambda: paraemter of poisson
# seed_len: the number of seed words per topic


## functions
########################################################################################

## generate regular topic word distribution
## beta_vec: a vector with length V
## topic_len: a integer K, the number of dirichlet draws you want to create
Gen_phiR <- function(topic_len, beta_vec){
	phiR <- MCMCpack::rdirichlet(topic_len, beta_vec)
	return(phiR)
}

## generate seeds topic word distribution
Gen_seedwords <- function(K, seeds_len, V){
	seed_words <- matrix(NA, K, seeds_len)
	for(k in 1:K){
		seed_words[k, ] <- paste0(as.character(k), "S", as.character(sample(1:V, seeds_len)))
	}
	return(seed_words)	
}

manual_phiS <- function(topic_len, prob){
	prob <- prob / sum(prob)
	tmp <- t( matrix(rep(prob, topic_len), length(prob), topic_len) )
	return(tmp)
}

Gen_phiS <- function(topic_len, beta_vec){
	phiS <- MCMCpack::rdirichlet(topic_len, beta_vec)
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

Gen_theta <- function(doc_len, alpha_vec){
	theta <- MCMCpack::rdirichlet(doc_len, alpha_vec)
	return(theta)
}

## generate topic (z) given topic-document distribution
Gen_z <- function(theta, d, doc_len){
	z <- LaplacesDemon::rcat(doc_len, theta[d,])
	return(z)
}

## generate indicator (x) whether draw from seeds or usual topic-word distribution
Gen_x <- function(topic_num, probX){
	x <- LaplacesDemon::rbern(1, probX[topic_num])
	return(x)
}

## generate word (w) given topic and topic-word distribution
Gen_w <- function(phi, topic){
	tmp <- LaplacesDemon::rcat(1, phi[topic,])
	return(tmp)
}

## generate word (w) given seeds, uniform 
Gen_w_seeds <- function(index, phiR, phiS, z, x, seed_words){
	topic <- z[index]
	indicator <- x[index]

	if (indicator == 0){
		# Regular words
		word <- Gen_w(phiR, topic)
		word <- paste0(as.character(topic), "R", as.character(word)) # no overlap of words between topics
	} else {
		# Seed words
		seed_index <- LaplacesDemon::rcat(1, phiS[topic, ])
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
Save_info <- function(topic_voc, doc_voc, doc_topic, phi, p, saveDir){
	## save topic-word assignment
	write.csv(topic_voc, paste0(saveDir,"topic_word.csv"), na="0", row.names=FALSE)
	write.csv(doc_voc, paste0(saveDir,"doc_voc.csv"), na="0", row.names=FALSE)
	write.csv(doc_topic, paste0(saveDir,"doc_topic.csv"), na="0", row.names=FALSE)
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
Save_param <- function(D, V, K, alpha, beta_r, beta_s, gamma1, gamma2, lambda, seeds_len, saveDir){
	txt <- paste( paste("D", D, sep = ":"),
		paste("K", K, sep = ":"),
		paste("V (used)", V, sep = ":"),
		paste("alpha", alpha, sep = ":"),
		paste("beta_r", beta_r, sep = ":"),
		paste("beta_s", beta_s, sep = ":"),
		paste("gamam1", gamma1, sep = ":"),
		paste("gamma2", gamma2, sep = ":"),
		paste("lambda", lambda, sep = ":"),
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
########################################################################################

## function
########################################################################################


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
		z_tmp <- Reduce("+", map(z_tmp, length))

		p[i + 1] <- sum_enum_vec/z_tmp
		print(p[i + 1])
	}
	return(p)
}
