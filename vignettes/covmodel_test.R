# Try Covariate Model
# install.packages(c("Rcpp", "tidyverse"), dependencies=T)
setwd("/Users/Shusei/Dropbox/Study/Project/ImaiText/topicdict")

remove.packages("topicdict")
devtools::install() ; devtools::document()

library(tidyverse)
library(quanteda)

#############
# Functions #
#############
create_model <- function(docs, seed_list, extra_k, covariates=NULL){
  set.seed(225)
  names(seed_list) <- 1:length(seed_list)
  dict <- quanteda::dictionary(seed_list)
  model <- topicdict_model(docs,
               dict = dict, 
							 covariates_data = covariates,
							 covariates_formula = ~.+0, # no intercept
							 extra_k = extra_k,
               remove_numbers = FALSE, 
               remove_punct = TRUE,
               remove_symbols = TRUE,
               remove_separators = TRUE)

  return(model)
}

assign_trueZ <- function(model, Z_folder){

	files <- list.files(Z_folder, pattern = "*.txt", full.names = TRUE)
	texts <- data.frame(text = unlist(lapply(files, function(x){ paste0(readLines(x, encoding = "UTF-8"),
                                                              collapse = "\n") })),
                     stringsAsFactors = FALSE)


	for(d in 1:length(model$Z)){
		model$Z[[d]] <- as.integer(strsplit(texts[d,], "\\s")[[1]]) - 1
	}

	return(model)

}

assign_trueX <- function(model, X_folder){

	files <- list.files(X_folder, pattern = "*.txt", full.names = TRUE)
	texts <- data.frame(text = unlist(lapply(files, function(x){ paste0(readLines(x, encoding = "UTF-8"),
                                                              collapse = "\n") })),
                     stringsAsFactors = FALSE)


	for(d in 1:length(model$X)){
		model$X[[d]] <- as.integer(strsplit(texts[d,], "\\s")[[1]])
	}

	return(model)

}
###############
###############

doc_folder <- "/Users/Shusei/Desktop/temp/SimulationData/SeededCov/W/"
Z_folder <- "/Users/Shusei/Desktop/temp/SimulationData/SeededCov/Z/"
X_folder <- "/Users/Shusei/Desktop/temp/SimulationData/SeededCov/X/"
docs <- list.files(doc_folder, pattern = "*.txt", full.names = TRUE)
seed_list <- list(
									c("W54t1 W229t1 W202t1 W212t1 W201t1"),
									c("W105t2 W244t2 W124t2 W248t2 W85t2"),
									c("W17t3 W226t3 W55t3 W74t3 W163t3"),
									c("W142t4 W242t4 W27t4 W135t4 W226t4")
									)
seed_list <- lapply(seed_list, function(x){strsplit(x, " ")[[1]]}) # need it if you input keywords as  single string
covariates <- read_csv("/Users/Shusei/Desktop/temp/SimulationData/SeededCov/doc_covariates.csv")
covariates %>% arrange(CName) %>% select(-CName) -> covariates
true_Lambda <- read_csv("/Users/Shusei/Desktop/temp/SimulationData/SeededCov/doc_Lambda.csv")
# covariates <- covariates[, 2:ncol(covariates)]
set.seed(125) ; model <- create_model(docs, seed_list, extra_k=0, covariates)

# model <- assign_trueZ(model, Z_folder)
# model <- assign_trueX(model, X_folder)
res <- topicdict_train_cov(model, iter=500)
res$sampling_info

# exp(as.matrix(covariates) %*% t(as.matrix(true_Lambda))) # Alpha

res_Lambda <- res$Lambda[200:500]
hist(unlist(map(res_Lambda, `[`, 1, 1))) ; mean(unlist(map(res_Lambda, `[`, 1, 1)))
hist(unlist(map(res_Lambda, `[`, 1, 2))) ; mean(unlist(map(res_Lambda, `[`, 1, 2)))
hist(unlist(map(res_Lambda, `[`, 1, 3))) ; mean(unlist(map(res_Lambda, `[`, 1, 3)))

hist(unlist(map(res_Lambda, `[`, 2, 1))) ; mean(unlist(map(res_Lambda, `[`, 2, 1)))
hist(unlist(map(res_Lambda, `[`, 3, 1))) ; mean(unlist(map(res_Lambda, `[`, 3, 1)))

# Reduce("+", res_Lambda)/ length(res_Lambda)
reduce(res_Lambda, `+`) / length(res_Lambda)
true_Lambda


visualize_sample <- function(true_Lambda, res_Lambda, dim1, dim2){
	truth <- as.numeric(true_Lambda[dim1, dim2])
	sampled <- tibble(lambda = unlist(map(res_Lambda, `[`, dim1, dim2)))

	p <- ggplot() +
		geom_histogram(data=sampled, aes(x=lambda), bins=20) +
		geom_vline(xintercept = truth, colour="red") +
		ggtitle(paste0("Estimated Lambda: ", "Dim: [", dim1, ",", dim2, "]")) +
		theme_bw() + theme(plot.title = element_text(hjust = 0.5))
	return(p)
}

for(i in 1:3){
	for(m in 1:3){
		print(visualize_sample(true_Lambda, res_Lambda, i, m))
		Sys.sleep(1.0)
	}
}


A <- 0.5
x <- seq(-10, 10, 0.2)
x <- seq(-5, 5, 0.5)
y <- 1.0 / (1.0 + exp( -x * A))
plot(x, y)


