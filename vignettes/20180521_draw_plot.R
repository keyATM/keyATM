## code to check results and create plots

## simulation studies 1: estimate with 5 seed topics, 5 keywords each, pick up the correct, strongest keywords
## simulation studies 2:estimate with 5 seed topics, {3,5,7} keywords each, pick up the correct, {strong, middle weak} keywords


## load packages
library(quanteda)
library(tibble)
library(ggplot2)
library(dplyr)
library(topicmodels)
library(topicdict)
library(tidytext)
library(tm)
library(gtools)
library(pforeach)
library(magrittr)
library(rlist)
library(Matrix)
library(readr)
library(purrr)


## ggplot themes
ng1 <- theme(panel.background = element_rect(fill = "white", colour = "white"), 
	panel.grid.major = element_line(colour = NA), 
	axis.line = element_line(size = 1, colour="black"), 
	axis.ticks=element_line(color="black"), 
	axis.text=element_text(color="black",size=10), 
	axis.title=element_text(color="black",size=15), 
	panel.grid.minor = element_line(colour = NA), 
	legend.text = element_text(size=12), 
	legend.key = element_rect(fill = "white"), 
	legend.title = element_blank(),
	legend.key.size = unit(1, "cm"))


## set working directory
setwd("/Users/tomoyasasaki/Documents/Projects/UN_bill/tmp/sLDA/20180521")
source("../sLDA_sim_function.R")


##################################
## read dictionary file to match with LDA results
##################################

dic <- read_csv("/Users/tomoyasasaki/Documents/Projects/UN_bill/tmp/simData/simData_k20_d2000_alpha1/dict.csv", col_names = F)
colnames(dic) <- c("word", "id")


##################################
## read LDA results (sim1 and sim2 both)
##################################
## results from LDA estimated used python
LDAdir <- "/Users/tomoyasasaki/Documents/Projects/UN_bill/tmp/simRes/Sim_LDA_K5/"
K5 <- list.files(LDAdir, pattern = "*beta.txt", recursive = T, full.names = TRUE)

K5.list <- list()
for (i in 1:length(K5)){
	K5.list[[i]] <- read_table(K5[i], col_names = F)
}

K5.list2 <- lapply(K5.list, function(x, y){cbind(y,x)}, dic)
K5.list2 <- lapply(K5.list2, function(x){x[,-2] })

K5.list2 <- lapply(K5.list2, function(x){t(x)})
K5.list2 <- lapply(K5.list2, function(x){colnames(x) = x[1, ]; x[-1, ] })
# K5.list3 <- lapply(K5.list2, function(x) {sapply(x, as.numeric)} )
K5.list3 <- list()
for (i in 1:length(K5.list2)){
	d <- matrix(nrow = nrow(K5.list2[[i]]), ncol = ncol(K5.list2[[i]]))
	for (j in 1:nrow(K5.list2[[i]])){
		d[j, ] <- as.numeric(K5.list2[[i]][j,])
	}
	K5.list3[[i]] <- d
}

####################################################################
####################################################################
################################## simulation 1
####################################################################
####################################################################


##################################
## read sLDA results sim 1
##################################

slda_sim1 <- list.files("/Users/tomoyasasaki/Documents/Projects/UN_bill/tmp/simRes/Sim_sLDA_K5_sim1/", full.names = TRUE)

files_sim1 <- list()
for (i in 1:length(slda_sim1)){
  # name1 <- paste0("../20180508/res/0508sim1_", i, ".Rdata")
  name1 <- slda_sim1[i]
  load(name1)
  files_sim1[[i]] <- post
}

mats <- list()
prm <- proc.time()
for (i in 1:length(files_sim1)){
	mats[[i]] <- topic_term_mat_slda(files_sim1[[i]])
	cat( proc.time() - prm, "\n"  )
}
proc.time() - prm

##################################
## read true words
##################################

true_word_sim1 <- read_true_word("/Users/tomoyasasaki/Documents/Projects/UN_bill/tmp/simData/simData_k20_d2000_alpha1/W/")
unique_term <- colnames(mats[[1]])
true_mat <- topic_term_mat_true(true_word_sim1, unique_term, 40)



##################################
## calculate KL divergence sLDA sim 1
##################################
seedK <- 5
KL_res <- data.frame(matrix(0, nrow = 10 * seedK, ncol = 3))
for (i in 0:(length(files_sim1) - 1) ){
	for (j in 1:seedK){
		KL_res[i * seedK + j, 1] <- paste0("Model", i)
		KL_res[i * seedK + j, 2] <- paste0("Topic", j)
		cat(sprintf("model %d topic %d \n", i, j))
		KL_res[i * seedK + j, 3] <- cal_KL(mats[[i+1]][j,], true_mat[j,])
	}
}

KL_res$X1 <- as.factor(KL_res$X1)

topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5")

# topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")
KL_res$X2 <- factor(KL_res$X2, level = topic_level)

model_level <- c("K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K", "K + 7*K","K + 8*K", "K + 9*K", "K + 10*K")

## BE CAREFUL with the order!!!
KL_res$X4 <- ifelse(KL_res$X1 == "Model0", model_level[1],
				ifelse(KL_res$X1 == "Model1", model_level[2],
					ifelse(KL_res$X1 == "Model2", model_level[3],
						ifelse(KL_res$X1 == "Model3",model_level[4] ,
							ifelse(KL_res$X1 == "Model4",model_level[5],
								ifelse(KL_res$X1 == "Model5",model_level[6],
									ifelse(KL_res$X1 == "Model6",model_level[7],
										ifelse(KL_res$X1 == "Model7", model_level[8],
											ifelse(KL_res$X1 == "Model8", model_level[9], model_level[10])))))))))

KL_res$X4 <- factor(KL_res$X4, level = model_level)
if ("plyr" %in% (.packages())){detach(package:plyr)}
# KL_res %>% group_by(X1) %>%
# 			summarise(mean = mean(X3)) %>% KL_res2

# select_model <- c("K + 1", "K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K")
select_model <- c("K + K", "K + 2*K", "K + 3*K", "K + 8*K", "K + 9*K", "K + 10*K")
KL_res %>% filter(X4 %in% select_model) -> KL_res_plot
			

KL_res_plot %>% select(c(X2, X3, X4)) -> KL_res_plot

# KL_res_plot$size <- factor(c(rep(1, 40), rep(2, 10)))
g <- ggplot(KL_res_plot, aes(x = X2, y = X3, group = X4)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X4, color = X4)) + 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = select_model, 
			values = c(.4, .8, 1.2, 0.6, 1.2, 1.8, 2.4)) +
		scale_color_manual(name = "# of topics", label = select_model, 
			values=c("black", "black", "black", "red", "red" , "red")) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") + 
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_slda.pdf", g,width = 8, height = 4)


##################################
## calculate KL divergence LDA sim1
##################################
KL_res_lda <- data.frame(matrix(0, nrow = 10 * seedK, ncol = 3))
for (i in 0:(length(files_sim1) - 1)){
	for (j in 1:seedK){
		KL_res_lda[i * seedK + j, 1] <- paste0("Model", i)
		KL_res_lda[i * seedK + j, 2] <- paste0("Topic", j)
		cat(sprintf("model %d topic %d \n", i, j))

		## calculate the topic which gives you the least KL divergence
		vec <- numeric(nrow(K5.list3[[i+1]]))
		for (l in 1:nrow(K5.list3[[i+1]]) ) {
			 vec[l] <- cal_KL(K5.list3[[i+1]][l,], true_mat[j,])
		}
		KL_res_lda[i * seedK + j, 3] <- min(vec)
	}
}

KL_res_lda$X1 <- as.factor(KL_res$X1)

topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")
KL_res_lda$X2 <- factor(KL_res_lda$X2, level = topic_level)

model_level <- c("K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K", "K + 7*K","K + 8*K", "K + 9*K", "K + 10*K")

KL_res_lda$X4 <- ifelse(KL_res_lda$X1 == "Model0", model_level[1],
				ifelse(KL_res_lda$X1 == "Model1", model_level[2],
					ifelse(KL_res_lda$X1 == "Model2", model_level[3],
						ifelse(KL_res_lda$X1 == "Model3",model_level[4] ,
							ifelse(KL_res_lda$X1 == "Model4",model_level[5],
								ifelse(KL_res_lda$X1 == "Model5",model_level[6],
									ifelse(KL_res_lda$X1 == "Model6",model_level[7],
										ifelse(KL_res_lda$X1 == "Model7", model_level[8],
											ifelse(KL_res_lda$X1 == "Model8", model_level[9], model_level[10])))))))))

KL_res_lda$X4 <- factor(KL_res_lda$X4, level = model_level)
if ("plyr" %in% (.packages())){detach(package:plyr)}
# KL_res %>% group_by(X1) %>%
# 			summarise(mean = mean(X3)) %>% KL_res2

select_model <- c("K + K", "K + 2*K", "K + 3*K", "K + 8*K", "K + 9*K", "K + 10*K")
KL_res_lda %>% filter(X4 %in% select_model) -> KL_res_lda_plot
			

KL_res_lda_plot %>% select(c(X2, X3, X4)) -> KL_res_lda_plot

# KL_res_lda_plot$size <- factor(c(rep(1, 40), rep(2, 10)))
g2 <- ggplot(KL_res_lda_plot, aes(x = X2, y = X3, group = X4)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X4, color = X4)) + 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = select_model, 
			values = c(.4, .8, 1.2, .6, 1.2, 2)) +
		scale_color_manual(name = "# of topics", label = select_model, 
			values=c("black", "black", "black", "blue", "blue", "blue")) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") +
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_lda.pdf", g2, width = 8, height = 4)

## Bind two data frame
binded <- rbind( cbind(KL_res_plot, type = "sLDA"),  cbind( KL_res_lda_plot, type="LDA"))

binded %>% filter(X4 %in%  c("K + 8*K", "K + 9*K", "K + 10*K")) -> binded_plot

binded_plot$X5 <- factor(paste(binded_plot$X4, binded_plot$type), 
					level = c("K + 8*K sLDA", "K + 9*K sLDA", "K + 10*K sLDA", "K + 8*K LDA", "K + 9*K LDA", "K + 10*K LDA"))

g3 <- ggplot(binded_plot, aes(x = X2, y = X3, group = X5)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X5, color = X5) )+ 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = c("K + 8*K sLDA", "K + 9*K sLDA", "K + 10*K sLDA", "K + 8*K LDA", "K + 9*K LDA", "K + 10*K LDA"),
			values = c(.6, 1.2, 2, .6, 1.2, 2)) +
		scale_color_manual(name = "# of topics", label = c("K + 8*K sLDA", "K + 9*K sLDA", "K + 10*K sLDA", "K + 8*K LDA", "K + 9*K LDA", "K + 10*K LDA"), 
			values=c("red", "red", "red", "blue", "blue", "blue")) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") +
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_both.pdf", g3, width = 8, height = 4)



#######################################
## Calculate KL divergence by removing seed words
## simulation 1
#######################################
## dictionary used in simulation 1
sim1_key <- files_sim1[[1]]$dict
sim1_key <- lapply(sim1_key, function(x){gsub("w", "W", x)} )

##################################
## calculate KL divergence sLDA without keywords
##################################
seedK <- 5
KL_res <- data.frame(matrix(0, nrow = 10 * seedK, ncol = 3))
for (i in 0:(length(files_sim1) - 1) ){
	for (j in 1:seedK){
		KL_res[i * seedK + j, 1] <- paste0("Model", i)
		KL_res[i * seedK + j, 2] <- paste0("Topic", j)
		cat(sprintf("model %d topic %d \n", i, j))

		mats2 <- mats[[i+1]][j,][-which(names(mats[[i+1]][j,]) %in% sim1_key[[j]])]
		true_mat2 <- true_mat[j,][-which(names(true_mat[j,]) %in% sim1_key[[j]])]

		KL_res[i * seedK + j, 3] <- cal_KL(mats2, true_mat2)
	}
}



KL_res$X1 <- as.factor(KL_res$X1)

topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5")

# topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")
KL_res$X2 <- factor(KL_res$X2, level = topic_level)

model_level <- c("K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K", "K + 7*K","K + 8*K", "K + 9*K", "K + 10*K")

## BE CAREFUL with the order!!!
KL_res$X4 <- ifelse(KL_res$X1 == "Model0", model_level[1],
				ifelse(KL_res$X1 == "Model1", model_level[2],
					ifelse(KL_res$X1 == "Model2", model_level[3],
						ifelse(KL_res$X1 == "Model3",model_level[4] ,
							ifelse(KL_res$X1 == "Model4",model_level[5],
								ifelse(KL_res$X1 == "Model5",model_level[6],
									ifelse(KL_res$X1 == "Model6",model_level[7],
										ifelse(KL_res$X1 == "Model7", model_level[8],
											ifelse(KL_res$X1 == "Model8", model_level[9], model_level[10])))))))))

KL_res$X4 <- factor(KL_res$X4, level = model_level)
if ("plyr" %in% (.packages())){detach(package:plyr)}
# KL_res %>% group_by(X1) %>%
# 			summarise(mean = mean(X3)) %>% KL_res2

# select_model <- c("K + 1", "K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K")
select_model <- c("K + K", "K + 2*K", "K + 3*K", "K + 8*K", "K + 9*K", "K + 10*K")
KL_res %>% filter(X4 %in% select_model) -> KL_res_plot
			

KL_res_plot %>% select(c(X2, X3, X4)) -> KL_res_plot

# KL_res_plot$size <- factor(c(rep(1, 40), rep(2, 10)))
g <- ggplot(KL_res_plot, aes(x = X2, y = X3, group = X4)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X4, color = X4)) + 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = select_model, 
			values = c(.4, .8, 1.2, 0.6, 1.2, 1.8, 2.4)) +
		scale_color_manual(name = "# of topics", label = select_model, 
			values=c("black", "black", "black", "red", "red" , "red", "red")) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") + 
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_slda_wo_key.pdf", g, width = 8, height = 4)


##################################
## calculate KL divergence LDA without keywords
##################################
K5.list3 <- lapply(K5.list3, function(x){data.frame(x)})
K5.list3 <- lapply(K5.list3, function(x){colnames(x) <- gsub("X", "W", colnames(x)); x})

KL_res_lda <- data.frame(matrix(0, nrow = 10 * seedK, ncol = 3))
for (i in 0:(length(files_sim1) - 1)){
	for (j in 1:seedK){
		KL_res_lda[i * seedK + j, 1] <- paste0("Model", i)
		KL_res_lda[i * seedK + j, 2] <- paste0("Topic", j)
		cat(sprintf("model %d topic %d \n", i, j))
		K5.list3_2 <- K5.list3[[i+1]][-which(names(K5.list3[[i+1]]) %in% sim1_key[[j]])]
		true_mat2 <- true_mat[j,][-which(names(true_mat[j,]) %in% sim1_key[[j]])]

		## calculate the topic which gives you the least KL divergence
		vec <- numeric(nrow(K5.list3_2))
		K5.list3_2mat <- as.matrix(K5.list3_2)
		for (l in 1:nrow(K5.list3_2) ) {
			 vec[l] <- cal_KL(K5.list3_2mat[l,], true_mat2)
		}
		KL_res_lda[i * seedK + j, 3] <- min(vec)
	}
}

KL_res_lda$X1 <- as.factor(KL_res$X1)

topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")
KL_res_lda$X2 <- factor(KL_res_lda$X2, level = topic_level)

model_level <- c("K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K", "K + 7*K","K + 8*K", "K + 9*K", "K + 10*K")

KL_res_lda$X4 <- ifelse(KL_res_lda$X1 == "Model0", model_level[1],
				ifelse(KL_res_lda$X1 == "Model1", model_level[2],
					ifelse(KL_res_lda$X1 == "Model2", model_level[3],
						ifelse(KL_res_lda$X1 == "Model3",model_level[4] ,
							ifelse(KL_res_lda$X1 == "Model4",model_level[5],
								ifelse(KL_res_lda$X1 == "Model5",model_level[6],
									ifelse(KL_res_lda$X1 == "Model6",model_level[7],
										ifelse(KL_res_lda$X1 == "Model7", model_level[8],
											ifelse(KL_res_lda$X1 == "Model8", model_level[9], model_level[10])))))))))

KL_res_lda$X4 <- factor(KL_res_lda$X4, level = model_level)
if ("plyr" %in% (.packages())){detach(package:plyr)}
# KL_res %>% group_by(X1) %>%
# 			summarise(mean = mean(X3)) %>% KL_res2

select_model <- c("K + K", "K + 2*K", "K + 3*K", "K + 8*K", "K + 9*K", "K + 10*K")
KL_res_lda %>% filter(X4 %in% select_model) -> KL_res_lda_plot
			

KL_res_lda_plot %>% select(c(X2, X3, X4)) -> KL_res_lda_plot

# KL_res_lda_plot$size <- factor(c(rep(1, 40), rep(2, 10)))
g2 <- ggplot(KL_res_lda_plot, aes(x = X2, y = X3, group = X4)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X4, color = X4)) + 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = select_model, 
			values = c(.4, .8, 1.2, .6, 1.2, 2)) +
		scale_color_manual(name = "# of topics", label = select_model, 
			values=c("black", "black", "black", "blue", "blue", "blue")) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") +
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_lda_wo_key.pdf", g2, width = 8, height = 4)

## Bind two data frame
binded <- rbind( cbind(KL_res_plot, type = "sLDA"),  cbind( KL_res_lda_plot, type="LDA"))

binded %>% filter(X4 %in%  c("K + 8*K", "K + 9*K", "K + 10*K")) -> binded_plot

binded_plot$X5 <- factor(paste(binded_plot$X4, binded_plot$type), 
					level = c("K + 8*K sLDA", "K + 9*K sLDA", "K + 10*K sLDA", "K + 8*K LDA", "K + 9*K LDA", "K + 10*K LDA"))

g3 <- ggplot(binded_plot, aes(x = X2, y = X3, group = X5)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X5, color = X5) )+ 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = c("K + 8*K sLDA", "K + 9*K sLDA", "K + 10*K sLDA", "K + 8*K LDA", "K + 9*K LDA", "K + 10*K LDA"),
			values = c(.6, 1.2, 2, .6, 1.2, 2)) +
		scale_color_manual(name = "# of topics", label = c("K + 8*K sLDA", "K + 9*K sLDA", "K + 10*K sLDA", "K + 8*K LDA", "K + 9*K LDA", "K + 10*K LDA"), 
			values=c("red", "red", "red", "blue", "blue", "blue")) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") +
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_both_wo_key.pdf", g3, width = 8, height = 4)




####################################################################
####################################################################
################################## simulation 2
####################################################################
####################################################################

##################################
## read sLDA results
##################################

slda_sim2 <- list.files("/Users/tomoyasasaki/Documents/Projects/UN_bill/tmp/simRes/Sim_sLDA_K5_sim2/", full.names = TRUE)

files_sim2 <- list()
for (i in 1:length(slda_sim2)){
  # name1 <- paste0("../20180508/res/0508sim1_", i, ".Rdata")
  name1 <- slda_sim2[i]
  load(name1)
  files_sim2[[i]] <- post
}

mats <- list()
prm <- proc.time()
for (i in 1:length(files_sim2)){
	mats[[i]] <- topic_term_mat_slda(files_sim2[[i]])
	cat( proc.time() - prm, "\n"  )
}
proc.time() - prm


##################################
## calculate KL divergence sLDA
##################################
seedK <- 5
KL_res <- data.frame(matrix(0, nrow = 12 * seedK, ncol = 3))
for (i in 0:(length(files_sim2) - 1) ){
	for (j in 1:seedK){
		KL_res[i * seedK + j, 1] <- paste0("Model", i)
		KL_res[i * seedK + j, 2] <- paste0("Topic", j)
		cat(sprintf("model %d topic %d \n", i, j))
		KL_res[i * seedK + j, 3] <- cal_KL(mats[[i+1]][j,], true_mat[j,])
	}
}

KL_res$X1 <- as.factor(KL_res$X1)

topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5")

# topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")
KL_res$X2 <- factor(KL_res$X2, level = topic_level)

model_level <- c("3 seeds, W", "3 seeds, M", "3 seeds, S", "5 seeds, W", "5 seeds, M", "5 seeds, S", "7 seeds, W","7 seeds, M", "7 seeds, S", "9 seeds, W", "9 seeds, M", "9 seeds, S")

## BE CAREFUL with the order!!!
KL_res$X4 <- ifelse(KL_res$X1 == "Model0", model_level[1],
				ifelse(KL_res$X1 == "Model1", model_level[2],
					ifelse(KL_res$X1 == "Model2", model_level[3],
						ifelse(KL_res$X1 == "Model3",model_level[4] ,
							ifelse(KL_res$X1 == "Model4",model_level[5],
								ifelse(KL_res$X1 == "Model5",model_level[6],
									ifelse(KL_res$X1 == "Model6",model_level[7],
										ifelse(KL_res$X1 == "Model7", model_level[8],
											ifelse(KL_res$X1 == "Model8", model_level[9], 
												ifelse(KL_res$X1 == "Model9", model_level[10],
													ifelse(KL_res$X1 == "Model10", model_level[11],model_level[12])))))))))))

KL_res$X4 <- factor(KL_res$X4, level = model_level)
if ("plyr" %in% (.packages())){detach(package:plyr)}
# KL_res %>% group_by(X1) %>%
# 			summarise(mean = mean(X3)) %>% KL_res2

# select_model <- c("K + 1", "K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K")
# select_model <- c("K + K", "K + 2*K", "K + 3*K", "K + 8*K", "K + 9*K", "K + 10*K")
# KL_res %>% filter(X4 %in% select_model) -> KL_res_plot
			

KL_res %>% select(c(X2, X3, X4)) -> KL_res
KL_res <- KL_res[c(1:45), ]

# KL_res_plot$size <- factor(c(rep(1, 40), rep(2, 10)))
g <- ggplot(KL_res, aes(x = X2, y = X3, group = X4)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X4, color = X4)) + 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = model_level, 
			values = rep(c(.8, 1.6, 2.4), 4) ) +
		scale_color_manual(name = "# of topics", label = model_level, 
			values = rep(c("black", "red", "blue" , "green"), each = 3)) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") + 
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_slda_sim2.pdf", g,width = 8, height = 4)


##################################
## calculate KL divergence sLDA without keywords
##################################
seedK <- 5
KL_res <- data.frame(matrix(0, nrow = 12 * seedK, ncol = 3))
for (i in 0:(length(files_sim2) - 1) ){
	sim2_key <- files_sim2[[i+1]]$dict
	sim2_key <- lapply(sim2_key, function(x){gsub("w", "W", x)} )
	for (j in 1:seedK){
		KL_res[i * seedK + j, 1] <- paste0("Model", i)
		KL_res[i * seedK + j, 2] <- paste0("Topic", j)
		cat(sprintf("model %d topic %d \n", i, j))

		mats2 <- mats[[i+1]][j,][-which(names(mats[[i+1]][j,]) %in% sim2_key[[j]])]
		true_mat2 <- true_mat[j,][-which(names(true_mat[j,]) %in% sim2_key[[j]])]

		KL_res[i * seedK + j, 3] <- cal_KL(mats2, true_mat2)
	}
}

KL_res$X1 <- as.factor(KL_res$X1)

topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5")

# topic_level <- c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")
KL_res$X2 <- factor(KL_res$X2, level = topic_level)

model_level <- c("3 seeds, W", "3 seeds, M", "3 seeds, S", "5 seeds, W", "5 seeds, M", "5 seeds, S", "7 seeds, W","7 seeds, M", "7 seeds, S", "9 seeds, W", "9 seeds, M", "9 seeds, S")

## BE CAREFUL with the order!!!
KL_res$X4 <- ifelse(KL_res$X1 == "Model0", model_level[1],
				ifelse(KL_res$X1 == "Model1", model_level[2],
					ifelse(KL_res$X1 == "Model2", model_level[3],
						ifelse(KL_res$X1 == "Model3",model_level[4] ,
							ifelse(KL_res$X1 == "Model4",model_level[5],
								ifelse(KL_res$X1 == "Model5",model_level[6],
									ifelse(KL_res$X1 == "Model6",model_level[7],
										ifelse(KL_res$X1 == "Model7", model_level[8],
											ifelse(KL_res$X1 == "Model8", model_level[9], 
												ifelse(KL_res$X1 == "Model9", model_level[10],
													ifelse(KL_res$X1 == "Model10", model_level[11],model_level[12])))))))))))

KL_res$X4 <- factor(KL_res$X4, level = model_level)
if ("plyr" %in% (.packages())){detach(package:plyr)}
# KL_res %>% group_by(X1) %>%
# 			summarise(mean = mean(X3)) %>% KL_res2

# select_model <- c("K + 1", "K + K", "K + 2*K", "K + 3*K", "K + 4*K", "K + 5*K", "K + 6*K")
# select_model <- c("K + K", "K + 2*K", "K + 3*K", "K + 8*K", "K + 9*K", "K + 10*K")
# KL_res %>% filter(X4 %in% select_model) -> KL_res_plot
			

KL_res %>% select(c(X2, X3, X4)) -> KL_res
KL_res <- KL_res[c(1:45), ]

# KL_res_plot$size <- factor(c(rep(1, 40), rep(2, 10)))
g <- ggplot(KL_res, aes(x = X2, y = X3, group = X4)) + 
		# geom_line(aes(group = X4, linetype = X4, size = size, color = size)) + 
		geom_line(aes(size = X4, color = X4)) + 
		# geom_point(aes(size = X4, color = X4)) + 
		scale_size_manual(name = "# of topics", label = model_level, 
			values = rep(c(.8, 1.6, 2.4), 4) ) +
		scale_color_manual(name = "# of topics", label = model_level, 
			values = rep(c("black", "red", "blue" , "green"), each = 3)) +
		# scale_linetype_manual(name = "# of topics", label = select_model,
			# values = c("dotted", "dotted", "dotted", "dotted", "solid")) +
		# scale_size_discrete(range=c(.5, 3), guide=FALSE) +
		# scale_color_manual(values=c("black", "red"), guide = FALSE) +
		# scale_color_manual(values=c("black", "black", "black", "black", "red", "black", "black"), guide = FALSE) +
		# scale_linetype_manual(values = c("K + 10*K" = "solid"), width = 2) +
		ng1 + 
		ylab("KL divergence") + 
		xlab("Topic index") + 
		scale_y_continuous(breaks=seq(0, 8, 2), limits = c(0, 10)) +
		NULL

ggsave("./kl_k5_slda_sim2_wo_key.pdf", g,width = 8, height = 4)


