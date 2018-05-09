library(topicdict)
library(quanteda)
library(dplyr)
library(ggplot2)
theme_set(theme_minimal())

data("corpus_uk_platforms")
post70s_corpus <- corpus_subset(corpus_uk_platforms, date > 1970)

dictfile <- system.file("extdata/laver-garry-ajps.ykd", package = "topicdict")
# We'll use just the economics sections of this dictionary
dict <- dictionary(file = dictfile)[["Laver and Garry"]][["State in Economy"]]
dict$Neut <- NULL # not interested in this category! Just Pro / Con
# load the corpus
data("corpus_uk_platforms")

new_dict <- preprocess_dictionary(dict, post70s_corpus, min.freq = 5,
                                  remove = stopwords(), remove_numbers = TRUE,
                                  remove_punct = TRUE, remove_symbols = TRUE,
                                  remove_separators = TRUE, remove_hyphens = FALSE)

partyabbrev <- c("Con", "Lib", "Lab", "LD", "LibSDP")
docvars(post70s_corpus, "main") <- docvars(post70s_corpus, "party") %in% partyabbrev

docvars(post70s_corpus, "text_number") <- paste0("text", 1:ndoc(post70s_corpus))

plot_with_k <- function(extra){
  set.seed(1234)
  mod <- topicdict_model(post70s_corpus, dict = new_dict,
                         stopwords = stopwords(), remove_numbers = TRUE,
                         remove_punct = TRUE, remove_symbols = TRUE,
                         remove_separators = TRUE, remove_hyphens = FALSE,
                         extra_k = extra)
  mod <- topicdict_train(mod, 340)
  post <- posterior(mod)

  th <- data.frame(post$theta)
  th$doc_id <- rownames(th) # for merging
  th$modPOS <- (th[,'Con'] - th[,'Pro']) / (th[,'Con'] + th[,'Pro'])
  th$logitPOS <- log(th[,'Con']) - log(th[,'Pro'])
  th$n <- summary(post70s_corpus)$Tokens
  together <- left_join(docvars(post70s_corpus), th, by = c("text_number" = "doc_id")) %>%
    filter(party %in% partyabbrev)
  together$se <- sqrt(1 / (together$Con * together$n) +
                        1 / (together$Pro * together$n))
  together$ci.lower <- together$logitPOS - 2 * together$se
  together$ci.upper <- together$logitPOS + 2 * together$se
  # together$prop <- together$Con + together$Pro > 0.6
  p <- ggplot(together, aes(date, logitPOS, color = party)) +
    geom_line() +
    geom_point() + #aes(size = prop)) +
    geom_vline(xintercept = c(1992, 1997), alpha = 0.1, size = 5 ) +
    geom_ribbon(aes(x = date, ymin = ci.lower, ymax = ci.upper,
                    fill = party, color = NA),
                alpha = 0.25) +
    scale_color_manual(values = c("blue", "red", "orange", "orange", "orange")) +
    scale_fill_manual(values = c("blue", "red", "orange", "orange", "orange")) +
    coord_flip() +
    ggtitle(paste(extra, "extra topics"))
  ggsave(paste0("plots/plot_extra_", extra, ".pdf"), p)
}

for (i in 1:25){
  plot_with_k(i)
}

