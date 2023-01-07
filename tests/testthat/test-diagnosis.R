library(quanteda)
library(magrittr)
data(data_corpus_inaugural, package = "quanteda")
data_corpus_inaugural <- head(data_corpus_inaugural, n = 58)
data_tokens <- tokens(data_corpus_inaugural, remove_numbers = TRUE,
                      remove_punct = TRUE, remove_symbols = TRUE,
                      remove_separators = TRUE, remove_url = TRUE) %>%
                 tokens_tolower() %>%
                 tokens_remove(c(stopwords("english"),
                               "may", "shall", "can",
                               "must", "upon", "with", "without")) %>%
                 tokens_select(min_nchar = 3)
data_dfm <- dfm(data_tokens) %>% dfm_trim(min_termfreq = 5, min_docfreq = 2)
keyATM_docs <- keyATM_read(data_dfm)
keywords <- list(Government     = c("laws", "law", "executive"),
                 Congress       = c("congress", "party"),
                 Peace          = c("peace", "world", "freedom"),
                 Constitution   = c("constitution", "rights"),
                 ForeignAffairs  = c("foreign", "war"))
out <- keyATM(docs = keyATM_docs, no_keyword_topics = 5, keywords = keywords,
  model = "base", options = list(seed = 250, iterations = 100))

test_that("Semantic coherence", {
  res <- semantic_coherence(out, keyATM_docs)
  skip_on_os(c("windows", "linux")) ; skip_on_cran()
  expect_equal(as.numeric(res[3]), -7.450149, tolerance = 0.0001)
})
