library(quanteda)
library(magrittr)
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
keyATM_docs <- keyATM_read(bills_dfm)

keywords <- list(Government     = c("laws", "law", "executive"),
                 Congress       = c("congress", "party"),
                 Peace          = c("world", "freedom"),
                 Constitution   = c("constitution"),
                 ForeignAffairs  = c("foreign", "war"))
out <- keyATM(docs = keyATM_docs, no_keyword_topics = 5, keywords = keywords,
  model = "base", options = list(seed = 250, iterations = 100))

test_that("Semantic coherence", {
  res <- semantic_coherence(out, keyATM_docs)
  skip_on_os(c("windows", "linux")) ; skip_on_cran()
  expect_equal(as.numeric(res[3]), -35.252358, tolerance = 0.0001)
})
