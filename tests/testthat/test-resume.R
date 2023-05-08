if (compareVersion(paste0(version$major, ".", version$minor), "3.6") < 0) {
  skip("Randomization algorithm has changed from R 3.6")
}

# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
keyATM_docs <- keyATM_read(bills_dfm)

# keyATM
test_that("keyATM base", {
  all <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "base",
    options = list(seed = 250, iterations = 19))

  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "base",
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "base",
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))
  fs::file_delete(paste0(tempdir(), "/resume.rds"))

  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$alpha_iter$alpha[5], resumed$values_iter$alpha_iter$alpha[6])
})

# Weighted LDA
test_that("weightedLDA base", {
  all <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "base",
    options = list(seed = 250, iterations = 19))

  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "base",
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "base",
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))

  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$alpha_iter$alpha[5], resumed$values_iter$alpha_iter$alpha[6])
})
