if (compareVersion(paste0(version$major, ".", version$minor), "3.6") < 0) {
  skip("Randomization algorithm has changed from R 3.6")
}

# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
labels_use <- keyATM_data_bills$labels
keyATM_docs <- keyATM_read(bills_dfm)

out <- keyATMvb(docs = keyATM_docs,  # text input
                no_keyword_topics = 3,  # number of regular topics
                keywords = bills_keywords,  # keywords
                model = "base",  # select the model
                options = list(seed = 250, iterations = 3),
                vb_options = list(convtol = 0.05))

test_that("keyATMvb MCMC", {
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(out$vb_options$Perplexity_VB$value[3], 9255.628, tolerance = 0.0001)
  expect_equal(top_words(out)[1, 1], "education [\U2713]")
  expect_equal(top_words(out)[2, 5], "technology")
  expect_equal(out$pi$Proportion[2], 6.181089, tolerance = 0.001)
})

out <- keyATMvb(docs = keyATM_docs,  # text input
                no_keyword_topics = 3,  # number of regular topics
                keywords = bills_keywords,  # keywords
                model = "base",  # select the model
                options = list(seed = 250, iterations = 3),
                vb_options = list(convtol = 0.05, init = "random"))

test_that("keyATMvb Random", {
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(out$vb_options$Perplexity_VB$value[3], 9757.913, tolerance = 0.00001)
  expect_equal(top_words(out)[1, 1], "education [\U2713]")
  expect_equal(top_words(out)[2, 5], "facility")
  expect_equal(out$pi$Proportion[2], 6.102534, tolerance = 0.0001)
})
