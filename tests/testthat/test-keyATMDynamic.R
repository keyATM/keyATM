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


# Dynamic
dyn <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "dynamic",
              model_settings = list(time_index = bills_time_index - 100,
                                    num_states = 5),
              options = list(seed = 250, verbose = FALSE, iterations = 15, thinning = 2))

test_that("keyATM dynamic", {
  expect_s3_class(plot_alpha(dyn, start = 10), "keyATM_fig")

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(dyn$model_fit$Perplexity[3], 2171.842, tolerance = 0.00001)
  expect_equal(top_words(dyn)[1, 1], "education [\U2713]")
  expect_equal(top_words(dyn)[2, 5], "security")
  expect_equal(dyn$pi$Proportion[2], 2.960897, tolerance = 0.00001)
})


