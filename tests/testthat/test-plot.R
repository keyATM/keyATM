# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
labels_use <- keyATM_data_bills$labels
keyATM_docs <- keyATM_read(bills_dfm)

# Base
base <- keyATM(docs = keyATM_docs,
               no_keyword_topics = 3,
               keywords = bills_keywords,
               model = "base",
               options = list(seed = 250, store_theta = TRUE, iterations = 12,
                              store_pi = 1, use_weights = 1, thinning = 1))


test_that("Plot modelfit", {
  p <- plot_modelfit(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(values_fig(p)$value[2], -9903822, tolerance = 0.001)
})


test_that("Plot alpha", {
  p <- plot_alpha(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(values_fig(p)$alpha[2], 3.553541, tolerance = 0.001)
})


test_that("Plot pi", {
  p <- plot_pi(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(values_fig(p)$uq[2], 0.03964407, tolerance = 0.001)
})

# Dynamic
dyn <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "dynamic",
              model_settings = list(time_index = bills_time_index - 100,
                                    num_states = 5),
              options = list(seed = 250, verbose = FALSE, store_theta = TRUE, iterations = 20, thinning = 1))

test_that("Time series: with intervals", {
  p <- plot_timetrend(dyn) ; expect_s3_class(p, "keyATM_fig")
  p <- plot_timetrend(dyn, time_index_label = bills_time_index) ; expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(values_fig(p)$Lower[2], 0.04844414, tolerance = 0.001)
})


dyn <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "dynamic",
              model_settings = list(time_index = bills_time_index - 100,
                                    num_states = 5),
              options = list(seed = 250, verbose = FALSE, store_theta = FALSE, iterations = 20, thinning = 1))

test_that("Time series: without intervals", {
  p <- plot_timetrend(dyn) ; expect_s3_class(p, "keyATM_fig")
  p <- plot_timetrend(dyn, time_index_label = bills_time_index) ; expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(values_fig(p)$Proportion[2], 0.05016963, tolerance = 0.001)
})
