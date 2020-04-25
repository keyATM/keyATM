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
               options = list(seed = 250, store_theta = T, iterations = 12,
                              store_pi = 1, use_weights = 1, thinning = 1))


test_that("Plot modelfit", {
  p <- plot_modelfit(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
})


test_that("Plot alpha", {
  p <- plot_alpha(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
})


test_that("Plot pi", {
  p <- plot_pi(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
})
