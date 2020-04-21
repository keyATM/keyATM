data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
keyATM_docs <- keyATM_read(bills_dfm)

test_that("Reading documents from quanteda dfm", {
  expect_identical(keyATM_docs[[1]][3], "one")
  expect_identical(keyATM_docs[[10]][10], "congress")
  expect_identical(keyATM_docs[[140]][100], "number")
})


test_that("Visualizing keywords", {
  p <- visualize_keywords(keyATM_docs, bills_keywords)
  expect_s3_class(p, "keyATM_viz")
  skip_on_cran()
  expect_message(save_fig(p, paste0(tempdir(), "/test.pdf")), "Saving 7 x 7 in image")
})


test_that("Parallel initialization", {
  skip_on_cran() ; skip_on_travis()
  out <- keyATM(docs = keyATM_docs,  # text input
                no_keyword_topics = 3,  # number of regular topics
                keywords = bills_keywords,  # keywords
                model = "base",  # select the model
                # model_settings = list(labels = labels_use),
                options = list(seed = 250, iterations = 0, parallel_init = TRUE))
  expect_identical(out$Z[[1]][3], 1L)
  expect_identical(out$Z[[140]][15], 6L)
})
