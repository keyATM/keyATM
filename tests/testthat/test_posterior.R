library(topicdict)
context("Posterior post-processing")

set.seed(1234)
dict_list <- list(mac = c("macavity", "mystery", "cat"),
                  victims = c("milk", "admiralty", "trellis", "drawings"))
folder <- system.file("extdata/macavity", package = "topicdict")
# docs <- list.files(folder, pattern = "macavity", full.names = TRUE)
ll <- topicdict_model(file.path(folder, "macavity*"),
                      dict = quanteda::dictionary(dict_list),
                      remove_numbers = TRUE,
                      remove_punct = TRUE,
                      remove_symbols = TRUE,
                      remove_separators = TRUE)

test_that("posterior function", {
  post <- posterior(ll)
  expect_named(post, c("seed_K", "extra_K", "V", "N", "theta", "beta", "topic_counts",
                       "word_counts", "doc_lens", "vocab"))
  expect_equal(dim(post$beta), c(3, 220))
  expect_equal(post$seed_K, 2)
  expect_equal(post$extra_K, 1)
  expect_equal(post$topic_counts, c(`0` = 149, `1` = 147, `2` = 143))
  expect_equal(dim(post$theta), c(7, 3))
  expect_equal(post$theta[7,] * post$doc_lens[7],
               c(`0` = 30, `1` = 26, `2` = 23))
})


