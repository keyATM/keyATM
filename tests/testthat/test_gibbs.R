library(topicdict)
context("Gibbs sampler")

set.seed(1234)
dict_list <- list(mac = c("macavity", "mystery", "cat"),
                  victims = c("milk", "admiralty", "trellis", "drawings"))
folder <- system.file("extdata/macavity", package = "topicdict")
docs <- list.files(folder, pattern = "macavity", full.names = TRUE)
ll <- init(docs,
           dict = quanteda::dictionary(dict_list),
           remove_numbers = TRUE,
           remove_punct = TRUE,
           remove_symbols = TRUE,
           remove_separators = TRUE)

test_that("posterior function", {
  mod <- train(ll, 2, 1, 50, iter = 0)
  expect_equal(mod, ll)
  mod <- train(ll, 2, 1, 50, iter = 1)
  expect_equal(mod$W, ll$W)
})
