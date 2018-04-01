library(topicdict)
context("Gibbs sampler")

set.seed(1234)
dict_list <- list(mac = c("macavity", "mystery", "cat"),
                  victims = c("milk", "admiralty", "trellis", "drawings"))
folder <- system.file("extdata/macavity", package = "topicdict")
docs <- list.files(folder, pattern = "macavity", full.names = TRUE)
ll <- topicdict_model(docs,
                      dict = quanteda::dictionary(dict_list),
                      remove_numbers = TRUE,
                      remove_punct = TRUE,
                      remove_symbols = TRUE,
                      remove_separators = TRUE)

test_that("posterior function", {
  mod <- topicdict_train(ll, iter = 0)
  expect_equal(mod, ll)
  mod <- topicdict_train(ll, iter = 20)
  expect_equal(mod$W, ll$W)
})
