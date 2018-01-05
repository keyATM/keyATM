library(topicdict)
context("Reading in texts and initializing model")

test_that("init function", {
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
  v1 <- data.frame(id = ll$W[[1]],
                   wd = ll$vocab[ll$W[[1]] + 1],
                   z = ll$Z[[1]],
                   x = ll$X[[1]], stringsAsFactors = FALSE)

  expect_equal(ll$vocab[220], "napoleon")
  expect_equal(ll$seeds$mac, c(30,2,3))
  expect_equal(ll$dict$mac, c("macavity", "mystery", "cat"))
  expect_equal(v1$wd[1:4], c("macavity's", "a", "mystery", "cat"))
  expect_equal(v1$id[1:6], 0:5)
  expect_equal(v1$z[1:6], c(0,1,0,0,1,1))
  expect_equal(v1$x[1:6], c(0,0,1,1,0,0))
})
