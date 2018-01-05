library(topicdoct)
context("Reading in texts")

test_that("init function", {
  set.seed(1234)
  dict_list <- list(mac = c("macavity", "mystery", "cat"),
                    victims = c("milk", "admiralty", "trellis", "drawings"))
  ll <- init(list.files("inst/extdata", pattern = "macavity", full.names = TRUE),
             dict = quanteda::dictionary(dict_list),
             remove_numbers = TRUE,
             remove_punct = TRUE,
             remove_symbols = TRUE,
             remove_separators = TRUE)
  v1 <- data.frame(id = ll$W[[1]],
                   wd = ll$vocab[ll$W[[1]] + 1],
                   z = ll$Z[[1]],
                   x = ll$X[[1]], stringsAsFactors = FALSE)

  expect_equal(ll$seeds$mac, c(30,2,3))
  expect_equal(ll$dict$mac, c("macavity", "mystery", "cat"))
  expect_equal(v1$wd, c("macavity's", "a", "mystery", "cat"))
  expect_equal(v1$vocab[220], "napoleon")
  expect_equal(v1$id, 0:5)
  expect_equal(v1$z, c(0,1,0,0,1,1))
  expect_equal(v1$x, c(0,0,1,1,0,0))
})
