test_that("get_doc_index", {
  # Single
  W_raw <- list(d1 = c("a", "b", "c"), d2 = c())
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw, check = TRUE), "Index to check: 2")
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw, check = FALSE), "Index to check: 2")

  # Multiple
  W_raw <- list(d1 = c("a", "b", "c"), d2 = c(), d3 = c())
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw, check = TRUE), "Indexes to check: 2 and 3")
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw, check = FALSE), "Indexes to check: 2 and 3")
})

test_that("check_keywords", {
  keywords <- list(t1 = c("a", "b", "c"), t2 = c("d", "e"))

  unique_words <- c("a", "b", "c", "d")
  expect_warning(tmp <- keyATM:::check_keywords(unique_words, keywords, prune = TRUE), "A keyword is pruned")
  expect_error(tmp <- keyATM:::check_keywords(unique_words, keywords, prune = FALSE), "A keyword is not found")

  unique_words <- c("a", "b", "d")
  expect_warning(tmp <- keyATM:::check_keywords(unique_words, keywords, prune = TRUE), "Keywords are pruned")
  expect_error(tmp <- keyATM:::check_keywords(unique_words, keywords, prune = FALSE), "Keywords are not found")
})

test_that("get_empty_index", {
  keywords <- list(t1 = c("a", "b", "c"), t2 = c("d", "e"), t3 = c("f"))

  unique_words <- c("a", "b", "c", "d")
  expect_warning(tmp <- keyATM:::get_empty_index(unique_words, keywords))
  expect_equal(tmp, 3)
  expect_error(tmp <- keyATM:::get_empty_index(unique_words, keywords, prune = FALSE))
})
