data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
keyATM_docs <- keyATM_read(bills_dfm)
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index

test_that("Reading documents from quanteda dfm", {
  expect_identical(keyATM_docs$W_raw[[1]][3], "one")
  expect_identical(keyATM_docs$W_raw[[10]][10], "congress")
  expect_identical(keyATM_docs$W_raw[[140]][100], "number")
})


test_that("Visualizing keywords", {
  p <- visualize_keywords(keyATM_docs, bills_keywords)
  expect_s3_class(p, "keyATM_viz")
  expect_identical(values_fig(p)$WordCount[2], 318L)
  skip_on_cran()
  expect_message(
    save_fig(p, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )
})


test_that("Parallel initialization", {
  skip_on_cran()
  skip_on_travis()
  skip_on_os("linux")
  out <- keyATM(
    docs = keyATM_docs, # text input
    no_keyword_topics = 3, # number of regular topics
    keywords = bills_keywords, # keywords
    model = "base", # select the model
    options = list(seed = 250, iterations = 0, parallel_init = TRUE)
  )
  expect_identical(out$Z[[1]][3], 1L)
  expect_identical(out$Z[[140]][15], 3L)
})


# Time index
test_that("keyATM Dynamic: Initialization (correct time index)", {
  expect_message(
    out <- keyATM(
      docs = keyATM_docs, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "dynamic",
      model_settings = list(
        time_index = bills_time_index - 100,
        num_states = 5
      ),
      options = list(seed = 250, iterations = 0)
    )
  )
})


test_that("keyATM Dynamic: Initialization (wrong time index)", {
  bills_time_index_wrong <- bills_time_index
  bills_time_index_wrong[11] <- 103
  expect_error(
    out <- keyATM(
      docs = keyATM_docs, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "dynamic",
      model_settings = list(
        time_index = bills_time_index_wrong - 100,
        num_states = 5
      ),
      options = list(seed = 250, iterations = 0)
    )
  )
})


# Keep document names
test_that("Keep document names", {
  docs <- keyATM_read(bills_dfm, keep_docnames = TRUE)
  expect_identical(docs$docnames[10], "101th-congress_senate-bill_1726")

  out <- keyATM(
    docs = docs, # text input
    no_keyword_topics = 3, # number of regular topics
    keywords = bills_keywords, # keywords
    model = "base",
    options = list(seed = 250, iterations = 3)
  )

  expect_identical(row.names(out$theta)[10], "101th-congress_senate-bill_1726")
})


# Documents with 0 length

bills_dfm_0 <- bills_dfm
bills_dfm_0[10, ] <- rep(0, ncol(bills_dfm_0))
bills_dfm_0[50, ] <- rep(0, ncol(bills_dfm_0))
class(bills_dfm_0) <- c("dfm")

test_that("Documents with length 0: base", {
  expect_warning(docs0 <- keyATM_read(bills_dfm_0))

  skip_on_cran()
  expect_warning(
    out <- keyATM(
      docs = docs0, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "base", # select the model
      options = list(seed = 250, iterations = 0)
    )
  )
  expect_identical(length(out$Z), 138L)

  expect_warning(
    out <- keyATM(
      docs = docs0, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "base", # select the model
      options = list(seed = 250, iterations = 1)
    )
  )
  expect_identical(length(out$kept_values$doc_index_used), 138L)
})


test_that("Documents with length 0: covariate", {
  expect_warning(docs0 <- keyATM_read(bills_dfm_0))

  skip_on_cran()
  expect_warning(
    out <- keyATM(
      docs = docs0, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "covariates",
      model_settings = list(
        covariates_data = bills_cov,
        standardize = "all",
        covariates_formula = ~.,
        covariates_model = "DirMulti"
      ),
      options = list(seed = 250, iterations = 0)
    )
  )

  expect_identical(length(out$Z), 138L)
  expect_identical(nrow(out$model_settings$covariates_data_use), 138L)

  expect_warning(
    out <- keyATM(
      docs = docs0, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "covariates",
      model_settings = list(
        covariates_data = bills_cov,
        standardize = "all",
        covariates_formula = ~.,
        covariates_model = "DirMulti"
      ),
      options = list(seed = 250, iterations = 1)
    )
  )
  expect_identical(length(out$kept_values$doc_index_used), 138L)
})


test_that("Documents with length 0: dynamic", {
  expect_warning(docs0 <- keyATM_read(bills_dfm_0))

  skip_on_cran()
  expect_warning(
    out <- keyATM(
      docs = docs0, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "dynamic",
      model_settings = list(
        time_index = bills_time_index - 100,
        num_states = 5
      ),
      options = list(seed = 250, iterations = 0)
    )
  )

  expect_identical(length(out$Z), 138L)
  expect_identical(length(out$model_settings$time_index), 138L)

  expect_warning(
    out <- keyATM(
      docs = docs0, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "dynamic",
      model_settings = list(
        time_index = bills_time_index - 100,
        num_states = 5
      ),
      options = list(seed = 250, iterations = 1)
    )
  )
  expect_identical(length(out$kept_values$doc_index_used), 138L)
  expect_identical(length(out$kept_values$model_settings$time_index), 138L)
})

# keyATM runs without setting seed
test_that("Initialize without setting seed", {
  skip_on_cran()
  expect_no_error(
    out <- keyATM(
      docs = keyATM_docs, # text input
      no_keyword_topics = 3, # number of regular topics
      keywords = bills_keywords, # keywords
      model = "base", # select the model
      options = list(iterations = 0)
    )
  )
})


# Check `get_doc_index()` function
test_that("`get_doc_index() function", {
  W_raw <- list(c("a", "b"), c(), c("a", "b", "c"))
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw, check = TRUE))
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw, check = FALSE))
  W_raw2 <- list(c("a", "b"), c(), c("a", "b", "c"), c())
  expect_warning(tmp <- keyATM:::get_doc_index(W_raw2, check = FALSE))
})
