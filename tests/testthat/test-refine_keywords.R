data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords


bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
keyATM_docs <- keyATM_read(bills_dfm)

test_that("refine_keywords", {
  ## add a videogame topic
  bills_keywords$Videogame <- c("metroid", "castlevania", "balatro")
  expect_true("Videogame" %in% names(bills_keywords))

  ## prune = TRUE (default)
  expect_warning(new_keywords <- refine_keywords(bills_keywords, keyATM_docs))
  expect_false("Videogame" %in% names(new_keywords))

  ## prune = TRUE (explicit)
  expect_warning(
    new_keywords_explicit <- refine_keywords(
      bills_keywords,
      keyATM_docs,
      prune = TRUE
    )
  )
  expect_equal(new_keywords, new_keywords_explicit)

  ## prune = FALSE
  expect_error(
    refine_keywords(bills_keywords, keyATM_docs, prune = FALSE),
    "not found in the documents"
  )

  ## prune = FALSE (correct keywords)
  valid_keywords <- keyATM_data_bills$keywords
  expect_silent(
    out <- refine_keywords(valid_keywords, keyATM_docs, prune = FALSE)
  )
  expect_equal(out, valid_keywords)

  ## All topics are pruned
  videogame_keywords <- list(
    "Videogame" = c("metroid", "castlevania", "balatro")
  )
  expect_error(suppressWarnings(refine_keywords(
    videogame_keywords,
    keyATM_docs
  )))
})

test_that("integration", {
  skip_on_os("linux")
  skip_on_cran()
  suppressWarnings(
    mod <- keyATM(
      docs = keyATM_docs,
      no_keyword_topics = 3,
      keywords = refine_keywords(bills_keywords, keyATM_docs),
      model = "base",
      options = list(
        seed = 250,
        store_theta = TRUE,
        iterations = 30,
        store_pi = 1,
        use_weights = 1
      )
    )
  )
  expect_false(any(grepl("Videogame", colnames(mod$theta))))
})
