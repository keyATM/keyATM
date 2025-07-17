## copy from keyATM tests

data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords

## add a videogame topic
bills_keywords$Videogame <- c("metroid", "castlevania", "balatro")

bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
keyATM_docs <- keyATM_read(bills_dfm)

test_that("default with error", {
  expect_error(suppressWarnings(
    mod <- keyATM(docs = keyATM_docs,
                  no_keyword_topics = 3,
                  keywords = bills_keywords,
                  model = "base",
                  options = list(seed = 250,
                                 store_theta = TRUE,
                                 iterations = 30,
                                 store_pi = 1,
                                 use_weights = 1))
  ))
})

test_that("drop but not compensate", {
  no_keyword_topics <- 3
  ## default compensate_empty_topics is FALSE
  expect_warning(mod <- keyATM(docs = keyATM_docs,
                             no_keyword_topics = no_keyword_topics,
                             keywords = bills_keywords,
                             model = "base",
                             options = list(seed = 250,
                                            store_theta = TRUE,
                                            iterations = 30,
                                            store_pi = 1,
                                            use_weights = 1,
                                            drop_empty_topics = TRUE))
                 )
  expect_false(any(grepl("Videogame", colnames(mod$theta))))
  expect_true(any(grepl("Drug", colnames(mod$theta))))
  expect_equal(ncol(mod$theta), length(bills_keywords) - 1 + no_keyword_topics)
})

test_that("drop and compensate", {
  no_keyword_topics <- 3
  ## default compensate_empty_topics is FALSE
  expect_warning(mod <- keyATM(docs = keyATM_docs,
                             no_keyword_topics = no_keyword_topics,
                             keywords = bills_keywords,
                             model = "base",
                             options = list(seed = 250,
                                            store_theta = TRUE,
                                            iterations = 30,
                                            store_pi = 1,
                                            use_weights = 1,
                                            drop_empty_topics = TRUE,
                                            compensate_empty_topics = TRUE))
                 )
  expect_false(any(grepl("Videogame", colnames(mod$theta))))
  expect_true(any(grepl("Drug", colnames(mod$theta))))
  expect_equal(ncol(mod$theta), length(bills_keywords) + no_keyword_topics)
  expect_true("Other_4" %in% colnames(mod$theta))
})

test_that("same thing but keyATMvb", {
  ## copy from keyATMvb tests
  no_keyword_topics <- 3
  expect_error(suppressWarnings(mod <- keyATMvb(docs = keyATM_docs,  # text input
                                                no_keyword_topics = no_keyword_topics,  # number of regular topics
                                                keywords = bills_keywords,  # keywords
                                                model = "base",  # select the model
                                                options = list(seed = 250, iterations = 3),
                                                vb_options = list(convtol = 0.05))
                                ))
  expect_warning(mod <- keyATMvb(docs = keyATM_docs,  # text input
                                 no_keyword_topics = no_keyword_topics,  # number of regular topics
                                 keywords = bills_keywords,  # keywords
                                 model = "base",  # select the model
                                 options = list(seed = 250, iterations = 3,
                                                drop_empty_topics = TRUE),
                                 vb_options = list(convtol = 0.05)))
  expect_false(any(grepl("Videogame", colnames(mod$theta))))
  expect_true(any(grepl("Drug", colnames(mod$theta))))
  expect_equal(ncol(mod$theta), length(bills_keywords) - 1 + no_keyword_topics)

  expect_warning(mod <- keyATMvb(docs = keyATM_docs,  # text input
                                 no_keyword_topics = no_keyword_topics,  # number of regular topics
                                 keywords = bills_keywords,  # keywords
                                 model = "base",  # select the model
                                 options = list(seed = 250, iterations = 3,
                                                drop_empty_topics = TRUE,
                                                compensate_empty_topics = TRUE),
                                 vb_options = list(convtol = 0.05)))
  expect_false(any(grepl("Videogame", colnames(mod$theta))))
  expect_true(any(grepl("Drug", colnames(mod$theta))))
  expect_equal(ncol(mod$theta), length(bills_keywords) + no_keyword_topics)
  expect_true("Other_4" %in% colnames(mod$theta))
})
