if (compareVersion(paste0(version$major, ".", version$minor), "3.6") < 0) {
  skip("Randomization algorithm has changed from R 3.6")
}


# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
labels_use <- keyATM_data_bills$labels
keyATM_docs <- keyATM_read(bills_dfm)


# Covariate
cov <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "covariates",
              model_settings = list(covariates_data = bills_cov, standardize = "all", 
                                  covariates_formula = ~.
                                  ),
              options = list(seed = 250, store_theta = TRUE, iterations = 20,
                             store_pi = 1, thinning = 5, verbose = FALSE),
              keep = c("Z", "S")
             )

test_that("keyATM covariate", {
  expect_output(covariates_info(cov))
  expect_type(covariates_get(cov), "double")
  expect_error(covariates_info(base))

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(cov$model_fit$Perplexity[3], 1874.663, tolerance = 0.001)
  expect_equal(top_words(cov)[1, 1], "education [\U2713]")
  expect_equal(top_words(cov)[3, 3], "care")
  expect_equal(cov$pi$Proportion[2], 4.836863, tolerance = 0.00001)
})


# Heterogeneity
test_that("keyATM Heterogeneity Doc-Topic", {
  strata_topic <- by_strata_DocTopic(cov, by_var = "RepParty", labels = c("Dem", "Rep"), parallel = FALSE, posterior_mean = FALSE)

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(summary(strata_topic, method = "eti")[[2]]$Lower[2], 0.13578, tolerance = 0.00001)

  p <- plot(strata_topic, show_topic = c(1,2,3,4), by = "covariate", method = "eti")
  expect_s3_class(p, "keyATM_fig")

  expect_message(suppressWarnings(save_fig(p, paste0(tempdir(), "/test.pdf"))), "Saving 7 x 7 in image")
})

test_that("keyATM Heterogeneity Topic-Word", {
  RepParty <- as.vector(bills_cov[, "RepParty"])  # the length should be the same as the number of documents
  strata_tw <- by_strata_TopicWord(cov, keyATM_docs, by = RepParty)
  
  RepParty_chr <- ifelse(bills_cov[, "RepParty"] == 0, "Democrat", "Republican")
  strata_tw_chr <- by_strata_TopicWord(cov, keyATM_docs, RepParty_chr)
  expect_equal(top_words(strata_tw_chr, n = 3)$Republican[1, 3], "public [\U2713]")
})


# Covariate settings: Standardize
bills_cov_modified <- as.data.frame(bills_cov)
bills_cov_modified$factor <- factor(rep(1:10, nrow(bills_cov_modified) / 10))
bills_cov_modified$numerical <- 1:nrow(bills_cov_modified)

cov <- suppressWarnings(keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "covariates",
              model_settings = list(covariates_data = bills_cov, standardize = "none", 
                                  covariates_formula = NULL
                                  ),
              options = list(seed = 250, store_theta = TRUE, iterations = 5,
                             store_pi = 1, thinning = 5, verbose = FALSE)
             ))

test_that("Covariates settings: Standardize - none, no formula", {
  expect_identical(cov$kept_values$model_settings$covariates_data_use[2, 1], 0L)
  expect_identical(ncol(cov$kept_values$model_settings$covariates_data_use), 1L)
  
  skip_on_os("linux") ; skip_on_cran()
  expect_error(predict(cov, bills_cov_modified))
  expect_equal(as.numeric(suppressWarnings(predict(cov, bills_cov, transform = TRUE))[3, 3]), 0.1734721, tolerance = 0.000001)

  bills_cov_copy <- bills_cov
  bills_cov_copy[, 1] <- 1
  expect_equal(as.numeric(suppressWarnings(predict(cov, bills_cov_copy, transform = TRUE))[3, 3]), 0.1926955, tolerance = 0.000001)
  expect_equal(as.numeric(suppressWarnings(predict(cov, bills_cov_copy))[3, 3]), 0.1926955, tolerance = 0.000001)
})


cov <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "covariates",
              model_settings = list(covariates_data = bills_cov_modified, standardize = "none", 
                                  covariates_formula = ~.
                                  ),
              options = list(seed = 250, store_theta = TRUE, iterations = 5,
                             store_pi = 1, thinning = 5, verbose = FALSE)
             )

test_that("Covariates settings: Standardize - none", {
  expect_identical(cov$kept_values$model_settings$covariates_data_use[2, 2], 0)
  expect_identical(ncol(cov$kept_values$model_settings$covariates_data_use), 12L)


  skip_on_os("linux") ; skip_on_cran()
  expect_error(predict(cov, bills_cov_modified))
  expect_equal(as.numeric(suppressMessages(predict(cov, bills_cov_modified, transform = TRUE))[3, 3]), 0.008454946, tolerance = 0.000001)
})


cov <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "covariates",
              model_settings = list(covariates_data = bills_cov_modified, standardize = "non-factor", 
                                  covariates_formula = ~.
                                  ),
              options = list(seed = 250, store_theta = TRUE, iterations = 5,
                             store_pi = 1, thinning = 5, verbose = FALSE)
             )

test_that("Covariates settings: Standardize - non-factor", {
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[2, 2]), -1.257464, tolerance = 0.00001)
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[5, 12]), -1.614947, tolerance = 0.00001)
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[8, 1]), 1, tolerance = 0.00001)
  expect_identical(ncol(cov$kept_values$model_settings$covariates_data_use), 12L)

  skip_on_os("linux") ; skip_on_cran()
  expect_error(predict(cov, bills_cov_modified))
  expect_equal(as.numeric(suppressMessages(predict(cov, bills_cov_modified, transform = TRUE))[2, 3]), 0.1179603, tolerance = 0.000001)
})



cov <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords,
              model = "covariates",
              model_settings = list(covariates_data = bills_cov_modified, standardize = "all", 
                                  covariates_formula = ~.
                                  ),
              options = list(seed = 250, store_theta = TRUE, iterations = 5,
                             store_pi = 1, thinning = 5, verbose = FALSE)
             )

test_that("Covariates settings: Standardize - all", {
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[2, 2]), -1.257464, tolerance = 0.00001)
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[3, 3]), -0.3321407, tolerance = 0.00001)
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[5, 12]), -1.614947, tolerance = 0.00001)
  expect_equal(as.numeric(cov$kept_values$model_settings$covariates_data_use[8, 1]), 1, tolerance = 0.00001)
  expect_identical(ncol(cov$kept_values$model_settings$covariates_data_use), 12L)

  skip_on_os("linux") ; skip_on_cran()
  expect_error(predict(cov, bills_cov_modified))
  expect_equal(as.numeric(suppressMessages(predict(cov, bills_cov_modified, transform = TRUE))[5, 2]), 0.086891, tolerance = 0.000001)
})

