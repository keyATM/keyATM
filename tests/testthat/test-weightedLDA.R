if (compareVersion(paste0(version$major, ".", version$minor), "3.6") < 0) {
  skip("Randomization algorithm has changed from R 3.6")
}
skip_on_cran()

# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
keyATM_docs <- keyATM_read(bills_dfm)

# weightedLDA base
base <- weightedLDA(docs = keyATM_docs,
                      number_of_topics = 5,
                      model="base",
                      options=list(seed = 100, iterations = 30))

test_that("weightedLDA base", {
  expect_error(plot_alpha(base, start = 10))
  expect_error(plot_pi(base))

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(base$model_fit$Perplexity[3], 1981.469, tolerance = 0.0001)
  expect_equal(top_words(base)[3, 1], "end")
  expect_equal(base$pi, NULL)
})


# weightedLDA cov
cov <- weightedLDA(docs = keyATM_docs,
                   number_of_topics = 3,
                   model="covariates",
                   model_settings = list(covariates_data = bills_cov, standardize = "none",
                                         covariates_formula = ~.),
                   options = list(seed = 100, iterations = 10))

test_that("weightedLDA cov", {
  expect_error(plot_alpha(cov, start = 10))
  expect_error(plot_pi(cov))

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(cov$model_fit$Perplexity[2], 2179.786, tolerance = 0.0001)
  expect_equal(top_words(cov)[3, 1], "research")
})


# weightedLDA HMM
dyn <- weightedLDA(docs = keyATM_docs,
                   number_of_topics = 3,
                   model = "dynamic",
                   model_settings = list(time_index = keyATM_data_bills$time_index - 100, 
                                         num_states = 5),
                   options = list(seed = 100, iterations = 10))

test_that("weightedLDA dynamic", {
  expect_error(plot_alpha(dyn, start = 10))
  expect_error(plot_pi(dyn))

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(dyn$model_fit$Perplexity[2], 2087.217, tolerance = 0.0001)
  expect_equal(top_words(dyn)[3, 1], "commission")
})
