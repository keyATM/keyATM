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

# Base
base <- keyATM(docs = keyATM_docs,
               no_keyword_topics = 3,
               keywords = bills_keywords,
               model = "base",
               options = list(seed = 250, store_theta = TRUE, iterations = 30,
                              store_pi = 1, use_weights = 1))

test_that("keyATM base", {
  expect_s3_class(plot_alpha(base, start = 10), "keyATM_fig")
  expect_s3_class(plot_pi(base, method = "eti"), "keyATM_fig")

  skip_on_os("linux") ; skip_on_cran()
  expect_equal(base$model_fit$Perplexity[3], 1861.29, tolerance = 0.00001)
  expect_equal(top_words(base)[1, 1], "education [\U2713]")
  expect_equal(top_words(base)[3, 1], "educational")
  expect_equal(base$pi$Proportion[3], 6.403216, tolerance = 0.00001)
})


# Only one keyword
out <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 3,
              keywords = bills_keywords[1],
              model = "base",
              options = list(seed = 250, iterations = 5))
test_that("keyATM onle one keyword topic", {
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(out$model_fit$Perplexity[2], 3064.172, tolerance = 0.0001)
})


# Overlapping keywords
bills_keywords2 <- list(
                        Education = c("education", "child", "student"),	
                              Law = c("court", "law", "attorney"),	
                           Health = c("public", "health", "student"),	
                             Drug = c("drug", "court")	
                       )

out <- keyATM(docs = keyATM_docs, 
              no_keyword_topics = 3,
              keywords = bills_keywords2,
              model = "base",
              options = list(seed = 250, store_theta = TRUE, iterations = 12,
                             thinning = 2, use_weights = 1))

test_that("keyATM overlapping keywords", {
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(out$model_fit$Perplexity[2], 2283.335, tolerance = 0.0001)
  expect_equal(top_words(out)[1, 1], "education [\U2713]")
  expect_equal(top_words(out)[2, 5], "commission")
  expect_equal(out$pi$Proportion[2], 4.750078, tolerance = 0.00001)
})

# Same keywords in multiple topics
bills_keywords_multiple <- bills_keywords
bills_keywords_multiple$New <- c("public", "drug", "health")
bills_keywords_multiple$New2 <- c("law", "public", "health")
bills_keywords_multiple$Education <- c("education", "child", "student", "law")

out <- keyATM(docs = keyATM_docs,
              no_keyword_topics = 0,
              keywords = bills_keywords_multiple,
              model = "base",
              options = list(seed = 250, store_theta = TRUE, iterations = 10))

test_that("keyATM same keywords in multiple topics", {
  skip_on_os("linux") ; skip_on_cran()
  expect_equal(out$model_fit$Perplexity[2], 2246.162, tolerance = 0.0001)
  expect_equal(top_words(out)[1, 1], "education [\U2713]")
  expect_equal(top_words(out)[2, 5], "follow")
  expect_equal(top_words(out)[9, 6], "law [\U2713]")
})
