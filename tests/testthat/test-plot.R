# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
keyATM_docs <- keyATM_read(bills_dfm)


#
# With intervals
#

# Base
base <- keyATM(
  docs = keyATM_docs,
  no_keyword_topics = 3,
  keywords = bills_keywords,
  model = "base",
  options = list(
    seed = 250,
    store_theta = TRUE,
    iterations = 12,
    store_pi = 1,
    use_weights = 1,
    thinning = 1
  )
)


test_that("Plot modelfit: the base model", {
  p <- plot_modelfit(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(
    save_fig(p, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )
  skip_on_os("linux")
  skip_on_cran()
  expect_equal(values_fig(p)$value[2], -9903822, tolerance = 0.00001)
})


test_that("Plot alpha: the base model", {
  p <- plot_alpha(base)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(
    save_fig(p, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )
  skip_on_os("linux")
  skip_on_cran()
  expect_equal(values_fig(p)$alpha[2], 3.553541, tolerance = 0.00001)
  expect_identical(values_fig(p)$Topic[3], "3_Health")

  p <- plot_alpha(base, show_topic = c(1, 4))
  expect_identical(values_fig(p)$Topic[3], "1_Education")
})


test_that("Plot pi: the base model", {
  p <- plot_pi(base, ci = 0.95, method = "eti")
  p2 <- plot_pi(base, method = "hdi")
  expect_s3_class(p, "keyATM_fig")

  skip_on_cran()
  expect_message(
    save_fig(p, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )

  skip_on_os("linux")
  skip_on_cran()
  expect_equal(
    as.numeric(values_fig(p)$Upper[2]),
    0.03964407,
    tolerance = 0.00001
  )
  expect_equal(
    as.numeric(values_fig(p2)$Lower[2]),
    0.03485987,
    tolerance = 0.00001
  )
})

test_that("Plot topic plot: the base model", {
  p <- plot_topicprop(base)
  expect_s3_class(p, "keyATM_fig")
  p <- plot_topicprop(base, order = "topicid")
  expect_s3_class(p, "keyATM_fig")
  p <- plot_topicprop(base, n = 5)
  expect_s3_class(p, "keyATM_fig")
  p <- plot_topicprop(base, show_topwords = FALSE)
  expect_s3_class(p, "keyATM_fig")
  p <- plot_topicprop(base, show_topic = 1:3, label_topic = paste0("T", 1:3))
  expect_s3_class(p, "keyATM_fig")
  p <- plot_topicprop(
    base,
    show_topic = c(1, 3, 5),
    label_topic = paste0("T", c(1, 3, 5)),
    order = "topicid"
  )
  expect_s3_class(p, "keyATM_fig")
  p <- plot_topicprop(base, show_topic = c(1, 4, 7))
  expect_s3_class(p, "keyATM_fig")

  expect_error(plot_topicprop(base, show_topic = -1))
  expect_error(plot_topicprop(base, show_topic = 5:100))
  expect_error(plot_topicprop(base, label_topic = 5:100))
})

# Dynamic
dyn <- keyATM(
  docs = keyATM_docs,
  no_keyword_topics = 3,
  keywords = bills_keywords,
  model = "dynamic",
  model_settings = list(time_index = bills_time_index - 100, num_states = 5),
  options = list(
    seed = 250,
    verbose = FALSE,
    store_theta = TRUE,
    iterations = 20,
    thinning = 1
  )
)

test_that("Plot time series: with intervals", {
  p1 <- plot_timetrend(dyn)
  expect_s3_class(p1, "keyATM_fig")
  p2 <- plot_timetrend(dyn, time_index_label = bills_time_index, method = "eti")
  p3 <- plot_timetrend(
    dyn,
    time_index_label = bills_time_index,
    method = "hdi",
    ci = 0.89
  )
  expect_s3_class(p2, "keyATM_fig")

  skip_on_cran()
  expect_message(
    save_fig(p2, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )

  skip_on_os("linux")
  skip_on_cran()
  expect_equal(values_fig(p2)$Lower[2], 0.04844414, tolerance = 0.00001)
  expect_equal(values_fig(p3)$Upper[2], 0.1310682, tolerance = 0.00001)

  expect_equal(values_fig(p3)$state_id[56], 5, tolerance = 0.00001)

  expect_equal(values_fig(p1)$time_index[36], 9, tolerance = 0.00001)
  expect_equal(values_fig(p3)$time_index[36], 109, tolerance = 0.00001)

  p <- plot_topicprop(dyn)
  expect_s3_class(p, "keyATM_fig")
})

test_that("Plot alpha: the dynamic model", {
  skip_on_cran()
  p <- plot_alpha(dyn)
  expect_equal(tail(values_fig(p), 1)$alpha[1], 1.341934, tolerance = 0.00001)
  expect_identical(values_fig(p)$Topic[10], "2_Law")

  p <- plot_alpha(dyn, show_topic = c(1, 4))
  expect_identical(values_fig(p)$Topic[3], "1_Education")
})


#
# Without intervals
#

base <- keyATM(
  docs = keyATM_docs,
  no_keyword_topics = 3,
  keywords = bills_keywords,
  model = "base",
  options = list(
    seed = 250,
    store_theta = TRUE,
    iterations = 12,
    store_pi = FALSE,
    use_weights = 1,
    thinning = 1
  )
)

test_that("Plot pi without interval", {
  p <- plot_pi(base)
  expect_s3_class(p, "keyATM_fig")

  skip_on_cran()
  expect_message(
    save_fig(p, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )

  skip_on_os("linux")
  skip_on_cran()
  expect_equal(
    as.numeric(values_fig(p)$Probability[3]),
    0.06074006,
    tolerance = 0.00001
  )
})


dyn <- keyATM(
  docs = keyATM_docs,
  no_keyword_topics = 3,
  keywords = bills_keywords,
  model = "dynamic",
  model_settings = list(time_index = bills_time_index - 100, num_states = 5),
  options = list(
    seed = 250,
    verbose = FALSE,
    store_theta = FALSE,
    iterations = 20,
    thinning = 1
  )
)

test_that("Time series: without intervals", {
  p <- plot_timetrend(dyn)
  expect_s3_class(p, "keyATM_fig")
  p <- plot_timetrend(dyn, time_index_label = bills_time_index)
  expect_s3_class(p, "keyATM_fig")
  skip_on_cran()
  expect_message(
    save_fig(p, paste0(tempdir(), "/test.pdf")),
    "Saving 7 x 7 in image"
  )
  skip_on_os("linux")
  skip_on_cran()
  expect_equal(values_fig(p)$Proportion[2], 0.05016963, tolerance = 0.00001)
})
