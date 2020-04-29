if (compareVersion(paste0(version$major, ".", version$minor), "3.6") < 0) {
  skip("Randomization algorithm has changed from R 3.6")
}
skip_on_cran()

# Prepare corpus
library(quanteda)
library(magrittr)

if (!require(dplyr)) {
  skip("dplyr is not installed")
}

data(data_corpus_inaugural, package = "quanteda")
data_tokens <- tokens(data_corpus_inaugural, remove_numbers = TRUE, 
                      remove_punct = TRUE, remove_symbols = TRUE,
                      remove_separators = TRUE, remove_url = TRUE) %>%
                 tokens_tolower() %>%
                 tokens_remove(c(stopwords("english"), 
                               "may", "shall", "can",
                               "must", "upon", "with", "without")) %>%
                 tokens_select(min_nchar = 3)
data_dfm <- dfm(data_tokens) %>%
               dfm_trim(min_termfreq = 5, min_docfreq = 2)

# Read texts into keyATM
keyATM_docs <- keyATM_read(data_dfm)
keywords <- list(Government     = c("laws", "law", "executive"),
                 Congress       = c("congress", "party"),
                 Peace          = c("peace", "world", "freedom"),
                 Constitution   = c("constitution", "rights"),
                 ForeignAffairs  = c("foreign", "war"))

# Create covariates
vars <- docvars(data_corpus_inaugural)
vars %>%
  as_tibble() %>%
  mutate(Period = case_when(Year <= 1899 ~ "18_19c",
                            TRUE ~ "20_21c")) %>%
  mutate(Party = case_when(Party == "Democratic" ~ "Democratic",
                           Party == "Republican" ~ "Republican",
                           TRUE ~ "Other")) %>%
  select(Party, Period) -> vars_selected

# Set the base line
vars_selected %>%
  mutate(Party  = factor(Party, 
                         levels = c("Other", "Republican", "Democratic")),
         Period = factor(Period, 
                         levels = c("18_19c", "20_21c"))) -> vars_selected

out <- keyATM(docs              = keyATM_docs,
              no_keyword_topics = 5,
              keywords          = keywords,
              model             = "covariates",
              model_settings    = list(covariates_data    = vars_selected, 
                                       covariates_formula = ~ Party + Period),
              options           = list(seed = 250, iterations = 20),
              keep              = c("Z", "S")
             )


test_that("Covariates info", {
  expect_output(mat <- covariates_info(out))
  expect_type(mat, "double")
})


test_that("Covariates get", {
  mat <- covariates_get(out)
  expect_type(mat, "double")
  expect_equivalent(dim(mat), c(58L, 4L))
})


test_that("Doc Topic", {
  strata_topic <- by_strata_DocTopic(out, by_var = "Period20_21c",
                                          labels = c("18_19c", "20_21c"))
  p <- plot(strata_topic, var_name = "Period", topics = c(1,2,3,4))
  expect_s3_class(p, "keyATM_fig")
  skip_on_os("linux")
  expect_equal(summary(strata_topic)[[2]]$Point[1], 0.1343383, tolerance = 0.001)
})


test_that("Topic Word", {
  skip_on_os("linux")
  strata_tw <- by_strata_TopicWord(out, keyATM_docs, by = as.vector(vars_selected$Party))
  top <- top_words(strata_tw, n = 3) 
  expect_equivalent(top$Democratic[1, 3], "world [\U2713]")
  expect_equivalent(top$Republican[1, 5], "war [\U2713]")
})



