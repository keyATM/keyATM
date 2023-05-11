if (compareVersion(paste0(version$major, ".", version$minor), "3.6") < 0) {
  skip("Randomization algorithm has changed from R 3.6")
}

# Read Data
data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
bills_keywords <- keyATM_data_bills$keywords
bills_cov <- keyATM_data_bills$cov
bills_time_index <- keyATM_data_bills$time_index
keyATM_docs <- keyATM_read(bills_dfm)


#
# keyATM
#
test_that("keyATM base, resume", {
  all <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "base",
    options = list(seed = 250, iterations = 19))

  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "base",
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "base",
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$alpha_iter$alpha[5], resumed$values_iter$alpha_iter$alpha[6])
})


test_that("keyATM cov, resume", {
  all <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "covariates",
    model_settings = list(covariates_data = bills_cov, standardize = "all",
                        covariates_formula = ~., covariates_model = "DirMulti"),
    options = list(seed = 250, iterations = 19))

  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "covariates",
    model_settings = list(covariates_data = bills_cov, standardize = "all",
                        covariates_formula = ~., covariates_model = "DirMulti"),
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "covariates",
    model_settings = list(covariates_data = bills_cov, standardize = "all",
                        covariates_formula = ~., covariates_model = "DirMulti"),
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$Lambda_iter[[5]][7, 1], resumed$values_iter$Lambda_iter[[6]][7, 1])
})


test_that("keyATM dynamic, resume (without storing the transition matrix)", {
  all <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100, num_states = 5),
    options = list(
      seed = 250, iterations = 19, store_transition_matrix = 0
    ))

  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100, num_states = 5),
    options = list(
      seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds"),
      store_transition_matrix = 0
    ))
  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100, num_states = 5),
    options = list(
      seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds"),
      store_transition_matrix = 0
    ))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$alpha_iter$alpha[175], resumed$values_iter$alpha_iter$alpha[210])
})


test_that("keyATM dynamic, resume (with storing the transition matrix)", {
  all <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100, num_states = 5),
    options = list(
      seed = 250, iterations = 19, store_transition_matrix = 1
    ))

  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100, num_states = 5),
    options = list(
      seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds"),
      store_transition_matrix = 1
    ))
  resumed <- keyATM(
    docs = keyATM_docs,
    no_keyword_topics = 3,
    keywords = bills_keywords,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100, num_states = 5),
    options = list(
      seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds"),
      store_transition_matrix = 1
    ))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$alpha_iter$alpha[175], resumed$values_iter$alpha_iter$alpha[210])
})


#
# Weighted LDA
#
test_that("weightedLDA base, resume", {
  all <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "base",
    options = list(seed = 250, iterations = 19))

  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "base",
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "base",
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
  expect_equal(all$values_iter$alpha_iter$alpha[5], resumed$values_iter$alpha_iter$alpha[6])
})


test_that("weightedLDA cov, resume", {
  all <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "covariates",
    model_settings = list(covariates_data = bills_cov, standardize = "all",
                        covariates_formula = ~., covariates_model = "DirMulti"),
    options = list(seed = 250, iterations = 19))

  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "covariates",
    model_settings = list(covariates_data = bills_cov, standardize = "all",
                        covariates_formula = ~., covariates_model = "DirMulti"),
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "covariates",
    model_settings = list(covariates_data = bills_cov, standardize = "all",
                        covariates_formula = ~., covariates_model = "DirMulti"),
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
})



test_that("weightedLDA dynamic, resume", {
  all <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100,
                          num_states = 5),
    options = list(seed = 250, iterations = 19))

  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100,
                          num_states = 5),
    options = list(seed = 250, iterations = 7, resume = paste0(tempdir(), "/resume.rds")))
  resumed <- weightedLDA(
    docs = keyATM_docs,
    number_of_topics = 5,
    model = "dynamic",
    model_settings = list(time_index = bills_time_index - 100,
                          num_states = 5),
    options = list(seed = 250, iterations = 12, resume = paste0(tempdir(), "/resume.rds")))

  fs::file_delete(paste0(tempdir(), "/resume.rds"))
  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
})


#
# Using the quanteda data
#

test_that("keyATM cov, resume with quanteda", {
  skip_on_cran()

  data(data_corpus_inaugural, package = "quanteda")
  data_corpus_inaugural <- head(data_corpus_inaugural, n = 58)
  data_tokens <- tokens(data_corpus_inaugural,remove_numbers = TRUE, remove_punct = TRUE, remove_symbols = TRUE, remove_separators = TRUE, remove_url = TRUE) %>%
                  tokens_tolower() %>%
                  tokens_remove(c(stopwords("english"),
                                "may", "shall", "can",
                                "must", "upon", "with", "without")) %>%
                  tokens_select(min_nchar = 3)
  data_dfm <- dfm(data_tokens) %>%
                dfm_trim(min_termfreq = 5, min_docfreq = 2)
  keyATM_docs <- keyATM_read(texts = data_dfm)
  keywords <- list(Government     = c("laws", "law", "executive"),
                  ForeignAffairs = c("foreign", "war"))

  vars <- docvars(data_corpus_inaugural)
  library(dplyr)
  vars %>%
    as_tibble() %>%
    mutate(Period = case_when(Year <= 1899 ~ "18_19c",
                              TRUE ~ "20_21c")) %>%
    mutate(Party = case_when(Party == "Democratic" ~ "Democratic",
                            Party == "Republican" ~ "Republican",
                            TRUE ~ "Other")) %>%
    select(Party, Period) -> vars_selected

  vars_selected %>%
    mutate(Party  = factor(Party,
                          levels = c("Other", "Republican", "Democratic")),
          Period = factor(Period,
                          levels = c("18_19c", "20_21c"))) -> vars_selected


  all <-  keyATM(
    docs              = keyATM_docs,    # text input
    no_keyword_topics = 5,              # number of topics without keywords
    keywords          = keywords,       # keywords
    model             = "covariates",
    model_settings    = list(covariates_data    = vars_selected,
                            covariates_formula = ~ Party + Period),
    options           = list(seed = 250, iterations = 11)
  )

  resumed <- keyATM(
    docs              = keyATM_docs,    # text input
    no_keyword_topics = 5,              # number of topics without keywords
    keywords          = keywords,       # keywords
    model             = "covariates",
    model_settings    = list(covariates_data    = vars_selected,
                            covariates_formula = ~ Party + Period),
    options           = list(seed = 250, iterations = 3, resume = paste0(tempdir(), "/resume.rds"))
  )
  resumed <- keyATM(
    docs              = keyATM_docs,    # text input
    no_keyword_topics = 5,              # number of topics without keywords
    keywords          = keywords,       # keywords
    model             = "covariates",
    model_settings    = list(covariates_data    = vars_selected,
                            covariates_formula = ~ Party + Period),
    options           = list(seed = 250, iterations = 8,  resume = paste0(tempdir(), "/resume.rds"))
  )
  fs::file_delete(paste0(tempdir(), "/resume.rds"))

  expect_equal(all$model_fit$Perplexity[3], resumed$model_fit$Perplexity[4])
})
