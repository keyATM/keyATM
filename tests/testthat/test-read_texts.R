data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
keyATM_docs <- keyATM_read(bills_dfm)
keyATM_docsSplit <- keyATM_read(bills_dfm, split = 0.3)

test_that("keyATM_read: split dfm", {
  expect_true("W_split" %in% names(keyATM_docsSplit))
  expect_identical(length(keyATM_docsSplit$W_split$W_raw), 140L)
  expect_identical(length(keyATM_docsSplit$W_raw), 140L)

  expect_error(keyATM_read(bills_dfm, split = -1))
  expect_error(keyATM_read(bills_dfm, split = 1))
  expect_error(keyATM_read(bills_dfm, split = 1.5))
})

test_that("keyATM_docs: check print function", {
  expect_message(print(keyATM_docs), "keyATM_docs object of 140 documents.")
  expect_message(print(keyATM_docsSplit), "keyATM_docs object of 140 documents.")
  expect_message(print(keyATM_docsSplit$W_split), "keyATM_docs object of 140 documents.")
})

test_that("keyATM_docs: check summary function", {
  expect_message(summary(keyATM_docs), "Number of unique words")
  expect_message(summary(keyATM_docsSplit), "Number of unique words")
  expect_message(summary(keyATM_docsSplit$W_split), "Number of unique words")
})

test_that("keyATM_docs: Docs names are kept", {
  keyATM_docs_names <- keyATM_read(bills_dfm, keep_docnames = TRUE)
  expect_identical(docnames(bills_dfm), keyATM_docs_names$docnames)
})


#
# Check alternative ways to load texts into keyATM
#
test_that("keyATM_read: reading from files", {
  skip_on_cran()

  docnames <- docnames(bills_dfm)
  for (i in seq_len(length(keyATM_docs$W_raw))) {
    filepath <- paste0(tempdir(), "/test_", docnames[i], ".txt")
    sentences <- paste(keyATM_docs$W_raw[[i]], collapse = " ")
    fileConn <- file(filepath)
    writeLines(sentences, fileConn)
    close(fileConn)
  }

  textfiles <- list.files(tempdir(), pattern = "test_.*\\.txt", full.names = TRUE)
  temp <- keyATM_read(textfiles)
  expect_identical(keyATM_docs, temp)
})

test_that("keyATM_read: reading from the data frame", {
  skip_on_cran()
  docs <- tibble::tibble(
    text = purrr::map_chr(keyATM_docs$W_raw, paste, collapse = " ")
  )
  temp <- keyATM_read(docs)
  expect_identical(keyATM_docs, temp)
})
