data(keyATM_data_bills)
bills_dfm <- keyATM_data_bills$doc_dfm
keyATM_docs <- keyATM_read(bills_dfm)
keyATM_docsSplit <- keyATM_read(bills_dfm, split = 0.3)

test_that("Split dfm", {
  expect_true("W_split" %in% names(keyATM_docsSplit))
  expect_identical(length(keyATM_docsSplit$W_split$W_raw), 140L)
  expect_identical(length(keyATM_docsSplit$W_raw), 140L)

  expect_error(keyATM_read(bills_dfm, split = -1))
  expect_error(keyATM_read(bills_dfm, split = 1))
  expect_error(keyATM_read(bills_dfm, split = 1.5))
})

test_that("Check summary function", {
  expect_output(print(keyATM_docs), "keyATM_docs object of 140 documents.")
  expect_output(print(keyATM_docsSplit), "keyATM_docs object of 140 documents.")
  expect_output(print(keyATM_docsSplit$W_split), "keyATM_docs object of 140 documents.")
})


