context("model-setup.R unit tests")
rm(list = ls())

txseq1 <- txseq(first = "erlotinib",
                second = c("osimertinib", "PBDC"),
                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
                second = c("osimertinib", "PBDC"),
                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
txseq3 <- txseq(first = "osimertinib",
                second = c("PBDC", "PBDC"),
                second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))  
txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)

test_that("model_structure", {
  struct <- model_structure(txseqs)
  expect_true(inherits(struct, "model_structure"))
  expect_error(model_structure(txseqs = 2))
  
  txseqs_v2 <- txseq_list(seq1 = txseq1, seq2 = txseq2, 
                          start_line = "second", mutation = "positive")
  expect_error(model_structure(txseqs = txseqs_v2, n_states = "four"))
})

test_that("create_states", {
  struct <- model_structure(txseqs)
  states <- create_states(struct)
  expect_equal(nrow(states), 4)
  
  struct <- model_structure(txseqs, n_states = "three")
  states <- create_states(struct)
  expect_equal(states$state_name, c("S1", "P1/S2", "D"))
  
  txseqs_v2 <- txseq_list(seq1 = txseq1, seq2 = txseq2, 
                          start_line = "second", mutation = "negative")
  struct <- model_structure(txseqs_v2, n_states = "three")
  states <- create_states(struct)
  expect_equal(states$state_name, c("P1/S2", "P2", "D"))
  
  expect_error(create_states(2, 2))
  expect_error(create_states(txseqs, 2))
})

test_that("create_trans_mat", {
  struct <- model_structure(txseqs, n_states = "four")
  tmat <- create_trans_mat(struct)
  expect_equal(nrow(tmat), 4)
  
  struct <- model_structure(txseqs, n_states = "three")
  tmat <- create_trans_mat(struct)
  expect_equal(nrow(tmat), 3)
  
  txseqs_v2 <- txseq_list(seq1 = txseq1, seq2 = txseq2, 
                          start_line = "second", mutation = "positive")
  struct <- model_structure(txseqs_v2, n_states = "three")
  tmat <- create_trans_mat(struct)
  expect_equal(nrow(tmat), 3)  
  
  expect_error(create_trans_mat("two"))
})

test_that("create_patients", {
  patients <- create_patients(n = 3)
  expect_true(inherits(patients, "data.table"))
  expect_equal(nrow(patients), 3)
})
