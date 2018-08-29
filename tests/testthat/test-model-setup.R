context("model-setup.R unit tests")
rm(list = ls())

txseq1 <- txseq(first = "erlotinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 

test_that("create_states", {
  states <- create_states(txseqs)
  expect_true(inherits(states, "states"))
  expect_equal(nrow(states), 5)
})

test_that("create_trans_mats", {
  tmats <- create_trans_mats(txseqs)
  expect_true(inherits(tmats, "trans_mats"))
  expect_equal(length(tmats), 3)
})