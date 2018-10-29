context("adverse-event.R unit tests")
library("data.table")
rm(list = ls())

# Treatment sequences
txseq1 <- txseq(first = "erlotinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))


test_that("ae_probs", {
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 
  struct <- model_structure(txseqs)
  ae_probs <- ae_probs(n = 4, struct = struct)
  expect_equal(ncol(ae_probs[[1]]), 2)
  expect_equal(attr(ae_probs, "tx_abb"),
               c("erl", "gef"))
  
  tidy_ae_probs <- tidy(ae_probs)
  expect_equal(tidy_ae_probs[strategy_id == 2 & ae_name == "Diarrhea", prob],
               ae_probs$diarrhea[, 2])
})