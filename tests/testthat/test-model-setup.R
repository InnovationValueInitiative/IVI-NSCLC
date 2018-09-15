context("stateval.R unit tests")
rm(list = ls())

txseq1 <- txseq(first = "erlotinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2_v1 <- txseq(first = "gefitinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
txseq2_v2 <- txseq(first = "osimertinib",
               second = c("PBDC", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
txseqs_v1 <- txseq_list(seq1 = txseq1, seq2 = txseq2_v1) 
txseqs_v2 <- txseq_list(seq1 = txseq1, seq2 = txseq2_v2)

test_that("create_states", {
  states <- create_states(txseqs_v1)
  expect_true(all(sapply(states, nrow) == 4))
  states <- create_states(txseqs_v2)
  expect_equal(nrow(states[[1]]), 4)
  expect_equal(nrow(states[[2]]), 3)
})

test_that("create_trans_mat", {
  tmat <- create_trans_mats(txseqs_v1)
  expect_true(all(sapply(tmat, nrow) == 4))
  expect_true(inherits(tmat[[1]], "matrix"))
  tmat <- create_trans_mats(txseqs_v2)
  expect_equal(nrow(tmat[[1]]), 4)
  expect_equal(nrow(tmat[[2]]), 3)
})

test_that("create_patients", {
  patients <- create_patients(n = 3)
  expect_true(inherits(patients, "data.table"))
  expect_equal(nrow(patients), 3)
})

# test_that("create_strategies", {
#   strategies <- create_strategies(txseqs_v1)
#   expect_equal(strategies$name_first, sapply(txseqs_v1, function (x) x$first), 
#                check.attributes = FALSE)
#   expect_equal(strategies$name_second_pos, sapply(txseqs_v1, function (x) x$second["pos"]), 
#                check.attributes = FALSE)  
#   
#   strategies <- create_strategies(txseqs_v2)
#   expect_equal(strategies$name_second_neg, sapply(txseqs_v1, function (x) x$second["neg"]), 
#                check.attributes = FALSE)  
#   expect_equal(strategies$name_second_plus_pos, sapply(txseqs_v1, function (x) x$second_plus["pos"]), 
#                check.attributes = FALSE)    
#   expect_equal(strategies$name_second_plus_neg, sapply(txseqs_v1, function (x) x$second_plus["neg"]), 
#                check.attributes = FALSE)  
# })
