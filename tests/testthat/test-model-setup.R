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
  # Line 1
  tmat <- create_trans_mat(line = "1")
  expect_equal(as.numeric(tmat[1, ]),
               c(NA, 1, NA, NA, 2))

  # Line 2
  tmat <- create_trans_mat(line = "2") 
  expect_equal(as.numeric(tmat[2, ]),
               c(NA, NA, 1, NA, 2))  

  # Line 3  
  tmat <- create_trans_mat(line = "2+") 
  expect_equal(as.numeric(tmat[3, ]),
               c(NA, NA, NA, 1, 2))
  expect_equal(as.numeric(tmat[4, ]),
               c(NA, NA, NA, NA, 3))  
})