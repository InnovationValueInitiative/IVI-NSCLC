context("txseq.R unit tests")
rm(list = ls())

test_that("tx_1L", {
  expect_error(tx_1L(), NA)
})

test_that("tx_2L", {
  expect_error(tx_2L(tx_1L()))
  
  # 1st/2nd generation TKIs
  second <- tx_2L("erlotinib")
  expect_equal(second$pos, "osimertinib")
  expect_equal(second$neg, c("PBDC", "PBDC + bevacizumab"))
  
  # Osimertinib
  second <- tx_2L("osimertinib")
  expect_equal(second$pos, second$neg)
})

test_that("tx_2LP", {
  expect_error(tx_2LP("osimertinib"))
  second_plus <- tx_2LP(c("osimertinib", "PBDC"))
  expect_equal(names(second_plus), c("pos", "neg"))
})

test_that("txseq", {
  txseq <- txseq(first = "erlotinib",
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
  expect_true(inherits(txseq, "txseq"))
  expect_error(txseq(first = 1))  
  
  # Incorrect selections
  expect_error(txseq(first = "PBDC",
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("osimertinib", "PBDC + bevacizumab")))   
  expect_error(txseq(first = "erlotinib",
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("osimertinib", "PBDC + bevacizumab"))) 
    
  # Incorrect types
  expect_error(txseq(first = c("erlotinib", "gefitinib"),
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")))
  expect_error(txseq(first = c("erlotinib"),
                  second = c("osimertinib"),
                  second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")))  
  expect_error(txseq(first = c("erlotinib"),
                  second = 2,
                  second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")))    
  expect_error(txseq(first = c("erlotinib"),
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("PBDC + bevacizumab")))
  expect_error(txseq(first = "erlotinib",
                  second = c("osimertinib", "PBDC"),
                  second_plus = 3)) 
})

