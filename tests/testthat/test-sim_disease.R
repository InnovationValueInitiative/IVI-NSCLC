context("sim_disease.R unit tests")
rm(list = ls())

txseq1 <- txseq(first = "erlotinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2) 

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

test_that("create_strategies", {
  # Line 1
  strategies <- create_strategies(txseqs, line = "1") 
  expect_equal(strategies$name,
                as.character(sapply(txseqs, function(x) x[[1]])))
  expect_equal(strategies$d_erl, c(1, 0))
  
  # Line 2
  strategies <- create_strategies(txseqs, line = "2", mutation = "pos") 
  expect_equal(strategies$name,
                as.character(sapply(txseqs, function(x) x[[2]])["pos", ])) 
  expect_equal(strategies$d_osi, c(1, 1))
  
  # Line 2+
  strategies <- create_strategies(txseqs, line = "2+", mutation = "neg")
  expect_equal(strategies$name,
                as.character(sapply(txseqs, function(x) x[[3]])["pos", ]))   
})