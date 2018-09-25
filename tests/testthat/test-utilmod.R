context("utilmod.R unit tests")
rm(list = ls())

# Model setup
## Treatment sequences
txseq1 <- txseq(first = "erlotinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))

# Patient population
pats <- create_patients(n = 4)


test_that("create_utilmod first line,  4 health states", {
  # Model structure
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull")
  params <- sample_params(n = 2)
  utilmod <- create_utilmod(params, struct, pats)
  expect_true(inherits(utilmod, "StateVals"))
})

test_that("create_utilmod first line,  3 health states", {
  # Model structure
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull",
                            n_states = "three")
  params <- sample_params(n = 2)
  utilmod <- create_utilmod(params, struct, pats)
  expect_true(inherits(utilmod, "StateVals"))
})

