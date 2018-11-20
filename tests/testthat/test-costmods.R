context("costmods.R unit tests")
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

test_that("create_costmods first line, 4 health states", {
  n_samples <- 5
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull")
  ae_probs <- ae_probs(n_samples, struct)
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats,
                              ae_probs = ae_probs)
  expect_true(inherits(costmods, "list"))
  expect_true(all(sapply(costmods, function(x) inherits(x, "StateVals"))))
})

test_that("create_costmods first line, 3 health states", {
  n_samples <- 5
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull",
                            n_states = "three")
  ae_probs <- ae_probs(n_samples, struct)
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats,
                              ae_probs = ae_probs)
  expect_true(inherits(costmods, "list"))
  expect_true(all(sapply(costmods, function(x) inherits(x, "StateVals"))))
})