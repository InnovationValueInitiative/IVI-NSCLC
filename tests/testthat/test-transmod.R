context("transmod.R unit tests")
rm(list = ls())

expect_all_equal <- function(object, value) {
  expect_true(length(object) >= 1)
  expect_true(all(object == value))
}

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


test_that("create_transmod_data: first line,  4 health states", {
  # Model structure
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull")
  tmat <- create_trans_mat(struct)
  
  # Create data
  transmod_data <- create_transmod_data(struct, tmat, pats)
  expect_true(inherits(transmod_data, "expanded_hesim_data"))
  expect_equal(attributes(transmod_data)$dist, "weibull")
  expect_equal(length(unique(transmod_data$transition_id)), 5)
  expect_equal(max(transmod_data$transition_id), 5)  
  
  sub_dt <- transmod_data[transition_id == 1]
  expect_all_equal(sub_dt$osi_s1p1_scale, 1)
  expect_all_equal(sub_dt$osi_s1p1_shape, 1)
  
  sub_dt <- transmod_data[transition_id == 1 & tx_abb == "erl"]
  expect_all_equal(sub_dt$d_erl_s1p1_scale, 1)
  expect_all_equal(sub_dt$d_erl_s1p1_shape, 1)
  expect_all_equal(sub_dt$d_erl_s1d_scale, 0)
  expect_all_equal(sub_dt$d_gef_s1p1_scale, 0)
  
  sub_dt <- transmod_data[transition_id == 2 & tx_abb == "gef"]
  expect_all_equal(sub_dt$d_gef_s1p1_scale, 0)
  expect_all_equal(sub_dt$d_gef_s1d_scale, 1)
  
  sub_dt <- transmod_data[transition_id == 3 & tx_abb == "osi"]  
  expect_all_equal(sub_dt$osi_s2p2_scale, 1)
})

test_that("create_transmod_data: first line,  3 health states", {
  # Model structure
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, n_states = "three", dist = "weibull")
  tmat <- create_trans_mat(struct)
  
  # Create data
  transmod_data <- create_transmod_data(struct, tmat, pats)  
  expect_true(inherits(transmod_data, "expanded_hesim_data"))
  expect_equal(length(unique(transmod_data$transition_id)), 3)
  expect_equal(max(transmod_data$transition_id), 3)
  
  sub_dt <- transmod_data[transition_id == 1]
  expect_all_equal(sub_dt$osi_s1p1_scale, 1)
  expect_all_equal(sub_dt$osi_s1p1_shape, 1) 
  
  sub_dt <- transmod_data[transition_id == 1 & tx_abb == "erl"]
  expect_all_equal(sub_dt$d_erl_s1p1_scale, 1)
  expect_all_equal(sub_dt$d_erl_s1p1_shape, 1)
  expect_all_equal(sub_dt$d_erl_s1d_scale, 0)
  expect_all_equal(sub_dt$d_gef_s1p1_scale, 0)  
})