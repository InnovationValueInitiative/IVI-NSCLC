context("mcda.R unit tests")
library("data.table")
library("flexsurv")
library("hesim")
rm(list = ls())

n_samples <- 5
strategies <- c("Strategy 1", "Strategy 2")
outcome1 <- c(rnorm(n_samples, mean = 10, sd = 5),
              rnorm(n_samples, mean = 8, sd = 4))
outcome2 <- c(rnorm(n_samples, mean = 1500, sd = 90),
              rnorm(n_samples, mean = 1000, sd = 100))
outcomes <- data.table(sample = rep(1:n_samples, length(strategies)),
                       strategy_id = rep(strategies, each = n_samples),
                       criteria1 = outcome1,
                       criteria2 = outcome2)

test_that("performance_matrix", {
  pmat <- performance_matrix(outcomes, 
                             strategy = "strategy_id", 
                             criteria = c("criteria1", "criteria2"),
                             rownames = c("Criteria 1", "Criteria 2"), 
                             colnames = strategies)
  expect_true(inherits(pmat, "matrix"))
  expect_equal(nrow(pmat), ncol(outcomes) - 2)
  expect_equal(ncol(pmat), length(strategies))
  
  pos <- gregexpr(pattern = "\\(", pmat[1,1])[[1]] - 2
  expect_equal(round(mean(outcomes[strategy_id == "Strategy 1", criteria1]), 2),
              as.numeric(substr(pmat[1,1],1, pos[[1]])))
})

test_that("mcda", {
  weights <- c(.7, .3)
  mcda <- mcda(outcomes, sample = "sample", strategy = "strategy_id",
             criteria = c("criteria1", "criteria2"),
             weights = weights,
             optimal = c("low", "high"))
  
  expect_true(inherits(mcda, "list"))
  expect_true(all(mcda$scores$criteria1 >= 0 &
                  mcda$scores$criteria1 <= 100.01)) # .01 because of strange floating point imprecision 
  expect_true(all(mcda$weighted_scores$weighted_score >= 0 & 
                mcda$weighted_scores$weighted_score <= 100.01))  
  expect_true(all(mcda$total_value$score >= 0 & 
                mcda$total_value$score <= 100.01))
  expect_true(all(mcda$prob_rank$prob <= 1 & 
                  mcda$prob_rank$prob >= 0))
})

test_that("txattr_performance", {
  # Simulate economic model
  ictstm <- example_IndivCtstm(n_strategies = 2, n_patients = 3)
  pats <- create_patients(n = 3)
  
  # Set up economic model
  txseq1 <- txseq(first = "erlotinib",
                 second = c("osimertinib", "PBDC"),
                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
  txseq2 <- txseq(first = "gefitinib",
                 second = c("osimertinib", "PBDC"),
                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab")) 
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull")
  
  # Test
  txattr <- txattr_performance(struct = struct, patients = pats, econmod = ictstm)
  expect_true(max(txattr$route) <= 1)
  expect_true(min(txattr$route) >= 0)
  
  ## Errors
  expect_error(txattr_performance(struct = 2, patients = patients, econmod = ictstm))
  expect_error(txattr_performance(struct = struct, patients = 2, econmod = ictstm))
  expect_error(txattr_performance(struct = struct, patients = patients, econmod = 2))
})