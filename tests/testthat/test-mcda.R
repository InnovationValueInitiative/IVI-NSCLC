context("mcda.R unit tests")
library("data.table")
library("flexsurv")
library("hesim")
rm(list = ls())

strategies <- c("Strategy 1", "Strategy 2")
outcome1 <- c(c(10, 18, 12, 15, 14), # Strategy 1
              c(9, 7, 6, 8, 5))
outcome2 <- c(c(1600, 1500, 1700, 1750, 1800),
                1100, 1150, 1050, 1250, 1300)
n_samples <- length(outcome1)/2
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
                  mcda$scores$criteria1 <= 100)) 
  expect_true(all(mcda$weighted_scores$weighted_score >= 0 & 
                mcda$weighted_scores$weighted_score <= 100))  
  expect_true(all(mcda$total_value$score >= 0 & 
                mcda$total_value$score <= 100))
  expect_true(all(mcda$prob_rank$prob <= 1 & 
                  mcda$prob_rank$prob >= 0))
  
  # With minimum criteria
  ## all criteria 1 less than maximum score (lower scores are better)
  mcda <- mcda(outcomes, sample = "sample", strategy = "strategy_id",
             criteria = c("criteria1", "criteria2"),
             criteria_min = c(max(outcome1) + 10, min(outcome2)),
             criteria_max =  c(max(outcome1) + 1, max(outcome2)),
             weights = weights,
             optimal = c("low", "high"))  
  expect_true(all(mcda$scores$criteria1 == 0))

  ## all criteria 2 greater than maximum score (higher scores are better)
  mcda <- mcda(outcomes, sample = "sample", strategy = "strategy_id",
             criteria = c("criteria1", "criteria2"),
             criteria_min = c(max(outcome1), min(outcome2) - 5),
             criteria_max = c(min(outcome1), min(outcome2) - 1),
             weights = weights,
             optimal = c("low", "high")) 
  expect_true(all(mcda$scores$criteria2 == 100))
  
  # Error
  expect_error(mcda(outcomes, sample = "sample", strategy = "strategy_id",
               criteria = c("criteria1", "criteria2"),
               criteria_min = c(0, 0),
               criteria_max = c(0, 0),
               weights = weights,
               optimal = c("low", "high")))
  
})

test_that("txattr_performance", {
  # Simulate economic model
  ictstm <- example_IndivCtstm(n_patients = 3)
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