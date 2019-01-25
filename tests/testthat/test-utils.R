context("utils.R unit tests")
library("data.table")
library("flexsurv")
library("hesim")
rm(list = ls())


test_that("summarize_outcomes", {
  econmod <- example_IndivCtstm(n_samples = 2, n_patients = 3)
  pats <- create_patients(n = 3)
  prodcosts <- sim_prod_costs(econmod, patients = pats)
  strategy_names <-  c("S1", "S2")
  
  # Normal behavior
  ## With productivity
  res <- summarize_outcomes(econmod = econmod, prod_costs = prodcosts, dr_qalys = .03,
                     dr_costs = .03, strategy_names = strategy_names)
  expect_true(inherits(res, "data.table"))
  expect_equal(colnames(res), c("Outcome", strategy_names))
  
  ## Without productivity
  res <- summarize_outcomes(econmod = econmod, prod_costs = NULL, dr_qalys = .03,
                     dr_costs = .03, strategy_names = strategy_names)
  expect_true(inherits(res, "data.table"))
  expect_false("Productivity costs" %in% res$Outcome)
  expect_false("Societal costs" %in% res$Outcome)
  
  # Errors
  expect_error(summarize_outcomes(econmod = 2, prod_costs = prodcosts, dr_qalys = .03,
                                  dr_costs = .03, strategy_names = strategy_names))
  expect_error(summarize_outcomes(econmod = econmod, prod_costs = 3, dr_qalys = .03,
                                  dr_costs = .03, strategy_names = strategy_names))  
  expect_error(summarize_outcomes(econmod = econmod, prod_costs = prodcosts, dr_qalys = .03,
                     dr_costs = .03, prob = 2))
  econmod$costs_ <- NULL
  expect_error(summarize_outcomes(econmod = econmod, prod_costs = prodcosts, dr_qalys = .03,
                     dr_costs = .03, strategy_names = strategy_names))
  econmod$qalys_ <- NULL
  expect_error(summarize_outcomes(econmod = econmod, prod_costs = prodcosts, dr_qalys = .03,
                     dr_costs = .03, strategy_names = strategy_names))
})