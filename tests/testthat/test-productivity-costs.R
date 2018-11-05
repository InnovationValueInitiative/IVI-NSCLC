context("productivity-costs.R unit tests")
library("data.table")
library("flexsurv")
library("hesim")
rm(list = ls())


test_that("prod_costs and add_prod_costs", {
  econmod <- example_IndivCtstm(n_samples = 2, n_patients = 3)
  ce <- econmod$summarize()
  pats <- create_patients(n = 3)
  
  # prod_costs
  ## Works
  prodcosts <- sim_prod_costs(econmod, patients = pats)
  expect_true(inherits(prodcosts, "prod_costs"))
  
  ## Errors
  expect_error(sim_prod_costs(econmod = 2, patients = pats))
  expect_error(sim_prod_costs(econmod, patients = 2))
  
  # add_prod_costs
  ## Works
  ce2 <- add_prod_costs(ce, prodcosts)
  expect_true(inherits(ce2, "ce"))
  
  # Errors
  expect_error(add_prod_costs(ce = 2, prodcosts))
  expect_error(add_prod_costs(ce, prod_costs = 2))
})