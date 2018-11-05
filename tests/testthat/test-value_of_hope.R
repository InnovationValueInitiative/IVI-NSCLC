context("value_of_hope.R unit tests")
library("data.table")
library("flexsurv")
library("hesim")
rm(list = ls())


test_that("value_of_hope", {
  econmod <- example_IndivCtstm(n_samples = 2, n_patients = 3)

  # Works
  voh <- value_of_hope(econmod, comparator = 1)
  expect_true(inherits(voh, "data.table"))
  
  # Errors
  expect_error(value_of_hope(econmod = 2, comparator = 1))

})