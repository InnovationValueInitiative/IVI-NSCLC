context("stateval.R unit tests")
rm(list = ls())

test_that("model_structure", {
  struct <- model_structure()
  expect_true(inherits(struct, "model_structure"))
  expect_error(model_structure(start_line = "second"))
})

test_that("create_states", {
  states <- create_states(model_structure())
  expect_equal(nrow(states), 4)
  
  states <- create_states(model_structure(n_states = "three"))
  expect_equal(states$state_name, c("S1", "P1", "D"))
  
  states <- create_states(model_structure(n_states = "three", start_line = "second"))
  expect_equal(states$state_name, c("S2", "P2", "D"))
  
  expect_error(create_states(2))
})

test_that("create_trans_mat", {
  struct <- model_structure(n_states = "four", start_line = "first")
  tmat <- create_trans_mat(struct)
  expect_equal(nrow(tmat), 4)
  
  struct <- model_structure(n_states = "three", start_line = "first")
  tmat <- create_trans_mat(struct)
  expect_equal(nrow(tmat), 3)
  
  struct <- model_structure(n_states = "three", start_line = "second")
  tmat <- create_trans_mat(struct)
  expect_equal(nrow(tmat), 3)  
  
  expect_error(create_trans_mat("two"))
})

test_that("create_patients", {
  patients <- create_patients(n = 3)
  expect_true(inherits(patients, "data.table"))
  expect_equal(nrow(patients), 3)
})
