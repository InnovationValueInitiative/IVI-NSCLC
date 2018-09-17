context("sample_params.R unit tests")
rm(list = ls())

test_that("sample_params, n > 1", {
  params <- sample_params(n = 2)
  expect_true(inherits(params, "sampled_params"))
  
  # mstate_nma
  expect_true(inherits(params$mstate_nma$weibull, "params_surv"))
  expect_equal(nrow(params$mstate_nma$weibull$coefs$scale), 2)
  expect_equal(params$mstate_nma$weibull$n_samples, 2)
})

test_that("sample_params, n == 1", {
  params <- sample_params(n = 1)
  
  # mstate_nma
  expect_true(is.matrix(params$mstate_nma$weibull$coefs$shape))
})

