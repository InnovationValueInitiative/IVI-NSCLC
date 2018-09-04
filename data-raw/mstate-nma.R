rm(list = ls())
library("data.table")
library("MASS")
library("hesim")

# For now, we will simulate efficacy data
n_samples <- 100
n_treatments <- nrow(iviNSCLC:::treatments)

## shape
shape <- matrix(runif(n_samples, .8, 1.2), ncol = 1)
colnames(shape) <- "Intercept"

## scale
scale_mean <- c(mu = 1,
                runif(n_treatments, -.03, 0))
scale_se <- matrix(0, nrow = length(scale_mean),
                   ncol = length(scale_mean))
diag(scale_se) <- c(1,
                    rep(.02, n_treatments))
scale <- MASS::mvrnorm(n_samples, mu = scale_mean, Sigma = scale_se)
colnames(scale) <- c("mu", paste0("d_", iviNSCLC:::treatments$abb))

## params_surv_list object
params_surv1 <- hesim::params_surv(coefs = list(scale = scale, shape = shape),
                                  dist = "weibullNMA")
params_mstate <- list(weibull = list(first = params_surv_list(sp = params_surv1,
                                                              pd = params_surv1),
                                    second = params_surv_list(sp = params_surv1,
                                                              pd = params_surv1),
                                    second_plus = params_surv_list(sp = params_surv1,
                                                                   pd = params_surv1)))
save(params_mstate, file = "../data/params_mstate.rda", compress = "bzip2")

