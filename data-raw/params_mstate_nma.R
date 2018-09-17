rm(list = ls())
library("data.table")
library("MASS")
library("hesim")

# For now, we will simulate efficacy data with a Weibull distribution
n_samples <- 100

# Treatment variables
## First line
vars_1L <- c("osi", "d_erl", "d_gef", "d_afa", "d_dac")
vars_1L <- c(paste0(vars_1L, "_s1p1"),
             paste0(vars_1L, "_s1d"),
             paste0(vars_1L, "_p1d"))

## Second line
vars_2L <- c("osi", "pbdc", "d_bev", "d_pbdc_bev", "d_erl", "d_gef", "d_afa", "d_dac")
vars_2L <- c(paste0(vars_2L, "_s2p2"),
             paste0(vars_2L, "_s2d"),
             paste0(vars_2L, "_p2d")) 

## Combine and reorder
vars <- c(vars_1L, vars_2L)
d_inds <- which(startsWith(vars, "d"))
mu_inds <- which(!startsWith(vars, "d"))
vars <- vars[c(mu_inds, d_inds)]
vars_scale <- paste0(vars, "_scale")
vars_shape <- paste0(vars, "_shape")
vars <- c(vars_scale, vars_shape)

# Weibull parameters
## Scale
scale_mean_mu <- rep(1, length(mu_inds)) 
scale_mean_d <- runif(length(d_inds), -.03, 0)
scale_mean <- c(scale_mean_mu, scale_mean_d)
scale_se <- matrix(0, 
                   nrow = length(scale_mean),
                   ncol = length(scale_mean))
diag(scale_se) <- c(rep(1, length(mu_inds)),
                    rep(.02, length(d_inds)))
scale <- MASS::mvrnorm(n_samples, mu = scale_mean, Sigma = scale_se)
colnames(scale) <- vars_scale

## Shape
shape_mean_mu <- rep(1, length(mu_inds)) 
shape_mean_d <- runif(length(d_inds), -.03, 0)
shape_mean <- c(shape_mean_mu, shape_mean_d)
shape_se <- matrix(0, 
                   nrow = length(shape_mean),
                   ncol = length(shape_mean))
diag(shape_se) <- c(rep(1, length(mu_inds)),
                    rep(.02, length(d_inds)))
shape <- MASS::mvrnorm(n_samples, mu = shape_mean, Sigma = shape_se)
colnames(shape) <- vars_shape

# Save object
params_mstate_nma <- list()
params_mstate_nma$weibull <- params_surv(coef = list(scale = scale, shape = shape),
                                         dist = "weibullNMA")
save(params_mstate_nma, file = "../data/params_mstate_nma.rda", compress = "bzip2")

