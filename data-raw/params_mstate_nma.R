rm(list = ls())
library("data.table")
library("MASS")
library("hesim")

# For now, we will simulate efficacy data with a Weibull distribution
n_samples <- 100

## First line
vars_1L <- c("osi", "d_erl", "d_gef", "d_afa", "d_dac")

coefs_1L <- function(trans = c("s1p1", "s1d", "p1d"),
                     est_mean = c(0, -1, -.5, -.3, -.1),
                     est_se = c(0, 0, 0, 0, 0)) {
  Sigma <- matrix(0, 
                   nrow = length(est_mean),
                   ncol = length(est_mean))
  diag(Sigma) <- est_se
  coefs <- MASS::mvrnorm(n_samples, mu = est_mean, Sigma = Sigma)
  colnames(coefs) <- paste0(vars_1L, "_", trans, "_scale")
  return(coefs)
}

### Scale
coefs_1L_s1p1_scale <- coefs_1L(trans = "s1p1", 
                                est_mean = c(2, -1, -.5, -.3, -.1))
coefs_1L_s1d_scale <- coefs_1L(trans = "s1d")
coefs_1L_p1d_scale <- coefs_1L(trans = "p1d")

### Shape
zeros <- rep(0, length(vars_1L))
coefs_1L_s1p1_shape <- coefs_1L(trans = "s1p1", est_mean = zeros)
coefs_1L_s1d_shape <- coefs_1L(trans = "s1d", est_mean = zeros)
coefs_1L_p1d_shape <- coefs_1L(trans = "p1d", est_mean = zeros)

## Second line
vars_2L <- c("osi", "pbdc", "d_bev", "d_pbdc_bev", "d_erl", "d_gef", "d_afa", "d_dac")

coefs_2L <- function(trans = c("s2p2", "s2d", "p2d"),
                     est_mean = c(0, 0, -1, -.5, -.3, -.1, .1, .3),
                     est_se = c(0, 0, 0, 0, 0, 0, 0, 0)) {
  Sigma <- matrix(0, 
                   nrow = length(est_mean),
                   ncol = length(est_mean))
  diag(Sigma) <- est_se
  coefs <- MASS::mvrnorm(n_samples, mu = est_mean, Sigma = Sigma)
  colnames(coefs) <- paste0(vars_2L, "_", trans, "_scale")
  return(coefs)
}

### Scale
coefs_2L_s2p2_scale <- coefs_2L(trans = "s2p2")
coefs_2L_s2d_scale <- coefs_2L(trans = "s2d",
                               est_mean = c(-.5, -.5, -1, -.5, -.3, -.1, .1, .3))
coefs_2L_p2d_scale <- coefs_2L(trans = "p2d")

### Shape
zeros <- rep(0, length(vars_2L))
coefs_2L_s2p2_shape <- coefs_2L(trans = "s2p2", est_mean = zeros)
coefs_2L_s2d_shape <- coefs_2L(trans = "s2d", est_mean = zeros)
coefs_2L_p2d_shape <- coefs_2L(trans = "p2d", est_mean = zeros)

# Save object
coefs_scale <- cbind(coefs_1L_s1p1_scale, coefs_1L_s1d_scale, coefs_1L_p1d_scale,
                     coefs_2L_s2p2_scale, coefs_2L_s2d_scale, coefs_2L_p2d_scale)
coefs_shape <- cbind(coefs_1L_s1p1_shape, coefs_1L_s1d_shape, coefs_1L_p1d_shape,
                     coefs_2L_s2p2_shape, coefs_2L_s2d_shape, coefs_2L_p2d_shape)

params_mstate_nma <- list()
params_mstate_nma$weibull <- params_surv(coef = list(scale = coefs_scale, 
                                                     shape = coefs_shape),
                                         dist = "weibull")
save(params_mstate_nma, file = "../data/params_mstate_nma.rda", compress = "bzip2")

