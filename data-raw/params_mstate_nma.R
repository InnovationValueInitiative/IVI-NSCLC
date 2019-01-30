rm(list = ls())
library("data.table")
library("readxl")
library("MASS")
library("hesim")
library("ggplot2")
treatments <- fread("treatments.csv")
max_months = 50 # Maximum months for plots reflective of follow-up time in RCTs

# Helpful functions ------------------------------------------------------------
d_col <- function(i, j){
  res <- paste0("d[", i, ",", j, "]")
  return(res)
}

mu_col <- function(i){
  res <- paste0("MU[", i, "]")
  return(res)
}

time_vec <- function(t, powers) {
  n_obs <- length(t)
  n_poly <- length(powers)
  X <- matrix(0, nrow = n_obs, ncol = n_poly)
  x1 <- ifelse(powers[1] != rep(0, n_obs), t^powers[1], log(t))
  X[, 1] <- x1
  if (n_poly >= 2) {
    for (i in 2:n_poly) {
      if (powers[i] == powers[(i - 1)]) 
        x2 <- log(t) * x1
      else x2 <- ifelse(powers[i] != rep(0, n_obs), t^powers[i], 
        log(t))
      X[, i] <- x2
      x1 <- x2
      }
  }
  return(cbind(1, X))
}

# Function to compute probability in stable and progressed states and PFS/OS
# by first computing transition specific hazard functions
state_probs <- function(gamma, powers, t, n_sims){
  p_s <- p_p <- haz_sp <- haz_sd <- haz_pd <- matrix(NA, nrow = length(t), ncol = n_sims)
  p_s[1, ] <- 1 # PFS is 1 in first time period
  haz_sp[1, ] <- haz_sd[1, ] <- haz_pd[1, ] <- 0 # Hazard is 0 in first period
  p_p[1, ] <- 0 # Proportion in the progressed state is 0 in first time period
  for (j in 2:length(t)){
    tvec <- time_vec(t[j], powers)
    haz_sp[j, ] <- exp(gamma[[1]] %*% t(tvec))
    haz_sd[j, ] <- exp(gamma[[2]] %*% t(tvec))
    haz_pd[j, ] <- exp(gamma[[3]] %*% t(tvec))
    p_s[j, ] <- p_s[j - 1, ] * exp(-(haz_sp[j, ] + haz_sd[j, ]))
    term2_num <- p_s[j - 1, ] * haz_sp[j, ] * 
                 (exp(-(haz_sp[j, ] + haz_sd[j, ])) - exp(-haz_pd[j, ]))
    term2_denom <- haz_pd[j, ] - haz_sp[j, ] - haz_sd[j, ]
    p_p[j, ] <- p_p[j - 1, ] * exp(-haz_pd[j, ]) + term2_num/term2_denom
  } # End loop over time periods
  pfs <- p_s
  os <- p_s + p_p  
  return(list(haz_sp = haz_sp, haz_sd = haz_sd, haz_pd = haz_pd,
              pfs = pfs, os = os))
}

# Function to compute hazard ratios
hazard_ratios <- function(dvec, powers, t, n_sims){
  hr_sp <- hr_sd <- hr_pd <- matrix(NA, nrow = length(t), ncol = n_sims)
  hr_sp[1, ] <- hr_sd[1, ] <- hr_pd[1, ] <- 1 # HR is 1 first period
  for (j in 2:length(t)){
    tvec <- time_vec(t[j], powers)
    hr_sp[j, ] <- exp(dvec[[1]] %*% t(tvec))
    hr_sd[j, ] <- exp(dvec[[2]] %*% t(tvec))
    hr_pd[j, ] <- exp(dvec[[3]] %*% t(tvec))
  } # End loop over time periods
  return(list(hr_sp = hr_sp, hr_sd = hr_sd, hr_pd = hr_pd))  
}

# Models types -----------------------------------------------------------------
mod_names <- c("weibull", "gompertz", "fracpoly1", "fracpoly2")
fp_powers <- list(weibull = c(0),
                  gompertz = c(1),
                  fracpoly1 = c(0, 0),
                  fracpoly2 = c(0, 1))
n_models <- length(mod_names)
mod_dists <- c("weibullNMA", "gompertz", "fracpoly", "fracpoly")
mod_aux <- list(weibull = NULL,
                gompertz = NULL,
                fracpoly1 = list(powers = fp_powers$fracpoly1,
                                 cumhaz_method = "riemann",
                                 step = .02),
                fracpoly2 = list(powers = fp_powers$fracpoly1,
                                 cumhaz_method = "riemann",
                                 step = .02))

model_lookup <- function(powers){
  powers <- paste(as.character(powers), collapse = ", ")
  model <- switch(powers,
                  "0" = "Weibull",
                  "1" = "Gompertz",
                  "0, 0" = "Fractional polynomial (0, 0)",
                  "0, 1" = "Fractional polynomial (0, 1)")
  return(model)
}

# Followup time ----------------------------------------------------------------
pfs_followup <- 36
os_followup <- 48

# NMA parameter estimates ------------------------------------------------------
sample_posterior <- function(x, n_sims, random = FALSE){
  for (i in 1:length(x)){ 
    if (nrow(x[[i]]) != n_sims){
      if (random){ # Randomly sample from posterior
        sampled_rows <- sample.int(nrow(x[[i]]), n_sims, replace = FALSE)
        x[[i]] <- x[[i]][sampled_rows, ]
      } else{ # Use last n_sims rows
        last_row <- nrow(x[[i]])
        first_row <- last_row - n_sims + 1
        x[[i]] <- x[[i]][first_row:last_row, ]
      }
    }
  }
  return (x)
}

select_vars <- function(patterns, x) {
  colnums <- grep(paste(patterns, collapse = "|"),
                  colnames(x))
  return(x[, colnums, with = FALSE])
}

# 1L
read_nma_1L <- function(filename){
  dat <- fread(filename)
  dat <- select_vars(patterns = c("d\\[", "PFS\\[", "OS\\["),
                                       dat)
  dat <- as.matrix(dat)
  return(dat)
}
nma_1L <- list(weibull = read_nma_1L("mstate-nma/nma-1L-fe-fp-p0.csv"),
               gompertz = read_nma_1L("mstate-nma/nma-1L-fe-fp-p1.csv"), 
               fracpoly1 = read_nma_1L("mstate-nma/nma-1L-fe-fp-p00.csv"), # (p1 = 0, p2 = 0)
               fracpoly2 = read_nma_1L("mstate-nma/nma-1L-fe-fp-p01.csv")) # (p1 = 0, p2 = 1)
ma_1L <- list(weibull = as.matrix(fread("mstate-nma/ma-1L-fe-gef-fp-p0.csv")),
              gompertz = as.matrix(fread("mstate-nma/ma-1L-fe-gef-fp-p1.csv")), 
              fracpoly1 = as.matrix(fread("mstate-nma/ma-1L-fe-gef-fp-p00.csv")), # (p1 = 0, p2 = 0)
              fracpoly2 = as.matrix(fread("mstate-nma/ma-1L-fe-gef-fp-p01.csv"))) # (p1 = 0, p2 = 1)
n_sims <- min(c(sapply(nma_1L, nrow), sapply(ma_1L, nrow)))
nma_1L <- sample_posterior(nma_1L, n_sims)
ma_1L <- sample_posterior(ma_1L, n_sims)

# 2L (osimertinib and T790M+)
ma_2L_t790m_osi <- list(weibull = as.matrix(fread("mstate-nma/ma-2L-fe-t790m-osi-fp-p0.csv")),
                        gompertz = as.matrix(fread("mstate-nma/ma-2L-fe-t790m-osi-fp-p1.csv")), 
                        fracpoly1 = as.matrix(fread("mstate-nma/ma-2L-fe-t790m-osi-fp-p00.csv")), # (p1 = 0, p2 = 0)
                        fracpoly2 = as.matrix(fread("mstate-nma/ma-2L-fe-t790m-osi-fp-p01.csv"))) # (p1 = 0, p2 = 1)
ma_2L_t790m_osi <- sample_posterior(ma_2L_t790m_osi, n_sims)

# 2L (PBDC)
ma_2L_pbdc <- list(weibull = as.matrix(fread("mstate-nma/ma-2L-fe-pbdc-fp-p0.csv")),
                   gompertz = as.matrix(fread("mstate-nma/ma-2L-fe-pbdc-fp-p1.csv")), 
                   fracpoly1 = as.matrix(fread("mstate-nma/ma-2L-fe-pbdc-fp-p00.csv")), # (p1 = 0, p2 = 0)
                   fracpoly2 = as.matrix(fread("mstate-nma/ma-2L-fe-pbdc-fp-p01.csv"))) # (p1 = 0, p2 = 1)
ma_2L_pbdc <- sample_posterior(ma_2L_pbdc, n_sims)

# First line parameter estimates -----------------------------------------------
# NMA parameter lookups
nma_params_lookup_1L <- vector(mode = "list", length = n_models)
names(nma_params_lookup_1L) <- mod_names
for (i in 1:n_models){ 
    nma_params_lookup_1L[[i]] <- data.table(read_excel("mstate-nma/params-lookup-1L.xlsx",
                                                       sheet = mod_names[i]))
} # End loop over models

# Fist line treatments
nma_tx_lookup_1L <- fread("mstate-nma/tx-lookup-1L.csv")
econmod_tx_1L <- tx_1L()
row <- match(econmod_tx_1L, nma_tx_lookup_1L$tx_name)
econmod_tx_lookup_1L <- data.table(name = econmod_tx_1L,
                                   num = row,
                                   abb = nma_tx_lookup_1L[row, tx_abb])

# Function to extract parameters
create_mstate_coefs_1L <- function(nma_post, ma_gef_post, econmod_tx_lookup, 
                                        nma_params_lookup){
  # Args:
  # nma_post: Posterior distribution of parameters from NMA for estimating
  #           relative treatment effects.
  # ma_gef_post: Posterior distribution of parameters from meta-analysis for
  #              estimating absolute effects for gefitinib (the reference arm).
  # econmod_tx_lookup: Numeric ID for each treatment in economic model in the NMA. 
  # nma_params_lookup: Numeric ID for "d" columns and "mu" columns in NMA (for
  #                    relative treatment effects) and MA (for absolute effects)
  #                    by transition and parameter.
  
  # Number of simulations
  if (nrow(nma_post) != nrow(ma_gef_post)){
    stop("The number of posterior simulations must be the same for 'nma_post' and 'ma_gef_post'.")
  } else{
    n_sims <- nrow(nma_post)
  }
    
  # Parameters
  nma_params_lookup <- nma_params_lookup[order(transition_id, param_id)]
  params <- unique(nma_params_lookup$param) # Note: must be sorted in same order as in params_surv objects in hesim
  n_params <- length(params)
  
  # Transitions
  trans <- unique(nma_params_lookup$transition)
  n_trans <- length(trans)
  
  # Store results
  coefs <- vector(mode = "list", length = n_params)
  names(coefs) <- params
  n_cols <- n_trans * nrow(econmod_tx_lookup)
  for (i in 1:n_params){
    coefs[[i]] <- matrix(0, nrow = n_sims, ncol = n_cols)
    colnames(coefs[[i]]) <- rep("tmp", n_cols)
    nma_params_lookup_i <- nma_params_lookup[param == params[i]]
    cntr <- 1
    for (j in 1:nrow(econmod_tx_lookup)){
      for (k in 1:n_trans){
        tx_abb <- econmod_tx_lookup$abb[j]
        if (tx_abb == "gef"){ ### The reference arm
          mu_num_k <- nma_params_lookup_i[k, mu_num]
          if (!is.na(mu_num_k)){
            coefs[[i]][, cntr] <- ma_gef_post[, mu_col(mu_num_k)]
          }
          colnames(coefs[[i]])[cntr] <- paste0(tx_abb, "_", trans[k], "_", params[i])
        } else{ ### The relative treatment effects
          d_num_k <- nma_params_lookup_i[k, d_num]
          if (!is.na(d_num_k)){
            col_name <- d_col(econmod_tx_lookup$num[j], d_num_k)
            coefs[[i]][, cntr] <- nma_post[, col_name]
          }
          colnames(coefs[[i]])[cntr] <- paste0("d_", tx_abb, "_", trans[k], "_", params[i])
        } 
        cntr <- cntr + 1          
        } # End look over transitions
    } # End loop over treatments
  } # End loop over parameters
  return(coefs)
}

mstate_coefs_1L <- list()
mstate_coefs_1L$weibull <- create_mstate_coefs_1L(nma_post = nma_1L$weibull,
                                                  ma_gef_post = ma_1L$weibull,
                                                  econmod_tx_lookup = econmod_tx_lookup_1L,
                                                  nma_params_lookup = nma_params_lookup_1L$weibull)
mstate_coefs_1L$gompertz <- create_mstate_coefs_1L(nma_post = nma_1L$gompertz,
                                                   ma_gef_post = ma_1L$gompertz,
                                                   econmod_tx_lookup = econmod_tx_lookup_1L,
                                                   nma_params_lookup = nma_params_lookup_1L$gompertz)
mstate_coefs_1L$fracpoly1 <- create_mstate_coefs_1L(nma_post = nma_1L$fracpoly1,
                                                    ma_gef_post = ma_1L$fracpoly1,
                                                    econmod_tx_lookup = econmod_tx_lookup_1L,
                                                    nma_params_lookup = nma_params_lookup_1L$fracpoly1)
mstate_coefs_1L$fracpoly2 <- create_mstate_coefs_1L(nma_post = nma_1L$fracpoly2,
                                                    ma_gef_post = ma_1L$fracpoly2,
                                                    econmod_tx_lookup = econmod_tx_lookup_1L,
                                                    nma_params_lookup = nma_params_lookup_1L$fracpoly2)

# Second line parameter estimates (T790M positive) -----------------------------
# MA parameter lookups
ma_params_lookup_2L_t790m_osi <- vector(mode = "list", length = n_models)
names(ma_params_lookup_2L_t790m_osi) <- mod_names
for (i in 1:n_models){ 
    ma_params_lookup_2L_t790m_osi[[i]] <- data.table(read_excel("mstate-nma/params-lookup-2L-t790m-osi.xlsx",
                                                       sheet = mod_names[i]))
} # End loop over models

create_mstate_coefs_2L_t790m_osi <- function(ma_post, ma_params_lookup){
  # Args:
  # ma_post: Posterior distribution of parameters from MA for estimating
  #          absolute effects for T790M+ patients at second line. Only
  #          relevant for osimertinib.
  # ma_params_lookup: Numeric ID for"mu" columns in MA (for absolute effects)
  #                    by transition and parameter.

  n_sims <- nrow(ma_post)
  ma_params_lookup <- ma_params_lookup[order(transition_id, param_id)]
    
  # Parameters
  params <- unique(ma_params_lookup$param)
  n_params <- length(params)

  # Transitions
  trans <- unique(ma_params_lookup$transition)
  n_trans <- length(trans)   
  
  # Store results
  coefs <- vector(mode = "list", length = n_params)
  names(coefs) <- params
  n_cols <- n_trans 
  for (i in 1:n_params){
    coefs[[i]] <- matrix(0, nrow = n_sims, ncol = n_cols)
    colnames(coefs[[i]]) <- paste0("osi", "_", trans, "_", params[i])
    ma_params_lookup_i <- ma_params_lookup[param == params[i]]
    for (k in 1:n_trans){
      mu_num_k <- ma_params_lookup_i[k, mu_num]
      if (!is.na(mu_num_k)){
       coefs[[i]][, k] <- ma_post[, mu_col(mu_num_k)] 
      }
    } # End loop over transitions
  } # End loop over parameters
  return(coefs)  
}

mstate_coefs_2L_t790m_osi <- list()
mstate_coefs_2L_t790m_osi$weibull <- create_mstate_coefs_2L_t790m_osi(ma_post = ma_2L_t790m_osi$weibull,
                                                                      ma_params_lookup = ma_params_lookup_2L_t790m_osi$weibull)
mstate_coefs_2L_t790m_osi$gompertz <- create_mstate_coefs_2L_t790m_osi(ma_post = ma_2L_t790m_osi$gompertz,
                                                                       ma_params_lookup = ma_params_lookup_2L_t790m_osi$gompertz)
mstate_coefs_2L_t790m_osi$fracpoly1 <- create_mstate_coefs_2L_t790m_osi(ma_post = ma_2L_t790m_osi$fracpoly1,
                                                                       ma_params_lookup = ma_params_lookup_2L_t790m_osi$fracpoly1)
mstate_coefs_2L_t790m_osi$fracpoly2 <- create_mstate_coefs_2L_t790m_osi(ma_post = ma_2L_t790m_osi$fracpoly2,
                                                                       ma_params_lookup = ma_params_lookup_2L_t790m_osi$fracpoly2)
                   
# Second line parameter estimates (T790M negative) -----------------------------
# MA parameter lookups
ma_params_lookup_2L_pbdc <- vector(mode = "list", length = n_models)
names(ma_params_lookup_2L_pbdc) <- mod_names
for (i in 1:n_models){ 
    ma_params_lookup_2L_pbdc[[i]] <- data.table(read_excel("mstate-nma/params-lookup-2L-pbdc.xlsx",
                                                           sheet = mod_names[i]))
} # End loop over models

# All possible 2nd line treatments other than osimertinib
econmod_tx_2L <- vector(mode = "list", length = length(econmod_tx_1L))
for (i in 1:length(econmod_tx_2L)){
  econmod_tx_2L[[i]] <- unlist(tx_2L(econmod_tx_1L[i]))
}
econmod_tx_2L <- unique(unlist(econmod_tx_2L))
tx_abb_2L <- treatments[match(econmod_tx_2L, treatments$tx_name), tx_abb]
tx_abb_2L <- tx_abb_2L[tx_abb_2L != "osi"]

create_mstate_coefs_2L_pbdc <- function(ma_post, ma_params_lookup, tx_abb){
  # Args:
  # ma_post: Posterior distribution of parameters from MA for estimating
  #          absolute effects for patients using PBDC at second line. 2nd line
  #          treatment effects are assumed to be the same for PBDC, 
  #          PBC + anti-VEGF therapy, and 1st/2nd generation TKIs.
  # ma_params_lookup: Numeric ID for"mu" columns in MA (for absolute effects)
  #                    by transition and parameter.
  # tx_abb: Abbreviations for 2nd line treatments other than osimertinib.

  n_sims <- nrow(ma_post)
  ma_params_lookup <- ma_params_lookup[order(transition_id, param_id)]
  
  # Parameters
  params <- unique(ma_params_lookup$param)
  n_params <- length(params)

  # Transitions
  trans <- unique(ma_params_lookup$transition)
  n_trans <- length(trans)  
  
  # Treatments
  n_tx <- length(tx_abb)
  
  # Store results
  coefs <- vector(mode = "list", length = n_params)
  names(coefs) <- params
  n_cols <- n_trans * n_tx 
  for (i in 1:n_params){
    coefs[[i]] <- matrix(0, nrow = n_sims, ncol = n_cols)
    colnames(coefs[[i]]) <- colnames(coefs[[i]]) <- rep("tmp", n_cols)
    ma_params_lookup_i <- ma_params_lookup[param == params[i]]
    cntr <- 1
    for (j in 1:length(tx_abb)){
      for (k in 1:n_trans){
        if (tx_abb[j] == "pbdc"){ # The reference arm
          mu_num_k <- ma_params_lookup_i[k, mu_num]
          if (!is.na(mu_num_k)){
            coefs[[i]][, cntr] <- ma_post[, mu_col(mu_num_k)] 
          }
          colnames(coefs[[i]])[cntr] <- paste0(tx_abb[j], "_", trans[k], "_", params[i])  
        } else{ # All relative treatment effects are 0
          colnames(coefs[[i]])[cntr] <- paste0("d_", tx_abb[j], "_", trans[k], "_", params[i]) 
        }
        cntr <- cntr + 1
    } # End loop over transitions
    } # End loop over treatments
  } # End loop over parameters
  return(coefs)  
}

mstate_coefs_2L_pbdc <- list()
mstate_coefs_2L_pbdc$weibull <- create_mstate_coefs_2L_pbdc(ma_post = ma_2L_pbdc$weibull,
                                                            ma_params_lookup = ma_params_lookup_2L_pbdc$weibull,
                                                            tx_abb = tx_abb_2L)
mstate_coefs_2L_pbdc$gompertz <- create_mstate_coefs_2L_pbdc(ma_post = ma_2L_pbdc$gompertz,
                                                            ma_params_lookup = ma_params_lookup_2L_pbdc$gompertz,
                                                            tx_abb = tx_abb_2L)
mstate_coefs_2L_pbdc$fracpoly1 <- create_mstate_coefs_2L_pbdc(ma_post = ma_2L_pbdc$fracpoly1,
                                                            ma_params_lookup = ma_params_lookup_2L_pbdc$fracpoly1,
                                                            tx_abb = tx_abb_2L)
mstate_coefs_2L_pbdc$fracpoly2 <- create_mstate_coefs_2L_pbdc(ma_post = ma_2L_pbdc$fracpoly2,
                                                            ma_params_lookup = ma_params_lookup_2L_pbdc$fracpoly2,
                                                            tx_abb = tx_abb_2L)

# Combine parameter estimates --------------------------------------------------
params_mstate_nma <- vector(mode = "list", length = n_models)
names(params_mstate_nma) <- mod_names
for (i in 1:n_models){
  n_params <- length(mstate_coefs_1L[[i]])
  mstate_coefs_i <- vector(mode = "list", length = n_params)
  names(mstate_coefs_i) <- names(mstate_coefs_1L[[i]])
  for (j in 1:n_params){
    mstate_coefs_i[[j]] <- cbind(mstate_coefs_1L[[i]][[j]],
                                 mstate_coefs_2L_t790m_osi[[i]][[j]],
                                 mstate_coefs_2L_pbdc[[i]][[j]]) 
  } # End loop over parameters
  params_mstate_nma[[i]] <- hesim::params_surv(coefs = mstate_coefs_i,
                                               dist = mod_dists[i],
                                               aux = mod_aux[[i]])
} # End loop over models

# Regression tables ------------------------------------------------------------
coef_summary <- function(x){
  probs <- c(.025, .25, .5, .75, .975)
  ests <- matrix(NA,  nrow = ncol(x), ncol = length(probs))
  for (i in 1:ncol(ests)){
    ests[, i] <- apply(x, 2, quantile, prob = probs[i])
  }
  colnames(ests) <- paste0("q", probs)
  return(ests)
}

nma_coef_table <- function(x, econmod_tx_lookup, params_lookup,
                          powers){
  
  # Extract relevant columns
  x <- x[, grep(paste0("d", "\\["), colnames(x))]
  
  # Compute quantiles of regression coefficients
  ests <- coef_summary(x)
  
  # Add ID variables
  comma_pos <- regexpr(",", colnames(x))
  tx_num <- coef_num <- rep(NA, ncol(x))
  for (i in 1:ncol(x)){
    tx_num[i] <- substr(colnames(x)[i], comma_pos[i] - 1, comma_pos[i] - 1)
    coef_num[i] <-substr(colnames(x)[i], comma_pos[i] + 1, comma_pos[i] + 1) 
  }
  tx_num <- as.integer(tx_num)
  coef_num <- as.integer(coef_num)
  tx_name <- econmod_tx_lookup$name[match(tx_num, econmod_tx_lookup$num)]
  transition <- params_lookup$transition[match(coef_num, params_lookup$d_num)]
  
  # Return
  tbl <- data.table(model = model_lookup(powers),
                    tx_name = tx_name, transition = transition,
                    coef_num = coef_num, ests)
  tbl <- tbl[!is.na(tx_name)]
  return(tbl)
}

ma_coef_table <- function(x, params_lookup, powers){
  
  # Extract relevant columns
  x <- x[, grep(paste0("MU", "\\["), colnames(x))]
  
  # Compute quantiles of regression coefficients
  ests <- coef_summary(x)
  
  # Add ID variables
  bracket_pos <- regexpr("\\[", colnames(x))
  tx_num <- coef_num <- rep(NA, ncol(x))
  for (i in 1:ncol(x)){
    coef_num[i] <-substr(colnames(x)[i], bracket_pos[i] + 1, bracket_pos[i] + 1) 
  }
  coef_num <- as.integer(coef_num)
  transition <- params_lookup$transition[match(coef_num, params_lookup$mu_num)]  
  
  # Return
  tbl <- data.table(model = model_lookup(powers),
                    transition = transition,
                    coef_num = coef_num, ests)
  return(tbl)  
}

# 1st line NMA
coef_nma_1L <- list()
coef_nma_1L$wei <- nma_coef_table(x = nma_1L$weibull, 
                                     econmod_tx_lookup = econmod_tx_lookup_1L,
                                     params_lookup = nma_params_lookup_1L$weibull,
                                     powers = 0)
coef_nma_1L$gomp <- nma_coef_table(x = nma_1L$gompertz, 
                                      econmod_tx_lookup = econmod_tx_lookup_1L,
                                      params_lookup = nma_params_lookup_1L$gompertz,
                                      powers = 1)
coef_nma_1L$fracpoly1 <- nma_coef_table(x = nma_1L$fracpoly1, 
                                            econmod_tx_lookup = econmod_tx_lookup_1L,
                                            params_lookup = nma_params_lookup_1L$fracpoly1,
                                            powers = c(0, 0))
coef_nma_1L$fracpoly2 <- nma_coef_table(x = nma_1L$fracpoly2, 
                                            econmod_tx_lookup = econmod_tx_lookup_1L,
                                            params_lookup = nma_params_lookup_1L$fracpoly2,
                                            powers = c(0, 1))
coef_nma_1L <- rbindlist(coef_nma_1L)
coef_nma_1L[, line := 1]
setcolorder(coef_nma_1L, "line")


# 1st line MA
coef_ma_1L <- list()
coef_ma_1L$wei <- ma_coef_table(x = ma_1L$weibull, 
                                    params_lookup = nma_params_lookup_1L$weibull,
                                    powers = 0)
coef_ma_1L$gomp <- ma_coef_table(x = ma_1L$gompertz, 
                                    params_lookup = nma_params_lookup_1L$gompertz,
                                    powers = 1)
coef_ma_1L$fracpoly1 <- ma_coef_table(x = ma_1L$fracpoly1, 
                                          params_lookup = nma_params_lookup_1L$fracpoly1,
                                          powers = c(0, 0))
coef_ma_1L$fracpoly2 <- ma_coef_table(x = ma_1L$fracpoly2, 
                                          params_lookup = nma_params_lookup_1L$fracpoly2,
                                          powers = c(0, 1))
coef_ma_1L <- rbindlist(coef_ma_1L)
coef_ma_1L[, line := 1]

# 2nd line MA (T790M positive - osimertinib)
coef_ma_2L_t790m_osi <- list()
coef_ma_2L_t790m_osi$wei <- ma_coef_table(x = ma_2L_t790m_osi$weibull, 
                                                         params_lookup = ma_params_lookup_2L_t790m_osi$weibull,
                                                         powers = 0)
coef_ma_2L_t790m_osi$gomp <- ma_coef_table(x = ma_2L_t790m_osi$gompertz, 
                                              params_lookup = ma_params_lookup_2L_t790m_osi$gompertz,
                                              powers = 1)
coef_ma_2L_t790m_osi$fracpoly1 <- ma_coef_table(x = ma_2L_t790m_osi$fracpoly1, 
                                                    params_lookup = ma_params_lookup_2L_t790m_osi$fracpoly1,
                                                    powers = c(0, 0))
coef_ma_2L_t790m_osi$fracpoly2 <- ma_coef_table(x = ma_2L_t790m_osi$fracpoly2, 
                                                    params_lookup = ma_params_lookup_2L_t790m_osi$fracpoly2,
                                                    powers = c(0, 1))
coef_ma_2L_t790m_osi <- rbindlist(coef_ma_2L_t790m_osi)
coef_ma_2L_t790m_osi[, line := 2]
coef_ma_2L_t790m_osi[, mutation := 1]

# 2nd line MA (PBDC)
coef_ma_2L_pbdc <- list()
coef_ma_2L_pbdc$wei <- ma_coef_table(x = ma_2L_pbdc$weibull, 
                                        params_lookup = ma_params_lookup_2L_pbdc$weibull,
                                        powers = 0)
coef_ma_2L_pbdc$gomp <- ma_coef_table(x = ma_2L_pbdc$gompertz, 
                                          params_lookup = ma_params_lookup_2L_pbdc$gompertz,
                                          powers = 1)
coef_ma_2L_pbdc$fracpoly1 <- ma_coef_table(x = ma_2L_pbdc$fracpoly1, 
                                               params_lookup = ma_params_lookup_2L_pbdc$fracpoly1,
                                               powers = c(0, 0))
coef_ma_2L_pbdc$fracpoly2 <- ma_coef_table(x = ma_2L_t790m_osi$fracpoly2, 
                                               params_lookup = ma_params_lookup_2L_pbdc$fracpoly2,
                                               powers = c(0, 1))
coef_ma_2L_pbdc <- rbindlist(coef_ma_2L_pbdc)
coef_ma_2L_pbdc[, line := 2]
coef_ma_2L_pbdc[, mutation := 0]

# Combine
mstate_nma_coef <- coef_nma_1L
mstate_ma_coef <- rbind(coef_ma_1L,
                     coef_ma_2L_pbdc,
                     coef_ma_2L_t790m_osi,
                     fill = TRUE)
setcolorder(mstate_ma_coef, c("line", "mutation"))

# 1st line PFS/OS --------------------------------------------------------------
surv_1L <- function(n_months, nma_post, ma_gef_post, econmod_tx_lookup,
                    nma_params_lookup, powers){
  
  n_params <- length(unique(nma_params_lookup$param))
  n_trans <- 3
  t <- seq(0, n_months)
  n_powers <- length(powers)
  
  # Number of simulations
  if (nrow(nma_post) != nrow(ma_gef_post)){
    stop("The number of posterior simulations must be the same for 'nma_post' and 'ma_gef_post'.")
  } else{
    n_sims <- nrow(nma_post)
  }  

  # Loop over treatments and months
  surv <- vector(mode = "list", length = nrow(econmod_tx_lookup))
  names(surv) <- econmod_tx_lookup$name
  
  for (i in 1:nrow(econmod_tx_lookup)){
    
    # Parameter estimates for treatment i (i.e., "gamma")
    tx_num <- econmod_tx_lookup$num[i]
    gamma <- dvec <- vector(mode = "list", length = n_trans)
    for (j in 1:n_trans){
      nma_params_lookup_j <- nma_params_lookup[transition_id == j]
      gamma[[j]] <- matrix(NA, nrow = n_sims, ncol = nrow(nma_params_lookup_j))
      dvec[[j]] <- matrix(0, nrow = n_sims, ncol = nrow(nma_params_lookup_j))
      for (k in 1:nrow(nma_params_lookup_j)){
        d_num <- nma_params_lookup_j$d_num[k]
        if (!is.na(d_num)){
          d <- nma_post[, d_col(tx_num, d_num)]
        } else{
          d <- rep(0, n_sims)
        }
        dvec[[j]][, k] <- d
        mu_num <- nma_params_lookup_j$mu_num[k]
        if (!is.na(mu_num)){
          mu <- ma_gef_post[, mu_col(mu_num)]
        } else{
          mu <- rep(0, n_sims)
        }
        gamma[[j]][, k] <- mu + d
      } # End loop over parameters
    } # End loop over transitions
    
    # Compute PFS/OS given parameter estimates
    hr <- hazard_ratios(dvec, powers, t, n_sims)
    stprobs <- state_probs(gamma, powers, t, n_sims)
    surv[[i]] <- data.table(line = 1,
                           mutation = NA,
                           model = model_lookup(powers),
                           tx_name = econmod_tx_lookup$name[i],
                           sample = rep(1:n_sims, each = length(t)),
                           month = rep(t, n_sims),
                           hr_sp = c(hr$hr_sp),
                           hr_sd = c(hr$hr_sd),
                           hr_pd = c(hr$hr_pd),
                           haz_sp = c(stprobs$haz_sp),
                           haz_sd = c(stprobs$haz_sd),
                           haz_pd = c(stprobs$haz_pd),
                           pfs = c(stprobs$pfs),
                           os = c(stprobs$os))
  } # End loop over treatments
  surv <- rbindlist(surv)
  return(surv)
}

surv_1L_est <- list()
surv_1L_est$wei <- surv_1L(n_months = max_months, 
                           nma_post = nma_1L$weibull, 
                           ma_gef_post = ma_1L$weibull,
                           econmod_tx_lookup = econmod_tx_lookup_1L,
                           nma_params_lookup = nma_params_lookup_1L$weibull,
                           powers = 0)
surv_1L_est$gomp <- surv_1L(n_months = max_months, 
                            nma_post = nma_1L$gompertz, 
                            ma_gef_post = ma_1L$gompertz,
                            econmod_tx_lookup = econmod_tx_lookup_1L,
                            nma_params_lookup = nma_params_lookup_1L$gompertz,
                            powers = 1)
surv_1L_est$fracpoly1 <- surv_1L(n_months = max_months, 
                                 nma_post = nma_1L$fracpoly1, 
                                 ma_gef_post = ma_1L$fracpoly1,
                                 econmod_tx_lookup = econmod_tx_lookup_1L,
                                 nma_params_lookup = nma_params_lookup_1L$fracpoly1,
                                 powers = c(0, 0))
surv_1L_est$fracpoly2 <- surv_1L(n_months = max_months, 
                                 nma_post = nma_1L$fracpoly2, 
                                 ma_gef_post = ma_1L$fracpoly2,
                                 econmod_tx_lookup = econmod_tx_lookup_1L,
                                 nma_params_lookup = nma_params_lookup_1L$fracpoly2,
                                 powers = c(0, 1))
surv_1L_est <- rbindlist(surv_1L_est)

# Second line PFS/OS -----------------------------------------------------------
surv_2L <- function(n_months, ma_post,
                    ma_params_lookup, powers, tx_name, mutation){
  
  n_params <- length(unique(ma_params_lookup$param))
  n_trans <- 3
  t <- seq(0, n_months)
  n_powers <- length(powers)
  n_sims <- nrow(ma_post)

  # Loop over treatments and months
  pfs <- os <- matrix(NA, nrow = length(t), ncol = n_sims)
  pfs[1, ] <- os[1, ] <- 1 # PFS and OS are 1 in f
  
  # Parameter estimates for treatment (i.e., "gamma")
  gamma <- vector(mode = "list", length = n_trans)
    
  for (j in 1:n_trans){
    ma_params_lookup_j <- ma_params_lookup[transition_id == j + 2]
    gamma[[j]] <- matrix(NA, nrow = n_sims, ncol = nrow(ma_params_lookup_j))
    for (k in 1:nrow(ma_params_lookup_j)){
      mu_num <- ma_params_lookup_j$mu_num[k]
      if (!is.na(mu_num)){
        mu <- ma_post[, mu_col(mu_num)]
      } else{
        mu <- rep(0, n_sims)
      }
      gamma[[j]][, k] <- mu
    } # End loop over parameters
  } # End loop over transitions
  
  # Compute PFS/OS give parameter estimates
  stprobs <- state_probs(gamma, powers, t, n_sims)
  surv <- data.table(line = 2,
                     mutation = mutation,
                     model = model_lookup(powers),
                     tx_name = tx_name,
                     sample = rep(1:n_sims, each = length(t)),
                     month = rep(t, n_sims),
                     haz_sp = c(stprobs$haz_sp),
                     haz_sd = c(stprobs$haz_sd),
                     haz_pd = c(stprobs$haz_pd),                     
                     pfs = c(stprobs$pfs),
                     os = c(stprobs$os))
  return(surv)
}

# T790M positive (osimertinib)
surv_2L_t790m_osi_est <- list()
surv_2L_t790m_osi_est$wei <- surv_2L(n_months = max_months, 
                                     ma_post = ma_2L_t790m_osi$weibull, 
                                     ma_params_lookup = ma_params_lookup_2L_t790m_osi$weibull,
                                     powers = 0,
                                     "osimertinib",
                                     mutation = 1)
surv_2L_t790m_osi_est$gomp <- surv_2L(n_months = max_months, 
                                      ma_post = ma_2L_t790m_osi$gompertz, 
                                      ma_params_lookup = ma_params_lookup_2L_t790m_osi$gompertz,
                                      powers = 1,
                                      "osimertinib", 
                                      mutation = 1)
surv_2L_t790m_osi_est$fracpoly1 <- surv_2L(n_months = max_months, 
                                           ma_post = ma_2L_t790m_osi$fracpoly1, 
                                           ma_params_lookup = ma_params_lookup_2L_t790m_osi$fracpoly1,
                                           powers = c(0, 0),
                                           "osimertinib",
                                           mutation = 1)
surv_2L_t790m_osi_est$fracpoly2 <- surv_2L(n_months = max_months, 
                                           ma_post = ma_2L_t790m_osi$fracpoly2, 
                                           ma_params_lookup = ma_params_lookup_2L_t790m_osi$fracpoly2,
                                           powers = c(0, 1),
                                           "osimertinib",
                                           mutation = 1)
surv_2L_t790m_osi_est <- rbindlist(surv_2L_t790m_osi_est)

# T790M negative (PBDC)
surv_2L_pbdc_est <- list()
surv_2L_pbdc_est$wei <- surv_2L(n_months = max_months, 
                                ma_post = ma_2L_pbdc$weibull, 
                                ma_params_lookup = ma_params_lookup_2L_pbdc$weibull,
                                powers = 0,
                                "PBDC",
                                 mutation = 0)
surv_2L_pbdc_est$gomp <- surv_2L(n_months = max_months, 
                                 ma_post = ma_2L_pbdc$gompertz, 
                                 ma_params_lookup = ma_params_lookup_2L_pbdc$gompertz,
                                 powers = 1,
                                 "PBDC",
                                 mutation = 0)
surv_2L_pbdc_est$fracpoly1 <- surv_2L(n_months = max_months, 
                                      ma_post = ma_2L_pbdc$fracpoly1, 
                                      ma_params_lookup = ma_params_lookup_2L_pbdc$fracpoly1,
                                      powers = c(0, 0),
                                      "PBDC",
                                      mutation = 0)
surv_2L_pbdc_est$fracpoly2 <- surv_2L(n_months = max_months, 
                                      ma_post = ma_2L_pbdc$fracpoly2, 
                                      ma_params_lookup = ma_params_lookup_2L_pbdc$fracpoly2,
                                      powers = c(0, 1),
                                      "PBDC",
                                      mutation = 0)
surv_2L_pbdc_est <- rbindlist(surv_2L_pbdc_est)

# Combine hazards, PFS, and OS -------------------------------------------------
surv_est <- rbind(surv_1L_est,
                  surv_2L_t790m_osi_est,
                  surv_2L_pbdc_est,
                  fill = TRUE)

# PFS
mstate_nma_pfs <- surv_est[, .(mean = mean(pfs),
                               median = median(pfs),
                               l95 = quantile(pfs, .025),
                               u95 = quantile(pfs, .975)),
                            by = c("line", "mutation", "model", "tx_name", "month")]

# OS
mstate_nma_os <- surv_est[, .(mean = mean(os),
                              median = median(os),
                              l95 = quantile(os, .025),
                              u95 = quantile(os, .975)),
                          by = c("line", "mutation", "model", "tx_name", "month")]

# Hazards
hazard_est <- melt(surv_est, 
                 id.vars = c("line", "mutation", "model", "tx_name", 
                             "sample", "month"),
                 measure.vars = c("haz_sp", "haz_sd", "haz_pd"),
                 variable.name = "transition",
                 value.name = "hazard")
hazard_est[, transition := factor(transition,
                                  levels = c("haz_sp", "haz_sd", "haz_pd"),
                                  labels = c("Stable to progression",
                                             "Stable to death",
                                             "Progression to death"))]
mstate_nma_hazard <- hazard_est[, .(mean = mean(hazard),
                                    median = median(hazard),
                                    l95 = quantile(hazard, .025),
                                    u95 = quantile(hazard, .975)),
                                 by = c("line", "mutation", "model", "tx_name",
                                        "month", "transition")]

# Hazard ratios
hr_est <- melt(surv_est, 
                id.vars = c("line", "mutation", "model", "tx_name", 
                             "sample", "month"),
                 measure.vars = c("hr_sp", "hr_sd", "hr_pd"),
                 variable.name = "transition",
                 value.name = "hr")
hr_est[, transition := factor(transition,
                              levels = c("hr_sp", "hr_sd", "hr_pd"),
                               labels = c("Stable to progression",
                                          "Stable to death",
                                           "Progression to death"))]
mstate_nma_hr <- hr_est[line == 1 & tx_name != "gefitinib",
                        .(mean = mean(hr),
                          median = median(hr),
                          l95 = quantile(hr, .025),
                          u95 = quantile(hr, .975)),
                        by = c("line", "mutation", "model", "tx_name",
                                "month", "transition")]
mstate_nma_hr[, mutation := NULL]

# Save -------------------------------------------------------------------------
save(params_mstate_nma, file = "../data/params_mstate_nma.rda", compress = "bzip2")
save(mstate_nma_pfs, file = "../data/mstate_nma_pfs.rda", compress = "bzip2")
save(mstate_nma_os, file = "../data/mstate_nma_os.rda", compress = "bzip2")
save(mstate_nma_hazard, file = "../data/mstate_nma_hazard.rda", compress = "bzip2")
save(mstate_nma_hr, file = "../data/mstate_nma_hr.rda", compress = "bzip2")
save(mstate_nma_coef, file = "../data/mstate_nma_coef.rda", compress = "bzip2")
save(mstate_ma_coef, file = "../data/mstate_ma_coef.rda", compress = "bzip2")

# Check PFS against JAGS code --------------------------------------------------
# The mu's
jags_v_R <- function(ma_post, line, tx_name, outcome = c("PFS", "OS")){
  
  outcome <- match.arg(outcome)
  
  # JAGS estimates
  surv_jags <- vector(mode = "list", length = n_models)
  for (i in 1:n_models){
    cols <- grep(outcome, colnames(ma_post[[i]]))
    n_months <- length(cols)
    surv_jags[[i]] <- data.table(computation = "JAGS",
                                model = model_lookup(fp_powers[[i]]),
                                sim = rep(1:n_sims, each = n_months),
                                month = rep(1:n_months, n_sims),
                                survival = c(t(ma_post[[i]][, cols]))) 
  }
  surv_jags <- rbindlist(surv_jags)
  surv_jags <- surv_jags[, .(mean = mean(survival)),
                          by = c("computation", "model", "month")]
  
  # R estimates
  tx_name_env <- tx_name
  line_env <- line
  surv_R <- switch(outcome,
                   "PFS" = mstate_nma_pfs[line == line_env & tx_name == tx_name_env],
                   "OS" = mstate_nma_os[line == line_env & tx_name == tx_name_env]
                   )
  surv_R[, computation := "R"]
  
  # Plot
  surv <- rbind(surv_jags, 
               surv_R[, colnames(surv_jags), with = FALSE])  
  y_lab <- switch(outcome,
                  "PFS" = "Progression-free survival",
                  "OS" = "Overall survival"
                  )
  p <- ggplot(surv, aes(x = month, y = mean, col = computation)) + 
        geom_line(position = position_jitter()) +
        facet_wrap(~model) +
        xlab("Month") + ylab(y_lab) +
        scale_color_discrete(name = "") + 
        theme(legend.position = "bottom")
  return(p)
}

## 1st Line
p <- jags_v_R(ma_1L, line = 1, tx_name = "gefitinib", outcome = "PFS")
ggsave("figs/pfs-1L-gef-check.pdf", p, width = 7, height = 7)

p <- jags_v_R(ma_1L, line = 1, tx_name = "gefitinib", outcome = "OS")
ggsave("figs/os-1L-gef-check.pdf", p, width = 7, height = 7)

## Second line (T790m positive, osimertinib)
p <- jags_v_R(ma_2L_t790m_osi, line = 2, tx_name = "osimertinib", outcome = "PFS")
ggsave("figs/pfs-2L-t790m-osi-check.pdf", p, width = 7, height = 7)

p <- jags_v_R(ma_2L_t790m_osi, line = 2, tx_name = "osimertinib", outcome = "OS")
ggsave("figs/os-2L-t790m-osi-check.pdf", p, width = 7, height = 7)

## Second line (PBDC)
p <- jags_v_R(ma_2L_pbdc, line = 2, tx_name = "PBDC", outcome = "PFS")
ggsave("figs/pfs-2L-pbdc-check.pdf", p, width = 7, height = 7)

p <- jags_v_R(ma_2L_pbdc, line = 2, tx_name = "PBDC", outcome = "OS")
ggsave("figs/os-2L-pbdc-check.pdf", p, width = 7, height = 7)

# The mu's and the d's
jags_v_R_d <- function(nma_post, ma_post, econmod_tx_lookup, outcome = c("PFS", "OS")){
  outcome <- match.arg(outcome)
  y_lab <- switch(outcome,
                  "PFS" = "Progression-free survival",
                  "OS" = "Overall survival"
                  )  
  
  # JAGS estimates
  n_tx <- nrow(econmod_tx_lookup)
  p <- vector(mode = "list", length = n_tx)
  names(p) <- econmod_tx_lookup$name
  for (j in 1:nrow(econmod_tx_lookup)){
    num <- econmod_tx_lookup[j, num]
    surv_jags <- vector(mode = "list", length = n_models)
    for (i in 1:n_models){
      cols <- grep(paste0(outcome, "\\[", num), colnames(nma_post[[i]]))
      n_months <- length(cols)
      surv_jags[[i]] <- data.table(computation = "JAGS",
                                      model = model_lookup(fp_powers[[i]]),
                                      tx_name = econmod_tx_lookup[j, name],
                                      sim = rep(1:n_sims, each = n_months),
                                      month = rep(1:n_months, n_sims),
                                      survival = c(t(nma_post[[i]][, cols])))
    } # End loop over models
    surv_jags <- rbindlist(surv_jags)
    surv_jags <- surv_jags[, .(mean = mean(survival)),
                          by = c("computation", "model", "month")]  
    
    # R estimates
    tx_name_env <- econmod_tx_lookup[j, name]
    surv_R <- switch(outcome,
                   "PFS" = mstate_nma_pfs[line == 1 & tx_name == tx_name_env],
                   "OS" = mstate_nma_os[line == 1 & tx_name == tx_name_env]
                   )  
    surv_R[, computation := "R"]
    
    # Plot
    surv <- rbind(surv_jags, 
               surv_R[, colnames(surv_jags), with = FALSE])
    p[[j]] <- ggplot(surv, aes(x = month, y = mean, col = computation)) + 
                    geom_line(position = position_jitter()) +
                    facet_wrap(~model) +
                    xlab("Month") + ylab(y_lab) +
                    scale_color_discrete(name = "") + 
                    theme(legend.position = "bottom",
                          plot.title = element_text(hjust = 0.5)) +
                    scale_y_continuous(breaks = seq(0, 1, by = .2)) +
                    ggtitle(tx_name_env) 
  } # End loop over treatments
  return(p)
}

## PFS
p <- jags_v_R_d(nma_1L, ma_1L, econmod_tx_lookup = econmod_tx_lookup_1L,
                outcome = "PFS")
pdf("figs/pfs-1L-check.pdf")
for (i in 1:length(p)) {
  print(p[[i]])
}
dev.off()

## OS
p <- jags_v_R_d(nma_1L, ma_1L, econmod_tx_lookup = econmod_tx_lookup_1L,
                outcome = "OS")
pdf("figs/os-1L-check.pdf")
for (i in 1:length(p)) {
  print(p[[i]])
}
dev.off()

# Check the coefficients -------------------------------------------------------
# 1L mu's
## Print the mean of the posterior distribution
check_mu_1L <- function(dist = "weibull"){
  cols <- grep("MU", colnames(ma_1L[[dist]]))
  print(apply(ma_1L[[dist]][, cols], 2, mean))
  print(nma_params_lookup_1L[[dist]])
}
check_mu_1L("weibull")
apply(params_mstate_nma$weibull$coefs$a0, 2, mean)

# 1L d's
## Print the mean of the posterior distribution
check_d <- function(dist = "weibull"){
  cols <- grep("d", colnames(nma_1L[[dist]]))
  print(apply(nma_1L[[dist]][, cols], 2, mean))
  print(econmod_tx_lookup_1L)
  print(nma_params_lookup_1L[[dist]])
}

# Weibull
check_d("weibull")
apply(params_mstate_nma$weibull$coefs$a0, 2, mean)
apply(params_mstate_nma$weibull$coefs$a1, 2, mean)




