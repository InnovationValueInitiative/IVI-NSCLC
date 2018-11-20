rm(list = ls())
library("data.table")
library("readxl")
library("MASS")
library("hesim")
library("iviNSCLC")

d_col <- function(i, j){
  res <- paste0("d[", i, ",", j, "]")
  return(res)
}

mu_col <- function(i, j){
  res <- paste0("mu[", i, ",", j, "]")
  return(res)
} 

# Create list of params_surv objects from NMA estimates ------------------------

# Fitted models
mod_names <- c("weibullNMA", "gompertz", "fracpoly1", "fracpoly2")
dists <- c("weibullNMA", "gompertz", "fracpoly", "fracpoly")
n_models <- length(mod_names)
params_mstate_nma <- vector(mode = "list", length = n_models)
names(params_mstate_nma) <- mod_names
param_names <- list(weibullNMA = c("a0", "a1"),
                    gompertz = c("shape", "rate"),
                    fracpoly1 = c("gamma1", "gamma2", "gamma3"),
                    fracpoly2 = c("gamma1", "gamma2", "gamma3"))
aux <- list(weibullNMA = list(powers = 0),
            gompertz = list(powers = 1),
            fracpoly1 = list(powers = c(0, 0)),
            fracpoly2 = list(powers = c(0, 1)))


## Read data
### Parameter lookup files
nma_params_lookup <- vector(mode = "list", length = 2) 
names(nma_params_lookup) <- line_names
for (i in 1:length(nma_params_lookup)){ # Start loop over lines
 filename <- paste0("mstate_nma/params_lookup_", i, "L.xlsx")
 for (j in 1:n_models){ # Start loop over models
    nma_params_lookup[[i]][[j]] <- data.table(read_excel(filename,
                                                  sheet = mod_names[j]))
  } # End loop over models
  names(nma_params_lookup[[i]]) <-  mod_names 
} # End loop over lines

### Treatment lookup files
nma_tx_lookup <- vector(mode = "list", length = 2) 
nma_tx_lookup[[1]] <- fread("mstate_nma/tx_lookup_1L.csv")
nma_tx_lookup[[2]] <- fread("mstate_nma/tx_lookup_2L.csv")

### Parameter estimates
fp <- vector(mode = "list", length = 2) 
fp[[1]] <- list(weibull = as.matrix(fread("mstate_nma/fp_p0_1L.csv")),
                gompertz = as.matrix(fread("mstate_nma/fp_p1_1L.csv")), 
                fracpoly1 = as.matrix(fread("mstate_nma/fp_p00_1L.csv")), # (p1 = 0, p2 = 0)
                fracpoly2 = as.matrix(fread("mstate_nma/fp_p01_1L.csv"))) # (p1 = 0, p2 = 1)
fp[[2]] <- list(weibull = as.matrix(fread("mstate_nma/fp_p0_2L.csv")),
                gompertz = as.matrix(fread("mstate_nma/fp_p1_2L.csv")), 
                fracpoly1 = as.matrix(fread("mstate_nma/fp_p00_2L.csv")), # (p1 = 0, p2 = 0)
                fracpoly2 = as.matrix(fread("mstate_nma/fp_p01_2L.csv"))) # (p1 = 0, p2 = 1)

# Create objects
# There is one for each fitted model
# Within a model, coefficients vary by treatment, transition, and parameter (e.g., scale/shape)

## Transitions
trans_names <- n_trans <- vector(mode = "list", length = 2)
trans_names[[1]] <- c("s1p1", "s1d", "p1d")
trans_names[[2]] <- c("s2p2", "s2d", "p2d")
n_trans[[1]] <- length(trans_names[[1]])
n_trans[[2]] <- length(trans_names[[2]])

## Treatments
tx <- tx_lookup <- vector(mode = "list", length = 2)

### 1L
tx[[1]] <- tx_1L()
row <- match(tx[[1]], nma_tx_lookup[[1]]$tx_name)
tx_lookup[[1]] <- data.table(name = tx[[1]],
                             num = row,
                             abb = nma_tx_lookup[[1]][row, tx_abb])

### 2L
#### All possible 2nd line treatments
tx[[2]] <- vector(mode = "list", length = length(tx[[1]]))
for (i in 1:length(tx_1L)){
  tx[[2]][[i]] <- unlist(tx_2L(tx[[1]][i]))
}
tx[[2]] <- unique(unlist(tx[[2]]))

#### Lookup table
row <- match(tx[[2]], nma_tx_lookup[[2]]$tx_name)
tx_lookup[[2]] <- data.table(name = tx[[2]],
                            num = row,
                            abb = nma_tx_lookup[[2]][row, tx_abb])

#### Container for NMA parameters
alpha <- vector(mode = "list", length = 2)
for (i in 1:2){
  alpha[[i]] <- vector(mode = "list", length = n_models)
  names(alpha[[i]]) <- mod_names
  for (j in 1:n_models){
    alpha[[i]][[j]] <- vector(mode = "list", length = n_trans[[i]])
    names(alpha[[i]][[j]]) <- trans_names[[i]]
    for (k in 1:n_trans[[i]]){
      alpha[[i]][[j]][[k]] <- vector(mode = "list", length = length(tx[[i]]))
      names(alpha[[i]][[j]][[k]]) <- tx[[i]]
      for (l in 1:length(tx[[i]])){
        alpha[[i]][[j]][[k]][[l]] <- matrix(NA, 
                                           nrow = nrow(fp[[i]][[j]]),
                                           ncol = length(param_names[[j]]))
      }
    }
  }
}

## Create object
for (i in 1:n_models){ # Start loop over models (need to expand to 2nd order fractional polynomial models)
  modname_i <- names(params_mstate_nma)[i]
  params_i <- param_names[[i]]
  n_params_i <- length(params_i)
  coefs_i <- vector(mode = "list", length = n_params_i)
  names(coefs_i) <- params_i  
  
  for (j in 1:n_params_i){ # Start loop over parameters of distribution for model i
    coefs_i[[j]] <- vector(mode = "list", length = 2) # For each line of therapy
    
    for (k in 1:2){ # Start loop over line of therapy
      params_lookup_ijk <- nma_params_lookup[[k]][[i]][param == params_i[j]]
      n_cols_k <- n_trans[[k]] * length(tx[[k]])
      mat_ijk <- matrix(0, ncol = n_cols_k, nrow = nrow(fp[[k]][[i]]))
      colnames(mat_ijk) <- rep("tmp", ncol(mat_ijk))

      cntr <- 1
      for (l in 1:n_trans[[k]]){ # Start loop over line specific transitions
        d_num <- params_lookup_ijk[transition == trans_names[[k]][l], d_num] 
        mu_num <-  params_lookup_ijk[transition == trans_names[[k]][l], mu_num] 
        
        for (m in 1:length(tx[[k]])) { # Start loop over treatments
          tx_num <- tx_lookup[[k]][m, num]
          
          if (tx_lookup[[k]][m, name] != "osimertinib"){
            if (!is.na(d_num)){ # d is only nonzero for transitions estimated by NMA
              mat_ijk[, cntr] <- fp[[k]][[i]][, d_col(tx_num, d_num)]
            } 
            var_name <- paste0("d_", tx_lookup[[k]][m, abb], "_",
                              trans_names[[k]][l], "_",
                              params_i[j])
           alpha[[k]][[i]][[l]][[m]][, j] <- fp[[k]][[i]][, mu_col(12, mu_num)] +
                                              mat_ijk[, cntr]
          } else{
            mat_ijk[, cntr] <- fp[[k]][[i]][, mu_col(12, mu_num)]
            var_name <- paste0(tx_lookup[[k]][m, abb], "_",
                              trans_names[[k]][l], "_",
                              params_i[j])     
           alpha[[k]][[i]][[l]][[m]][, j] <- fp[[k]][[i]][, mu_col(12, mu_num)]
          }           
          colnames(mat_ijk)[cntr] <- var_name
          cntr <- cntr + 1          
        } # End loop over treatments
        
      } # End loop over transitions
      coefs_i[[j]][[k]] <- mat_ijk
    
    } # End loop over line of therapy  
    coefs_i[[j]] <- do.call("cbind", coefs_i[[j]])
    
  } # End loop over parameters
  params_mstate_nma[[i]] <- params_surv(coefs = coefs_i, dist = dists[i],
                                        aux = aux[[i]])
} # End loop over models

# Compute PFS/OS ---------------------------------------------------------------
tdata <- function(time, powers){
  basis_power <- function(x, p){
    if (p == 0){
      return (log(x))
    } else{
      return (x^p)
    }    
  }
  
  var <- matrix(NA, nrow = length(time), ncol = length(powers) + 1)
  var[, 1] <- 1
  var[, 2] <- basis_power(time, powers[1])
  if (length(powers) > 1){
    xp_old = var[, 2];
    for (i in 2:(length(powers))){
      if (powers[i] == powers[i - 1]){
        xp_new <- log(time) * xp_old
      } else{
        xp_new <- basis_power(time, powers[i])
      }
      var[, i + 1] <- xp_new
      xp_old <- xp_new
    }
  }
  return(var)
}

month_time <- seq(1, 4 * 12)
pfs <- vector(mode = "list", length = 2)
pfs[[1]] <- vector(mode = "list", length = length(tx[[1]]))
names(pfs[[1]]) <- tx[[1]]
pfs[[2]] <- vector(mode = "list", length = length(tx[[2]]))
names(pfs[[2]]) <- tx[[2]]

for (i in 1:length(alpha)){ # Start loop over lines of therapy
  pfs[[i]] <- vector(mode = "list", length = n_models)  
  names(pfs[[i]]) <- mod_names
  for (j in 1:n_models){ # Start loop over models
    pfs[[i]][[j]] <- vector(mode = "list", length = length(tx[[i]])) 
    names(pfs[[i]][[j]]) <-tx[[i]]
    tdata_j <- tdata(time = month_time, aux[[j]]$powers)
    for (k in 1:length(tx[[i]])){
      hazard_sp <- exp(alpha[[i]][[j]][[1]][[k]] %*% t(tdata_j)) # Rows are simulations, columns are time periods
      hazard_sd <- exp(alpha[[i]][[j]][[2]][[k]] %*% t(tdata_j)) 
      hazard_pd <- exp(alpha[[i]][[j]][[3]][[k]] %*% t(tdata_j)) 
      cumhazard_sp <- t(apply(hazard_sp, 1, cumsum)) 
      cumhazard_sd <- t(apply(hazard_sd, 1, cumsum))   
      pfs[[i]][[j]][[k]] <- exp(-(cumhazard_sp + cumhazard_sd))
      pfs[[i]][[j]][[k]] <- cbind(1, pfs[[i]][[j]][[k]])
      # Let's take the discrete cumulative hazard here
    }
  } # End loop over models
} # End loop over lines of therapy

post <- apply(fp[[1]]$weibull, 2, median)
beta1 <- post[mu_col(12, 1)] + post[d_col(2, 1)]
beta2 <- post[mu_col(12, 2)] + post[d_col(2, 2)]
beta3 <- post[mu_col(12, 3)] 
beta4 <- post[mu_col(12, 4)] 
beta5 <- post[mu_col(12, 5)] + post[d_col(2, 3)]
beta6 <- post[mu_col(12, 6)]

tmp <- rweibullNMA(1000, a0 = post[mu_col(12, 1)], a1 =  post[mu_col(12, 2)])
summary(tmp)

t <- seq(1, 10)
haz_sp <- exp(beta1 + beta2 * log(t))
haz_sd <- exp(beta3 + beta4 * log(t))
haz_pd <- exp(beta5 + beta6 * log(t))

pfs <- rep(NA, length(t) + 1)
pfs[1] <- 1
for (i in 1:length(t)){
  pfs[i + 1] <- pfs[i] * exp(-(haz_sp[i] + haz_sd[i]))
}

# Superseded -------------------------------------------------------------------

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

