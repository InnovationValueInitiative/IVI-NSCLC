context("costmods.R unit tests")
rm(list = ls())

as_data_table <- function(x){
  x2 <- data.table(mu = x$params$mu[, 1],
                   strategy_id = x$input_mats$strategy_id,
                   patient_id = x$input_mats$patient_id,
                   state_id = x$input_mats$state_id,
                   time_id = x$input_mats$time_id)
  return(x2)
}

get_tx_costs <- function(x, agents, 
                         type = c("acquisition_costs", "administration_costs")){
  type <- match.arg(type)
  n_agents <- length(agents)
  costs <- rep(NA, n_agents)
  for (i in 1:n_agents){
    indx <- match(agents[i], x$agent_name) 
    costs[i] <- unlist(c(x[indx, type, with = FALSE]))
  }
  return(sum(costs))
}

get_agents <- function(lookup, tx_name){
  tx_name_env <- tx_name
  colnums <- grep("agent_name", colnames(lookup))
  agents <- unlist(c(lookup[tx_name == tx_name_env, colnums, with = FALSE]))
  return(agents[!is.na(agents)])
}

# Function to check treatment costs assuming discount rates are set to 0. 
check_tx_costs <- function(costmods, 
                           type = c("acquisition_costs", "administration_costs")){
  type <- match.arg(type)
  if (type == "acquisition_costs"){
    tx_dt <- as_data_table(costmods$tx_ac)
    expect_true(all(apply(costmods$tx_ac$params$mu, 1, sd) == 0)) # Since discount rates are set to 0
  } else{
    tx_dt <- as_data_table(costmods$tx_admin)
    expect_true(all(apply(costmods$tx_admin$params$mu, 1, sd) == 0)) # Since discount rates are set to 0
  }
    
  ## Costs in state 1
  costs <- get_tx_costs(ann_costs, txseq1$first, type)
  expect_true(all(tx_dt[strategy_id == 1 & state_id == 1, mu] == costs))
  
  ## Costs in state 2 (T790M+)
  costs <- get_tx_costs(ann_costs, txseq1$second["pos"], type)
  expect_true(all(tx_dt[patient_id %in% pat_mut1$patient_id & strategy_id == 1 &
                        state_id == 2]$mu == costs))
  
  ## Costs in state 2 (T790M-)
  agents <- get_agents(params_costs_tx2$lookup, txseq2$second["neg"])
  costs <- get_tx_costs(ann_costs, agents, type)
  expect_true(all(tx_dt[patient_id %in% pat_mut0$patient_id & strategy_id == 2 &
                          state_id == 2 & 
                        time_id == 1]$mu == costs))
  
  ## Costs in state 3 (T790M+)
  agents <- get_agents(params_costs_tx2$lookup, txseq2$second_plus["pos"])
  costs <- get_tx_costs(ann_costs, agents, type)
  pat_mut0 <- pats[mutation == 1]
  expect_true(all(tx_dt[patient_id %in% pat_mut0$patient_id & strategy_id == 2 &
                        state_id == 3 & 
                        time_id == 1]$mu == costs))    
}

# Model setup
## Treatment sequences
txseq1 <- txseq(first = "erlotinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
txseq2 <- txseq(first = "gefitinib",
               second = c("osimertinib", "PBDC"),
               second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))

# Patient population
pats <- create_patients(n = 4)
pat_mut1 <- pats[mutation == 1]
pat_mut0 <- pats[mutation == 0]

# Treatment costs with no discounts
params_costs_tx2 <- copy(iviNSCLC::params_costs_tx)
params_costs_tx2$discounts[, discount_lower := 0]
params_costs_tx2$discounts[, discount_upper := 0]
ann_costs <- annualized_tx_costs(params_costs_tx2)

test_that("create_costmods first line, 4 health states", {
  n_samples <- 5
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull")
  ae_probs <- ae_probs(n_samples, struct)
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats,
                              ae_probs = ae_probs)
  expect_true(inherits(costmods, "list"))
  expect_true(all(sapply(costmods, function(x) inherits(x, "StateVals"))))
  
  # Check treatment costs
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats,
                              ae_probs = ae_probs, params_costs_tx = params_costs_tx2)
  check_tx_costs(costmods, "acquisition_costs")
  check_tx_costs(costmods, "administration_costs")
})

test_that("create_costmods first line, 3 health states", {
  n_samples <- 5
  txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
  struct <- model_structure(txseqs, dist = "weibull",
                            n_states = "three")
  ae_probs <- ae_probs(n_samples, struct)
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats,
                              ae_probs = ae_probs)
  expect_true(inherits(costmods, "list"))
  expect_true(all(sapply(costmods, function(x) inherits(x, "StateVals"))))
  
  # Check treatment costs
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats,
                              ae_probs = ae_probs, params_costs_tx = params_costs_tx2) 
  check_tx_costs(costmods, "acquisition_costs")
  check_tx_costs(costmods, "administration_costs")  
})