example_IndivCtstm <- function(n_samples = 2, n_patients = 3){
  
  # Treatment strategies, target population, and model structure
  n_strategies <- 2
  strategies <- data.frame(strategy_id = 1:n_strategies)
  patients <- create_patients(n = n_patients)
  states <- data.frame(state_id = c(1, 2))
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients,
                          states = states)
  tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
  
  # Parameterization
  ## Multi-state model
  surv_dat <- data.table(ctstm3_exdata$transitions)
  fits <- vector(length = max(tmat, na.rm = TRUE), mode = "list")  
  for (i in 1:length(fits)){
    fits[[i]] <- flexsurvreg(Surv(years, status) ~ factor(strategy_id), 
                             data = surv_dat,
                             subset = (trans == i),
                             dist = "weibull")
  }  
  
  ## Utility
  utility_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       mean = ctstm3_exdata$utility$mean,
                                       se = ctstm3_exdata$utility$se),
                            dist = "beta",
                            hesim_data = hesim_dat)
  ## Costs
  drugcost_tbl <- stateval_tbl(data.frame(strategy_id = strategies$strategy_id,
                                         est = ctstm3_exdata$costs$drugs$costs),
                              dist = "fixed",
                              hesim_data = hesim_dat) 
  medcost_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                         mean = ctstm3_exdata$costs$medical$mean,
                                         se = ctstm3_exdata$costs$medical$se),
                              dist = "gamma",
                              hesim_data = hesim_dat)    
  
  # Simulation
  ## Construct economic model
  ### Transitions
  transmod_data <- expand(hesim_dat)  
  transmod <- create_IndivCtstmTrans(flexsurvreg_list(fits), input_data = transmod_data,
                                     trans_mat = tmat,
                                     n = n_samples)  
  
  ### Utility 
  utilitymod <- create_StateVals(utility_tbl, n = n_samples)
  
  ### Costs
  drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
  medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
  costmods <- list(drugs = drugcostmod,
                   medical = medcostmod)
  
  ### Combine
  ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilitymod,
                         cost_models = costmods)


  ## Simulate outcomes
  head(ictstm$sim_disease()$disprog_)
  head(ictstm$sim_stateprobs(t = c(0, 5, 10))$stateprobs_[t == 5])
  ictstm$sim_qalys(dr = .03)
  ictstm$sim_costs(dr = .03)
  return(ictstm)
}