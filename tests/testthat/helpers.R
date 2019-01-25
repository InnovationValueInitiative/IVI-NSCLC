example_IndivCtstm <- function(n_samples = 2, n_patients = 3,
                               dr_qalys = .03, dr_costs = .03){
  
  # Sample patients
  pats <- create_patients(n = n_patients)
  
  # Create treatment sequences
  txseq1 <- txseq(first = "gefitinib",
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
  txseq2 <- txseq(first = "erlotinib",
                  second = c("osimertinib", "PBDC"),
                  second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
  txseqs <- txseq_list("Sequence 1" = txseq1, "Sequence 2" = txseq2, start_line = "first",
                       mutation = "unknown")  
  
  # Model structure
  struct <- model_structure(txseqs, dist = "weibull", n_states = "three")
  states <- create_states(struct)
  tmat <- create_trans_mat(struct)
  
  # Construct model
  ## Transition model
  transmod_data <- create_transmod_data(struct, tmat, pats)
  transmod_params <- create_transmod_params(n = n_samples, data = transmod_data)
  transmod <- create_transmod(params = transmod_params, data = transmod_data)
  
  ## Adverse events
  ae_probs <- ae_probs(n_samples, struct)
  
  ## Utility model
  utilmod <- create_utilmod(n = n_samples, struct = struct, patients = pats,
                          ae_probs = ae_probs)
  
  ## Cost models
  costmods <- create_costmods(n = n_samples, struct = struct, patients = pats, ae_probs = ae_probs)
  
  ## Create economic model
  econmod <- hesim::IndivCtstm$new(trans_model = transmod,
                                 utility_model = utilmod,
                                 cost_models = costmods)
  
  ## Simulate economic model
  max_yrs <- 20
  max_age <- 100
  econmod$sim_disease(max_t = max_yrs * 12, max_age = max_age * 12)
  econmod$disprog_[, ':=' (time_start = time_start/12, time_stop = time_stop/12)]
  econmod$sim_qalys(dr = dr_qalys)
  econmod$sim_costs(dr = dr_costs)
  return(econmod)
}