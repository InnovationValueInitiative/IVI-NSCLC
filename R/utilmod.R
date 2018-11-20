#' Create utility model
#' 
#' Create a model for health state utility given values of utility by
#'  health state, treatment, and time sampled from a probability distribution. 
#' @param n The number of random observations of the parameters to draw.
#' @param struct A \code{\link{model_structure}} object.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @param ae_probs An "ae_probs" object as returned by \code{\link{ae_probs}}.
#' @param params_utility Parameter estimates for health state utilities and
#' adverse event disutilities in the same format as \code{\link{params_utility}}.
#' @return An object of class "StateVals" from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package.
#' @examples
#' # Treatment sequences
#' txseq1 <- txseq(first = "erlotinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseq2 <- txseq(first = "gefitinib",
#'                 second = c("osimertinib", "PBDC"),
#'                 second_plus = c("PBDC + bevacizumab", "PBDC + bevacizumab"))
#' txseqs <- txseq_list(seq1 = txseq1, seq2 = txseq2)
#'
#' # Patient population
#' pats <- create_patients(n = 2)
#'
#' ## Model structure
#' struct <- model_structure(txseqs, dist = "weibull")
#'
#' ## Utility model
#' n_samples <- 2
#' ae_probs <- ae_probs(n = n_samples, struct = struct)
#' utilmod <- create_utilmod(n = n_samples, struct = struct, patients = pats,
#'                           ae_probs = ae_probs)
#' @export
create_utilmod <- function(n = 100, struct, patients,
                           ae_probs,
                           params_utility = iviNSCLC::params_utility){
  
  # hesim data
  strategies <- data.table(strategy_id = 1:length(struct$txseqs))
  hesim_dat <- hesim::hesim_data(strategies = strategies,
                                 patients = patients)
  # Utility by health state
  states <- create_states(struct)[get("state_name") != "D"]
  state_utility_tbl <- merge(states, params_utility$state_utility, by = "state_name")
  setorderv(state_utility_tbl, "state_id")
  state_utility_dist <- matrix(stats::rnorm(n * nrow(state_utility_tbl),
                                state_utility_tbl$mean, 
                                state_utility_tbl$se), 
                                nrow = n, byrow = TRUE)
  
  # Disutilities from adverse events
  disutility_ae_dist <- matrix(stats::rnorm(n * nrow(params_utility$ae_disutility),
                                params_utility$ae_disutility$mean, 
                                params_utility$ae_disutility$se), 
                                nrow = n, byrow = TRUE)
  colnames(disutility_ae_dist) <- params_utility$ae_disutility$ae_abb
  
  ## Weight by adverse events
  indices <- match(params_utility$ae_disutility$ae_abb, names(ae_probs))
  if (any(is.na(indices))){
    stop(paste0("The adverse event abbreviations in 'params_utility$ae_disutility' do not match ",
                "the names of the adverse events in 'ae_probs'."))
  }
  expected_disutility <- vector(mode = "list", length = length(ae_probs))
  names(expected_disutility) <- names(ae_probs)
  for (i in 1:length(ae_probs)){
    name_i <- names(ae_probs)[i]
    expected_disutility[[i]] <- ae_probs[[i]] * disutility_ae_dist[, name_i]
  }
  expected_disutility <- Reduce('+', expected_disutility)  
  
  # Create utility model
  tbl1 <- tbl2 <- vector(mode = "list", length = nrow(strategies))
  for (i in 1:nrow(strategies)){
    tbl1[[i]] <- data.table(strategy_id = i,
                            state_id = rep(states$state_id, each = n),
                            sample = rep(1:n, times = nrow(states)),
                            value = c(state_utility_dist) - expected_disutility[, i],
                            time_start = 0)
    tbl2[[i]] <- data.table(strategy_id = i,
                            state_id = rep(states$state_id, each = n),
                            sample = rep(1:n, times = nrow(states)),
                            value = c(state_utility_dist),
                            time_start = 1/12)    
  }
  tbl1 <- rbindlist(tbl1)
  tbl2 <- rbindlist(tbl2)
  tbl <- rbind(tbl1, tbl2)
  tbl <- hesim::stateval_tbl(tbl, dist = "custom", hesim_data = hesim_dat)
  utilmod <- hesim::create_StateVals(tbl, n = n)  
  return(utilmod)
}