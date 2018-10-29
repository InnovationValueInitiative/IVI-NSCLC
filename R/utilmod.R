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
  strategies <- data.table(strategy_id = 1:length(struct$txseqs))
  hesim_dat <- hesim::hesim_data(strategies = strategies,
                                 patients = patients)
  
  states <- create_states(struct)[get("state_name") != "D"]
  utility_tbl <- merge(states, params_utility$state_utility, by = "state_name")
  utility_tbl <- hesim::stateval_tbl(utility_tbl, dist = "beta", hesim_data = hesim_dat)
  
  utilmod <- hesim::create_StateVals(utility_tbl, n = n)

  return(utilmod)
}