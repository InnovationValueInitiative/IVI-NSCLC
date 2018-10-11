#' Create cost models
#' 
#' Create a cost model for four cost categories: (i) treatment costs (i.e., drug
#' acquisition and administration costs), (ii) inpatient medical costs, 
#' (iii) outpatient medical costs, and (iv) costs due to adverse events.
#' @param n The number of random observations of the parameters to draw.
#' @param struct A \code{\link{model_structure}} object.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @return An object of class "StateVals" from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package.
#' @param params_costs_op Parameter estimates for outpatient medical costs
#' in the same format as \code{\link{params_costs_op}}.
#' @param params_costs_inpt Parameter estimates for inpatient medical costs
#' in the same format as \code{\link{params_costs_inpt}}.
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
#' ## Cost models
#' costmods <- create_costmods(n = 2, struct = struct, patients = pats)
#' @export
create_costmods <- function(n = 100, struct, patients,
                            params_costs_op = iviNSCLC::params_costs_op,
                            params_costs_inpt = iviNSCLC::params_costs_inpt){
  costmods <- list()
  costmods$op <- create_costmods_default(n, struct, patients,
                                        params_costs_op)
  costmods$inpt <- create_costmods_default(n, struct, patients,
                                          params_costs_inpt)
  return(costmods)
}

create_costmods_default <- function(n = 100, 
                                    struct, patients, 
                                    params){
  
  strategies <- data.table(strategy_id = 1:length(struct$txseqs))
  hesim_dat <- hesim::hesim_data(strategies = strategies,
                                  patients = patients)
  
  states <- create_states(struct)[get("state_name") != "D"]
  tbl <- merge(states, params, by = "state_name")  
  tbl <- hesim::stateval_tbl(tbl, dist = "gamma", hesim_data = hesim_dat)
  
  mod <- hesim::create_StateVals(tbl, n = n)
  
  return(mod)  
}