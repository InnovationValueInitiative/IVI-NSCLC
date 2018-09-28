#' Create cost models
#' 
#' Create a cost model for four cost categories: (i) treatment costs (i.e., drug
#' acquisition and administration costs), (ii) inpatient medical costs, 
#' (iii) outpatient medical costs, and (iv) costs due to adverse events.
#' @param params A "sampled_params" object returned from \code{\link{sample_params}}.
#' @param struct A \code{\link{model_structure}} object.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @return An object of class "StateVals" from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package.
#' @examples
#' # Parameters
#' params <- sample_params(n = 2)
#' 
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
#' costmods <- create_costmods(params, struct, pats)
#' @export
create_costmods <- function(params, struct, patients){
  costmods <- list()
  costmods$op <- create_costmods_default(params, struct, patients,
                                          "costs_op")
  costmods$inpt <- create_costmods_default(params, struct, patients,
                                          "costs_inpt")
  return(costmods)
}

create_costmods_default <- function(params, struct, patients, 
                                    category = c("costs_op", "costs_inpt")){
  category <- match.arg(category)
  
  # Select health states based on model structure
  if (struct$n_states == "four"){
    costs <- params[[category]]
  } else if (struct$n_states == "three"){
    if (attributes(struct$txseqs)$start_line == "first"){
      costs <- params[[category]][, 1:2]
    } else {
      costs <- params[[category]][, 2:3]
    }
  }
  
  # Create "stateval_means" object
  n_strategies <- length(struct$txseqs)
  cost_means <- hesim::stateval_means(values = costs,
                                         strategy_id = 1:n_strategies,
                                         patient_id = patients$patient_id)
  
  # Create "StateVals" object
  costmod <- hesim::create_StateVals(cost_means)
  
  return(costmod)  
}