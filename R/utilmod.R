#' Create utility model
#' 
#' Create a model for health state utility given
#' @param params A "sampled_params" object returned from \code{\link{sample_params}}.
#' @param struct A \code{\link{model_structure}} object.
#' @param patients A data table returned from \code{\link{create_patients}}.
#' @return An object of class "StateVals" from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package.
#' @examples
#' ## Treatment sequences
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
#' utilmod <- create_utilmod(params, struct, pats)
#' @export
create_utilmod <- function(params, struct, patients){
  # Select health states based on model structure
  if (struct$n_states == "four"){
    utility <- params$utility
  } else if (struct$n_states == "three"){
    if (attributes(struct$txseqs)$start_line == "first"){
      utility <- params$utility[, 1:2]
    } else {
      utility <- params$utility[, 2:3]
    }
  }
  
  # Create "stateval_means" object
  n_strategies <- length(struct$txseqs)
  utility_means <- hesim::stateval_means(values = utility,
                                         strategy_id = 1:n_strategies,
                                         patient_id = patients$patient_id)
  
  # Create "StateVals" object
  utilmod <- create_StateVals(utility_means)
  
  return(utilmod)
}