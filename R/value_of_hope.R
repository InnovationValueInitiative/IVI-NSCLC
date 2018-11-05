#' Value of hope
#' 
#' Compute the quality-adjusted life-years (QALYs) that a patient would need to obtain
#' to be indifferent between treatment sequences relative to a reference treatment 
#' sequence (i.e., compute the "certainty equivalent"), give differences in the
#' distribution of QALYs.
#' @param econmod An economic model of class \code{"IndivCtstm"}. Disease progression
#' must have been previously simulated (i.e., \code{$disprog_} cannot be \code{NULL}.)
#' @param comparator The \code{strategy_id} from \code{econmod} to use as the comparator.
#' @param crra Constant relative risk aversion parameter. 
#' @param dr Discount rate.
#' @return A \code{data.table} with columns:
#' \describe{
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{comparator}{Equal to 1 if the treatment strategy is the comparator and 0 otherwise.}
#' \item{qalys}{Mean QALYs.}
#' \item{iqalys}{Incremental mean QALYs; that is, mean QALYs relative to the comparator.}
#' \item{ce}{The certainty equivalent.}
#' \item{voh}{The 'value of hope'. See 'Details'.}
#' }
#' @details The value of hope is the difference between the certainty equivalent and 
#' differences in mean QALYs between a given treatment strategy and a comparator. That is, 
#' if \eqn{\alpha} is the certainty equivalent, then the value of hope for treatment strategy
#' 2 relative to treatment strategy 1 is
#' \deqn{\alpha - (E[qalys_2] - E[qalys_1])}.
#' @seealso See the example in the \href{https://innovationvalueinitiative.github.io/IVI-NSCLC/articles/tutorial.html}{tutorial}.
#' @export
value_of_hope <- function(econmod, comparator, crra = .39, dr = .03){
  
  # Checks
  check_is_class(econmod, name = "econmod", class = "IndivCtstm")
  
  # Prevent CRAN warnings
  prob <- strategy_id <- qalys <- iqalys <- ce <- qalys_comparator <- voh <- NULL
  
  # Functions
  ## Utility function
  utility_fun <- function(x, r, alpha = 0){
    m <- x - alpha
    res <- ifelse(m <= 0, 0, (m)^r)
    return(res)
  }
  
  ## Certainty equivalent function
  certainty_equivalent_fun <- function(alpha, r, qalys, prob, eu1){
    return(sum(utility_fun(qalys, r, alpha) * prob) - eu1)
  }
  
  # Simulate QALYs by patient
  dr_env <- dr
  econmod2 <- econmod$clone(deep = TRUE)
  econmod2$sim_qalys(dr = dr_env, by_patient = TRUE, lys = FALSE)
  sim <- econmod2$qalys_[, lapply(.SD, sum),
                           .SDcols = "qalys",
                           by = c("sample", "strategy_id", "patient_id", "dr")]  
  
  # Compute value of hope
  strategy_ids <- sort(unique(sim$strategy_id))
  n_strategies <- length(strategy_ids)
  sim[, prob := 1/.N, by = "strategy_id"] # Empirical PDFs
  
  ## Expected QALYs 
  res <- sim[, lapply(.SD, mean),
                      .SDcols = "qalys",
                       by = c("strategy_id")]
  
  ## Expected utility for comparator
  sim_comparator <- sim[strategy_id == comparator]
  eu1 <- sum(sim_comparator$prob * utility_fun(sim_comparator$qalys, r = 1 + crra))
  
  ## Certainty equivalent
  sim_treatments <- sim[strategy_id != comparator]
  certainty_equivalent <- rep(NA, n_strategies)
  for (i in 1:n_strategies){
    if (strategy_ids[i] == comparator){
      certainty_equivalent[i] <- 0
    } else{
      sim_treatments_i <- sim_treatments[strategy_id == strategy_ids[i]]
      root_list <- stats::uniroot(certainty_equivalent_fun, 
                                c(-10, 10),
                                extendInt = "yes",
                                r = 1 + crra,
                                qalys = sim_treatments_i$qalys,
                                prob = sim_treatments_i$prob,
                                eu1 = eu1)  
      certainty_equivalent[i] <- root_list$root
    }
  }
  
  ## Value of hope
  mean_qalys_comparator <- res[strategy_id == comparator, qalys]
  res[, qalys_comparator := mean_qalys_comparator]
  res[, iqalys := qalys - qalys_comparator]
  res$ce <-  certainty_equivalent
  res[, voh := ce - (qalys - qalys_comparator)]
  res[, qalys_comparator := NULL]
  res[, ("comparator") := ifelse(strategy_id == comparator, 1, 0)]
  
  
  # Return
  setcolorder(res, c("strategy_id", "comparator", "qalys", "iqalys", "ce", "voh"))
  return(res[,])
  
}