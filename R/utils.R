#' Input validation for class objects
#' 
#' \code{check} is a generic function for validating the inputs of class objects.
#' @param object object to check.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return If validation is successful, returns the object in question; otherwise,
#' informs the user that an error has occurred.  
#' @keywords internal
check <- function (object, ...) {
  UseMethod("check")
}

check_is_class <- function(object, name, class){
  if (!inherits(object, class)){
      stop(paste0("'", name, "' must be of class '", class, "'"),
         call. = FALSE)
  }  
}

#' Summarize model outcomes
#' 
#' Summarize clinical and economic outcomes from the simulation.
#' @param econmod An economic model of class "IndivCtstm" from the
#' \href{https://hesim-dev.github.io/hesim/}{hesim} package. 
#'  QALYs and costs must have been previously simulated (i.e., $qalys_ and 
#'  $costs_ cannot be NULL).
#' @param prod_costs An object of class "prod_costs" simulated using 
#' \code{\link{sim_prod_costs}}.
#' @param dr_qalys Discount rate for QALYs.
#' @param dr_costs Discount rate for costs.
#' @param wtp Willingness to pay for a QALY.
#' @param digits_qalys Number of digits to use to report QALYs.
#' @param digits_costs Number of digits to use to report costs.
#' @param prob A numeric scalar in the interval \code{(0,1)} giving the credible interval.
#' Default is 0.95 for a 95 percent credible interval. 
#' @param strategy_names A character vector denoting names of treatment strategies.
#' 
#' @return A \code{data.table} summarizing model outcomes. 
#' @seealso See the example in the \href{https://innovationvalueinitiative.github.io/IVI-NSCLC/articles/tutorial.html}{tutorial}.
#' @export
summarize_outcomes <- function (econmod, prod_costs = NULL,
                                dr_qalys, dr_costs, 
                                wtp = 150000, digits_qalys = 2, digits_costs = 0,
                                prob = 0.95, strategy_names = NULL) {
  
  value <- dr <- lys <- qalys <- outcome <- costs <-  NULL
 
  # Checks
  check_is_class(econmod, "econmod", "IndivCtstm")
  if (is.null(econmod$qalys_)){
    stop("You must first simulate QALYs with 'econmod' using '$sim_qalys()'")
    if (!"lys" %in% colnames(econmod$qalys_)){
      stop("You must simulate life-years with '$sim_qalys()'.")
    }
  }
  if (is.null(econmod$costs_)){
    stop("You must first simulate costs with 'econmod' using '$sim_costs()'")
  }  
  if (!is.null(prod_costs)){
    check_is_class(prod_costs, "prod_costs", "prod_costs")
  }
  if (prob > 1 | prob < 0){
    stop("'prob' must be in the interval (0,1)",
         call. = FALSE)
  }

  # Summarize outcomes
  ## Health outcomes
  qalys_summary <- econmod$qalys_[dr == dr_qalys, 
                                   list(lys = sum(lys),
                                    qalys = sum(qalys)),
                                   by = c("sample", "strategy_id")]
  health_summary <- melt(qalys_summary,
                         id.vars = c("sample", "strategy_id"),
                         variable.name = "outcome")
  health_summary[, outcome := factor(outcome, 
                                     levels = c("lys", "qalys"),
                                     labels = c("Life-years", "QALYs"))]
  
  ## Health care sector costs
  cost_summary <- econmod$costs_[dr == dr_costs, 
                                  list(value = sum(costs)),
                                   by = c("category", "sample", "strategy_id")] 
  setnames(cost_summary, "category", "outcome")
  cost_summary[, outcome := factor(outcome, 
                                   levels = c("tx_ac", "tx_admin", "op", 
                                              "inpt", "ae"),
                                   labels = c("Drug acquisition costs",
                                              "Drug administration costs",
                                              "Outpatient medical costs",
                                              "Inpatient medical costs",
                                              "Adverse event costs"))]  
  cost_summary_total <- cost_summary[, list(value = sum(value)),
                                     by = c("sample", "strategy_id")]
  cost_summary_total[, outcome := "Health care sector costs"]
  cost_summary <- rbind(cost_summary, cost_summary_total)
  
  ## Productivity costs 
  if (!is.null(prod_costs)){
    prod_costs_summary <- prod_costs[, c("sample", "strategy_id", "costs", "category")]
    setnames(prod_costs_summary, c("category", "costs"), c("outcome","value"))
    prod_costs_summary[, outcome := "Productivity costs"] 
    cost_summary <- rbind(cost_summary, prod_costs_summary) 
    cost_societal_summary <- cost_summary[outcome == "Health care sector costs" |
                                          outcome == "Productivity costs",
                                          list(value = sum(value)),
                                          by = c("sample", "strategy_id")]
   cost_societal_summary[, outcome := "Societal costs"]
   cost_summary <- rbind(cost_summary, cost_societal_summary)  
  }

  ## Net monetary benefit
  if (!is.null(prod_costs)){
    nmb <- qalys_summary$qalys * wtp - cost_summary[outcome == "Societal costs"]$value
  } else {
    nmb <- qalys_summary$qalys * wtp - cost_summary[outcome == "Health care sector costs"]$value
  }
  nmb <- cbind(qalys_summary[, c("sample", "strategy_id")], nmb)
  nmb[, outcome := "Net monetary benefit"]
  setnames(nmb, "nmb", "value")
  cost_summary <- rbind(cost_summary, nmb)
  
  # Table
  prob_lower <- (1 - prob)/2
  prob_upper <- 1 - prob_lower
  
  format_costs <- function(x, digits){
    formatC(x, format = "f", digits = digits, big.mark = ",")
  }
  
  format_qalys <- function(x, digits){
    formatC(x, format = "f", digits = digits)
  }
  
  format_cri <- function(est, lower, upper, costs = TRUE, digits){
    if (costs){
      est <- format_costs(est, digits = digits)
      lower <- format_costs(lower, digits = digits)
      upper <- format_costs(upper, digits = digits)
    } else{
      est <- format_qalys(est, digits = digits)
      lower <- format_qalys(lower, digits = digits)
      upper <- format_qalys(upper, digits = digits)
    }
    paste0(est, " (",lower, ", ", upper, ")")
  }  
  
  format_tbl <- function(x, costs = TRUE, digits){
    res <- x[, list(mean = mean(value),
                 lower = stats::quantile(value, prob_lower),
                 upper = stats::quantile(value, prob_upper)),
             by = c("outcome", "strategy_id")]
    res$value <- format_cri(res$mean, res$lower, res$upper,
                       costs = costs, digits = digits)
    return(dcast(res, outcome ~strategy_id, value.var = "value"))
  }
  health_summary <- format_tbl(health_summary, costs = FALSE, digits = digits_qalys)
  cost_summary <- format_tbl(cost_summary, costs = TRUE, digits = digits_costs)
  summary_tbl <- rbind(health_summary, cost_summary)
  colnames(summary_tbl)[1] <- "Outcome"
  if (!is.null(strategy_names)){
   colnames(summary_tbl)[-1] <- strategy_names 
  }
  return(summary_tbl)
}

#' Tidy data
#' 
#' A generic function for creating tidy data.
#' @param object An object to create tidy data from.
#' @param ... Further arguments passed to or from other methods. 
#' @export
#' @keywords internal
#' @seealso \code{\link{tidy.ae_probs}}
tidy <- function(object, ...){
  UseMethod("tidy", object)
}

pv <- function(z, r, t1, t2){
  if (r == 0){
    return (z * (t2 - t1))
  } else{
    return (z * ((exp(-r * t1) - exp(-r * t2))/r))
  }  
}

tx_by_state <- function(struct){
  line <- state_id <- NULL
  
  n_txseqs <- length(struct$txseqs)
  txseq_dt <- vector(mode = "list", length = n_txseqs)
  for (i in 1:n_txseqs){
    txseq_i <- unlist(struct$txseqs[[i]])
    txseq_dt[[i]] <- data.table(strategy_id = i,
                                tx_name = c(txseq_i[1], txseq_i),
                                line = rep(c("first", "second", "second_plus"), 
                                           each = 2),
                                mutation = rep(c(1, 0), 3),
                                grp_id = rep(c(1, 2), 3)) 
  }
  txseq_dt <- rbindlist(txseq_dt)
  if (attr(struct$txseqs, "start_line") == "second"){
    txseq_dt <- txseq_dt[line != "first"]
  }
  txseq_dt[, state_id := .GRP, by = "line"]
  return(txseq_dt[,])
}