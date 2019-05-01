#' Compute productivity costs
#' 
#' Compute simulated productivity costs given simulated disease progression including 
#' costs from temporary disability, permanent disability, and premature mortality. 
#' @param econmod An economic model of class \code{"IndivCtstm"}. Disease progression
#' must have been previously simulated (i.e., \code{$disprog_} cannot be \code{NULL}.)
#' @param patients A \code{data.table} returned from \code{\link{create_patients}}.
#' @param params An object of the same format as \code{\link{params_costs_prod}}.
#' @param dr Discount rate.
#' @param retirement_age A numeric scalar denoting age at retirement.
#' @param method Method used to compute productivity losses. Currently only supports
#' the "human capital approach" (\code{hca}).
#' @return An object of class "prod_costs", which is a \code{data.table} 
#' containing productivity costs. Columns are:
#' \describe{
#' \item{category}{The cost category.}
#' \item{dr}{The discount rate.}
#' \item{sample}{\code{sample} from \code{econmod$disprog_}}
#' \item{strategy_id}{\code{strategy_id} from \code{econmod$disprog_}}
#' \item{cost}{Productivity costs.}
#' }
#' @seealso See the example in the \href{https://innovationvalueinitiative.github.io/IVI-NSCLC/articles/tutorial.html}{tutorial}.
#' @export
sim_prod_costs <- function(econmod, patients, params = iviNSCLC::params_costs_prod, 
                       dr = .03, retirement_age = 65, method = "hca"){
  
  # Checks
  check_is_class(econmod, name = "econmod", class = "IndivCtstm")
  check_is_class(patients, name = "patients", class = "patients") 
  
  # Prevent CRAN notes
  female <- employment_status <- gender <- weekly_wage <- prop <- hourly_wage <- NULL
  final <- to <- start_age <- end_age <- time_stop <- NULL 
  prodcost <- prodcost_mortality <- prodcost_permanent <- prodcost_temp <- NULL
  missed_yrs <- category <- NULL
  
  
  # Setting up
  method <- match.arg(method)
  patients <- copy(patients)
  sim <- econmod$disprog_[final == 1]
  n_samples <- max(sim$sample)
  death_state <- econmod$trans_model$death_state
  
  # Wages
  wages <- params$wages
  wages[, female := ifelse(gender == "female", 1, 0)]
  weekly_wage_dt <- wages[, list(weekly_wage = stats::weighted.mean(weekly_wage, w = prop)), by = "female"]
  wages[, hourly_wage := ifelse(employment_status == "full", 
                                weekly_wage/35,
                                ifelse(employment_status == "part", 
                                       weekly_wage/17.5,
                                       0))]
  hourly_wage_dt <- wages[, list(hourly_wage = stats::weighted.mean(hourly_wage, w = prop)), by = "female"]
  patients[, ("weekly_wage") := weekly_wage_dt$weekly_wage[match(patients$female, weekly_wage_dt$female)]]
  patients[, ("hourly_wage") := hourly_wage_dt$hourly_wage[match(patients$female, hourly_wage_dt$female)]]
  
  # Add patient specific data to econmod$disprog_
  indx <- match(sim$patient_id, patients$patient_id)
  sim[, weekly_wage := patients$weekly_wage[indx]]
  sim[, hourly_wage := patients$hourly_wage[indx]]
  sim[, start_age := patients$age[indx]]
  sim[, end_age := start_age + time_stop]
  
  # Productivity costs from premature mortality
  sim[, prodcost_mortality := ifelse(end_age <= 65 & to == death_state, 
                                     pv(z = weekly_wage * 52, r = dr, t1 = time_stop,
                                        t2 = time_stop + (retirement_age - end_age)),
                                     0)]
  
  # Productivity costs from temporary disability
  missed_days <- stats::runif(n_samples,
                              min = params$temporary_disability["missed_days_lower"],
                              max = params$temporary_disability["missed_days_upper"])
  sim[, missed_yrs := missed_days[sample]/365.25]
  sim[, missed_yrs := pmin(missed_yrs, time_stop)]
  sim[, prodcost_temp := pv(weekly_wage * 52, r = dr,
                            t1 = 0, t2 = missed_yrs)]
  
  # Productivity costs from permanent disability
  hours_reduction <- stats::runif(n_samples,
                                  min = params$permanent_disability["hours_reduction_lower"],
                                  max = params$permanent_disability["hours_reduction_upper"])
  sim[, hours_reduction := hours_reduction[sample]]
  sim[, prodcost_permanent := pv(hourly_wage * hours_reduction, r = dr,
                            t1 = missed_yrs, t2 = time_stop)]
  
  # Aggregate
  sim[, prodcost := prodcost_mortality + prodcost_temp + prodcost_permanent]
  object <- sim[, list(costs = mean(prodcost)), by = c("sample", "strategy_id")]
  
  # Return
  object[, category := "prod"]
  dr_env <- dr
  object[, dr := dr_env]
  setorderv(object, c("category", "dr", "sample", "strategy_id", "costs"))
  setattr(object, "class", c("prod_costs", "data.table", "data.frame"))
  return(object[, ])
}

#' Add productivity costs to a cost-effectiveness object
#' 
#' Add productivity costs to a 
#' \href{https://hesim-dev.github.io/hesim/reference/ce.html}{cost-effectiveness object}
#' from the \href{https://hesim-dev.github.io/hesim/index.html}{hesim} pacakge.
#' @param ce An object of class "ce". 
#' @param prod_costs An object of class "prod_costs". 
#' @return Unknown
#' @seealso See the example in the \href{https://innovationvalueinitiative.github.io/IVI-NSCLC/articles/tutorial.html}{tutorial}.
#' @export
add_prod_costs <- function(ce, prod_costs){
  # Checks
  check_is_class(ce, name = "ce", class = "ce")
  check_is_class(prod_costs, name = "prod_costs", class = "prod_costs")
  
  # Prevent CRAN warnings
  costs <- category <- NULL
  
  # Compute
  ce$costs <- ce$costs[category != "total"]
  ce$costs <- rbind(ce$costs, prod_costs)
  costs_total <- ce$costs[, list(costs = sum(costs)), 
                          by = c("dr", "sample", "strategy_id")]
  costs_total[, category := "total"]
  ce$costs <- rbind(ce$costs, costs_total)  
  return(ce)
}