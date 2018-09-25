#' Treatments
#'
#' Dataset of treatments for EGFR positive NSCLC.
#'
#' @format A data.table object with 1 row for each treatment and w columns:
#' \describe{
#'   \item{name}{Name of treatment.}
#'   \item{abb}{Abbreviation for the treatment.}
#' }
"treatments"

#' Multi-state NMA parameters
#'
#' The regression coefficients of the continuous time state transition models
#' estimated using the multi-state network meta-analysis. 
#'
#' @format A list of \code{\link[hesim]{params_surv}} objects from the 
#' \href{https://innovationvalueinitiative.github.io/hesim/}{hesim} package
#' package. The list contains regression coefficient estimates from the Weibull,
#' Gompertz, and 2nd order fractional polynomial survival distributions. The 
#' parameter (e.g., scale and shape for the Weibull distribution) of each survival 
#' distribution are predicted as a function of treatments, health state transitions, and
#' treatment history (for 2L treatments).
"params_mstate_nma"

#' Utility parameters
#'
#' Utility estimates by health state and disutilities by adverse event.
#'
#' @format A list containing the following elements:
#' \itemize{
#' \item{state_utility}{ Utility estimates by health state.}
#' \item{as_disutility}{ Disutility by health state.}
#' }
#' 
#' @section State utility:
#' The \code{state_utility} element is a data table with the following columns:
#' \describe{
#' \item{state_name}{Name of the health state.}
#' \item{mean}{Mean utility.}
#' \item{se}{Standard error of utility.}
#' \item{ref}{BibTeX reference for estimate.}
#' }
#' 
#' @section Adverse event disutilities:
#' The \code{as_disutility} element is a data table with the following columns:
#' \describe{
#' \item{ae_name}{Name of the adverse event.}
#' \item{ae_abb}{Abbreviation for the adverse event.}
#' \item{ae_mean}{Mean disutility.}
#' \item{ae_se}{Standard error of disutility.}
#' \item{ae_ref}{BibTeX reference for estimate.}
#' }
"params_utility"

#' Outpatient cost parameters
#'
#' Outpatient cost estimates by health state.
#'
#' @format A data table with the following columns:
#' \describe{
#' \item{state_name}{Name of the health state.}
#' \item{mean}{Mean outpatient costs.}
#' \item{se}{Standard error of outpatient costs.}
#' \item{ref}{BibTeX reference for estimate.}
#' }
#' @examples 
#' print(params_costs_op)
"params_costs_op"